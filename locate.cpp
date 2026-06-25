/**
 * locate.cpp  –  K-mer Read Locator
 *
 * Reads a reference genome FASTA and a reads FASTA (produced by simulate.cpp).
 * For each read it calls locate(read) which returns every genome position at
 * which any k-mer of the read occurs.  Each hit is then classified:
 *
 *   in_interval  – hit position falls within [orig_pos, orig_pos + len - 1]
 *                  (the interval the read was originally sampled from)
 *   out_interval – hit position is outside that window
 *
 * Build:  g++ -std=c++17 -O2 -o locate locate.cpp
 *
 * Usage:  ./locate -r reads.fasta -g genome.fasta -k 15
 *         ./locate -r reads.fasta -g genome.fasta -k 15 --verbose
 */

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
// Argument parsing
// ---------------------------------------------------------------------------

struct Args {
    std::string reads;
    std::string genome;
    int k        = -1;
    bool verbose = false;
};

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " -r <reads.fasta> -g <genome.fasta> -k <kmer_size> [--verbose]\n";
    std::exit(1);
}

static Args parse_args(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string key = argv[i];
        auto need = [&](const char* name) -> std::string {
            if (++i >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                std::exit(1);
            }
            return argv[i];
        };
        if      (key == "-r")        a.reads   = need("-r");
        else if (key == "-g")        a.genome  = need("-g");
        else if (key == "-k")        a.k       = std::stoi(need("-k"));
        else if (key == "--verbose") a.verbose = true;
        else if (key == "-h" || key == "--help") usage(argv[0]);
        else { std::cerr << "Unknown argument: " << key << "\n"; usage(argv[0]); }
    }
    if (a.reads.empty() || a.genome.empty() || a.k < 1) {
        std::cerr << "Arguments -r, -g, and -k are required (k >= 1).\n";
        usage(argv[0]);
    }
    return a;
}

// ---------------------------------------------------------------------------
// FASTA parsing
// ---------------------------------------------------------------------------

struct FastaRecord {
    std::string header;   // without leading '>'
    std::string sequence;
};

static std::vector<FastaRecord> parse_fasta(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open '" + path + "'");

    std::vector<FastaRecord> records;
    std::string line;
    FastaRecord cur;

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur.header.empty())
                records.push_back(std::move(cur));
            cur = FastaRecord{line.substr(1), ""};
        } else {
            cur.sequence += line;
        }
    }
    if (!cur.header.empty())
        records.push_back(std::move(cur));

    return records;
}

// ---------------------------------------------------------------------------
// Header metadata parsing
// ---------------------------------------------------------------------------

static bool parse_int_tag(const std::string& header,
                           const std::string& tag, int& out) {
    auto pos = header.find(tag);
    if (pos == std::string::npos) return false;
    pos += tag.size();
    try {
        std::size_t len = 0;
        out = std::stoi(header.substr(pos), &len);
        return len > 0;
    } catch (...) { return false; }
}

// ---------------------------------------------------------------------------
// K-mer index
// ---------------------------------------------------------------------------

using KmerIndex = std::unordered_map<std::string, std::vector<int>>;

static KmerIndex build_kmer_index(const std::string& genome, int k) {
    KmerIndex idx;
    int n = static_cast<int>(genome.size());
    idx.reserve(n);
    for (int i = 0; i <= n - k; ++i)
        idx[genome.substr(i, k)].push_back(i);
    return idx;
}

// ---------------------------------------------------------------------------
// locate
// ---------------------------------------------------------------------------

/**
 * Returns the sorted, deduplicated list of genome positions at which any
 * k-mer of `read` occurs according to `kmer_index`.
 */
static std::vector<int> locate(const std::string& read,
                                int k,
                                const KmerIndex& kmer_index) {
    std::vector<int> hits;
    int n = static_cast<int>(read.size());
    for (int i = 0; i <= n - k; ++i) {
        auto it = kmer_index.find(read.substr(i, k));
        if (it != kmer_index.end())
            hits.insert(hits.end(), it->second.begin(), it->second.end());
    }
    // Deduplicate
    std::sort(hits.begin(), hits.end());
    hits.erase(std::unique(hits.begin(), hits.end()), hits.end());
    return hits;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    Args a = parse_args(argc, argv);

    // Load genome
    auto genome_records = parse_fasta(a.genome);
    if (genome_records.empty()) {
        std::cerr << "ERROR: no sequences found in '" << a.genome << "'\n";
        return 1;
    }
    const std::string& genome = genome_records[0].sequence;
    std::cout << "Genome length : " << genome.size() << "\n";

    if (a.k > static_cast<int>(genome.size())) {
        std::cerr << "ERROR: k=" << a.k << " is larger than the genome ("
                  << genome.size() << ")\n";
        return 1;
    }

    // Build k-mer index
    std::cout << "Building " << a.k << "-mer index …\n";
    KmerIndex kmer_index = build_kmer_index(genome, a.k);
    std::cout << "Distinct " << a.k << "-mers in genome: "
              << kmer_index.size() << "\n";

    // Load reads
    auto reads = parse_fasta(a.reads);
    std::cout << "Reads to process: " << reads.size() << "\n\n";

    long long total_in      = 0;
    long long total_out     = 0;
    int reads_no_origin     = 0;
    int reads_with_origin   = 0;
    int intervals_with_hit  = 0;   // reads where at least one k-mer hit lands in the interval

    for (const auto& rec : reads) {
        int orig_pos = -1, read_len = -1;
        bool has_origin = parse_int_tag(rec.header, "pos=",  orig_pos)
                       && parse_int_tag(rec.header, "len=",  read_len);

        std::vector<int> hits = locate(rec.sequence, a.k, kmer_index);

        if (!has_origin) {
            ++reads_no_origin;
            total_out += static_cast<long long>(hits.size());
            if (a.verbose)
                std::cout << "  [" << rec.header << "]  no origin info — "
                          << hits.size() << " hits (all out)\n";
            continue;
        }

        ++reads_with_origin;
        int interval_end = orig_pos + read_len - 1;
        long long r_in  = 0;
        for (int p : hits)
            if (p >= orig_pos && p <= interval_end) ++r_in;
        long long r_out = static_cast<long long>(hits.size()) - r_in;

        total_in  += r_in;
        total_out += r_out;
        if (r_in > 0) ++intervals_with_hit;

        if (a.verbose) {
            std::cout << "  [" << rec.header << "]"
                      << "  interval=[" << orig_pos << "," << interval_end << "]"
                      << "  hits="  << hits.size()
                      << "  in="    << r_in
                      << "  out="   << r_out << "\n";
        }
    }

    // Summary
    long long total = total_in + total_out;
    auto pct_hits = [&](long long x) -> std::string {
        if (total == 0) return "";
        std::ostringstream os;
        os.precision(1);
        os << std::fixed << "  (" << (100.0 * x / total) << "%)";
        return os.str();
    };
    auto pct_reads = [&](int x) -> std::string {
        if (reads_with_origin == 0) return "";
        std::ostringstream os;
        os.precision(1);
        os << std::fixed << "  (" << (100.0 * x / reads_with_origin) << "%)";
        return os.str();
    };

    std::string sep(58, '=');
    std::cout << sep << "\n"
              << "k-mer size                              : " << a.k << "\n"
              << "\n"
              << "K-mer hit statistics\n"
              << "  Total hit positions examined          : " << total << "\n"
              << "  Within original read interval         : " << total_in  << pct_hits(total_in)  << "\n"
              << "  Outside original read interval        : " << total_out << pct_hits(total_out) << "\n"
              << "\n"
              << "Interval coverage\n"
              << "  Reads with known origin               : " << reads_with_origin << "\n"
              << "  Intervals with >= 1 k-mer hit         : " << intervals_with_hit
                                                              << pct_reads(intervals_with_hit) << "\n"
              << "  Intervals with no k-mer hit           : " << (reads_with_origin - intervals_with_hit)
                                                              << pct_reads(reads_with_origin - intervals_with_hit) << "\n";
    if (reads_no_origin)
        std::cout << "  (reads without origin info            : " << reads_no_origin << ")\n";
    std::cout << sep << "\n";

    return 0;
}
