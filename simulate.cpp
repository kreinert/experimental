/**
 * simulate.cpp  –  DNA Genome Simulator
 *
 * Generates a random DNA genome of length s, then samples n reads of length l
 * at evenly-spaced positions (step = ceil(s/n)), inserts e substitution errors
 * into each read, and writes the genome and reads as FASTA files.
 *
 * Build:  g++ -std=c++17 -O2 -o simulate simulate.cpp
 *
 * Usage:  ./simulate -s 10000 -n 100 -l 100 -e 2
 *         ./simulate -s 10000 -n 100 -l 100 -e 2 -o reads.fasta --genome genome.fasta --seed 42
 */

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Argument parsing
// ---------------------------------------------------------------------------

struct Args {
    int s = -1;               // genome length
    int n = -1;               // number of reads
    int l = -1;               // read length
    int e = -1;               // errors per read
    std::string output  = "reads.fasta";
    std::string genome  = "genome.fasta";
    unsigned seed       = std::random_device{}();
    bool seed_set       = false;
};

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog
        << " -s <genome_len> -n <num_reads> -l <read_len> -e <errors>"
           " [-o reads.fasta] [--genome genome.fasta] [--seed <uint>]\n";
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
        if      (key == "-s")       a.s      = std::stoi(need("-s"));
        else if (key == "-n")       a.n      = std::stoi(need("-n"));
        else if (key == "-l")       a.l      = std::stoi(need("-l"));
        else if (key == "-e")       a.e      = std::stoi(need("-e"));
        else if (key == "-o")       a.output = need("-o");
        else if (key == "--genome") a.genome = need("--genome");
        else if (key == "--seed") { a.seed   = std::stoul(need("--seed")); a.seed_set = true; }
        else if (key == "-h" || key == "--help") usage(argv[0]);
        else { std::cerr << "Unknown argument: " << key << "\n"; usage(argv[0]); }
    }
    if (a.s < 1 || a.n < 1 || a.l < 1 || a.e < 0) {
        std::cerr << "All of -s, -n, -l, -e are required and must be positive.\n";
        usage(argv[0]);
    }
    if (a.l > a.s) {
        std::cerr << "Read length l=" << a.l << " exceeds genome length s=" << a.s << "\n";
        std::exit(1);
    }
    return a;
}

// ---------------------------------------------------------------------------
// Core helpers
// ---------------------------------------------------------------------------

static const char BASES[4] = {'A', 'C', 'G', 'T'};

static std::string generate_genome(int s, std::mt19937& rng) {
    std::uniform_int_distribution<int> dist(0, 3);
    std::string genome(s, ' ');
    for (char& c : genome)
        c = BASES[dist(rng)];
    return genome;
}

static std::string insert_errors(const std::string& read, int e, std::mt19937& rng) {
    if (e == 0) return read;
    std::string out = read;
    int len = static_cast<int>(out.size());
    int actual_e = std::min(e, len);

    // Choose 'actual_e' distinct positions (partial Fisher-Yates)
    std::vector<int> indices(len);
    std::iota(indices.begin(), indices.end(), 0);
    for (int i = 0; i < actual_e; ++i) {
        std::uniform_int_distribution<int> pick(i, len - 1);
        std::swap(indices[i], indices[pick(rng)]);
    }

    std::uniform_int_distribution<int> alt(0, 2);  // pick from 3 alternatives
    for (int i = 0; i < actual_e; ++i) {
        int pos = indices[i];
        // Find index of original base and pick a different one
        int orig_idx = 0;
        for (int b = 0; b < 4; ++b) if (BASES[b] == out[pos]) { orig_idx = b; break; }
        int shift = alt(rng) + 1;              // 1, 2, or 3
        out[pos] = BASES[(orig_idx + shift) % 4];
    }
    return out;
}

static void write_fasta(const std::string& path,
                        const std::vector<std::pair<std::string, std::string>>& seqs) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open '" + path + "' for writing");
    for (const auto& [hdr, seq] : seqs)
        f << '>' << hdr << '\n' << seq << '\n';
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    Args a = parse_args(argc, argv);

    std::mt19937 rng(a.seed);

    // Generate and save genome
    std::string genome = generate_genome(a.s, rng);
    write_fasta(a.genome, {{"genome length=" + std::to_string(a.s), genome}});
    std::cout << "Genome of length " << a.s << " written to '" << a.genome << "'\n";

    // Sample reads
    int step = static_cast<int>(std::ceil(static_cast<double>(a.s) / a.n));

    std::vector<std::pair<std::string, std::string>> reads;
    reads.reserve(a.n);

    for (int i = 0; i < a.n; ++i) {
        int pos = i * step;
        // Clamp so read stays within genome
        if (pos + a.l > a.s)
            pos = a.s - a.l;

        std::string raw = genome.substr(pos, a.l);
        std::string read = insert_errors(raw, a.e, rng);

        std::string hdr = "read_" + std::to_string(i)
                        + " pos="    + std::to_string(pos)
                        + " len="    + std::to_string(a.l)
                        + " errors=" + std::to_string(a.e);
        reads.emplace_back(std::move(hdr), std::move(read));
    }

    write_fasta(a.output, reads);
    std::cout << a.n << " reads of length " << a.l
              << " (step=" << step << ", errors=" << a.e << ")"
              << " written to '" << a.output << "'\n";

    return 0;
}
