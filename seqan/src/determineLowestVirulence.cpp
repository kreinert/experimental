// ==========================================================================
//                          determineLowestVirulence
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@fu-berlin.de
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/translation.h>

#include <seqan/arg_parse.h>

#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

//#include <lemon/list_graph.h>

#include <gurobi_c++.h>
using namespace std;

using namespace seqan;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

// This struct stores the options from the command line.
//
// You might want to rename this to reflect the name of your app.

struct AppOptions
{
    // Verbosity level.  0 -- quiet, 1 -- normal, 2 -- verbose, 3 -- very verbose.
    int verbosity;

    // The first (and only) argument of the program is stored here.
    String<char> codonTable;
    String<char> genomeFile;
    String<char> probFile;
    String<char> outFileGenome;
    String<char> outFileObjective;

    bool objective; // minimize = 0

    AppOptions() :
        verbosity(1),
        objective(false)
    {}
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function parseCommandLine()
// --------------------------------------------------------------------------

seqan::ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("determineLowestVirulence");
    // Set short description, version, and date.
    setShortDescription(parser, "This program shuffles codons but keeps the amino acid sequence");
    setVersion(parser, "0.1");
    setDate(parser, "Mai 2014");

    // We require one argument.
   // addArgument(parser, seqan::ArgParseArgument(seqan::ArgParseArgument::STRING, "TEXT"));

    addOption(parser, seqan::ArgParseOption("q", "quiet", "Set verbosity to a minimum."));
    addOption(parser, seqan::ArgParseOption("v", "verbose", "Enable verbose output."));
    addOption(parser, seqan::ArgParseOption("vv", "very-verbose", "Enable very verbose output."));

    addOption(parser, seqan::ArgParseOption("ct", "codonTable", "Name of the file the codon table is stored in.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("g", "genomeTable", "Name of the genome file.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("pf", "probFile", "Name of the probability file.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("max", "maximize", "Objective is to maximize"));
    addOption(parser, seqan::ArgParseOption("og", "outFileGenome", "Name of the genome output file.", ArgParseArgument::STRING, "TEXT"));
    addOption(parser, seqan::ArgParseOption("oo", "outFileObjective", "Name of the objective output file.", ArgParseArgument::STRING, "TEXT"));

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    if (isSet(parser, "quiet"))
        options.verbosity = 0;
    if (isSet(parser, "verbose"))
        options.verbosity = 2;
    if (isSet(parser, "very-verbose"))
        options.verbosity = 3;
   // seqan::getArgumentValue(options.text, parser, 0);

    getOptionValue(options.codonTable, parser, "ct");
    getOptionValue(options.genomeFile, parser, "g");
    getOptionValue(options.probFile, parser, "pf");
    getOptionValue(options.outFileGenome, parser, "og");
    getOptionValue(options.outFileObjective, parser, "oo");
    if (isSet(parser, "max"))
        options.objective = true;

    return seqan::ArgumentParser::PARSE_OK;
}

struct AminoAcidToDna_
{
    String<StringSet<String<char> > > table;
    AminoAcidToDna_(AppOptions const & options);
};

AminoAcidToDna_::AminoAcidToDna_(AppOptions const & options)
{
    String<char, MMap<> > codonString;
    open(codonString, toCString(options.codonTable));
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(codonString);

    String<char> codon;
    AminoAcid aminoAcid;
    resize(table, 24);
    while (!atEnd(iter))
    {
        if (value(iter) == '#')
        {
            for (; value(iter) != '\n'; ++iter) ;
            ++iter;
        }
        else
        {
            clear(codon);
            for (; value(iter) != '\t'; ++iter) 
            {
                appendValue(codon, value(iter));
            }
            ++iter;

            aminoAcid = value(iter);
            ++iter;
            ++iter;

            appendValue(table[ordValue(aminoAcid)], codon);
        }
    }

    /*
    for (unsigned i = 0; i < length(table); ++i)
    {
        std::cout << AminoAcid(i) << '\t';
        for (unsigned j = 0; j < length(table[i]); ++j)
            std::cout << toCString(table[i][j]) << '\t';
        std::cout << std::endl;
    }
    */
}

template <typename TSequence>
bool isSameAS(TSequence const & codon1, TSequence const & codon2)
{
   return TranslateTableDnaToAminoAcid_<GeneticCode<Canonical> >::VALUE[ordValue(Dna(codon1[0]))][ordValue(Dna(codon1[1]))][ordValue(Dna(codon1[2]))] == 
          TranslateTableDnaToAminoAcid_<GeneticCode<Canonical> >::VALUE[ordValue(Dna(codon2[0]))][ordValue(Dna(codon2[1]))][ordValue(Dna(codon2[2]))];
}

template <typename TSequence>
unsigned hashCodon(TSequence const & codon)
{
    return ordValue(Dna(codon[2])) + 4 * ordValue(Dna(codon[1])) + 16 * ordValue(Dna(codon[0]));
}

bool initTransitionProb(double (&transitionProb)[64][64], AppOptions const & options)
{
    String<char, MMap<> > inString;
    if (!open(inString, toCString(options.probFile)))
        return false;
    Iterator<String<char, MMap<> >, Rooted>::Type iter = begin(inString);

    String<char> temp;
    String<Dna> tempDna;
    StringSet<String<Dna> > codonOrder;
    while (!atEnd(iter))
    {
        clear(temp);
        clear(tempDna);
        if (value(iter) == '#')
        {
            for (; value(iter) != '\n'; ++iter) ;
            ++iter;
        }
        else if (value(iter) == ';')
        {
            ++iter;
            for (; value(iter) != '\n'; ++iter)
            {
                if(value(iter) == ';')
                {
                    appendValue(codonOrder, tempDna);
                    clear(tempDna);
                }
                else

            appendValue(codonOrder, tempDna);
            ++iter;
        }
        else
        {
            appendValue(tempDna, value(iter));
            ++iter;
            appendValue(tempDna, value(iter));
            ++iter;
            appendValue(tempDna, value(iter));
            ++iter;

            // ;
            ++iter;
            unsigned codonPos = hashCodon(tempDna);
            for (unsigned i = 0; i < 64; ++i)
            {
                clear(temp);
                for (; value(iter) != ';' && value(iter) != '\n'; ++iter)
                {
                    appendValue(temp, value(iter));
                }
                ++iter;
                lexicalCast2(transitionProb[codonPos][hashCodon(codonOrder[i])], temp);
            }
        }
    }
    return true;
}

template<typename TValue>
GRBLinExpr comp(TValue x, TValue y, TValue z)
{
    return x + y + 2 * z;
}

template<typename TValue>
GRBLinExpr coloumnsConstraints(String<TValue> const & columns)
{
    GRBLinExpr temp = 0.0;
    for (unsigned i = 0; i < length(columns); ++i)
        temp += columns[i];
    return temp;
}

template<typename TModel, typename TValue>
void edgeConstraints(TModel & model, StringSet<String<TValue> > const & columns, StringSet<String<TValue> > const & edges, unsigned column, bool left)
{
    GRBQuadExpr temp = 0.0;
    if (left)
    {
        for (unsigned i = 0; i < length(columns[column]); ++i)
            for (unsigned j = 0; j < length(columns[column + 1]); ++j)
                model.addConstr(edges[column][length(columns[column + 1]) * i + j] <= columns[column][i]);

//                temp += columns[column][i] * edges[column][length(columns[column + 1]) * i + j];
    }
    else 
    {
        for (unsigned j = 0; j < length(columns[column]); ++j)
            for (unsigned i = 0; i < length(columns[column - 1]); ++i)
                model.addConstr(columns[column][j] >= edges[column - 1][length(columns[column]) * i + j]);
                //temp += columns[column][j] * edges[column - 1][length(columns[column]) * i + j];
    }
    //return temp;
}


template<typename TValue, typename TValue2, typename TSequence>
GRBLinExpr objective(double (&transitionProb)[64][64], StringSet<String<TValue> > const & edges, StringSet<String<TValue2> > const & ids, TSequence const & seq)
{
    GRBLinExpr temp = 0.0;
    for (unsigned i = 1; i < length(ids); ++i)
        for (unsigned j = 0; j < length(ids[i - 1]); ++j)
            for (unsigned k = 0; k < length(ids[i]); ++k)
            {
                temp += edges[i - 1][length(ids[i]) * j + k] * transitionProb[hashCodon(infix(seq, ids[i - 1][j], ids[i - 1][j] + 3))][hashCodon(infix(seq, ids[i][k], ids[i][k] + 3))];
            }
    return temp;
}


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
    // Parse the command line.
    seqan::ArgumentParser parser;
    AppOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::cout << "determineLowestVirulence\n"
              << "========================\n\n";

    // reading and creating amino acid table
    AminoAcidToDna_ _table(options);

    typedef double TCargo;
    typedef Graph<Directed<TCargo> > TGraph;
    typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef String<String<char> > TProperties;


    
    // reading codon transition probabilities
    double transitionProb[64][64];
    initTransitionProb(transitionProb, options);

    String<char, MMap<> > outString;
    open(outString, toCString(options.outFileGenome));
    
    String<char, MMap<> > outStringObjective;
    open(outStringObjective, toCString(options.outFileObjective));

    SequenceStream seqIO(toCString(options.genomeFile));
    while (!atEnd(seqIO))
    {
        // Reading genome sequence.
        CharString id;
        Dna5String seq;

        int resRead = readRecord(id, seq, seqIO);
        if (resRead != 0)
        {
            std::cerr << "ERROR: Could not read record!\n";
            return 1;
        }

        if (length(seq) > 3000)
        {
	    
            /*
	    This is the exact solution obtained by a graph search.

            TProperties codons;
            String<AminoAcid> text;
            translate(text, seq);

            StringSet<String<TVertexDescriptor> > vertices;
            resize(vertices, length(text) + 2);
            resize(vertices[0], 1);
            resize(vertices[length(vertices) - 1], 1);

            TGraph g;

            vertices[0][0] = addVertex(g);
            vertices[length(vertices) - 1][0] = addVertex(g);

            for (unsigned i = 0; i < length(text); ++i)
            {
                resize(vertices[i + 1], length(_table.table[ordValue(text[i])]));
                for (unsigned j = 0; j < length(vertices[i + 1]); ++j)
                    vertices[i + 1][j] = addVertex(g);
            }


            resizeVertexMap(g, codons);
            assignProperty(codons, vertices[0][0], "start");
            assignProperty(codons, vertices[length(vertices) - 1][0], "end");
            for (unsigned i = 0; i < length(text); ++i)
            {
                for (unsigned j = 0; j < length(vertices[i + 1]); ++j)
                {
                    assignProperty(codons, vertices[i + 1][j], _table.table[ordValue(text[i])][j]);
                    for (unsigned k = 0; k < length(vertices[i]); ++k)
                    {
                        if (i == 0)
                            addEdge(g, vertices[i][k], vertices[i+1][j], 0.0);
                        else
                        {
                            addEdge(g, vertices[i][k], vertices[i+1][j], transitionProb[hashCodon(getProperty(codons, vertices[i][k]))][hashCodon(getProperty(codons,vertices[i+1][j]))]);
                        }
                    }
                }
            }

            for (unsigned i = 0; i < length(vertices[length(text)]) ; ++i)
                addEdge(g, vertices[length(text)][i], vertices[length(text) + 1][0], 0.0);

            // Run Dijkstra's algorithm from vertex 0.
            String<unsigned> predMap;
            String<double> distMap;
            InternalMap<TCargo> cargoMap;
            bellmanFordAlgorithm(g, 0, cargoMap, predMap, distMap);

            // Print results to stdout.
            std::cout << "The computed sequence is: \n";
            _printPath(g,predMap,(TVertexDescriptor) 0, vertices[length(vertices) - 1][0], codons);
            std::cout << std::endl;

            std::cout << "The score is: " << distMap[1] / (double)(length(text) - 1) << std::endl;
        */

	    // ILP
            try {

                std::cout << "====================== Starting with: ";
                std::cout << id;
                std::cout << " ======================" << std::endl;

                GRBEnv env = GRBEnv();

                GRBModel model = GRBModel(env);
                //setIntFeasTol(0.000000001);
        //        model.getEnv().set(GRB_DoubleParam_IntFeasTol, 1e-9); // C++

                // Create variables

                StringSet<String<GRBVar> > grbVarsColumns;
                StringSet<String<GRBVar> > grbVarsRows;
                StringSet<String<unsigned> > ids;

                resize(grbVarsColumns, length(seq) / 3);
                resize(ids, length(seq) / 3);
                resize(grbVarsRows,length(seq) / 3);
                GRBVar temp;
                for (unsigned i = 0; i < length(seq); i += 3)
                {
                    for (unsigned j = 0; j < length(seq); j += 3)
                    {
                        if ( ( (i == j) || (infix(seq, i, i + 3) != infix(seq, j, j + 3) ) ) && isSameAS(infix(seq, i, i + 3), infix(seq, j, j+ 3)))
                        {
                            if (i == j)
                                temp = model.addVar(0.0, 1.0, 1.0, GRB_BINARY);
                            else
                                temp = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                            appendValue(grbVarsColumns[i/3], temp);
                            appendValue(ids[i/3], j);
                            appendValue(grbVarsRows[j / 3], temp);
                        }
                    }
                }

                StringSet<String<GRBVar> > grbVarsEdges;
                resize(grbVarsEdges, length(ids) - 1);
                for (unsigned i = 0; i < length(grbVarsEdges); ++i)
                {
                    for (unsigned j = 0; j < length(ids[i]); ++j)
                        for (unsigned k = 0; k < length(ids[i + 1]); ++k)
                            appendValue(grbVarsEdges[i], model.addVar(0.0, 1.0, 0.0, GRB_BINARY));
                }

                // Integrate new variables
                //model.getEnv().set("IntFeasTol", 1e-9);
                model.update();

                if (!options.objective)
                    model.setObjective(objective(transitionProb,grbVarsEdges, ids, seq), GRB_MINIMIZE);
                else 
                    model.setObjective(objective(transitionProb,grbVarsEdges, ids, seq), GRB_MAXIMIZE);

                // Add constraint
                for (unsigned i = 0; i < length(grbVarsColumns); ++i)
                    model.addConstr(coloumnsConstraints(grbVarsColumns[i]) == 1);

                for (unsigned i = 0; i < length(grbVarsRows); ++i)
                    model.addConstr(coloumnsConstraints(grbVarsRows[i]) == 1);

                for (unsigned i = 0; i < length(grbVarsEdges); ++i)
                    edgeConstraints(model, grbVarsColumns, grbVarsEdges, i, true);

                for (unsigned i = 0; i < length(grbVarsEdges); ++i)
                    edgeConstraints(model, grbVarsColumns, grbVarsEdges, i + 1, false);

                for (unsigned i = 0; i < length(grbVarsEdges); ++i)
                    model.addConstr(coloumnsConstraints(grbVarsEdges[i]) == 1);


                // Optimize model
                model.optimize();

                std::ostringstream sstream;
                sstream << model.get(GRB_DoubleAttr_ObjVal) / ((double)length(seq) / 3.0 - 1.0);

                append(outStringObjective, id);
                append(outStringObjective, "\t:\t");
                append(outStringObjective, sstream.str());
                append(outStringObjective, "\n");

                cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) / ((double)length(seq) / 3.0 - 1.0) << endl;

                appendValue(outString, '>');
                append(outString, id);
                appendValue(outString, '\n');

                // output
                for (unsigned i = 0; i < length(grbVarsColumns) - 1; ++i)
                    for (unsigned j = 0; j < length(grbVarsColumns[i]); ++j)
                        {
                            if (grbVarsColumns[i][j].get(GRB_DoubleAttr_X) != 0)
                            {
                                append(outString, infix(seq, ids[i][j], ids[i][j] + 3));
                            }
                        }

                // output of last column
                for (unsigned j = 0; j < length(grbVarsColumns[length(grbVarsColumns) - 1]); ++j)
                {
                    if (grbVarsColumns[length(grbVarsColumns) - 1][j].get(GRB_DoubleAttr_X) != 0)
                    {
                        append(outString, infix(seq, ids[length(grbVarsColumns) - 1][j], ids[length(grbVarsColumns) - 1][j] + 3));
                    }
                }
                appendValue(outString, '\n');

            } catch(GRBException e) 
            {
                    cout << "Error code = " << e.getErrorCode() << endl;
                    cout << e.getMessage() << endl;
            } catch(...) 
            {
                cout << "Exception during optimization" << endl;
            }
        }
    }

    return 0;
}

