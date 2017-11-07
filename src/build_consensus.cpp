#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: " << argv[0] << " reads.fa\n";
        return 1;
    }

    SeqFileIn seqFileIn;
    StringSet<CharString> ids;
    StringSet<DnaString> seqs;

    if (!open(seqFileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not load reads " << argv[1] << std::endl;
        return 1;
    }

    readRecords(ids, seqs, seqFileIn);

    if (length(seqs) == 0)
    {
        std::cerr << "ERROR: No reads to process. " << std::endl;
        return 1;
    }


    Align<DnaString> align;
    resize(rows(align), length(seqs));
    for (unsigned i = 0; i < length(seqs); ++i)
        assignSource(row(align, i), seqs[i]);

    globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));

    // create the profile string

    String<ProfileChar<Dna> > profile;
    resize(profile, length(row(align, 0)));
    for (unsigned rowNo = 0; rowNo < 4u; ++rowNo)
        for (unsigned i = 0; i < length(row(align, rowNo)); ++i)
            profile[i].count[ordValue(getValue(row(align, rowNo), i))] += 1;

    for (auto p : profile[0])
        std::cout << p;
    std::cout << std::endl;

    // call consensus from this string
    DnaString consensus;
    for (unsigned i = 0; i < length(profile); ++i)
    {
        int idx = getMaxIndex(profile[i]);
        if (idx < 4)  // is not gap
            appendValue(consensus, Dna(getMaxIndex(profile[i])));
    }

    std::cout << "> consensus\n"
              << consensus << "\n";

    return 0;
}
