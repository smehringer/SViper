#include <iostream>
#include <sstream>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    if (argc != 4)
    {
        std::cerr << "USAGE: " << argv[0] << " alignment.fa REF START_POS\n";
        return 1;
    }

    SeqFileIn seqFileIn;
    CharString read_id;
    CharString ref_id;
    CharString read_seq;
    CharString ref_seq;
    int mapping_pos = atoi(argv[3]);

    if (!open(seqFileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not load reads " << argv[1] << std::endl;
        return 1;
    }

    if (atEnd(seqFileIn))
    {
        std::cerr << "ERROR: No alignments to process. " << std::endl;
        return 1;
    }
    readRecord(read_id, read_seq, seqFileIn);

    if (atEnd(seqFileIn))
    {
        std::cerr << "ERROR: Only one sequence in input. " << std::endl;
        return 1;
    }
    readRecord(ref_id, ref_seq, seqFileIn);

    if (length(read_seq) != length(ref_seq))
    {
        std::cerr << "ERROR: aligned sequences must have the same length. " << std::endl;
        return 1;
    }

    unsigned i = 0;

    while (read_seq[i] == '-')
    {
        ++mapping_pos;
        ++i;
    }

    unsigned len{0};

    // build cigar
    unsigned qual{60};
    unsigned count{1};
    char operation{'M'};
    char tmp_operation{'M'};

    std::ostringstream cigar;
    std::ostringstream seq;
    seq << read_seq[i];
    ++i;
    while (true)
    {
        if (read_seq[i] == '-')
        {
            tmp_operation = 'D';
        }
        else if (ref_seq[i] == '-')
        {
            tmp_operation = 'I';
            seq << read_seq[i];
        }
        else
        {
            tmp_operation = 'M';
            seq << read_seq[i];
        }

        if (tmp_operation != operation)
        {
            if (operation != 'D')
                len += count;
            cigar << count;
            cigar << operation;
            count = 0;
        }

        ++count;

        if (i == length(read_seq) - 1)
        {
            if (tmp_operation != 'D')
            {
                len += count;
                cigar << count;
                cigar << tmp_operation;
            }
            break;
        }

        ++i;
        operation = tmp_operation;
    }

    if ((seq.str()).size() != len)
        std::cerr << "[ERROR] CIGAR and sequence length don't match: "
                  << (seq.str()).size() << "(seq length) != "
                  << len
                  << "(cigar length)." << std::endl;

    std::cout << read_id                    // 1 QNAME
              << "\t" << 0                  // 2 FLAG
              << "\t" << argv[2]            // 3 RNAME
              << "\t" << mapping_pos        // 4 POS
              << "\t" << qual               // 5 MAPQ
              << "\t" << cigar.str()        // 6 CIGAR
              << "\t" << "*"                // 7 RNEXT
              << "\t" << 0                  // 8 PNEXT
              << "\t" << (seq.str()).size() // 9 TLEN
              << "\t" << seq.str()          // 10 SEQ
              << "\t" << "*"                // 11 QUAL
              << std::endl;

    return 0;
}
