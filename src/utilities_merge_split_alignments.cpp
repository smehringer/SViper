#include <iostream>
#include <seqan/bam_io.h>

#include <merge_split_alignments.h>

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: file.sortedByName.sam";
        return 1;
    }

    seqan::BamFileIn bamfileIn;
    seqan::BamHeader header;

    if (!open(bamfileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << std::endl;
        return 1;
    }

    seqan::readHeader(header, bamfileIn);
    seqan::BamFileOut bamfileOut(seqan::context(bamfileIn));
    std::string out_filename = (std::string(argv[1]) + ".merged.sam");

    if (!seqan::open(bamfileOut, out_filename.c_str()))
    {
        std::cerr << "ERROR: Could not open " << out_filename << std::endl;
        return 1;
    }

    seqan::writeHeader(bamfileOut, header);

    // empty file must be check here otherwise the first read record will fail
    if (seqan::atEnd(bamfileIn))
        return 0;

    seqan::BamAlignmentRecord record;
    std::vector<seqan::BamAlignmentRecord> record_group;

    while (!atEnd(bamfileIn))
    {
        seqan::readRecord(record, bamfileIn);

        if (!record_group.empty() &&
            (record_group[record_group.size() - 1]).qName != record.qName)
        {
            seqan::BamAlignmentRecord merged_record = seqan::merge_record_group(record_group);
            seqan::writeRecord(bamfileOut, merged_record);
            record_group.clear();
            record_group.push_back(record);
        }
        else
        {
            record_group.push_back(record);
        }
    }
    //process last group
    seqan::BamAlignmentRecord merged_record = seqan::merge_record_group(record_group);
    seqan::writeRecord(bamfileOut, merged_record);

    return 0;
}
