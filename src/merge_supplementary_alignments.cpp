#include <iostream>
#include <seqan/bam_io.h>

#include <helper_functions.h>
#include <merge_split_alignments.h>

using namespace std;
using namespace seqan;

int main(int argc, char ** argv)
{
    if (argc != 2)
    {
        std::cerr << "USAGE: file.sortedByName.sam";
        return 1;
    }

    BamFileIn bamfileIn;
    BamHeader header;

    if (!open(bamfileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << std::endl;
        return 1;
    }

    readHeader(header, bamfileIn);
    BamFileOut bamfileOut(context(bamfileIn));
    std::string out_filename = (std::string(argv[1]) + ".merged.bam");

    if (!open(bamfileOut, out_filename.c_str()))
    {
        std::cerr << "ERROR: Could not open " << out_filename << std::endl;
        return 1;
    }

    writeHeader(bamfileOut, header);

    // empty file must be check here otherwise the first read record will fail
    if (atEnd(bamfileIn))
        return 0;

    BamAlignmentRecord record;
    vector<BamAlignmentRecord> record_group;

    while (!atEnd(bamfileIn))
    {
        readRecord(record, bamfileIn);

        if (!record_group.empty() &&
            (record_group[record_group.size() - 1]).qName != record.qName)
        {
            BamAlignmentRecord merged_record = merge_record_group(record_group);
            writeRecord(bamfileOut, merged_record);
            record_group.clear();
            record_group.push_back(record);
        }
        else
        {
            record_group.push_back(record);
        }
    }
    //process last group
    BamAlignmentRecord merged_record = merge_record_group(record_group);
    writeRecord(bamfileOut, merged_record);

    return 0;
}
