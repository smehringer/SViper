#include <iostream>
#include <seqan/bam_io.h>

#include <helper_functions.h>

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
    BamFileOut bamfileOut(context(bamfileIn));
    std::string out_filename = (std::string(argv[1]) + ".merged.sam");

    if (!open(bamfileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << std::endl;
        return 1;
    }

    if (!open(bamfileOut, out_filename.c_str()))
    {
        std::cerr << "ERROR: Could not open " << out_filename << std::endl;
        return 1;
    }

    // empty file must be check here otherwise the first read record will fail
    if (atEnd(bamfileIn))
        return 0;

    readHeader(header, bamfileIn);
    writeHeader(bamfileOut, header);

    BamAlignmentRecord record;
    vector<BamAlignmentRecord> record_group;

    while (!atEnd(bamfileIn))
    {
        readRecord(record, bamfileIn);

        if (!record_group.empty() &&
            (record_group[record_group.size() - 1]).qName != record.qName)
        {
            BamAlignmentRecord merged_record = process_record_group(record_group);
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
    BamAlignmentRecord merged_record = process_record_group(record_group);
    writeRecord(bamfileOut, merged_record);

    return 0;
}
