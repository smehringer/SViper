#include <iostream>
#include <fstream>

#include <basics.h>
#include <variant.h>

int main(int /*argc*/, const char ** argv)
{
    seqan::BamFileIn bam;
    seqan::open(bam, argv[1]);

    std::ifstream input_vcf(argv[2]);

    std::string line{"#"}; // TODO:: use seqan vcf parser instead (needs to be extended)
    while (line[0] == '#' && getline(input_vcf, line)) {}

    Variant var(line);

    while (!seqan::atEnd(bam))
    {
        seqan::BamAlignmentRecord record;
        seqan::readRecord(record, bam);
        std::cout << record.qName << " is supporting: " << record_supports_variant(record, var) << std::endl;
    }
}
