#include <iostream>
#include <tuple>
#include <seqan/bam_io.h>

#include <helper_functions.h>

int main(int argc, char ** argv)
{
	if (argc != 5)
	{
		std::cerr << "USAGE: bamfile";
		return 1;
	}

    seqan::BamFileIn bamfileIn;

    if (!seqan::open(bamfileIn, argv[1]))
    {
        std::cerr << "ERROR: Could not open " << argv[1] << std::endl;
        return 1;
    }

    seqan::BamAlignmentRecord record;

    if (!seqan::atEnd(bamfileIn))
        seqan::readRecord(record, bamfileIn);

    unsigned x1 = 0;
    unsigned x2 = 300;

    auto tup = get_read_region_boundaries(record, x1, x2);

	std::cout << std::get<0>(tup) << " " << std::get<1>(tup) << endl;
	return 0;
}
