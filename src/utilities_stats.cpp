#include <iostream>

#include <stats.h>

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cerr << "[ERROR] Please provide bam file name and a contig name to use (e.g. chr1)" << std::endl;
        return -1;
    }

    auto stat = stats_insert_size(std::string(argv[1]), std::string(argv[2]));

    std::cout << "Mean insert size: " << std::get<0>(stat) << std::endl
              << "Std  insert size: " << std::get<1>(stat) << std::endl;
}
