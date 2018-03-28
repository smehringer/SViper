#include <iostream>

#include <stats.h>

int main(int argc, char ** argv)
{
    if (argc != 2)
        std::cerr << "[ERROR] Please provide bam file name" << std::endl;

    auto stat = stats_insert_size(std::string(argv[1]));

    std::cout << "Mean insert size: " << get<0>(stat) << std::endl
              << "Std  insert size: " << get<1>(stat) << std::endl;
}
