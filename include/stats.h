#pragma once
#include <numeric>
#include <experimental/random>

#include <seqan/bam_io.h>

#include <basics.h>

std::pair<double, double> stats_insert_size(std::string const & bam_file_name, std::string contig_name)
{
    seqan::BamFileIn      bam_file;
    seqan::BamIndex<Bai>  bam_index;

    if (!open_file_success(bam_file, bam_file_name.c_str()))
        return {-1.0, -1.0};
    if (!open_file_success(bam_index, (bam_file_name + ".bai").c_str()))
        return {-1.0, -1.0};

    seqan::BamHeader bam_header;
    seqan::readHeader(bam_header, bam_file);

    std::vector<double> insert_sizes;
    seqan::BamAlignmentRecord record;
    bool hasAlignments = false;

    int rID = 0;
    if (!getIdByName(rID, contigNamesCache(context(bam_file)), contig_name))
        return {-1.0, -1.0};

    if (!jumpToRegion(bam_file, hasAlignments, rID, 7100001, 12100002, bam_index))
        return {-1.0, -1.0};
    if (!hasAlignments)
        return {-1.0, -1.0};

    for (unsigned i = 7100001; i < 117200000; ++i)
    {
        //int random_pos = std::experimental::randint(7100001, 117200000); // chr1 [p36.2-p12)

        seqan::readRecord(record, bam_file);

        if (seqan::hasFlagAllProper(record) && !seqan::hasFlagRC(record)) // only count pairs F1R2 since they are easy to compute
        {
            double sz = record.pNext - (record.beginPos /*+ seqan::getAlignmentLengthInRef(record) - 1*/);
            if (sz > 0 && sz < 1300)
                insert_sizes.push_back(sz);
        }
    }
    std::cout << std::endl;
    double sum = std::accumulate(insert_sizes.begin(), insert_sizes.end(), 0.0);
    double mean = sum / insert_sizes.size();

    std::vector<double> diff(insert_sizes.size());
    std::transform(insert_sizes.begin(), insert_sizes.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / insert_sizes.size());

    std::sort(insert_sizes.begin(), insert_sizes.end());

    return {insert_sizes[insert_sizes.size()/2], stdev};
}
