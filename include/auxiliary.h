#pragma once
#include <basics.h>
#include <config.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <iostream>
#include <memory>
#include <thread>

struct Auxiliary{
    std::ofstream                                  log_file;
    CmdOptions                                     options;
    std::vector<std::unique_ptr<seqan::BamFileIn>> long_read_file_handles;
    seqan::BamHeader                               long_read_header;   // The bam header object needed to fill bam context
    seqan::BamIndex<seqan::Bai>                    long_read_bai;
    std::vector<std::unique_ptr<seqan::BamFileIn>> short_read_file_handles;
    seqan::BamHeader                               short_read_header;  // The bam header object needed to fill bam context
    seqan::BamIndex<seqan::Bai>                    short_read_bai;
    std::vector<std::unique_ptr<seqan::FaiIndex>>  faidx_file_handles;
    std::vector<seqan::BamAlignmentRecord>         polished_reads; // stores records in case info.options.output-polished-bam is true

    Auxiliary(CmdOptions & options_) : options(options_) {}

    Auxiliary() = default;
    Auxiliary(const Auxiliary&) = default;
    Auxiliary(Auxiliary&&) = default;
    Auxiliary& operator=(const Auxiliary&) = default;
    Auxiliary& operator=(Auxiliary&&) = default;
};

bool prep_file_handles(Auxiliary & info)
{
    unsigned num_threads{info.options.threads};
    omp_set_num_threads(num_threads);

    info.long_read_file_handles.resize(num_threads);
    info.short_read_file_handles.resize(num_threads);
    info.faidx_file_handles.resize(num_threads);

    for (unsigned t = 0; t < num_threads; ++t)
    {
        try
        {
            info.short_read_file_handles[t] = std::make_unique<seqan::BamFileIn>(info.options.short_read_file_name.c_str());
            seqan::readHeader(info.short_read_header, *(info.short_read_file_handles[t]));
        }
        catch (seqan::Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted bam file " << info.options.short_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return false;
        }

        try
        {
            info.long_read_file_handles[t] = std::make_unique<seqan::BamFileIn>(info.options.long_read_file_name.c_str());
            seqan::readHeader(info.long_read_header, *(info.long_read_file_handles[t]));
        }
        catch (seqan::Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted bam file " << info.options.long_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return false;
        }

        try
        {
            info.faidx_file_handles[t] = std::make_unique<seqan::FaiIndex>();
            if (!open_file_success(*(info.faidx_file_handles[t]), info.options.reference_file_name.c_str()))
                return false;
        }
        catch (seqan::Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted faidx file " << info.options.long_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return false;
        }
    }

    return true;
}
