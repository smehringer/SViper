#include <sviper.h>

using namespace std;
using namespace seqan;

// Global file info struct.
file_info * info{};

int main(int argc, char const ** argv)
{
    // Parse Command Line Arguments
    // -------------------------------------------------------------------------
    // CmdOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(info->options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check files
    if (!open_file_success(info->long_read_bai, (info->options.long_read_file_name + ".bai").c_str()) ||
        !open_file_success(info->short_read_bai, (info->options.short_read_file_name + ".bai").c_str()) ||
        !open_file_success(info->log_file, (info->options.output_prefix + ".log").c_str()))
        return 1;
    // -------------------------------------------------------------------------
    // Read variants into container // TODO:: use seqan vcf parser instead (needs to be extended)
    // -------------------------------------------------------------------------
    std::vector<std::string> vcf_header;
    std::vector<Variant> variants;
    if (!read_vcf(variants, vcf_header, info))
        return 1;

    // Prepare file hangles for parallel computing
    // -------------------------------------------------------------------------
    unsigned num_threads{info->options.threads};
    omp_set_num_threads(num_threads);

    // std::vector<std::unique_ptr<BamFileIn>> long_read_file_handles;
    // std::vector<std::unique_ptr<BamFileIn>> short_read_file_handles;
    // std::vector<std::unique_ptr<FaiIndex>> faidx_file_handles;

    info->long_read_file_handles.resize(num_threads);
    info->short_read_file_handles.resize(num_threads);
    info->faidx_file_handles.resize(num_threads);

    for (unsigned t = 0; t < num_threads; ++t)
    {
        try
        {
            info->short_read_file_handles[t] = make_unique<BamFileIn>(info->options.short_read_file_name.c_str());
            readHeader(info->short_read_header, *(info->short_read_file_handles[t]));
        }
        catch (Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted bam file " << info->options.short_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return 1;
        }

        try
        {
            info->long_read_file_handles[t] = make_unique<BamFileIn>(info->options.long_read_file_name.c_str());
            readHeader(info->long_read_header, *(info->long_read_file_handles[t]));
        }
        catch (Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted bam file " << info->options.long_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return 1;
        }

        try
        {
            info->faidx_file_handles[t] = make_unique<FaiIndex>();
            if (!open_file_success(*(info->faidx_file_handles[t]), info->options.reference_file_name.c_str()))
                return 1;
        }
        catch (Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted faidx file " << info->options.long_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }

    // Polish variants
    // -------------------------------------------------------------------------
    info->log_file  << "======================================================================" << std::endl
                    << "START polishing variants in of file " << info->options.candidate_file_name << std::endl
                    << "======================================================================" << std::endl;

    // std::vector<seqan::BamAlignmentRecord> polished_reads; // stores records in case info->options.output-polished-bam is true

    polish_init(variants, info);

    // Write refined variants to output file
    // -------------------------------------------------------------------------
    if (!write_vcf(variants, vcf_header, info))
        return 1;

    info->log_file  << "======================================================================" << std::endl
              << "                                 DONE"  << std::endl
              << "======================================================================" << std::endl;

    return 0;
}
