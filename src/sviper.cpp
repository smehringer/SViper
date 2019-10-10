#include <sviper.h>

using namespace std;
using namespace seqan;

int main(int argc, char const ** argv)
{
    // Struct holding auxilary information.
    file_info info{};
    // Parse Command Line Arguments
    // -------------------------------------------------------------------------
    // CmdOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(info.options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check files
    if (!open_file_success(info.long_read_bai, (info.options.long_read_file_name + ".bai").c_str()) ||
        !open_file_success(info.short_read_bai, (info.options.short_read_file_name + ".bai").c_str()) ||
        !open_file_success(info.log_file, (info.options.output_prefix + ".log").c_str()))
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
    if (!prep_file_handles(info))
        return 1;

    // Polish variants
    // -------------------------------------------------------------------------
    info.log_file  << "======================================================================" << std::endl
                    << "START polishing variants in of file " << info.options.candidate_file_name << std::endl
                    << "======================================================================" << std::endl;

    // std::vector<seqan::BamAlignmentRecord> polished_reads; // stores records in case info.options.output-polished-bam is true
    #pragma omp parallel for schedule(guided)
    for (unsigned vidx = 0; vidx < variants.size(); ++vidx)
    {
        polish_init(variants[vidx], info);
    } // parallel for loop

    // Write refined variants to output file
    // -------------------------------------------------------------------------
    if (!write_vcf(variants, vcf_header, info))
        return 1;

    info.log_file  << "======================================================================" << std::endl
              << "                                 DONE"  << std::endl
              << "======================================================================" << std::endl;

    return 0;
}
