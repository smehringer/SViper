#include <sviper/io.h>
#include <sviper/sviper.h>

using namespace sviper;
int main(int argc, char const ** argv)
{
    // Struct holding input_output_information information.
    input_output_information info{};
    // Parse Command Line Arguments
    // -------------------------------------------------------------------------
    seqan::ArgumentParser::ParseResult res = parseCommandLine(info.cmd_options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check files
    if (!open_file_success(info.long_read_bai, (info.cmd_options.long_read_file_name + ".bai").c_str()) ||
        !open_file_success(info.short_read_bai, (info.cmd_options.short_read_file_name + ".bai").c_str()) ||
        !open_file_success(info.log_file, (info.cmd_options.output_prefix + ".log").c_str()))
        return 1;

    // Print all input options to log file. Useful for debugging when using SViper as a library and passing in options manually.
    print_log_header(info.cmd_options, info.log_file);
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
                   << "START polishing variants in of file " << info.cmd_options.candidate_file_name << std::endl
                   << "======================================================================" << std::endl;

    #pragma omp parallel for schedule(guided)
    for (unsigned vidx = 0; vidx < variants.size(); ++vidx)
    {
        polish_variant(variants[vidx], info);
    } // parallel for loop

    // Write refined variants to output file
    // -------------------------------------------------------------------------
    if (!write_vcf(variants, vcf_header, info))
        return 1;

    info.log_file << "======================================================================" << std::endl
                  << "                                 DONE" << std::endl
                  << "======================================================================" << std::endl;

    return 0;
}
