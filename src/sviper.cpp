#include <sviper.h>

using namespace std;
using namespace seqan;

// Global file info struct.
file_info * info{};

ArgumentParser::ParseResult
parseCommandLine(CmdOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("SViper");
    setVersion(parser, "2.0.0");

    addOption(parser, seqan::ArgParseOption(
        "c", "candidate-vcf",
        "A structural variant vcf file (with e.g. <DEL> tags), containing the potential variant sites to be looked at.",
        seqan::ArgParseArgument::INPUT_FILE, "VCF_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "s", "short-read-bam",
        "The indexed bam file containing short used for polishing at variant sites.",
        seqan::ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "l", "long-read-bam",
        "The indexed bam file containing long reads to be polished at variant sites.",
        seqan::ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "r", "reference",
        "The indexed (fai) reference file.",
        seqan::ArgParseArgument::INPUT_FILE, "FA_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "t", "threads",
        "The threads to use.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "k", "flanking-region",
        "The flanking region in bp's around a breakpoint to be considered for polishing",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "x", "coverage-short-reads",
        "The original short read mean coverage. This value is used to restrict short read coverage on extraction to avoid mapping bias",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "", "median-ins-size-short-reads",
        "The median of the short read insert size (end of read1 until beginning of read2). "
        "This value is used to compute a threshold for error correction.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "", "stdev-ins-size-short-reads",
        "The median of the short read insert size (end of read1 until beginning of read2). "
        "This value is used to compute a threshold for error correction..",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "o", "output-prefix",
        "A name for the output files. The current output is a log file and vcf file, that contains the "
        "polished sequences for each variant.",
        seqan::ArgParseArgument::INPUT_FILE, "PREFIX"));

    addOption(parser, seqan::ArgParseOption(
        "v", "verbose",
        "Turn on detailed information about the process."));

    addOption(parser, seqan::ArgParseOption(
        "", "output-polished-bam", "For debugging or manual inspection the polished reads can be written to a file."));

    setRequired(parser, "c");
    setRequired(parser, "l");
    setRequired(parser, "s");
    setRequired(parser, "r");

    setMinValue(parser, "k", "50");
    setMaxValue(parser, "k", "1000");
    setDefaultValue(parser, "k", "400");

    setDefaultValue(parser, "t", std::thread::hardware_concurrency());

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.candidate_file_name, parser, "candidate-vcf");
    getOptionValue(options.long_read_file_name, parser, "long-read-bam");
    getOptionValue(options.short_read_file_name, parser, "short-read-bam");
    getOptionValue(options.reference_file_name, parser, "reference");
    getOptionValue(options.output_prefix, parser, "output-prefix");
    getOptionValue(options.flanking_region, parser, "flanking-region");
    getOptionValue(options.mean_coverage_of_short_reads, parser, "coverage-short-reads");
    getOptionValue(options.mean_insert_size_of_short_reads, parser, "median-ins-size-short-reads");
    getOptionValue(options.stdev_insert_size_of_short_reads, parser, "stdev-ins-size-short-reads");
    getOptionValue(options.threads, parser, "threads");
    options.verbose = isSet(parser, "verbose");
    options.output_polished_bam = isSet(parser, "output-polished-bam");

    if (options.output_prefix.empty())
        options.output_prefix = options.candidate_file_name + "_polished";

    return seqan::ArgumentParser::PARSE_OK;
}

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
