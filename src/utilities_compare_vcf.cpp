#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include <algorithm>

#include <seqan/arg_parse.h>

#include <sviper/basics.h>
#include <sviper/variant.h>

using namespace sviper;
// global variable of the allowed deviation until SV's are still considered to be
// the same
const double ALLOWED_POS_DEVIATION = 0.95; // percent of the sv_length
const double ALLOWED_LENGTH_DEVIATION = 0.8;
const double SCORE_DEVIATION = 0.05;
//const double QUALITY_CUTOFF = 65.0;

struct CmdOptionsCompareVcf
{
    std::string vcf_file_name;
    std::string golden_vcf_file_name;
    unsigned score_cutoff;

    CmdOptionsCompareVcf() :
        score_cutoff(0)
    {}
};

seqan::ArgumentParser::ParseResult
parseCommandLine(CmdOptionsCompareVcf & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("comparing_vcf_files");

    // We require one argument.
    seqan::addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "INPUT_VCF"));

    seqan::addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "INPUT_GOLDEN_VCF"));

    // Define Options
    seqan::addOption(parser, seqan::ArgParseOption(
        "c", "cutoff", "Quality score cutoff, at which variants are regarded as false positives.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::setDefaultValue(parser, "c", 0);

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    seqan::getOptionValue(options.score_cutoff, parser, "c");
    seqan::getArgumentValue(options.vcf_file_name, parser, 0);
    seqan::getArgumentValue(options.golden_vcf_file_name, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;
}

// IMPORTANT: Those operators are unusual and can e.g. not be used for sorting
bool operator ==(const Variant& lhs, const Variant& rhs)
{
    int length_dev = std::max(80.0, (ALLOWED_LENGTH_DEVIATION * ((lhs.sv_length + rhs.sv_length)/2)));
    int pos_dev = std::max(100.0, (ALLOWED_POS_DEVIATION * ((lhs.sv_length + rhs.sv_length)/2)));

    return (lhs.ref_chrom == rhs.ref_chrom) &&
           (lhs.sv_type == rhs.sv_type) &&
           (std::abs(lhs.sv_length - rhs.sv_length) < length_dev) &&
           (std::abs(lhs.ref_pos - rhs.ref_pos) < pos_dev);
}

bool operator <(const Variant& lhs, const Variant& rhs)
{
    if (lhs.ref_chrom == rhs.ref_chrom)
        return lhs.ref_pos < rhs.ref_pos;
    else
        return lhs.ref_chrom < rhs.ref_chrom;
}

bool operator >(const Variant& lhs, const Variant& rhs)
{
    if (lhs.ref_chrom == rhs.ref_chrom)
        return lhs.ref_pos > rhs.ref_pos;
    else
        return lhs.ref_chrom > rhs.ref_chrom;
}

unsigned sum(std::vector<unsigned> v)
{
    unsigned sum{0};
    for (auto n : v)
        sum+=n;
    return sum;
}

int main(int argc, const char ** argv)
{
    CmdOptionsCompareVcf options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    std::ifstream vcf_file(options.vcf_file_name.c_str());  // must be sorted by position grouped by chromosome
    std::ifstream vcf_file_golden(options.golden_vcf_file_name.c_str());
    std::ofstream log_file((options.vcf_file_name + ".log").c_str());

    if (!vcf_file.is_open())
    {
        std::cerr << "[ERROR] Could not open file " << options.vcf_file_name << std::endl;
        return 1;
    }

    if (!vcf_file_golden.is_open())
    {
        std::cerr << "[ERROR] Could not open file " << options.golden_vcf_file_name << std::endl;
        vcf_file.close();
        return 1;
    }

    Variant curr_var;  // stores information of the current vcf line in file vcf_after
    Variant curr_var_golden; // stores information of the current vcf line in file vcf_golden

    bool golden{true};

    // Initialization:
    // -------------------------------------------------------------------------
    // Skip all headers and init variant with first line in file
    std::string dummy_line{"#"}; // used throughout processing to read into

    while (dummy_line[0] == '#')
        getline(vcf_file, dummy_line);
    curr_var = Variant(dummy_line);

    if (golden)
    {
        dummy_line = "#";
        while (dummy_line[0] == '#')
            getline(vcf_file_golden, dummy_line);
        curr_var_golden = Variant(dummy_line);
    }

    // Iteration:
    // -------------------------------------------------------------------------
    // Now go through all vcf files trying to pair the variants
    // stayed the same = total - impr - wors - notScored
    std::vector<std::vector<unsigned>> bp_deviation_after_to_golden{{}, {}, {}, {}, {}, {}}; // for golden
    std::vector<unsigned> variants_unique_to_golden{0, 0, 0, 0, 0, 0};

    log_file << "CLASS"
             << "\t" << "ID" /*which is the golden id if matched.*/
             << "\t" << "SV_TYPE"
             << "\t" << "CHROM"
             << "\t" << "POS"
             << "\t" << "POS_DIFF_TO_TRUTH"
             << "\t" << "POS_END"
             << "\t" << "POS_END_DIFF_TO_TRUTH"
             << "\t" << "SV_LENGTH"
             << "\t" << "SV_LENGTH_DIFF_TO_TRUTH"
             << "\t" << "POLISHING_SCORE"
             << "\t" << "ERROR_CODE"
             << std::endl;

    while (true) // will be broken if one of the vcf file is at end
    {
        while (golden &&
               !(curr_var_golden == curr_var) && // != is not equal to < so both must be tested
               (curr_var_golden < curr_var))
        {
            log_file << "FN"
                     << "\tGOLDEN_" << curr_var_golden.id
                     << "\t" << curr_var_golden.sv_type
                     << "\t" << curr_var_golden.ref_chrom
                     << "\t" << curr_var_golden.ref_pos
                     << "\tNA"
                     << "\t" << (curr_var_golden.ref_pos + curr_var_golden.sv_length)
                     << "\tNA"
                     << "\t" << curr_var_golden.sv_length
                     << "\tNA"
                     << "\tNA"
                     << "\tNA"
                     << std::endl;

            if (!getline(vcf_file_golden, dummy_line))
            {
                golden = false;
            }
            else
            {
                curr_var_golden = Variant(dummy_line);
            }

        }

        if (curr_var.filter == "PASS" && // variant passed polishing
            curr_var.quality >= options.score_cutoff)
        {
            if (golden && curr_var_golden == curr_var) // TP
            {
                assert(curr_var.ref_chrom == curr_var_golden.ref_chrom);
                assert(curr_var.sv_type == curr_var_golden.sv_type);

                log_file << "TP"
                         << "\tGOLDEN_" << curr_var_golden.id
                         << "\t" << curr_var.sv_type
                         << "\t" << curr_var.ref_chrom
                         << "\t" << curr_var.ref_pos
                         << "\t" << (curr_var.ref_pos - curr_var_golden.ref_pos)
                         << "\t" << (curr_var.ref_pos + curr_var.sv_length)
                         << "\t" << ((curr_var.ref_pos + curr_var.sv_length) - (curr_var_golden.ref_pos + curr_var_golden.sv_length))
                         << "\t" << curr_var.sv_length
                         << "\t" << (curr_var.sv_length - curr_var_golden.sv_length)
                         << "\t" << curr_var.quality
                         << "\t" << curr_var.filter
                         << std::endl;

                if (!getline(vcf_file_golden, dummy_line))
                    golden = false;
                else
                    curr_var_golden = Variant(dummy_line);
            }
            else // FP
            {
                log_file << "FP"
                         << "\t" << curr_var.id
                         << "\t" << curr_var.sv_type
                         << "\t" << curr_var.ref_chrom
                         << "\t" << curr_var.ref_pos
                         << "\tNA"
                         << "\t" << (curr_var.ref_pos + curr_var.sv_length)
                         << "\tNA"
                         << "\t" << curr_var.sv_length
                         << "\tNA"
                         << "\t" << curr_var.quality
                         << "\t" << curr_var.filter
                         << std::endl;
            }
        }
        else // variant failed polishing
        {
            if (golden && curr_var_golden == curr_var) // FN
            {
                assert(curr_var.ref_chrom == curr_var_golden.ref_chrom);
                assert(curr_var.sv_type == curr_var_golden.sv_type);

                log_file << "FN"
                         << "\tGOLDEN_" << curr_var_golden.id
                         << "\t" << curr_var.sv_type
                         << "\t" << curr_var.ref_chrom
                         << "\t" << curr_var.ref_pos
                         << "\t" << (curr_var.ref_pos - curr_var_golden.ref_pos)
                         << "\t" << (curr_var.ref_pos + curr_var.sv_length)
                         << "\t" << ((curr_var.ref_pos + curr_var.sv_length) - (curr_var_golden.ref_pos + curr_var_golden.sv_length))
                         << "\t" << curr_var.sv_length
                         << "\t" << (curr_var.sv_length - curr_var_golden.sv_length)
                         << "\t" << curr_var.quality
                         << "\t" << curr_var.filter
                         << std::endl;

                if (!getline(vcf_file_golden, dummy_line))
                    golden = false;
                else
                    curr_var_golden = Variant(dummy_line);
            }
            else // TN
            {
                log_file << "TN"
                         << "\t" << curr_var.id
                         << "\t" << curr_var.sv_type
                         << "\t" << curr_var.ref_chrom
                         << "\t" << curr_var.ref_pos
                         << "\tNA"
                         << "\t" << (curr_var.ref_pos + curr_var.sv_length)
                         << "\tNA"
                         << "\t" << curr_var.sv_length
                         << "\tNA"
                         << "\t" << curr_var.quality
                         << "\t" << curr_var.filter
                         << std::endl;
            }
        }

        // update variant
        if (!getline(vcf_file, dummy_line))
            break;
        curr_var = Variant(dummy_line);
    }

    // Finalizing:
    // -------------------------------------------------------------------------

	vcf_file.close();
    vcf_file_golden.close();

    return 0;
}
