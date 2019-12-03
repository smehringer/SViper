#pragma once

#include <vector>
#include <string>
#include <fstream>

#include <sviper/auxiliary.h>
#include <sviper/variant.h>

namespace sviper
{
bool read_vcf(std::vector<Variant> & variants, std::vector<std::string> & vcf_header, input_output_information & info)
{
    std::ifstream input_vcf{};           // The candidate variants to polish

    if (!open_file_success(input_vcf, info.cmd_options.candidate_file_name.c_str()))
        return false;

    std::string line{'#'};
    while (getline(input_vcf, line)) // skip header
    {
        if (line[0] != '#')
            break;
        else
            vcf_header.push_back(line); // copy header for output file
    }

    do variants.push_back(Variant(line)); while (getline(input_vcf, line)); // read in variants

    return true;
}

bool write_vcf(std::vector<Variant> & variants, std::vector<std::string> & vcf_header, input_output_information & info)
{
    std::ofstream output_vcf{};          // The polished variant as output
    if (!open_file_success(output_vcf, (info.cmd_options.output_prefix + ".vcf").c_str()))
        return false;

    for (size_t hl = 0; hl < vcf_header.size(); ++hl)
    {
        if (vcf_header[hl].substr(0, 6) == "##INFO")
        {
            bool seen_field_SEQ{false};

            while (vcf_header[hl].substr(0, 6) == "##INFO")
            {
                if (vcf_header[hl].substr(0, 14) == "##INFO=<ID=SEQ")
                    seen_field_SEQ = true;
                output_vcf << vcf_header[hl] << std::endl;
                ++hl;
            }

            // write out SEQ info field only if not already present in the header
            if (!seen_field_SEQ)
                output_vcf << "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"The alternative sequence.\">"
                           << std::endl;
        }
        else if (vcf_header[hl].substr(0, 8) == "##FILTER")
        {
            while (vcf_header[hl].substr(0, 8) == "##FILTER") // write out all existing filters
            {
                output_vcf << vcf_header[hl] << std::endl;
                ++hl;
            }

            // write out custom filters
            output_vcf << "##FILTER=<ID=FAIL0,Description=\"The fasta index has no entry for the given "
                       << "reference name of the variant.\">" << std::endl;
            output_vcf << "##FILTER=<ID=FAIL1,Description=\"No long reads in variant region.\">" << std::endl;
            output_vcf << "##FILTER=<ID=FAIL2,Description=\"No long reads support the variant.\">" << std::endl;
            output_vcf << "##FILTER=<ID=FAIL3,Description=\"The long read regions do not fit.\">" << std::endl;
            output_vcf << "##FILTER=<ID=FAIL4,Description=\"Not enough short reads.\">" << std::endl;
            output_vcf << "##FILTER=<ID=FAIL5,Description=\"The variant was polished away." << std::endl;
            output_vcf << "##FILTER=<ID=FAIL6,Description=\"The variant reference name does not exist in the " <<
                          "short read BAM file." << std::endl;
            output_vcf << "##FILTER=<ID=FAIL7,Description=\"The variant reference name does not exist in the " <<
                          "long read BAM file." << std::endl;
        }

        output_vcf << vcf_header[hl] << std::endl;
    }

    for (auto & var : variants)
        var.write(output_vcf);

    // Write polished reads if specified to output file
    // -------------------------------------------------------------------------
    if (info.cmd_options.output_polished_bam)
    {
        seqan::BamFileOut result_bam(seqan::context(*(info.long_read_file_handles[0]))); // The optional output bam file for polished reads

        if (!open_file_success(result_bam, seqan::toCString(info.cmd_options.output_prefix + "_polished_reads.bam")))
        {
            std::cerr << "Did not write resulting bam file." << std::endl;
        }
        else
        {
            seqan::writeHeader(result_bam, info.long_read_header);
            for (auto const & rec : info.polished_reads)
                seqan::writeRecord(result_bam, rec);
        }
    }
    return true;
}

void print_log_header(CmdOptions const & options, std::ofstream & log_file)
{
    log_file << "Long read file: " << options.long_read_file_name << std::endl
             << "Short read file: " << options.short_read_file_name << std::endl
             << "VCF file: " << options.candidate_file_name << std::endl
             << "Reference file: " << options.reference_file_name << std::endl
             << "Threads set: " << std::to_string(options.threads) << std::endl
             << "Flanking region: " << std::to_string(options.flanking_region) << std::endl
             << "Short read mean coverage: " << std::to_string(options.mean_coverage_of_short_reads) << std::endl
             << "Insert length mean, std dev: " << std::to_string(options.mean_insert_size_of_short_reads) << ", "
                                                << std::to_string(options.stdev_insert_size_of_short_reads) << std::endl
             << "Short read length: " << std::to_string(options.length_of_short_reads) << std::endl;
}
} // namespace sviper
