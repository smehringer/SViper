#pragma once
#include <cmath>
#include <chrono>
#include <limits>
#include <sstream>
#include <thread>
#include <vector>

#include <sviper/auxiliary.h>
#include <sviper/basics.h>
#include <sviper/config.h>
#include <sviper/evaluate_final_mapping.h>
#include <sviper/merge_split_alignments.h>
#include <sviper/polishing.h>
#include <sviper/variant.h>

#include <seqan/align.h>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/graph_msa.h>

namespace sviper
{
seqan::ArgumentParser::ParseResult parseCommandLine(CmdOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("SViper");
    setVersion(parser, "2.1.0");

    seqan::addOption(parser, seqan::ArgParseOption(
        "c", "candidate-vcf",
        "A structural variant vcf file (with e.g. <DEL> tags), containing the potential variant sites to be looked at.",
        seqan::ArgParseArgument::INPUT_FILE, "VCF_FILE"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "s", "short-read-bam",
        "The indexed bam file containing short used for polishing at variant sites.",
        seqan::ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "l", "long-read-bam",
        "The indexed bam file containing long reads to be polished at variant sites.",
        seqan::ArgParseArgument::INPUT_FILE, "BAM_FILE"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "r", "reference",
        "The indexed (fai) reference file.",
        seqan::ArgParseArgument::INPUT_FILE, "FA_FILE"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "t", "threads",
        "The threads to use.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "k", "flanking-region",
        "The flanking region in bp's around a breakpoint to be considered for polishing",
        seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "x", "coverage-short-reads",
        "The original short read mean coverage. This value is used to restrict short read coverage on extraction to avoid mapping bias",
        seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "", "median-ins-size-short-reads",
        "The median of the short read insert size (end of read1 until beginning of read2). "
        "This value is used to compute a threshold for error correction.",
        seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "", "stdev-ins-size-short-reads",
        "The median of the short read insert size (end of read1 until beginning of read2). "
        "This value is used to compute a threshold for error correction..",
        seqan::ArgParseArgument::INTEGER, "INT"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "o", "output-prefix",
        "A name for the output files. The current output is a log file and vcf file, that contains the "
        "polished sequences for each variant.",
        seqan::ArgParseArgument::INPUT_FILE, "PREFIX"));

    seqan::addOption(parser, seqan::ArgParseOption(
        "v", "verbose",
        "Turn on detailed information about the process."));

    seqan::addOption(parser, seqan::ArgParseOption(
        "", "output-polished-bam", "For debugging or manual inspection the polished reads can be written to a file."));

    seqan::setRequired(parser, "c");
    seqan::setRequired(parser, "l");
    seqan::setRequired(parser, "s");
    seqan::setRequired(parser, "r");

    seqan::setMinValue(parser, "k", "50");
    seqan::setMaxValue(parser, "k", "1000");
    seqan::setDefaultValue(parser, "k", "400");

    seqan::setDefaultValue(parser, "t", std::thread::hardware_concurrency());

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    seqan::getOptionValue(options.candidate_file_name, parser, "candidate-vcf");
    seqan::getOptionValue(options.long_read_file_name, parser, "long-read-bam");
    seqan::getOptionValue(options.short_read_file_name, parser, "short-read-bam");
    seqan::getOptionValue(options.reference_file_name, parser, "reference");
    seqan::getOptionValue(options.output_prefix, parser, "output-prefix");
    seqan::getOptionValue(options.flanking_region, parser, "flanking-region");
    seqan::getOptionValue(options.mean_coverage_of_short_reads, parser, "coverage-short-reads");
    seqan::getOptionValue(options.mean_insert_size_of_short_reads, parser, "median-ins-size-short-reads");
    seqan::getOptionValue(options.stdev_insert_size_of_short_reads, parser, "stdev-ins-size-short-reads");
    seqan::getOptionValue(options.threads, parser, "threads");
    options.verbose = isSet(parser, "verbose");
    options.output_polished_bam = isSet(parser, "output-polished-bam");

    if (options.output_prefix.empty())
        options.output_prefix = options.candidate_file_name + "_polished";

    return seqan::ArgumentParser::PARSE_OK;
}

bool polish_variant(Variant & var, input_output_information & info)
{
    std::stringstream localLog{};
    seqan::BamFileIn & short_read_bam = *(info.short_read_file_handles[omp_get_thread_num()]);
    seqan::BamFileIn & long_read_bam = *(info.long_read_file_handles[omp_get_thread_num()]);
    seqan::FaiIndex & faiIndex = *(info.faidx_file_handles[omp_get_thread_num()]);

    if (var.sv_type != SV_TYPE::DEL && var.sv_type != SV_TYPE::INS)
    {
        #pragma omp critical
        info.log_file << "----------------------------------------------------------------------" << std::endl
                      << " SKIP Variant " << var.id << " at " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
                      << "----------------------------------------------------------------------" << std::endl;
        var.filter = "SKIP";
        return false;
    }
    if (var.sv_length > 1000000)
    {
        #pragma omp critical
        info.log_file << "----------------------------------------------------------------------" << std::endl
                      << " SKIP too long Variant " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
                      << "----------------------------------------------------------------------" << std::endl;
        var.filter = "SKIP";
        return false;
    }

    localLog << "----------------------------------------------------------------------" << std::endl
             << " PROCESS Variant " << var.id << " at " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
             << "----------------------------------------------------------------------" << std::endl;

    // Compute reference length, start and end position of region of interest
    // ---------------------------------------------------------------------
    unsigned ref_fai_idx = 0;
    if (!seqan::getIdByName(ref_fai_idx, faiIndex, var.ref_chrom))
    {
        localLog << "[ ERROR ]: FAI index has no entry for reference name"
                 << var.ref_chrom << std::endl;
        #pragma omp critical
        info.log_file << localLog.str() << std::endl;
        var.filter = "FAIL0";
        return false;
    }

    // cash variables to avoid recomputing
    // Note that the positions are one based/ since the VCF format is one based
    int const ref_length       = seqan::sequenceLength(faiIndex, ref_fai_idx);
    int const ref_region_start = std::max(1, var.ref_pos - info.cmd_options.flanking_region);
    int const ref_region_end   = std::min(ref_length, var.ref_pos_end + info.cmd_options.flanking_region);

    int const var_ref_pos_add50     = std::min(ref_length, var.ref_pos + 50);
    int const var_ref_pos_sub50     = std::max(1, var.ref_pos - 50);
    int const var_ref_pos_end_add50 = std::min(ref_length, var.ref_pos_end + 50);
    int const var_ref_pos_end_sub50 = std::max(1, var.ref_pos_end - 50);

    localLog << "--- Reference region " << var.ref_chrom << ":"
             << ref_region_start << "-" << ref_region_end << std::endl;

    SEQAN_ASSERT_LEQ(ref_region_start, ref_region_end);

    // get reference id in bam File
    unsigned rID_short{};
    unsigned rID_long{};

    if (!seqan::getIdByName(rID_short, seqan::contigNamesCache(seqan::context(short_read_bam)), var.ref_chrom))
    {
        localLog << "[ ERROR ]: No reference sequence named "
                 << var.ref_chrom << " in short read bam file." << std::endl;
        var.filter = "FAIL6";
        #pragma omp critical
        info.log_file << localLog.str() << std::endl;
        return false;
    }

    if (!seqan::getIdByName(rID_long, seqan::contigNamesCache(seqan::context(long_read_bam)), var.ref_chrom))
    {
        localLog << "[ ERROR ]: No reference sequence named "
                 << var.ref_chrom << " in long read bam file." << std::endl;
        var.filter = "FAIL7";
        #pragma omp critical
        info.log_file << localLog.str() << std::endl;
        return false;
    }

    // Extract long reads
    // ---------------------------------------------------------------------
    std::vector<seqan::BamAlignmentRecord> supporting_records{};

    {
        std::vector<seqan::BamAlignmentRecord> long_reads{};

        // extract overlapping the start breakpoint +-50 bp's
        seqan::viewRecords(long_reads, long_read_bam, info.long_read_bai, rID_long, var_ref_pos_sub50, var_ref_pos_add50);
        // extract overlapping the end breakpoint +-50 bp's
        seqan::viewRecords(long_reads, long_read_bam, info.long_read_bai, rID_long, var_ref_pos_end_sub50, var_ref_pos_end_add50);

        if (long_reads.size() == 0)
        {
            localLog << "ERROR1: No long reads in reference region "
                     << var.ref_chrom << ":" << var_ref_pos_sub50 << "-" << var_ref_pos_add50 << " or "
                     << var.ref_chrom << ":" << var_ref_pos_end_sub50 << "-" << var_ref_pos_end_add50 << std::endl;

            var.filter = "FAIL1";
            #pragma omp critical
            info.log_file << localLog.str() << std::endl;
            return false;
        }

        localLog << "--- Extracted " << long_reads.size() << " long read(s). May include duplicates. " << std::endl;

        // Search for supporting reads
        // ---------------------------------------------------------------------
        // TODO check if var is not empty!

        localLog << "--- Searching in (reference) region ["
                 << (int)(var.ref_pos - DEV_POS * var.sv_length) << "-"
                 << (int)(var.ref_pos + var.sv_length + DEV_POS * var.sv_length) << "]"
                 << " for a variant of type " << var.alt_seq
                 << " of length " << (int)(var.sv_length - DEV_SIZE * var.sv_length) << "-"
                 << (int)(var.sv_length + DEV_SIZE * var.sv_length) << " bp's" << std::endl;

        for (auto const & rec : long_reads)
            if (record_supports_variant(rec, var) && length(rec.seq) > 0 /*sequence information is given*/)
                supporting_records.push_back(rec);

        if (supporting_records.size() == 0)
        {
            localLog << "--- No supporting reads that span the variant, start merging..." << std::endl;

            // Merge supplementary alignments to primary
            // ---------------------------------------------------------------------
            std::sort(long_reads.begin(), long_reads.end(), bamRecordNameLess());
            long_reads = merge_alignments(long_reads); // next to merging this will also get rid of duplicated reads

            localLog << "--- After merging " << long_reads.size() << " read(s) remain(s)." << std::endl;

            for (auto const & rec : long_reads)
                if (record_supports_variant(rec, var) && length(rec.seq) > 0 /*sequence information is given*/)
                    supporting_records.push_back(rec);

            if (supporting_records.size() == 0) // there are none at all
            {
                localLog << "ERROR2: No supporting long reads for a " << var.alt_seq
                         << " in region " << var.ref_chrom << ":"
                         << var_ref_pos_sub50 << "-" << var_ref_pos_end_add50
                         << std::endl;

                var.filter = "FAIL2";
                #pragma omp critical
                info.log_file << localLog.str() << std::endl;
                return false;
            }
        }
        else
        {
            // remove duplicates
            std::sort(supporting_records.begin(), supporting_records.end(), bamRecordNameLess());
            auto last = std::unique(supporting_records.begin(), supporting_records.end(), bamRecordEqual());
            supporting_records.erase(last, supporting_records.end());
        }

        localLog << "--- After searching for variant " << supporting_records.size()
                 << " supporting read(s) remain." << std::endl;
    } // scope of long_reads ends

    // Crop fasta sequence of each supporting read for consensus
    // ---------------------------------------------------------------------
    localLog << "--- Cropping long reads with a buffer of +-" << info.cmd_options.flanking_region << " around variants." << std::endl;

    seqan::StringSet<seqan::Dna5String> supporting_sequences{};
    std::vector<seqan::BamAlignmentRecord>::size_type maximum_long_reads = 5;

    // sort records such that the highest quality ones are chosen first
    std::sort(supporting_records.begin(), supporting_records.end(), bamRecordMapQGreater());

    for (unsigned i = 0; i < std::min(maximum_long_reads, supporting_records.size()); ++i)
    {
        auto region = get_read_region_boundaries(supporting_records[i], ref_region_start, ref_region_end);

        assert(std::get<0>(region) >= 0);
        assert(std::get<1>(region) >= 0);
        assert(std::get<0>(region) <= std::get<1>(region));
        assert(std::get<1>(region) - std::get<0>(region) <= static_cast<int>(length(supporting_records[i].seq)));

        seqan::Dna5String reg = seqan::infix(supporting_records[i].seq, std::get<0>(region), std::get<1>(region));

        // For deletions, the expected size of the subsequence is that of
        // the flanking region, since the rest is deleted. For insertions it
        // is that of the flanking region + the insertion length.
        int32_t expected_length{2*info.cmd_options.flanking_region};
        if (var.sv_type == SV_TYPE::INS)
            expected_length += var.sv_length;

        if (abs(static_cast<int32_t>(seqan::length(reg)) - expected_length) > info.cmd_options.flanking_region)
        {
            localLog << "------ Skip Read - Length:" << seqan::length(reg) << " Qual:" << supporting_records[i].mapQ
                     << " Name: "<< supporting_records[i].qName << std::endl;
            ++maximum_long_reads;
            return false; // do not use under or oversized region
        }

        seqan::appendValue(supporting_sequences, reg);

        localLog << "------ Region: [" << std::get<0>(region) << "-" << std::get<1>(region)
                 << "] Length:" << seqan::length(reg) << " Qual:" << supporting_records[i].mapQ
                 << " Name: "<< supporting_records[i].qName << std::endl;
    }

    if (seqan::length(supporting_sequences) == 0)
    {
        localLog << "ERROR3: No fitting regions for a " << var.alt_seq
                 << " in region " << var.ref_chrom << ":"
                 << var_ref_pos_sub50 << "-" << var_ref_pos_end_add50
                 << std::endl;

        var.filter = "FAIL3";
        #pragma omp critical
        info.log_file << localLog.str() << std::endl;
        return false;
    }

    // Build consensus of supporting read regions
    // ---------------------------------------------------------------------
    std::vector<double> mapping_qualities{};
    mapping_qualities.resize(supporting_records.size());
    for (unsigned i = 0; i < supporting_records.size(); ++i)
        mapping_qualities[i] = (supporting_records[i]).mapQ;

    seqan::Dna5String cns = build_consensus(supporting_sequences, mapping_qualities);

    localLog << "--- Built a consensus with a MSA of length " << seqan::length(cns) << "." << std::endl;

    seqan::Dna5String polished_ref{};
    SViperConfig config{info.cmd_options};
    config.ref_flank_length = 500;

    {
        seqan::StringSet<seqan::Dna5QString> short_reads_1{}; // reads (first in pair)
        seqan::StringSet<seqan::Dna5QString> short_reads_2{}; // mates (second in pair)

        {
            // Extract short reads in region
            // ---------------------------------------------------------------------
            std::vector<seqan::BamAlignmentRecord> short_reads{};
            // If the breakpoints are farther apart then illumina-read-length + 2 * flanking-region,
            // then extract reads for each break point separately.
            if (ref_region_end - ref_region_start > info.cmd_options.flanking_region * 2 + info.cmd_options.length_of_short_reads)
            {
                // extract reads left of the start of the variant [start-flanking_region, start+flanking_region]
                unsigned e = std::min(ref_length, var.ref_pos + info.cmd_options.flanking_region);
                seqan::viewRecords(short_reads, short_read_bam, info.short_read_bai, rID_short, ref_region_start, e);
                cut_down_high_coverage(short_reads, info.cmd_options.mean_coverage_of_short_reads);

                // and right of the end of the variant [end-flanking_region, end+flanking_region]
                std::vector<seqan::BamAlignmentRecord> tmp_short_reads;
                unsigned s = std::max(1, var.ref_pos_end - info.cmd_options.flanking_region);
                seqan::viewRecords(tmp_short_reads, short_read_bam, info.short_read_bai, rID_short, s, ref_region_end);
                cut_down_high_coverage(tmp_short_reads, info.cmd_options.mean_coverage_of_short_reads);
                append(short_reads, tmp_short_reads);
            }
            else
            {
                // extract reads left of the start of the variant [start-flanking_region, start]
                seqan::viewRecords(short_reads, short_read_bam, info.short_read_bai, rID_short, ref_region_start, ref_region_end);
                cut_down_high_coverage(short_reads, info.cmd_options.mean_coverage_of_short_reads);
            }

            if (short_reads.size() < 20)
            {
                localLog << "ERROR4: Not enough short reads (only " << short_reads.size()
                         << ") for variant of type " << var.alt_seq
                         << " in region " << var.ref_chrom << ":" << ref_region_start
                         << "-" << ref_region_end << std::endl;

                var.filter = "FAIL4";
                #pragma omp critical
                info.log_file << localLog.str() << std::endl;
                return false;
            }

            records_to_read_pairs(short_reads_1, short_reads_2, short_reads, short_read_bam, info.short_read_bai);

            localLog << "--- Extracted " << seqan::length(short_reads_1) << " pairs (proper or dummy pairs)." << std::endl;
        } // scope of short reads ends

        // Flank consensus sequence
        // ---------------------------------------------------------------------
        // Before polishing, append a reference flank to the conesnsus such that
        // the reads find a high quality anchor for mapping and pairs are correctly
        // identified.
        seqan::Dna5String flanked_consensus = append_ref_flanks(cns, faiIndex,
                                                                ref_fai_idx, ref_length,
                                                                ref_region_start, ref_region_end,
                                                                config.ref_flank_length);

        // Polish flanked consensus sequence with short reads
        // ---------------------------------------------------------------------
        compute_baseQ_stats(config, short_reads_1, short_reads_2); //TODO:: return qualities and assign to config outside

        localLog << "--- Short read base qualities: avg=" << config.baseQ_mean
                 << " stdev=" << config.baseQ_std << "." << std::endl;

        polished_ref = polish_to_perfection(short_reads_1, short_reads_2, flanked_consensus, config);

        localLog << "DONE POLISHING: Total of "
                 << config.substituted_bases << " substituted, "
                 << config.deleted_bases     << " deleted and "
                 << config.inserted_bases    << " inserted bases. "
                 << config.rounds            << " rounds."
                 << std::endl;
    } // scope of short_reads1 and short_reads2 ends

    for (unsigned i = config.ref_flank_length; i < seqan::length(config.cov_profile) - config.ref_flank_length; ++i)
        localLog << config.cov_profile[i] << " ";
    localLog << std::endl;

    // Align polished sequence to reference
    // ---------------------------------------------------------------------
    seqan::Dna5String ref_part{};
    seqan::readRegion(ref_part, faiIndex, ref_fai_idx,
                      std::max(1u, ref_region_start - config.ref_flank_length),
                      std::min(ref_region_end + config.ref_flank_length, static_cast<unsigned>(ref_length)));

    typedef seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> TGapsRead;
    typedef seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> TGapsRef;
    TGapsRef gapsRef(ref_part);
    TGapsRead gapsSeq(polished_ref);

    seqan::globalAlignment(gapsRef, gapsSeq,
                           seqan::Score<double, seqan::Simple>(config.MM, config.MX, config.GE, config.GO),
                           seqan::AlignConfig<false, false, false, false>(),
                           seqan::ConvexGaps());

    seqan::BamAlignmentRecord final_record{};
    final_record.beginPos = std::max(0u, ref_region_start - config.ref_flank_length);
    final_record.seq = polished_ref;
    seqan::getIdByName(final_record.rID, seqan::contigNamesCache(seqan::context(long_read_bam)), var.ref_chrom);
    seqan::getCigarString(final_record.cigar, gapsRef, gapsSeq, std::numeric_limits<unsigned>::max());
    // std::numeric_limits<unsigned>::max() because in function getCigarString the value is compared to type unsigned
    // And we never want replace a Deletions D with N (short read exon identification)

    // Evaluate Alignment
    // ---------------------------------------------------------------------
    // If refine_variant fails, the variant is not supported anymore but remains unchanged.
    // If not, the variant now might have different start/end positions and other information

    assign_quality(final_record, var, config); // assigns a score to record.mapQ and var.quality

    if (!refine_variant(final_record, var))
    {
        var.filter = "FAIL5";
        //assign_quality(record, var, false);
        localLog << "ERROR5: \"Polished away\" variant " << var.id << " at " << var.ref_chrom << ":"
                 << var.ref_pos << "\t" << var.alt_seq << "\tScore:\t" << var.quality
                 << std::endl;
    }
    else
    {
        //assign_quality(record, var, true);
        localLog << "SUCCESS: Polished variant " << var.id << " at " << var.ref_chrom << ":"
                 << var.ref_pos << "\t" << var.alt_seq << "\tScore:\t" << var.quality
                 << std::endl;
    }

    if (info.cmd_options.output_polished_bam)
    {
        std::string read_identifier = (std::string("polished_var") +
                                       ":" + var.ref_chrom +
                                       ":" + std::to_string(var.ref_pos) +
                                       ":" + std::to_string(var.ref_pos_end) +
                                       ":" + var.id);
        final_record.qName = read_identifier;

        #pragma omp critical
        info.polished_reads.push_back(final_record);
    }

    #pragma omp critical
    info.log_file << localLog.str() << std::endl;

    return true;
}
} // namespace sviper
