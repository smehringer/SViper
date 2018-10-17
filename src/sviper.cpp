#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <memory>
#include <thread>
#include <limits>

#include <basics.h>
#include <config.h>
#include <evaluate_final_mapping.h>
#include <merge_split_alignments.h>
#include <polishing.h>
#include <variant.h>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>


using namespace std;
using namespace seqan;

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
    CmdOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // Check files
    // -------------------------------------------------------------------------
    ifstream input_vcf;           // The candidate variants to polish
    BamHeader long_read_header;   // The bam header object needed to fill bam context
    BamHeader short_read_header;  // The bam header object needed to fill bam context
    BamIndex<Bai> long_read_bai;  // The bam index to the long read bam file
    BamIndex<Bai> short_read_bai; // The bam index to the short read bam file
    std::ofstream log_file;

    if (!open_file_success(input_vcf, options.candidate_file_name.c_str()) ||
        !open_file_success(long_read_bai, (options.long_read_file_name + ".bai").c_str()) ||
        !open_file_success(short_read_bai, (options.short_read_file_name + ".bai").c_str()) ||
        !open_file_success(log_file, (options.output_prefix + ".log").c_str()))
        return 1;

    // Read variants into container // TODO:: use seqan vcf parser instead (needs to be extended)
    // -------------------------------------------------------------------------
    std::vector<std::string> vcf_header;
    std::vector<Variant> variants;

    std::string line{'#'};
    while (getline(input_vcf, line)) // skip header
    {
        if (line[0] != '#')
            break;
        else
            vcf_header.push_back(line); // copy header for output file
    }

    do variants.push_back(Variant(line)); while (getline(input_vcf, line)); // read in variants

    // Prepare file hangles for parallel computing
    // -------------------------------------------------------------------------
    unsigned num_threads{options.threads};
    omp_set_num_threads(num_threads);

    std::vector<std::unique_ptr<BamFileIn>> long_read_file_handles;
    std::vector<std::unique_ptr<BamFileIn>> short_read_file_handles;
    std::vector<std::unique_ptr<FaiIndex>> faidx_file_handles;

    long_read_file_handles.resize(num_threads);
    short_read_file_handles.resize(num_threads);
    faidx_file_handles.resize(num_threads);

    for (unsigned t = 0; t < num_threads; ++t)
    {
        try
        {
            short_read_file_handles[t] = make_unique<BamFileIn>(options.short_read_file_name.c_str());
            readHeader(short_read_header, *(short_read_file_handles[t]));
        }
        catch (Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted bam file " << options.short_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return 1;
        }

        try
        {
            long_read_file_handles[t] = make_unique<BamFileIn>(options.long_read_file_name.c_str());
            readHeader(long_read_header, *(long_read_file_handles[t]));
        }
        catch (Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted bam file " << options.long_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return 1;
        }

        try
        {
            faidx_file_handles[t] = make_unique<FaiIndex>();
            if (!open_file_success(*(faidx_file_handles[t]), options.reference_file_name.c_str()))
                return 1;
        }
        catch (Exception & e)
        {
            std::cerr << "[ ERROR ] Corrupted faidx file " << options.long_read_file_name << std::endl;
            std::cerr << e.what() << std::endl;
            return 1;
        }
    }

    // Polish variants
    // -------------------------------------------------------------------------
    log_file  << "======================================================================" << std::endl
              << "START polishing variants in of file " << options.candidate_file_name << std::endl
              << "======================================================================" << std::endl;

    std::vector<seqan::BamAlignmentRecord> polished_reads; // stores records in case options.output-polished-bam is true

    #pragma omp parallel for schedule(guided)
    for (unsigned vidx = 0; vidx < variants.size(); ++vidx)
    {
        Variant & var = variants[vidx];
        std::stringstream localLog;
        seqan::BamFileIn & short_read_bam = *(short_read_file_handles[omp_get_thread_num()]);
        seqan::BamFileIn & long_read_bam = *(long_read_file_handles[omp_get_thread_num()]);
        seqan::FaiIndex & faiIndex = *(faidx_file_handles[omp_get_thread_num()]);

        if (var.alt_seq != "<DEL>" && var.alt_seq != "<INS>")
        {
            #pragma omp critical
            log_file << "----------------------------------------------------------------------" << std::endl
                     << " SKIP Variant " << var.id << " at " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
                     << "----------------------------------------------------------------------" << std::endl;
            var.filter = "SKIP";
            continue;
        }
        if (var.sv_length > 10000000)
        {
            #pragma omp critical
            log_file << "----------------------------------------------------------------------" << std::endl
                     << " SKIP too long Variant " << var.ref_chrom << ":" << var.ref_pos << " " << var.alt_seq << " L:" << var.sv_length << std::endl
                     << "----------------------------------------------------------------------" << std::endl;
            var.filter = "SKIP";
            continue;
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
            log_file << localLog.str() << std::endl;
            var.filter = "FAIL0";
            continue;
        }

        // cash variables to avoid recomputing
        // Note that the positions are one based/ since the VCF format is one based
        int const ref_length       = seqan::sequenceLength(faiIndex, ref_fai_idx);
        int const ref_region_start = std::max(1, var.ref_pos - options.flanking_region);
        int const ref_region_end   = std::min(ref_length, var.ref_pos_end + options.flanking_region);

        int const var_ref_pos_add50     = std::min(ref_length, var.ref_pos + 50);
        int const var_ref_pos_sub50     = std::max(1, var.ref_pos - 50);
        int const var_ref_pos_end_add50 = std::min(ref_length, var.ref_pos_end + 50);
        int const var_ref_pos_end_sub50 = std::max(1, var.ref_pos_end - 50);

        localLog << "--- Reference region " << var.ref_chrom << ":"
                 << ref_region_start << "-" << ref_region_end << std::endl;

        SEQAN_ASSERT_LEQ(ref_region_start, ref_region_end);

        // get reference id in bam File
        unsigned rID_short;
        unsigned rID_long;

        if (!seqan::getIdByName(rID_short, seqan::contigNamesCache(seqan::context(short_read_bam)), var.ref_chrom))
        {
            localLog << "[ ERROR ]: No reference sequence named "
                     << var.ref_chrom << " in short read bam file." << std::endl;
            var.filter = "FAIL4";
            #pragma omp critical
            log_file << localLog.str() << std::endl;
            continue;
        }

        if (!seqan::getIdByName(rID_long, seqan::contigNamesCache(seqan::context(long_read_bam)), var.ref_chrom))
        {
            localLog << "[ ERROR ]: No reference sequence named "
                     << var.ref_chrom << " in long read bam file." << std::endl;
            var.filter = "FAIL1";
            #pragma omp critical
            log_file << localLog.str() << std::endl;
            continue;
        }

        // Extract long reads
        // ---------------------------------------------------------------------
        vector<BamAlignmentRecord> supporting_records;

        {
            std::vector<seqan::BamAlignmentRecord> long_reads;

            // extract overlapping the start breakpoint +-50 bp's
            viewRecords(long_reads, long_read_bam, long_read_bai, rID_long, var_ref_pos_sub50, var_ref_pos_add50);
            // extract overlapping the end breakpoint +-50 bp's
            viewRecords(long_reads, long_read_bam, long_read_bai, rID_long, var_ref_pos_end_sub50, var_ref_pos_end_add50);

            if (long_reads.size() == 0)
            {
                localLog << "ERROR1: No long reads in reference region "
                         << var.ref_chrom << ":" << var_ref_pos_sub50 << "-" << var_ref_pos_add50 << " or "
                         << var.ref_chrom << ":" << var_ref_pos_end_sub50 << "-" << var_ref_pos_end_add50 << std::endl;

                var.filter = "FAIL1";
                #pragma omp critical
                log_file << localLog.str() << std::endl;
                continue;
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
                if (record_supports_variant(rec, var))
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
                    if (record_supports_variant(rec, var))
                        supporting_records.push_back(rec);

                if (supporting_records.size() == 0) // there are none at all
                {
                    localLog << "ERROR2: No supporting long reads for a " << var.alt_seq
                             << " in region " << var.ref_chrom << ":"
                             << var_ref_pos_sub50 << "-" << var_ref_pos_end_add50
                             << std::endl;

                    var.filter = "FAIL2";
                    #pragma omp critical
                    log_file << localLog.str() << std::endl;
                    continue;
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
        localLog << "--- Cropping long reads with a buffer of +-" << options.flanking_region << " around variants." << endl;

        StringSet<Dna5String> supporting_sequences;
        std::vector<seqan::BamAlignmentRecord>::size_type maximum_long_reads = 5;

        // sort records such that the highest quality ones are chosen first
        std::sort(supporting_records.begin(), supporting_records.end(), bamRecordMapQGreater());

        for (unsigned i = 0; i < std::min(maximum_long_reads, supporting_records.size()); ++i)
        {
            auto region = get_read_region_boundaries(supporting_records[i], ref_region_start, ref_region_end);
            Dna5String reg = seqan::infix(supporting_records[i].seq, get<0>(region), get<1>(region));

            // For deletions, the expected size of the subsequence is that of
            // the flanking region, since the rest is deleted. For insertions it
            // is that of the flanking region + the insertion length.
            int32_t expected_length{2*options.flanking_region};
            if (var.sv_type == SV_TYPE::INS)
                expected_length += var.sv_length;

            if (abs(static_cast<int32_t>(length(reg)) - expected_length) > options.flanking_region)
            {
                localLog << "------ Skip Read - Length:" << length(reg) << " Qual:" << supporting_records[i].mapQ
                         << " Name: "<< supporting_records[i].qName << endl;
                ++maximum_long_reads;
                continue; // do not use under or oversized region
            }

            appendValue(supporting_sequences, reg);

            localLog << "------ Region: [" << get<0>(region) << "-" << get<1>(region)
                     << "] Length:" << length(reg) << " Qual:" << supporting_records[i].mapQ
                     << " Name: "<< supporting_records[i].qName << endl;
        }

        if (length(supporting_sequences) == 0)
        {
            localLog << "ERROR3: No fitting regions for a " << var.alt_seq
                     << " in region " << var.ref_chrom << ":"
                     << var_ref_pos_sub50 << "-" << var_ref_pos_end_add50
                     << std::endl;

            var.filter = "FAIL3";
            #pragma omp critical
            log_file << localLog.str() << std::endl;
            continue;
        }

        // Build consensus of supporting read regions
        // ---------------------------------------------------------------------
        vector<double> mapping_qualities;
        mapping_qualities.resize(supporting_records.size());
        for (unsigned i = 0; i < supporting_records.size(); ++i)
            mapping_qualities[i] = (supporting_records[i]).mapQ;

        Dna5String cns = build_consensus(supporting_sequences, mapping_qualities);

        localLog << "--- Built a consensus with a MSA of length " << length(cns) << "." << endl;

        // ~supporting_records();   // not used any more
        // ~supporting_sequences(); // not used any more

        Dna5String polished_ref;
        SViperConfig config{options};
        config.ref_flank_length = 500;

        {
            StringSet<Dna5QString> short_reads_1; // reads (first in pair)
            StringSet<Dna5QString> short_reads_2; // mates (second in pair)

            {
                // Extract short reads in region
                // ---------------------------------------------------------------------
                vector<BamAlignmentRecord> short_reads;
                // If the breakpoints are farther apart then illumina-read-length + 2 * flanking-region,
                // then extract reads for each break point separately.
                if (ref_region_end - ref_region_start > options.flanking_region * 2 + options.length_of_short_reads)
                {
                    // extract reads left of the start of the variant [start-flanking_region, start+flanking_region]
                    unsigned e = std::min(ref_length, var.ref_pos + options.flanking_region);
                    viewRecords(short_reads, short_read_bam, short_read_bai, rID_short, ref_region_start, e);
                    cut_down_high_coverage(short_reads, options.mean_coverage_of_short_reads);

                    // and right of the end of the variant [end-flanking_region, end+flanking_region]
                    vector<BamAlignmentRecord> tmp_short_reads;
                    unsigned s = std::max(1, var.ref_pos_end - options.flanking_region);
                    viewRecords(tmp_short_reads, short_read_bam, short_read_bai, rID_short, s, ref_region_end);
                    cut_down_high_coverage(tmp_short_reads, options.mean_coverage_of_short_reads);
                    append(short_reads, tmp_short_reads);
                }
                else
                {
                    // extract reads left of the start of the variant [start-flanking_region, start]
                    viewRecords(short_reads, short_read_bam, short_read_bai, rID_short, ref_region_start, ref_region_end);
                    cut_down_high_coverage(short_reads, options.mean_coverage_of_short_reads);
                }

                if (short_reads.size() < 20)
                {
                    localLog << "ERROR4: Not enough short reads (only " << short_reads.size()
                             << ") for variant of type " << var.alt_seq
                             << " in region " << var.ref_chrom << ":" << ref_region_start
                             << "-" << ref_region_end << std::endl;

                    var.filter = "FAIL4";
                    #pragma omp critical
                    log_file << localLog.str() << std::endl;
                    continue;
                }

                records_to_read_pairs(short_reads_1, short_reads_2, short_reads, short_read_bam, short_read_bai);

                localLog << "--- Extracted " << length(short_reads_1) << " pairs (proper or dummy pairs)." << std::endl;
            } // scope of short reads ends

            // Flank consensus sequence
            // ---------------------------------------------------------------------
            // Before polishing, append a reference flank to the conesnsus such that
            // the reads find a high quality anchor for mapping and pairs are correctly
            // identified.
            Dna5String flanked_consensus = append_ref_flanks(cns, faiIndex,
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

        for (unsigned i = config.ref_flank_length; i < length(config.cov_profile) - config.ref_flank_length; ++i)
            localLog << config.cov_profile[i] << " ";
        localLog << std::endl;

        // Align polished sequence to reference
        // ---------------------------------------------------------------------
        Dna5String ref_part;
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

        BamAlignmentRecord final_record{};
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

        if (options.output_polished_bam)
        {
            std::string read_identifier = (string("polished_var") +
                                           ":" + var.ref_chrom +
                                           ":" + to_string(var.ref_pos) +
                                           ":" + to_string(var.ref_pos_end) +
                                           ":" + var.id);
            final_record.qName = read_identifier;

            #pragma omp critical
            polished_reads.push_back(final_record);
        }

        #pragma omp critical
        log_file << localLog.str() << std::endl;
    } // parallel for loop

    // Write refined variants to output file
    // -------------------------------------------------------------------------
    ofstream output_vcf;          // The polished variant as output
    if (!open_file_success(output_vcf, (options.output_prefix + ".vcf").c_str()))
        return 1;

    for (auto const & header_line : vcf_header)
        output_vcf << header_line << std::endl;

    for (auto & var : variants)
        var.write(output_vcf);


    // Write polished reads if specified to output file
    // -------------------------------------------------------------------------
    if (options.output_polished_bam)
    {
        BamFileOut result_bam(context(*(long_read_file_handles[0]))); // The optional output bam file for polished reads

        if (!open_file_success(result_bam, toCString(options.output_prefix + "_polished_reads.bam")))
        {
            std::cerr << "Did not write resulting bam file." << std::endl;
        }
        else
        {
            writeHeader(result_bam, long_read_header);
            for (auto const & rec : polished_reads)
                writeRecord(result_bam, rec);
        }
    }

    log_file  << "======================================================================" << std::endl
              << "                                 DONE"  << std::endl
              << "======================================================================" << std::endl;

    return 0;
}
