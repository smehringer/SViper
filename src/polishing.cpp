#include <iostream>
#include <fstream>
#include <cmath>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

#include <basics.h>
#include <merge_split_alignments.h>
#include <helper_functions.h>
#include <polishing.h>

using namespace std;
using namespace seqan;

struct CmdOptions
{
    bool verbose{false};
    bool veryVerbose{false};
    int flanking_region{400}; // size of flanking region for breakpoints
    string long_read_file_name;
    string short_read_file_name;
    string candidate_file_name;
    string out_vcf_file_name;
    string reference_file_name;
    string log_file_name;
};

ArgumentParser::ParseResult
parseCommandLine(CmdOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("Build Consensus");
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
        "k", "flanking-region",
        "The flanking region in bp's around a breakpoint to be considered for polishing",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "o", "output-vcf",
        "A filename for the polished vcf file.",
        seqan::ArgParseArgument::INPUT_FILE, "VCF_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "g", "log-file",
        "A filename for the log file.",
        seqan::ArgParseArgument::INPUT_FILE, "TXT_FILE"));

    addOption(parser, seqan::ArgParseOption(
        "v", "verbose",
        "Turn on detailed information about the process."));

    // addOption(parser, seqan::ArgParseOption(
    //     "vv", "very-verbose",
    //     "Turn on detailed information about the process and write out intermediate results."));

    setRequired(parser, "c");
    setRequired(parser, "l");
    setRequired(parser, "s");
    setRequired(parser, "r");

    setMinValue(parser, "k", "50");
    setMaxValue(parser, "k", "1000");
    setDefaultValue(parser, "k", "400");

    setDefaultValue(parser, "g", "polishing.log");

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
    getOptionValue(options.out_vcf_file_name, parser, "output-vcf");
    getOptionValue(options.log_file_name, parser, "log-file");
    getOptionValue(options.flanking_region, parser, "flanking-region");
    options.verbose = isSet(parser, "verbose");
    // options.veryVerbose = isSet(parser, "very-verbose");

    //if (options.veryVerbose)
    //    options.verbose = true;

    if (options.out_vcf_file_name.empty())
        options.out_vcf_file_name = options.candidate_file_name + "polished.vcf";

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
    ifstream input_vcf;              // The candidate variants to polish
    ofstream output_vcf;             // The polished variant as output
    FaiIndex faiIndex;               // The reference sequence fasta index
    BamHeader long_read_header;      // The bam header object needed to fill bam context
    BamHeader short_read_header;     // The bam header object needed to fill bam context
    BamFileIn long_read_bam;         // The long read bam file
    BamFileIn short_read_bam;        // THe short read bam file
    BamIndex<Bai> long_read_bai;     // The bam index to the long read bam file
    BamIndex<Bai> short_read_bai;    // The bam index to the short read bam file
    SeqFileOut final_fa("final.fa"); // The current output fasta file to be mapped
    ofstream log_file;               // The file to log all information

    if (!open_file_success(input_vcf, options.candidate_file_name.c_str()) ||
        !open_file_success(output_vcf, options.out_vcf_file_name.c_str()) ||
        !open_file_success(long_read_bam, options.long_read_file_name.c_str()) ||
        !open_file_success(short_read_bam, options.short_read_file_name.c_str()) ||
        !open_file_success(long_read_bai, (options.long_read_file_name + ".bai").c_str()) ||
        !open_file_success(short_read_bai, (options.short_read_file_name + ".bai").c_str()) ||
        !open_file_success(faiIndex, options.reference_file_name.c_str()) ||
        !open_file_success(log_file, options.log_file_name.c_str()))
        return 1;

    std::cout << "----------------------------------------------------------------------" << std::endl
              << "START polishing variants in of file" << options.candidate_file_name << std::endl
              << "----------------------------------------------------------------------" << std::endl;

    // Read bam file header
    // -------------------------------------------------------------------------
    // This is neccessary for the bam file context (containing reference ids and
    // more) is correctly initialized.
    readHeader(long_read_header, long_read_bam);
    readHeader(short_read_header, short_read_bam);

    // Scip VCF header
    // -------------------------------------------------------------------------
    std::string line{"#"};
    while (line[0] == '#' && getline(input_vcf, line)) {}

    do // per Variant
    {

        Variant var{line};

        if (var.alt_seq != "<DEL>" && var.alt_seq != "<INS>")
        {
            std::cout << "[ SKIP ] Variant " << var.ref_chrom << ":"
                      << var.ref_pos << " of tpye " << var.alt_seq
                      << " is currently not supported." << std::endl;
            continue;
        }

        // Compute reference length, start and end position of the variant
        // ---------------------------------------------------------------------
        unsigned ref_fai_idx = 0;
        if (!getIdByName(ref_fai_idx, faiIndex, var.ref_chrom))
        {
            log_file << "[ERROR]: FAI index has no entry for reference name"
                     << var.ref_chrom << std::endl;
            return 1;
        }

        int ref_length = sequenceLength(faiIndex, ref_fai_idx);
        int ref_region_start = max(0, var.ref_pos - options.flanking_region);
        int ref_region_end   = min(ref_length, var.ref_pos_end + options.flanking_region);

        // Extract long reads
        // ---------------------------------------------------------------------
        vector<BamAlignmentRecord> ont_reads = view_bam(long_read_bam,
                                                        long_read_bai,
                                                        var.ref_chrom,
                                                        var.ref_pos - 50,
                                                        var.ref_pos_end + 50,
                                                        true);
        if (ont_reads.size() == 0)
        {
            log_file << "[ ERROR1 ] No long reads in reference region "
                     << var.ref_chrom << " [" << ref_region_start << "-"
                     << ref_region_end << "]" << std::endl;
            continue;
        }

        log_file << "--- Found " << ont_reads.size() << " read(s) in reference region ["
                 << ref_region_start << "-" << ref_region_end << std::endl;


        // Merge supplementary alignments to primary
        // ---------------------------------------------------------------------
        ont_reads = merge_alignments(ont_reads);

        log_file << "--- After merging " << ont_reads.size() << " read(s) remain(s)." << endl;

        // Search for supporting reads
        // ---------------------------------------------------------------------
        vector<BamAlignmentRecord> supporting_records;
        // TODO check if var is not empty!

        log_file << "--- Searchin in reference region ["
                 << (int)(var.ref_pos - DEV_POS * var.sv_length) << "-"
                 << (int)(var.ref_pos + var.sv_length + DEV_POS * var.sv_length) << "]"
                 << " for a variant of type " << var.sv_type
                 << " of length " << (int)(var.sv_length - DEV_SIZE * var.sv_length) << "-"
                 << (int)(var.sv_length + DEV_SIZE * var.sv_length) << " bp's" << endl;

        for (auto const & rec : ont_reads)
            if (is_supporting(rec, var))
                supporting_records.push_back(rec);

        if (supporting_records.size() == 0)
        {
            std::cout << "[ ERROR2 ] No supporting reads in reference region "
                      << var.ref_chrom << " [" << ref_region_start << "-"
                      << ref_region_end << "]" << std::endl;
            continue;
        }

        log_file << "--- After searching " << supporting_records.size()
                 << " read(s) remain." << std::endl;

        // Crop fasta sequence of each supporting read for consensus
        // ---------------------------------------------------------------------
        log_file << "--- Cropping read regions +-" << options.flanking_region << " around variants." << endl;
        log_file << "--- Referene Region: [" << ref_region_start << "-" << ref_region_end << "]" << endl;
        log_file << "--- Corresponding Regions & Lengths in long reads:" << endl;


        StringSet<DnaString> supporting_sequences;

        for (auto rec : supporting_records)
        {
            auto region = get_read_region_boundaries(rec, ref_region_start, ref_region_end);
            DnaString reg = infix(rec.seq, get<0>(region), get<1>(region));
            appendValue(supporting_sequences, reg);

            log_file << "------ [" << get<0>(region) << "-" << get<1>(region) << "] L:" << length(reg) << endl;
        }

        // Build consensus of supporting read regions
        // ---------------------------------------------------------------------
        log_file << "--- Building a consensus with a multiple sequence alignment." << endl;

        vector<double> mapping_qualities;
        mapping_qualities.resize(supporting_records.size());
        for (unsigned i = 0; i < supporting_records.size(); ++i)
        {
            mapping_qualities[i] = (supporting_records[i]).mapQ;
        }

        DnaString consensus_sequence = build_consensus(supporting_sequences, mapping_qualities);

        // Extract short reads in region
        // ---------------------------------------------------------------------
        StringSet<Dna5QString> short_reads_1; // first in pair
        StringSet<Dna5QString> short_reads_2; // second in pair
        StringSet<CharString> short_ids_1; // first in pair
        StringSet<CharString> short_ids_2; // first in pair

        vector<BamAlignmentRecord> short_reads = view_bam(short_read_bam,
                                                          short_read_bai,
                                                          var.ref_chrom,
                                                          ref_region_start,
                                                          ref_region_end,
                                                          false);

        records_to_read_pairs(short_reads_1, short_ids_1,
                              short_reads_2, short_ids_2,
                              short_reads, short_read_bam, short_read_bai);

        log_file << "--- Extracted " << short_reads.size() << " short reads in"
                 << length(short_reads_1) << " pairs (proper or dummy pairs)."
                 << std::endl;

        // Flank consensus sequence
        // ---------------------------------------------------------------------
        // Before polishing, append a reference flank of length 150 (illumina
        // read length) to the conesnsus such that the reads find a high quality
        // anchor for mapping.
        String<char> id = "consensus";
        Dna5String ref = append_ref_flanks(consensus_sequence, faiIndex,
                                           ref_fai_idx, ref_length,
                                           ref_region_start, ref_region_end, 150);

        // Polish flanked consensus sequence with short reads
        // ---------------------------------------------------------------------
        ConsensusConfig config{}; // default
        config.verbose = options.verbose;
        compute_baseQ_stats(config, short_reads_1, short_reads_2);

        Dna5String polished_ref = polish_to_perfection(short_reads_1, short_ids_1,
                                                       short_reads_2, short_ids_2,
                                                       ref, id, config);

        log_file << "DONE POLISHING: Total of "
                 << config.substituted_bases << " substituted, "
                 << config.deleted_bases     << " deleted and "
                 << config.inserted_bases    << " inserted bases."
                 << std::endl;

        // Flank polished sequence
        // ---------------------------------------------------------------------
        // Append large reference flank for better mapping results of long reads.
        Dna5String final_sequence = append_ref_flanks(polished_ref,
                                                      faiIndex, ref_fai_idx,
                                                      ref_length,
                                                      ref_region_start - 150,
                                                      ref_region_end + 150,
                                                      4850);

        // Print polished and flanked sequence to file
        // ---------------------------------------------------------------------
        std::string read_identifier = (string("final") +
                                       "_" + var.ref_chrom +
                                       "_" + to_string(var.ref_pos) +
                                       "_" + var.id); // this name is important for evaluating
        writeRecord(final_fa, read_identifier, final_sequence);

        std::cout << "[ SUCCESS ] Polished variant " << var.ref_chrom << ":"
                  << var.ref_pos << "\t\t" << var.alt_seq << std::endl;
    }
    while (getline(input_vcf, line));

    return 0;
}
