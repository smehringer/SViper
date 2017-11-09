#include <iostream>
#include <fstream>
#include <cmath>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

#include <helper_functions.h>
#include <polishing.h>

using namespace std;
using namespace seqan;

struct CmdOptions
{
    bool verbose{false};
    int flanking_region{400}; // size of flanking region for breakpoints
    string bam_name;
    string vcf_name;
    string reference_name;
};

ArgumentParser::ParseResult
parseCommandLine(CmdOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("Build Consensus");
    setVersion(parser, "1.0.0");

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "REFERENCE"));

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "BAMFILE"));

    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INPUT_FILE, "VCFFILE"));

    // Define Options
    addOption(parser, seqan::ArgParseOption(
        "l", "flanking-region",
        "The flanking region in bp's around a breakpoint to be considered for polishing",
        seqan::ArgParseArgument::INTEGER, "INT"));

    addOption(parser, seqan::ArgParseOption(
        "v", "verbose",
        "Turn on detailed information about the process"));

    setMinValue(parser, "l", "50");
    setMaxValue(parser, "l", "1000");
    setDefaultValue(parser, "l", "400");

    // Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getOptionValue(options.flanking_region, parser, "flanking-region");
    options.verbose = isSet(parser, "verbose");

    seqan::getArgumentValue(options.reference_name, parser, 0);
    seqan::getArgumentValue(options.bam_name, parser, 1);
    seqan::getArgumentValue(options.vcf_name, parser, 2);

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
    // -------------------------------------------------------------------------
    // Parse Command Line Arguments
    // -------------------------------------------------------------------------
    CmdOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    // -------------------------------------------------------------------------
    // Set up
    // -------------------------------------------------------------------------


    // Check files
    // -------------------------------------------------------------------------
    ifstream vcf_file(options.vcf_name.c_str());
    BamFileIn bamfileIn;
    FaiIndex faiIndex;

    if (!vcf_file.is_open())
    {
        std::cerr << "ERROR: Could not open " << options.vcf_name << std::endl;
        return 1;
    }

    if (!open(bamfileIn, options.bam_name.c_str()))
    {
        std::cerr << "ERROR: Could not open " << options.bam_name << std::endl;
        return 1;
    }

    if (!open(faiIndex, options.reference_name.c_str()))
    {
        std::cout << "ERROR: Could not load FAI index " << options.reference_name << ".fai\n";
        return 1;
    }

    // Read Variant
    // -------------------------------------------------------------------------
    std::string line;
    Variant var;
    while (getline(vcf_file, line)) // get first variant in file
    {
        if (line[0] != '#') // no header line
        {
            var = Variant(line); // read from file
            break;
        }
    }

    // Compute some parameters
    // -------------------------------------------------------------------------
    unsigned ref_fai_idx = 0;
    if (!getIdByName(ref_fai_idx, faiIndex, var.ref_chrom))
    {
        std::cout << "ERROR: FAI index has no entry for " << var.ref_chrom<< std::endl;
        return 1;
    }
    int ref_length = sequenceLength(faiIndex, ref_fai_idx);
    int ref_region_start = max(0, var.ref_pos - options.flanking_region);
    int ref_region_end   = min(ref_length, var.ref_pos + var.sv_length + options.flanking_region); // TODO remove dirty cast

    // -------------------------------------------------------------------------
    // Merge supplementary alignments to primary
    // -------------------------------------------------------------------------
    vector<BamAlignmentRecord> merged_records;
    BamAlignmentRecord record;
    vector<BamAlignmentRecord> record_group;

    // empty file must be check here otherwise the first read record will fail
    if (atEnd(bamfileIn))
        return 0;

    while (!atEnd(bamfileIn))
    {
        readRecord(record, bamfileIn);

        if (!record_group.empty() &&
            (record_group[record_group.size() - 1]).qName != record.qName)
        {
            BamAlignmentRecord merged_record = process_record_group(record_group);
            merged_records.push_back(merged_record);
            record_group.clear();
            record_group.push_back(record);
        }
        else
        {
            record_group.push_back(record);
        }
    }
    //process last group
    BamAlignmentRecord merged_record = process_record_group(record_group);
    merged_records.push_back(merged_record);
    if (options.verbose)
        cout << "--- After merging, " << merged_records.size() << " record(s) remain(s)." << endl;

    // -------------------------------------------------------------------------
    // Search for supporting reads
    // -------------------------------------------------------------------------
    vector<BamAlignmentRecord> supporting_records;
    // TODO check if var is not empty!
    if (options.verbose)
        cout << "Searchin in reference region ["
             << (int)(var.ref_pos - DEV_POS * var.sv_length) << "-"
             << (int)(var.ref_pos + var.sv_length + DEV_POS * var.sv_length) << "]"
             << " for a variant of type " << var.sv_type
             << " of length " << (int)(var.sv_length - DEV_SIZE * var.sv_length) << "-"
             << (int)(var.sv_length + DEV_SIZE * var.sv_length) << " bp's" << endl;

    for (auto rec : merged_records)
        if (is_supporting(rec, var))
            supporting_records.push_back(rec);

    if (options.verbose)
        cout << "--- After searching for variant <"
             << var.sv_type << ":" << var.ref_chrom << ":" << var.ref_pos << ":" << var.sv_length
             << "> " << supporting_records.size()
             << ", record(s) remain(s)." << endl;

    // -------------------------------------------------------------------------
    // Crop fasta sequence of each supporting read for consensus
    // -------------------------------------------------------------------------
    if (options.verbose)
    {
        cout << "--- Cropping read regions +-" << options.flanking_region << " around variants." << endl;
        cout << "--- Regions & Lengths:" << endl;
    }

    StringSet<DnaString> supporting_sequences;

    for (auto rec : supporting_records)
    {
        auto region = get_read_region_boundaries(rec, ref_region_start, ref_region_end);
        DnaString reg = infix(rec.seq, get<0>(region), get<1>(region));
        appendValue(supporting_sequences, reg);
        if (options.verbose)
        {
            cout << "------ [" << get<0>(region) << "-" << get<1>(region) << "] L:" << length(reg) << endl;
        }
    }

    // -------------------------------------------------------------------------
    // Build consensus of supporting read regions
    // -------------------------------------------------------------------------
    if (options.verbose)
        cout << "--- Building a consensus with a multiple sequence alignment." << endl;

    if (supporting_records.size() == 0)
    {
        cerr << "[ERROR] No supporting reads found." << endl;
        return 1;
    }

    vector<double> mapping_qualities;
    mapping_qualities.resize(supporting_records.size());
    for (unsigned i = 0; i < supporting_records.size(); ++i)
    {
        mapping_qualities[i] = (supporting_records[i]).mapQ;
    }

    DnaString consensus_sequence = build_consensus(supporting_sequences, mapping_qualities);

    // -------------------------------------------------------------------------
    // Write consensus sequence to output file
    // -------------------------------------------------------------------------
    if (options.verbose)
    {
        cout << "--- Writing consensus sequence to file consensus.fa ." << endl;
        SeqFileOut seqFileOut("consensus.fa");
        writeRecord(seqFileOut, "consensus", consensus_sequence);

    }

    // -------------------------------------------------------------------------
    // Polish conensus
    // -------------------------------------------------------------------------

    // before polishing, append a reference flank of length 150 (illumina read length)
    // to the conesnsus such that the reads find an anchor for mapping.
    Dna5String leftFlank;
    Dna5String rightFlank;
    readRegion(leftFlank, faiIndex, ref_fai_idx, max(0, ref_region_start - 150), ref_region_start);
    readRegion(rightFlank, faiIndex, ref_fai_idx, ref_region_end, min(ref_region_end + 150, ref_length));

    String<char> id = "consensus";
    Dna5String ref = leftFlank;
    append(ref, consensus_sequence);
    append(ref, rightFlank);

    // now read in illumina files used for polishing
    SeqFileIn illumina_file_pair1("illumina.paired1.fastq");
    SeqFileIn illumina_file_pair2("illumina.paired2.fastq");

    StringSet<String<char>> ids1;
    StringSet<Dna5String> reads1;
    StringSet<String<char>> quals1;

    StringSet<String<char>> ids2;
    StringSet<Dna5String> reads2;
    StringSet<String<char>> quals2;

    readRecords(ids1, reads1, quals1, illumina_file_pair1);
    readRecords(ids2, reads2, quals2, illumina_file_pair2);

    Dna5String old_ref;
    ConsensusConfig config{}; // default

    unsigned round{1};

    while (ref != old_ref)
    {
        cout << "-------------------------------- SNV ROUND " << round
             << "--------------------------------" << endl;

        old_ref = ref; // store prior result

        ref = polish(reads1, ids1, quals1, reads2, ids2, quals2, ref, id, config);

        ++round;
        // break;
    }

    // after all substitutions has been corrected, correct insertions/deletions a few times with only proper pairs
    config.fix_indels = true;
    for (unsigned i= 0; i < 3; ++i)
    {
        cout << "-------------------------------- INDEL ROUND " << i
             << "--------------------------------" << endl;

        ref = polish(reads1, ids1, quals1, reads2, ids2, quals2, ref, id, config);

        ++round;
        // break;
    }

    config.only_proper_pairs = false;
    old_ref = ""; // reset
    while (ref != old_ref)
    {
        cout << "-------------------------------- INDEL NO PAIRS ROUND " << round
             << "--------------------------------" << endl;

        old_ref = ref; // store prior result

        ref = polish(reads1, ids1, quals1, reads2, ids2, quals2, ref, id, config);

        ++round;
        // break;
    }

    // append large reference flank for better mapping results of long reads
    Dna5String leftLargeFlank;
    Dna5String rightLargeFlank;
    readRegion(leftLargeFlank, faiIndex, ref_fai_idx,
               max(0, ref_region_start - 5000), max(0, ref_region_start - 150));
    readRegion(rightLargeFlank, faiIndex, ref_fai_idx,
               min(ref_region_end + 150, ref_length), min(ref_region_end + 5000, ref_length));

    Dna5String final_sequence = leftLargeFlank;
    append(final_sequence, ref);
    append(final_sequence, rightLargeFlank);

    SeqFileOut final("consensus-polished.fa");
    writeRecord(final, "polished_consensus", final_sequence);

    return 0;
}
