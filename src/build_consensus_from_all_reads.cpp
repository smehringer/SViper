#include <iostream>
#include <fstream>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

#include <helper_functions.h>

using namespace std;
using namespace seqan;

struct CmdOptions
{
    int flanking_region{400}; // size of flanking region for breakpoints
    string bam_name;
    string vcf_name;
};

ArgumentParser::ParseResult
parseCommandLine(CmdOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    seqan::ArgumentParser parser("Build Consensus");
    setVersion(parser, "1.0.0");

    // We require one argument.
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
        "Turn on detailed information about the process",
        seqan::ArgParseArgument::BOLEAN));

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

    seqan::getArgumentValue(options.bam_name, parser, 0);
    seqan::getArgumentValue(options.vcf_name, parser, 1);

    return seqan::ArgumentParser::PARSE_OK;
}

int main(int argc, char const ** argv)
{
    CmdOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;

    ifstream vcf_file(options.vcf_name.c_str());
    BamFileIn bamfileIn;

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

    // empty file must be check here otherwise the first read record will fail
    if (atEnd(bamfileIn))
        return 0;

    // -------------------------------------------------------------------------
    // Merge supplementary alignments to primary
    // -------------------------------------------------------------------------
    vector<BamAlignmentRecord> merged_records;
    BamAlignmentRecord record;
    vector<BamAlignmentRecord> record_group;

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
    if (verbose)
        cout << "--- After merging, " << merged_records.size() << " record(s) remain(s)." << endl;

    // -------------------------------------------------------------------------
    // Search for supporting reads
    // -------------------------------------------------------------------------
    vector<BamAlignmentRecord> supporting_records;
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
    // TODO check if var is not empty!
    if (verbose)
        cout << "Searchin in reference region ["
             << (int)(var.ref_pos - DEV_POS * var.sv_length) << "-"
             << (int)(var.ref_pos + var.sv_length + DEV_POS * var.sv_length) << "]"
             << " for a variant of type " << var.sv_type
             << " of length " << (int)(var.sv_length - DEV_SIZE * var.sv_length) << "-"
             << (int)(var.sv_length + DEV_SIZE * var.sv_length) << " bp's" << endl;

    for (auto rec : merged_records)
        if (is_supporting(rec, var))
            supporting_records.push_back(rec);

    if (verbose)
        cout << "--- After searching for variant <"
             << var.sv_type << ":" << var.ref_chrom << ":" << var.ref_pos << ":" << var.sv_length
             << "> " << supporting_records.size()
             << ", record(s) remain(s)." << endl;

    // -------------------------------------------------------------------------
    // Crop fasta sequence of each supporting read for consensus
    // -------------------------------------------------------------------------
    if (verbose)
    {
        cout << "--- Cropping read regions +-" << options.flanking_region << " around variants." << endl;
        cout << "--- Regions & Lengths:" << endl;
    }

    StringSet<DnaString> supporting_sequences;

    for (auto rec : supporting_records)
    {
        auto region = get_read_region_boundaries(rec, var.ref_pos - options.flanking_region, var.ref_pos + var.sv_length + options.flanking_region);
        DnaString reg = infix(rec.seq, get<0>(region), get<1>(region));
        appendValue(supporting_sequences, reg);
        if (verbose)
        {
            cout << "------ [" << get<0>(region) << "-" << get<1>(region) << "] L:" << length(reg) << endl;
        }
    }

    // -------------------------------------------------------------------------
    // Build consensus of supporting read regions
    // -------------------------------------------------------------------------
    if (verbose)
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
    if (verbose)
       cout << "--- Writing consensus sequence to file consensus.fa ." << endl;

    SeqFileOut seqFileOut("consensus.fa");
    writeRecord(seqFileOut, "conssensus", consensus_sequence);

    return 0;
}
