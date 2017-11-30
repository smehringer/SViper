#pragma once

#include <sstream>
#include <regex>
#include <string>

#include <seqan/bam_io.h>

#include <variant.h>

double compute_identity_measure(seqan::BamAlignmentRecord const & record)
{
    seqan::BamTagsDict tagsDict(record.tags);
    int id{-1};
    seqan::findTagKey(id, tagsDict, "NM");

    if (id == -1)
        return -1;

    double nm{0};
    seqan::extractTagValue(nm, tagsDict, id);

    return (1 - (nm/(seqan::length(record.seq) - 10000))) * 100;
}

double compute_fuzzyness_measure(seqan::BamAlignmentRecord const & record)
{
    unsigned count{0};

    // go over the whole cigar and add up the number of insertions and deletions
    unsigned cigar_pos{0};
    int read_pos{-1};                // -1 because ref position is 0 based
    int ref_pos{record.beginPos -1}; // -1 because ref position is 0 based

    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&] (int /**/, int /**/)
                     {
                        if ((record.cigar[cigar_pos]).operation != 'M')
                            ++count; // count insertions and deletions
                        return false; // never stop. Go over whole cigar
                     });

    return (1.0 - ((double)count / ((length(record.seq) - 10000) * 0.3 / 1.3) )) * 100.0;
}

void assign_quality(seqan::BamAlignmentRecord const & record,
                    Variant & variant,
                    bool true_variant)
{
    double identity = compute_identity_measure(record);
    double fuzzyness = compute_fuzzyness_measure(record);

    // subtract variant length from score so to not penalize the variant
    // that one actually expects (e.g. deletion of 1000bp's otherwise screws up
    // the edit distance)
    if (true_variant)
        identity += ((double)variant.sv_length / (length(record.seq) - 10000)) * 100;

    // write seperate scores into info field of the variant
    std::ostringstream ss;
    ss << ";IDENTITIY_SCORE=" << identity
       << ";FUZZYNESS_SCORE=" << fuzzyness;
    variant.info.append(ss.str());

    variant.quality = (identity + fuzzyness) / 2;
}

Variant evaluate_alignment(seqan::BamAlignmentRecord const & record,
                           std::map<std::string, Variant> const & variant_map)
{
    // in the pipeline the final.fa sequnces will have the chr position and if
    // of the variant they correspond to in their name plit by '_'.
    // This way they can be easily matched to the former variant.
    std::stringstream ss;
    ss << record.qName;

    // extract read information
    std::string name;
    std::string chrom;
    std::string pos;
    std::string id;

    std::getline(ss, name, ':');
    std::getline(ss, chrom, ':');
    std::getline(ss, pos, ':');
    std::getline(ss, id, ':');

    // now find the corresponding var
    auto it = variant_map.find(id);

    if (it == variant_map.end())
    {
        std::cerr << "Variant with id '" << id
                  << "' could not be found in vcf file." << std::endl;
        return Variant(string(chrom + "\t0\t" + id + "\tN\tN\t.\tNOTFOUND\tEND=0;\t*\t*"));
    }

    Variant var = it->second;

    if (var.ref_chrom != chrom || std::to_string(var.ref_pos) != pos)
    {
        std::cerr << "Variant with id '" << id << "' does not match with the read"
                  << " name information: "
                  << var.ref_chrom << ":" << var.ref_pos << " "
                  << "vs. " << chrom << ":" << pos << "." << std::endl;
    }

    if (!refine_variant(record, var)) // if this fails, the variant is not supported anymore
    {
        var.filter = "FAIL2";
        assign_quality(record, var, false);
    }
    else
    {
        assign_quality(record, var, true);
    }

    return var;
}
