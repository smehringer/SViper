#pragma once

#include <sstream>
#include <regex>
#include <string>

#include <seqan/bam_io.h>

#include <variant.h>

/*! Computes the identity score for evaluation.
 * The identity score shall represent how similar the polished sequence is to
 * the reference. The score therefore is computed by taking the edit
 * distance (NM tag) and dividing it by the length of the sequence.
 * This actual identity score exculdes the one (true) variant found in the
 * polished sequence, since this is expected and the correpsonging edit distance
 * will be added to the measurement afterwards.
 * @param record The alignment record containing the alignment of the polished
 *               sequence to the reference.
 * @return Returns measurement as double in range [0,100].
 */
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


/*! Computes the fuzzyness score for evaluation.
 * Since the edit distance of one big event (DEL/INS) can be the same as several
 * small ones, a second measurement is introduced for the evaluation of polishing.
 * The fuzzyness score represents the amount of events found in the polished
 * sequence alignment to the reference. It is computed by dividing the number of
 * events (DEL/INS) by the length of the alignment, subtracted off 1 (we want the
 * higher the better, analogous to the identity score).
 * NOTE: The higher the score the less fuzzy is the alignment and the better did
 *       did polishing work.
 * @param record The alignment record containing the alignment of the polished
 *               sequence to the reference.
 * @return Returns measurement as double in range [0,100].
 */
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

/*! Assigns a variant the average of identity and fuzzyness score.
 * Based on the polsihed sequence alignment given in record
 * @param record The alignment record containing the alignment of the polished
 *               sequence to the reference.
 * @param variant The variant to be scored (Quality member will be replaced).
 */
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
