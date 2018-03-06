#pragma once

#include <sstream>

#include <seqan/bam_io.h>

#include <variant.h>

/*! Assigns a variant the average of identity and fuzzyness score.
 * Based on the polsihed sequence alignment given in record
 * @param record The alignment record containing the alignment of the polished
 *               sequence to the reference.
 * @param variant The variant to be scored (Quality member will be replaced).
 */
void assign_quality(seqan::BamAlignmentRecord & record,
                    Variant & variant,
                    double align_score,
                    SViperConfig const & config)
{
    unsigned count{0};
    std::pair<int, unsigned> best_fitting_event{length(record.seq), 0};

    // go over the whole cigar and add up the number of insertions and deletions
    // and look for the best fitting variant
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
                        {
                            ++count; // count insertions and deletions
                            if (is_same_sv_type((record.cigar[cigar_pos]).operation, variant.sv_type))
                            {
                                if ((record.cigar[cigar_pos]).count > 30)
                                {
                                    if (get<0>(best_fitting_event) > std::abs(ref_pos - variant.ref_pos))
                                        best_fitting_event = {std::abs(ref_pos - variant.ref_pos),
                                                              (record.cigar[cigar_pos]).count};
                                }
                            }
                        }
                        return false; // never stop. Go over whole cigar
                     });

    // Helper computations
    //--------------------
    // maximum number of bases that can be matched
    double max_num_match = std::min(static_cast<unsigned>(length(record.seq)), seqan::getAlignmentLengthInRef(record));
    // The reference flanks should always match and thus can be subtracted the score to avoid bias
    double ref_flank_score = 2*config.ref_flank_length * config.MM;
    // The penalty of an event of length get<1>(best_fitting_event) == convex cgst cost model
    double event_penalty{0.0};
    if (get<1>(best_fitting_event) != 0) // only if a fitting event was found
        event_penalty = config.GO + config.GO/get<1>(best_fitting_event);
    // The maximum alignment score that an alignment could reach,
    // => all possible matches, but with expected event penalty
    double max_align_score = (max_num_match * config.MM) + event_penalty;

    // Score computation
    // -----------------
    double fuzzyness;
    double identity;

    if (ref_flank_score >= max_align_score) // should not happen since we flanked the variant region
    {
        SEQAN_FAIL("[ERROR] ref_flank_score %d is <= max_align_score %d.", ref_flank_score, max_align_score);
        identity = 0;
    }
    else
    {
        // Identity score: compute the percentage that the alignment score reached
        // regarding the maximum possible score, based on an all match alignment.
        // Subtract the match score for the reference flanks because those will always match.
        // Subtract the length of the best fitting event from the maximum match alignment
        // since this is expected and should not weaken the score.
        identity = std::max(0.0, (align_score - ref_flank_score) * 100.0 / (max_align_score - ref_flank_score));
    }

    if ((length(record.seq) - 2*config.ref_flank_length) == 0) // should not happen but avoid dividing by zero just in case
    {
        SEQAN_FAIL("[ERROR] length(record.seq) %d is == 2*config.ref_flank_length) %d.",
                   length(record.seq), (2*config.ref_flank_length));
        fuzzyness = 0;
    }
    else
    {
        fuzzyness = (1.0 - ((double)count / ((length(record.seq) - 2*config.ref_flank_length) * 0.06) )) * 100.0;
    }

    // write separate scores into info field of the variant
    std::ostringstream ss;
    ss << ";IDENTITIY_SCORE=" << identity
       << ";FUZZYNESS_SCORE=" << fuzzyness;
    variant.info.append(ss.str());

    variant.quality = (identity + fuzzyness) / 2;
    record.mapQ = variant.quality;
}
