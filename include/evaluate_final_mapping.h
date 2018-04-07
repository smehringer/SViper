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
                    SViperConfig const & config)
{
    // Score computation
    // -----------------
    double error_rate = ((double)length(record.cigar) - 1.0)/ (config.flanking_region * 2.0);
    double fuzzyness = (1.0 - error_rate/0.15) * 100.0;

    variant.quality = fuzzyness;
    record.mapQ = variant.quality;
}
