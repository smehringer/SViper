// This header only file is adapted from the mateHunter program written by
// Snædís <snaedisk@internal.decode.is> . Thank You :-)
#pragma once

#include <iostream>
#include <vector>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>

#include <config.h>

/*! Find the mate of a read by searching in the BAM file.
 *
 */
seqan::BamAlignmentRecord mateHunt(seqan::BamAlignmentRecord const & record,
                                   seqan::BamFileIn & bam_file,
                                   seqan::BamIndex<seqan::Bai> const & bam_index)
{
    seqan::BamAlignmentRecord mate{}; // default/empty record

    bool hasAlignments = false;
    if (!seqan::jumpToRegion(bam_file, hasAlignments, record.rNextId, record.pNext+1, record.pNext+2, bam_index))
    {
#ifndef NDEBUG
        std::cerr << "[ERROR] mateHunt - Could not jump to region "
                  << seqan::contigNames(seqan::context(bam_file))[record.rNextId] << ":" << record.pNext
                  << " to find a records mate." << std::endl;
#endif
        return mate;
    }

    if (!hasAlignments)
    {
#ifndef NDEBUG
        std::cerr << "[ERROR] mateHunt - Could not find any reads in region "
                  << seqan::contigNames(seqan::context(bam_file))[record.rNextId] << ":" << record.pNext
                  << " to find a records mate." << std::endl;
#endif
        return mate;
    }

    while (!seqan::atEnd(bam_file))
    {
        seqan::readRecord(mate, bam_file);

        if (mate.qName == record.qName)
            return mate;

        if (mate.beginPos > record.pNext)
            break;
    }

    // If I arrive here, no mate was found (return statement in while loop),
    // therefore clear mate entry before returning some last read record.
    mate = seqan::BamAlignmentRecord();
    return mate;
}
