// This header only file is adapted from the mateHunter program written by
// Snædís <snaedisk@internal.decode.is> . Thank You :-)

#include <iostream>
#include <vector>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

BamAlignmentRecord mateHunt(BamAlignmentRecord const & record,
                             BamFileIn & bam_file,
                             BamIndex<Bai> const & bam_index)
{
    BamAlignmentRecord mate{}; // default/empty record

    bool hasAlignments = false;
    if (!jumpToRegion(bam_file, hasAlignments, record.rNextId, record.pNext, record.pNext+1, bam_index))
    {
        std::cerr << "ERROR: Could not jump to region "
                  << record.rNextId << ":" << record.pNext
                  << " to find a records mate." << std::endl;
        return mate;
    }

    if (!hasAlignments)
    {
        std::cerr << "ERROR: Could not find any reads in region "
                  << record.rNextId << ":" << record.pNext
                  << " to find a records mate." << std::endl;
        return mate;
    }

    while (!atEnd(bam_file))
    {
        readRecord(mate, bam_file);

        if (mate.qName == record.qName)
            return mate;

        if (mate.beginPos > record.pNext)
            break;
    }

    // If I arrive here, no mate was found (return statement in while loop),
    // therefore clear mate entry before returning some last read record.
    mate = BamAlignmentRecord();
    return mate;
}
