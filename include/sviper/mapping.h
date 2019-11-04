#pragma once

#include <seqan/align.h>

namespace sviper
{
/*! Stores mapping information for one read pair.
 * This struct stores mapping information in the form of two pairwise alignments.
 * A pairwise alignment in seqan is represented as two gapped strings (Gaps object).
 * Therefore, there are 4 gaps bjects: A gapped read and reference sequence for
 * the alignment of the first read and the same again for the second read (mate).
 * In addition, the booleans read/mate_is_rc indicate weather the read was
 * mapped reverse complementedand and mapQread/mate stores the mapping quality.
 */
struct Mapping_object
{
    typedef seqan::Dna5QString TRead;
    typedef seqan::Dna5String  TRef;
    typedef seqan::Gaps<TRead, seqan::ArrayGaps> TGapsRead;
    typedef seqan::Gaps<TRef,  seqan::ArrayGaps> TGapsRef;

    TRead read;
    TGapsRead gapsRead;
    TGapsRef gapsRef;
    bool read_is_rc{false}; // rc = reverse-complemented
    double mapQRead{0.0};

    TRead mate;
    TGapsRead gapsMate;
    TGapsRef gapsRefMate;
    bool mate_is_rc{false}; // rc = reverse-complemented
    double mapQMate{0.0};

    bool proper_pair{false};

    Mapping_object(TRead r, TRead m, TRef ref):
        read{r}, mate{m}
    {
        seqan::assignSource(gapsRead, read);
        seqan::assignSource(gapsMate, mate);
        seqan::assignSource(gapsRef, ref);
        seqan::assignSource(gapsRefMate, ref);
    }

    bool hasMate() const
    {
        return (length(gapsMate) != 0);
    }

    double mapQ() const
    {
        if (this->hasMate())
            return (mapQRead + mapQMate)/2;
        return mapQRead;
    }
};

//! Computes the position of the last non-gap-character in a gaps object.
template <typename string_type>
inline int gapsEndPos(seqan::Gaps<string_type, seqan::ArrayGaps> const & gaps)
{
    return length(gaps) - seqan::countTrailingGaps(gaps);
}

//! Computes the position of the first non-gap-character in a gaps object.
template <typename string_type>
inline int gapsBeginPos(seqan::Gaps<string_type, seqan::ArrayGaps> const & gaps)
{
    return seqan::countLeadingGaps(gaps);
}

/*! Maps a single sequence (read) to a reference (ref).
 * This function takes a read sequence as input (`read`) and computes two
 * pairwise alignments against the reference (`ref`). One with the original
 * sequence and one with the reverse complemented sequence. The better alignment
 * is chosen and stored in the two output parameters gapsRead/RefOut.
 * @param gapsReadOut The outout parameter in which the gapped read sequence,
 *                    representing the alignment to the ref, is stored.
 * @param gapsRefOut  The outout parameter in which the gapped ref sequence,
 *                    representing the alignment to the read, is stored.
 * @param read        The read sequence to be alignmed.
 * @param ref         The ref sequence to align `read` to.
 */
inline bool map_single_read(seqan::Gaps<seqan::Dna5QString, seqan::ArrayGaps> & gapsReadOut,
                            seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> & gapsRefOut,
                            seqan::Dna5QString & read,
                            seqan::Dna5String & ref,
                            double & mapQ)
{
    typedef seqan::Gaps<seqan::Dna5QString, seqan::ArrayGaps> TGapsRead;
    typedef seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> TGapsRef;

    // BWA MEM 0.7.17 default parameters
    int const MATCH = 1;
    int const MISMATCH = -4;
    int const GAP_OPEN = -6;
    int const GAP_EXT = -1;

    if (length(read) == 0)
        return false;

    seqan::Dna5QString readRC = read; // copy for reverse complementing in place

    seqan::reverseComplement(readRC);

    TGapsRef gapsRef(ref);
    TGapsRef gapsRefRC(ref);
    TGapsRead gapsRead(read);
    TGapsRead gapsReadRC(readRC);

    int score = seqan::globalAlignment(gapsRef, gapsRead,
                                       seqan::Score<int, seqan::Simple>(MATCH, MISMATCH, GAP_EXT, GAP_OPEN),
                                       seqan::AlignConfig<true, false, false, true>(), // semi-global
                                       seqan::AffineGaps());

    int scoreRC = seqan::globalAlignment(gapsRefRC, gapsReadRC,
                                         seqan::Score<int, seqan::Simple>(MATCH, MISMATCH, GAP_EXT, GAP_OPEN),
                                         seqan::AlignConfig<true, false, false, true>(), // semi-global
                                         seqan::AffineGaps());

    if (score > scoreRC) // forward strand
    {
        seqan::copyGaps(gapsReadOut, gapsRead);
        seqan::copyGaps(gapsRefOut, gapsRef);
        mapQ = std::max(0, score); // avoid negative values in pileup
        return false;
    }
    else // reverse strand
    {
        read = readRC;
        seqan::assignSource(gapsReadOut, read); // reassign gaps source to rc read
        seqan::copyGaps(gapsReadOut, gapsReadRC);
        seqan::copyGaps(gapsRefOut, gapsRefRC);
        mapQ = std::max(0, scoreRC); // avoid negative values in pileup.
        return true;
    }
}

/*! Maps read pairs (in reads1 & reads2) to a reference (ref).
 * This function takes two vectors of read sequences as input and expects those,
 * to refer to read pairs (e.g. reads2[i] is the mate of reads1[i]). It stores
 * the result of each mapped pair of reads in a `MappingObject` and returns a
 * list of such.
 * @param reads1 The read sequences of all "first-in-pair" reads.
 * @param reads2 The read sequences of all "second-in-pair" reads.
 * @param ref         The ref sequence to align the read sequences to.
 */
std::vector<Mapping_object> mapping(seqan::StringSet<seqan::Dna5QString> const & reads1,
                                    seqan::StringSet<seqan::Dna5QString> const & reads2,
                                    seqan::Dna5String & ref,
                                    SViperConfig const & config)
{
    SEQAN_ASSERT_EQ(seqan::length(reads1), seqan::length(reads2));

    std::vector<Mapping_object> mobs;
    //mobs.resize(length(reads1)); TODO

    for (unsigned ridx = 0; ridx < seqan::length(reads1); ++ridx) // for every read (pair)
    {
        Mapping_object mob(reads1[ridx], reads2[ridx], ref);

        mob.read_is_rc = map_single_read(mob.gapsRead, mob.gapsRef, mob.read, ref, mob.mapQRead);

        if (mob.hasMate())
        {
            // check if proper pair
            mob.mate_is_rc = map_single_read(mob.gapsMate, mob.gapsRefMate, mob.mate, ref, mob.mapQMate); // TODO:: move to inner if clause
            if ((!mob.read_is_rc && mob.mate_is_rc &&
                    config.insert_size_in_range(gapsEndPos(mob.gapsRead), gapsBeginPos(mob.gapsMate))) ||
                (mob.read_is_rc && !mob.mate_is_rc &&
                    config.insert_size_in_range(gapsEndPos(mob.gapsMate), gapsBeginPos(mob.gapsRead))) )
                mob.proper_pair = true;
        }

        mobs.push_back(mob);
    }

    return mobs;
}
} // namespace sviper
