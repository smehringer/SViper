#pragma once

#include <seqan/align.h>

/*! Keeps mapping information for one read pair.
 *
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
    bool read_is_rc{false};
    double mapQRead{0.0};

    TRead mate;
    TGapsRead gapsMate;
    TGapsRef gapsRefMate;
    bool mate_is_rc{false};
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


template <typename string_type>
inline int gapsEndPos(seqan::Gaps<string_type, seqan::ArrayGaps> const & gaps)
{
    return length(gaps) - seqan::countTrailingGaps(gaps);
}

template <typename string_type>
inline int gapsBeginPos(seqan::Gaps<string_type, seqan::ArrayGaps> const & gaps)
{
    return seqan::countLeadingGaps(gaps);
}

double compute_mappQ(seqan::Gaps<Dna5QString, seqan::ArrayGaps> & gapsRead,
                     seqan::Gaps<Dna5String, seqan::ArrayGaps>  & gapsRef)
{
    double edit_distance = seqan::countGaps(gapsRead, seqan::countLeadingGaps(gapsRead))
                           - seqan::countTrailingGaps(gapsRead) +
                           seqan::countGaps(gapsRef, seqan::countLeadingGaps(gapsRef))
                           - seqan::countTrailingGaps(gapsRef);
    return (-10 * std::log10((edit_distance + 1)/(2*length(seqan::source(gapsRead))))); // +1 pseusocount
}

inline bool map_single_read(seqan::Gaps<Dna5QString, seqan::ArrayGaps> & gapsReadOut,
                            seqan::Gaps<Dna5String, seqan::ArrayGaps> & gapsRefOut,
                            seqan::Dna5QString & read,
                            seqan::Dna5String & ref)
{
    typedef seqan::Gaps<seqan::Dna5QString, seqan::ArrayGaps> TGapsRead;
    typedef seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> TGapsRef;

    int const MATCH = 3;
    int const MISMATCH = -2;
    int const GAP_OPEN = -3;
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
        return false;
    }
    else // reverse strand
    {
        read = readRC;
        seqan::assignSource(gapsReadOut, read); // reassign gaps source to rc read
        seqan::copyGaps(gapsReadOut, gapsReadRC);
        seqan::copyGaps(gapsRefOut, gapsRefRC);
        return true;
    }
}

vector<Mapping_object> mapping(seqan::StringSet<seqan::Dna5QString> const & reads1,
                               seqan::StringSet<seqan::Dna5QString> const & reads2,
                               seqan::Dna5String & ref)
{
    SEQAN_ASSERT_EQ(seqan::length(reads1), seqan::length(reads2));

    vector<Mapping_object> mobs;
    //mobs.resize(length(reads1)); TODO

    #pragma omp parallel for
    for (unsigned ridx = 0; ridx < seqan::length(reads1); ++ridx) // for every read (pair)
    {
        Mapping_object mob(reads1[ridx], reads2[ridx], ref);

        mob.read_is_rc = map_single_read(mob.gapsRead, mob.gapsRef, mob.read, ref);
        mob.mate_is_rc = map_single_read(mob.gapsMate, mob.gapsRefMate, mob.mate, ref); // TODO:: move to inner if clause
        mob.mapQRead = compute_mappQ(mob.gapsRead, mob.gapsRef);

        if (mob.hasMate())
        {
            mob.mapQMate = compute_mappQ(mob.gapsMate, mob.gapsRefMate);
            // check if proper pair
            if ((!mob.read_is_rc && mob.mate_is_rc && gapsEndPos(mob.gapsRead) < gapsBeginPos(mob.gapsMate)) ||
                (mob.read_is_rc && !mob.mate_is_rc && gapsEndPos(mob.gapsMate) < gapsBeginPos(mob.gapsRead)) ) // TODO check insert size
                mob.proper_pair = true;
        }

        #pragma omp critical
        mobs.push_back(mob);
    }

    return mobs;
}
