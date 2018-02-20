#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <utility> // pair

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/graph_msa.h>

#include <config.h>
#include <basics.h>
#include <mapping.h>

using namespace std;
using namespace seqan;

inline void add_base_to_profile(String<ProfileChar<Dna5, double> > & profile,
                                unsigned pos,
                                Dna5Q base,
                                double mapQ,
                                bool paired,
                                SViperConfig const & config)
{
    if (length(profile) <= pos)
        resize(profile, pos + 1);

    profile[pos].count[ordValue(base)] += config.add_to_profile(getQualityValue(base), mapQ, paired);
}

inline void add_gap_to_profile(String<ProfileChar<Dna5, double> > & profile,
                               unsigned pos,
                               double mapQ,
                               bool paired,
                               SViperConfig const & config)
{
    if (length(profile) <= pos)
        resize(profile, pos + 1);

    profile[pos].count[5] += config.add_to_profile(config.baseQ_mean, mapQ, paired);
}

inline void add_read_to_profile(String<ProfileChar<Dna5, double> > & profile,
                                vector<String<ProfileChar<Dna5, double> >> & ins_profiles,
                                Gaps<Dna5QString, ArrayGaps> const & gapsRead,
                                Gaps<Dna5String, ArrayGaps> const & gapsRef,
                                bool paired,
                                double mapQ,
                                SViperConfig const & config)
{
    SEQAN_ASSERT_EQ(length(gapsRead), length(gapsRef));

    unsigned begin = max(gapsBeginPos(gapsRef), gapsBeginPos(gapsRead));
    unsigned end = min(gapsEndPos(gapsRef), gapsEndPos(gapsRead));
    // append read bases to profile for every position in ref
    for (unsigned idx = begin; idx < end; ++idx)
    {
        unsigned pos_in_ref = toSourcePosition(gapsRef, idx);

        if (isGap(gapsRef, idx)) // this means the read features an insertion
        {
            SEQAN_ASSERT(!isGap(gapsRef, idx - 1)); // always handle full insertion

            unsigned ins_pos{0};

            while(isGap(gapsRef, idx))
            {
                SEQAN_ASSERT(!isGap(gapsRead, idx)); // always handle full insertion
                Dna5Q base = (source(gapsRead))[toSourcePosition(gapsRead, idx)];
                add_base_to_profile(ins_profiles[pos_in_ref], ins_pos, base, mapQ, paired, config);
                ++idx; ++ins_pos;
            }
        }
        else
        {
            if (!isGap(gapsRef, idx + 1)) // no insertion happened between idx and idx + 1
                add_gap_to_profile(ins_profiles[pos_in_ref], 0, mapQ, paired, config);
        }

        if (!isGap(gapsRead, idx))
        {
            Dna5Q base = (source(gapsRead))[toSourcePosition(gapsRead, idx)];
            add_base_to_profile(profile, pos_in_ref, base, mapQ, paired, config);
        }
        else
        {
            add_gap_to_profile(profile, pos_in_ref, mapQ, paired, config);
        }
    }
}

void fill_profile(String<ProfileChar<Dna5, double> > & profile,
                  Gaps<Dna5String, ArrayGaps> & gapsRef,
                  vector<Mapping_object> const & mobs,
                  SViperConfig const & config)
{
    // the reference sequence (consensus sequence) will have a profile over
    // every position. This will account for deletions and substitutions.
    // Insertions are first stored seperately because otherwise the positions
    // in read-ref-alignments are not "synchronized" and hard to handle.
    vector<String<ProfileChar<Dna5, double> >> insertion_profiles;
    insertion_profiles.resize(length(profile)); // there can be an insertion after every position

    for (auto const & mob : mobs) // for every read (pair)
    {
        if (!config.only_proper_pairs || mob.proper_pair)
        {
            add_read_to_profile(profile, insertion_profiles, mob.gapsRead, mob.gapsRef, mob.proper_pair, mob.mapQRead, config);

            if (mob.hasMate())
                add_read_to_profile(profile, insertion_profiles, mob.gapsMate, mob.gapsRefMate, mob.proper_pair, mob.mapQMate, config);
        }
    }

    // add insertions to profile
    unsigned pos{0}; // tracks the view position in gaps space of the reference
    for (auto & pro : insertion_profiles)
    {
        if (length(pro) != 0) // if there is an insertion
        {
            for (auto & ins : pro) // update gap qualitities stored at the first position
                ins.count[5] = pro[0].count[5];
            insert(profile, pos, pro);
            insertGaps(gapsRef, pos, length(pro));
            pos += length(pro) + 1; // +1 because we are also advancing one position anyway
        }
        else
        {
            ++pos;
        }
    }
}

void compute_mappQ_stats(SViperConfig & config,
                         vector<Mapping_object> const & mobs)
{
    if (length(mobs) == 0)
        return;

    double avg{0.0};
    double std{0.0};

    // compute average
    for (auto const & m : mobs)
        avg += m.mapQ();
    avg = avg / length(mobs);

    // compute standard deviation
    for (auto const & m : mobs)
        std += std::pow((m.mapQ() - avg), 2);
    std = std::sqrt(std / length(mobs));

    config.mappQ_mean = avg;
    config.mappQ_std = std;
}

void compute_baseQ_stats(SViperConfig & config,
                         StringSet<Dna5QString> const & reads1,
                         StringSet<Dna5QString> const & reads2)
{
    if (length(reads1) == 0)
        return;

    double avg{0.0};
    double std{0.0};

    // compute average
    for (auto const & read : reads1)
        for (auto const & c : read)
            avg += getQualityValue(c);

    for (auto const & read : reads2)
        for (auto const & c : read)
            avg += getQualityValue(c);

    // Note: this function assume all reads to be of equal length
    avg = avg / (length(reads1) * length(reads1[0]) + length(reads2) * length(reads2[0]));

    // compute standard deviation with avg
    for (auto const & read : reads1)
        for (auto const & c : read)
            std += std::pow((getQualityValue(c) - avg), 2);

    for (auto const & read : reads2)
        for (auto const & c : read)
            std += std::pow((getQualityValue(c) - avg), 2);

    // Note: this function assume all reads to be of equal length
    // and that reads1 and reads2 are of equal length
    std = std::sqrt(std / (2 * length(reads1) * length(reads1[0])));

    config.baseQ_mean = avg;
    config.baseQ_std = std;
}

inline Dna5String consensus_from_profile(String<ProfileChar<Dna5, double> > const & profile,
                                         Gaps<Dna5String, ArrayGaps> const & contigGaps,
                                         SViperConfig & config) // ebuffer at beginning and end
{
    SEQAN_ASSERT(length(profile) == length(contigGaps));

    Dna5String cns;
    unsigned begin = toViewPosition(contigGaps, config.buffer);
    unsigned end   = toViewPosition(contigGaps, length(source(contigGaps)) - config.buffer);

    // first append unpolished bases in the beginning (before begin)
    append(cns, prefix(source(contigGaps), config.buffer));

    // now fix conesnsus
    for (unsigned i = begin; i < end; ++i)
    {
        if (!config.fix_indels)       // if only substitutions are allowed
            if (isGap(contigGaps, i)) // do not insert bases
                continue;

        int idx = getMaxIndex(profile[i]);

        if (!config.score_passes_threshold(profile[i].count[idx]))
        {
            if (!isGap(contigGaps, i))
                appendValue(cns, source(contigGaps)[toSourcePosition(contigGaps, i)]);
            continue;
        }

        if (idx < 5)  // is not gap TODO replace by seqan alphabet size or something like that
        {
            appendValue(cns, Dna5(idx));

            if (isGap(contigGaps, i))
                ++config.inserted_bases;
            else if (source(contigGaps)[toSourcePosition(contigGaps, i)] != Dna5(getMaxIndex(profile[i])))
                ++config.substituted_bases;
        }
        else if (!config.fix_indels) // if idx < 5 but indels should not be fixed
        {
            if (!isGap(contigGaps, i)) // gap to gap alignment can happen because of insertion profiles.
                appendValue(cns, source(contigGaps)[toSourcePosition(contigGaps, i)]);
        }
        else // idx < 5 and fix_indels is true
        {
            if (!isGap(contigGaps, i))
                ++config.deleted_bases;
        }
    }

    // at last, append config.buffer at the end
    append(cns, suffix(source(contigGaps),  length(source(contigGaps)) - config.buffer));

    if (!config.fix_indels)
        SEQAN_ASSERT_EQ(length(source(contigGaps)), length(cns));

    return cns;
}

Dna5String polish(seqan::StringSet<seqan::Dna5QString> const & reads1,
                  seqan::StringSet<seqan::Dna5QString> const & reads2,
                  seqan::Dna5String ref, // TODO:: reference?
                  SViperConfig & config)
{
    seqan::Gaps<Dna5String, ArrayGaps> contigGaps(ref);

    std::vector<Mapping_object> mobs = mapping(reads1, reads2, ref);

    compute_mappQ_stats(config, mobs);
    config.alpha = config.mappQ_mean/config.baseQ_mean;

    seqan::String<seqan::ProfileChar<seqan::Dna5, double> > profile;
    seqan::resize(profile, seqan::length(ref));

    fill_profile(profile, contigGaps, mobs, config);

    // if (config.verbose)
    //     cout << "### MappQ avg: " << config.mappQ_mean << " and std: " << config.mappQ_std
    //          << ". Cov: " << config.mean_coverage << endl;

    // Do not correct bases to the left and right of the consensus sequences otherwise
    // we extend the alignment with low coverage.
    // Also skip 50 bp of the beginning because reads will map poorly there
    ref = consensus_from_profile(profile, contigGaps, config);

    return ref;
}

Dna5String polish_to_perfection(StringSet<Dna5QString> const & reads1,
                                StringSet<Dna5QString> const & reads2,
                                Dna5String ref,
                                SViperConfig & config)
{
    Dna5String old_ref;

    // first, do only fix mismatches
    while (ref != old_ref && config.rounds < 20)
    {
        old_ref = ref; // store prior result
        ref = polish(reads1, reads2, ref, config);
        ++config.rounds;
        // break;
    }

    // after all substitutions has been corrected, correct insertions/deletions
    // a few times with only proper pairs
    config.fix_indels = true;
    for (unsigned i= 0; i < 10; ++i)
    {
        ref = polish(reads1, reads2, ref, config);
        ++config.rounds;
        // break;
    }

    return ref;
}
