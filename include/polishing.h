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

inline void add_base_to_profile(seqan::String<seqan::ProfileChar<seqan::Dna5, double> > & profile,
                                std::vector<unsigned> & cov_profile,
                                unsigned pos,
                                seqan::Dna5Q base,
                                double mapQ,
                                bool paired,
                                SViperConfig const & config)
{
    if (seqan::length(profile) <= pos)
    {
        seqan::resize(profile, pos + 1);
        seqan::resize(cov_profile, pos + 1);
    }

    profile[pos].count[ordValue(base)] += config.add_to_profile(seqan::getQualityValue(base), mapQ, paired);
    ++cov_profile[pos];
}

inline void add_gap_to_profile(seqan::String<seqan::ProfileChar<seqan::Dna5, double> > & profile,
                               std::vector<unsigned> & cov_profile,
                               unsigned pos,
                               double mapQ,
                               bool paired,
                               SViperConfig const & config)
{
    if (seqan::length(profile) <= pos)
    {
        seqan::resize(profile, pos + 1);
        seqan::resize(cov_profile, pos + 1);
    }

    profile[pos].count[5] += config.add_to_profile(config.baseQ_mean, mapQ, paired);
    ++cov_profile[pos];
}

inline void add_read_to_profile(seqan::String<seqan::ProfileChar<seqan::Dna5, double> > & profile,
                                std::vector<unsigned> & cov_profile,
                                vector<seqan::String<seqan::ProfileChar<seqan::Dna5, double> >> & ins_profiles,
                                std::vector<std::vector<unsigned>> & ins_cov_profile,
                                std::vector<unsigned> & no_ins_cov_profile,
                                seqan::Gaps<seqan::Dna5QString, seqan::ArrayGaps> const & gapsRead,
                                seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> const & gapsRef,
                                bool paired,
                                double mapQ,
                                SViperConfig const & config)
{
    SEQAN_ASSERT_EQ(seqan::length(gapsRead), seqan::length(gapsRef));

    unsigned begin = max(gapsBeginPos(gapsRef), gapsBeginPos(gapsRead));
    unsigned end = min(gapsEndPos(gapsRef), gapsEndPos(gapsRead));
    // append read bases to profile for every position in ref
    for (unsigned idx = begin; idx < end; ++idx)
    {
        unsigned pos_in_ref = seqan::toSourcePosition(gapsRef, idx);

        if (seqan::isGap(gapsRef, idx)) // this means the read features an insertion
        {
            SEQAN_ASSERT(!seqan::isGap(gapsRef, idx - 1)); // always handle full insertion

            unsigned ins_pos{0};

            while(seqan::isGap(gapsRef, idx))
            {
                SEQAN_ASSERT(!seqan::isGap(gapsRead, idx)); // always handle full insertion
                seqan::Dna5Q base = (seqan::source(gapsRead))[seqan::toSourcePosition(gapsRead, idx)];
                add_base_to_profile(ins_profiles[pos_in_ref], ins_cov_profile[pos_in_ref], ins_pos, base, mapQ, paired, config);
                ++idx; ++ins_pos;
            }
        }
        else
        {
            if (!seqan::isGap(gapsRef, idx + 1)) // no insertion happened between idx and idx + 1
            {
                add_gap_to_profile(ins_profiles[pos_in_ref], ins_cov_profile[pos_in_ref], 0, mapQ, paired, config);
                ++no_ins_cov_profile[pos_in_ref];
            }
        }

        if (!seqan::isGap(gapsRead, idx))
        {
            seqan::Dna5Q base = (seqan::source(gapsRead))[seqan::toSourcePosition(gapsRead, idx)];
            add_base_to_profile(profile, cov_profile, pos_in_ref, base, mapQ, paired, config);
        }
        else
        {
            add_gap_to_profile(profile, cov_profile, pos_in_ref, mapQ, paired, config);
        }
    }
}

void fill_profile(seqan::String<seqan::ProfileChar<seqan::Dna5, double> > & profile,
                  std::vector<unsigned> & cov_profile,
                  seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> & gapsRef,
                  std::vector<Mapping_object> const & mobs,
                  SViperConfig const & config)
{
    // the reference sequence (consensus sequence) will have a profile over
    // every position. This will account for deletions and substitutions.
    // Insertions are first stored seperately because otherwise the positions
    // in read-ref-alignments are not "synchronized" and hard to handle.
    std::vector<seqan::String<seqan::ProfileChar<seqan::Dna5, double> >> insertion_profiles;
    std::vector<std::vector<unsigned>> insertion_cov_profiles;
    insertion_profiles.resize(seqan::length(profile)); // there can be an insertion after every position
    insertion_cov_profiles.resize(seqan::length(profile));

    // the ins_cov_profile count the coverage of insertions for every little profile
    // but additionally, one needs to count the number of reads that span the position
    // and do NOT display an insertion. Those are counted with non_insertion_coverage_count
    // for each read and the insertion profile is updated accordingly.
    std::vector<unsigned> non_insertion_coverage_count;
    non_insertion_coverage_count.resize(seqan::length(profile));
    for (auto & cov_count : non_insertion_coverage_count)
        cov_count = 0;

    for (auto const & mob : mobs) // for every read (pair)
    {
        if (mob.proper_pair)
        {
            add_read_to_profile(profile, cov_profile, insertion_profiles,
                                insertion_cov_profiles, non_insertion_coverage_count,
                                mob.gapsRead, mob.gapsRef,
                                mob.proper_pair, mob.mapQRead, config);

            if (mob.hasMate())
                add_read_to_profile(profile, cov_profile, insertion_profiles,
                                    insertion_cov_profiles, non_insertion_coverage_count,
                                    mob.gapsMate, mob.gapsRefMate, mob.proper_pair,
                                    mob.mapQMate, config);
        }
    }

    // add insertions to profile
    unsigned pos{0}; // tracks the view position in gaps space of the reference
    for (unsigned pdx = 0; pdx < seqan::length(insertion_profiles); ++pdx)
    {
        auto & pro = insertion_profiles[pdx];
        auto & cov_pro = insertion_cov_profiles[pdx];
        if (seqan::length(pro) != 0) // if there is an insertion
        {
            for (auto & ins : pro) // update gap qualitities stored at the first position
                ins.count[5] = pro[0].count[5];
            for (unsigned cdx = 1; cdx < seqan::length(cov_pro); ++cdx) // update gap coverage. skip first because it is already updated
                cov_pro[cdx] += non_insertion_coverage_count[pdx];
            seqan::insert(profile, pos, pro);
            seqan::insert(cov_profile, pos, cov_pro);
            seqan::insertGaps(gapsRef, pos, seqan::length(pro));
            pos += seqan::length(pro) + 1; // +1 because we are also advancing one position anyway
        }
        else
        {
            ++pos;
        }
    }
}

void compute_mappQ_stats(SViperConfig & config,
                         std::vector<Mapping_object> const & mobs)
{
    if (seqan::length(mobs) == 0)
        return;

    double avg{0.0};
    double std{0.0};

    // compute average
    for (auto const & m : mobs)
        avg += m.mapQ();
    avg = avg / seqan::length(mobs);

    // compute standard deviation
    for (auto const & m : mobs)
        std += std::pow((m.mapQ() - avg), 2);
    std = std::sqrt(std / seqan::length(mobs));

    config.mappQ_mean = avg;
    config.mappQ_std = std;
}

void compute_baseQ_stats(SViperConfig & config,
                         seqan::StringSet<seqan::Dna5QString> const & reads1,
                         seqan::StringSet<seqan::Dna5QString> const & reads2)
{
    if (seqan::length(reads1) == 0)
        return;

    double avg{0.0};
    double std{0.0};

    // compute average
    for (auto const & read : reads1)
        for (auto const & c : read)
            avg += seqan::getQualityValue(c);

    for (auto const & read : reads2)
        for (auto const & c : read)
            avg += seqan::getQualityValue(c);

    // Note: this function assume all reads to be of equal length
    avg = avg / (seqan::length(reads1) * seqan::length(reads1[0]) + seqan::length(reads2) * seqan::length(reads2[0]));

    // compute standard deviation with avg
    for (auto const & read : reads1)
        for (auto const & c : read)
            std += std::pow((seqan::getQualityValue(c) - avg), 2);

    for (auto const & read : reads2)
        for (auto const & c : read)
            std += std::pow((seqan::getQualityValue(c) - avg), 2);

    // Note: this function assume all reads to be of equal length
    // and that reads1 and reads2 are of equal length
    std = std::sqrt(std / (2 * seqan::length(reads1) * seqan::length(reads1[0])));

    config.baseQ_mean = avg;
    config.baseQ_std = std;
}

inline seqan::Dna5String consensus_from_profile(seqan::String<seqan::ProfileChar<seqan::Dna5, double> > const & profile,
                                         seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> const & contigGaps,
                                         SViperConfig & config) // ebuffer at beginning and end
{
    SEQAN_ASSERT(seqan::length(profile) == seqan::length(contigGaps));

    seqan::Dna5String cns;
    unsigned begin = seqan::toViewPosition(contigGaps, config.ref_flank_length);
    unsigned end   = seqan::toViewPosition(contigGaps, seqan::length(seqan::source(contigGaps)) - config.ref_flank_length);

    // first append unpolished bases in the beginning (before begin)
    seqan::append(cns, seqan::prefix(seqan::source(contigGaps), config.ref_flank_length));

    // now fix consensus
    for (unsigned i = begin; i < end; ++i)
    {
        if (!config.fix_indels)       // if only substitutions are allowed
            if (seqan::isGap(contigGaps, i)) // do not insert bases
                continue;

        int idx = seqan::getMaxIndex(profile[i]);

        if (!config.score_passes_threshold(profile[i].count[idx], i))
        {
            if (!seqan::isGap(contigGaps, i))
                seqan::appendValue(cns, seqan::source(contigGaps)[seqan::toSourcePosition(contigGaps, i)]);
            continue;
        }

        if (idx < 5)  // is not gap TODO replace by seqan alphabet size or something like that
        {
            seqan::appendValue(cns, seqan::Dna5(idx));

            if (seqan::isGap(contigGaps, i))
                ++config.inserted_bases;
            else if (seqan::source(contigGaps)[seqan::toSourcePosition(contigGaps, i)] != seqan::Dna5(seqan::getMaxIndex(profile[i])))
                ++config.substituted_bases;
        }
        else if (!config.fix_indels) // if idx < 5 but indels should not be fixed
        {
            if (!seqan::isGap(contigGaps, i)) // gap to gap alignment can happen because of insertion profiles.
                seqan::appendValue(cns, seqan::source(contigGaps)[seqan::toSourcePosition(contigGaps, i)]);
        }
        else // idx < 5 and fix_indels is true
        {
            if (!seqan::isGap(contigGaps, i))
                ++config.deleted_bases;
        }
    }

    // at last, append config.ref_flank_length at the end
    seqan::append(cns, seqan::suffix(seqan::source(contigGaps),  seqan::length(seqan::source(contigGaps)) - config.ref_flank_length));

    if (!config.fix_indels)
        SEQAN_ASSERT_EQ(seqan::length(seqan::source(contigGaps)), seqan::length(cns));

    return cns;
}

seqan::Dna5String polish(seqan::StringSet<seqan::Dna5QString> const & reads1,
                  seqan::StringSet<seqan::Dna5QString> const & reads2,
                  seqan::Dna5String ref, // TODO:: reference?
                  SViperConfig & config)
{
    seqan::Gaps<seqan::Dna5String, seqan::ArrayGaps> contigGaps(ref);

    std::vector<Mapping_object> mobs = mapping(reads1, reads2, ref, config);

    compute_mappQ_stats(config, mobs);
    config.alpha = config.mappQ_mean/config.baseQ_mean;

    seqan::String<seqan::ProfileChar<seqan::Dna5, double> > profile;
    seqan::resize(profile, seqan::length(ref));
    config.cov_profile.clear();
    seqan::resize(config.cov_profile, seqan::length(ref));

    fill_profile(profile, config.cov_profile, contigGaps, mobs, config);

    ref = consensus_from_profile(profile, contigGaps, config);

    return ref;
}

seqan::Dna5String polish_to_perfection(seqan::StringSet<seqan::Dna5QString> const & reads1,
                                seqan::StringSet<seqan::Dna5QString> const & reads2,
                                seqan::Dna5String ref,
                                SViperConfig & config)
{
    seqan::Dna5String old_ref;

    // first, do only fix mismatches
    while (ref != old_ref && config.rounds < 20)
    {
        old_ref = ref; // store prior result
        ref = polish(reads1, reads2, ref, config);
        ++config.rounds;
    }

    // after all substitutions has been corrected, correct insertions/deletions
    // a few times with only proper pairs
    config.fix_indels = true;
    seqan::clear(old_ref);
    while (ref != old_ref && config.rounds < 30)
    {
        old_ref = ref; // store prior result
        ref = polish(reads1, reads2, ref, config);
        ++config.rounds;
    }

    return ref;
}
