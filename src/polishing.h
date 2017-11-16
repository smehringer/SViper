#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <utility> // pair

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/store.h>
#include <seqan/align.h>

#include <helper_functions.h>

using namespace std;
using namespace seqan;

CharString get_quality_string(Dna5QString const & read)
{
    CharString qual;
    resize(qual, length(read));
    for (unsigned i = 0; i < length(read); ++i)
    {
        char c;
        convertQuality(c, getQualityValue(read[i]));
        qual[i] = c;
    }
    return qual;
}

BamAlignmentRecord map_single_read(Dna5QString read,       // copy for gaps object
                                   CharString const & id,
                                   Dna5String ref)        // copy for gaps object ?
{
    typedef String<CigarElement<char, unsigned>> TCigarString;
    typedef Gaps<Dna5QString, ArrayGaps> TGapsRead;
    typedef Gaps<Dna5String, ArrayGaps> TGapsRef;

    int const MATCH = 3;
    int const MISMATCH = -2;
    int const GAP_OPEN = -3;
    int const GAP_EXT = -1;

    if (length(read) == 0)
    {
        // return unmapped read
        BamAlignmentRecord record;
        record.qName = id;
        record.rID = BamAlignmentRecord::INVALID_REFID;
        record.beginPos = BamAlignmentRecord::INVALID_POS;
        record.tLen = BamAlignmentRecord::INVALID_LEN;
        return record;
    }

    Dna5QString readRC = read; // copy for reverse complementing in place
    Dna5String refRC  = ref;   // copy for assigning to gaps as ref

    reverseComplement(readRC);

    TGapsRef gapsRef(ref);
    TGapsRef gapsRefRC(refRC);
    TGapsRead gapsRead(read);
    TGapsRead gapsReadRC(readRC);

    int score = globalAlignment(gapsRef, gapsRead,
                                Score<int, Simple>(MATCH, MISMATCH, GAP_EXT, GAP_OPEN),
                                AlignConfig<true, false, false, true>(),
                                AffineGaps());

    int scoreRC = globalAlignment(gapsRefRC, gapsReadRC,
                                  Score<int, Simple>(MATCH, MISMATCH, GAP_EXT, GAP_OPEN),
                                  AlignConfig<true, false, false, true>(),
                                  AffineGaps());

    BamAlignmentRecord record;
    record.rID = 0;
    record.qName = id;

    TCigarString cigar;

    if (score > scoreRC) // forward strand
    {
        getCigarString(cigar, gapsRef, gapsRead);

        double mapP = (((double)score / (double)(length(read)*MATCH)));
        unsigned short mapQ = mapP * 100;

        record.beginPos = countLeadingGaps(gapsRead);
        record.flag = 0;
        record.seq = read;
        record.qual = get_quality_string(read);
        record.mapQ = mapQ;
    }
    else // reverse strand
    {
        getCigarString(cigar, gapsRefRC, gapsReadRC);

        double mapP = (((double)scoreRC / (double)(length(readRC)*MATCH)));
        unsigned short mapQ = mapP * 100;

        record.beginPos = countLeadingGaps(gapsReadRC);
        record.flag = 16;
        record.seq = readRC;
        record.qual = get_quality_string(readRC);
        record.mapQ = mapQ;
    }

    // fix cigar
    if ((cigar[0]).operation == 'N' || (cigar[0]).operation == 'D')
        erase(cigar, 0);
    else if ((cigar[0]).operation == 'I')
        (cigar[0]).operation = 'S';
    if ((cigar[length(cigar) - 1]).operation == 'N' || (cigar[length(cigar) - 1]).operation == 'D')
        erase(cigar, length(cigar) - 1);
    else if ((cigar[0]).operation == 'I')
        (cigar[0]).operation = 'S';

    record.cigar = cigar;

    return record;
}

inline void paired_mapping(StringSet<Dna5QString> const & reads,
                           StringSet<String<char>> const & ids,
                           StringSet<Dna5QString> const & reads2,
                           StringSet<String<char>> const & ids2,
                           Dna5String ref,
                           CharString id) // ref name
{
    // prepare bam context
    StringSet<CharString> contigNameStore;
    appendValue(contigNameStore, id);
    NameStoreCache<StringSet<CharString> > contigNameStoreCache(contigNameStore);
    BamIOContext<StringSet<CharString> > bamIOContext(contigNameStore, contigNameStoreCache);
    appendValue(contigLengths(bamIOContext), length(ref));

    // BamFileOut bamFileOut(refer);
    ofstream pseudo_bamfile("pseudo.sam");

    #pragma omp parallel for
    for (unsigned i = 0; i < length(reads); ++i)
    {
        if (ids[i] != ids2[i])
            throw std::ios_base::failure("[ERROR] mapping - paired end data are not correct");

        BamAlignmentRecord record1 = map_single_read(reads[i], ids[i], ref);

        if (length(reads2[i]) == 0) // if mate is a dummy mate
        {
            #pragma omp critical
            write(pseudo_bamfile, record1, bamIOContext, Sam());
            continue;
        }

        BamAlignmentRecord record2 = map_single_read(reads2[i], ids2[i], ref);

        record1.pNext = record2.beginPos;
        record2.pNext = record1.beginPos;
        record1.rNextId = record2.rID;
        record2.rNextId = record1.rID;

        // right now, the records only have a flag for reverse set
        if (hasFlagRC(record2))
            record1.flag |= 0x20; // mate is reversed
        if (hasFlagRC(record1))
            record2.flag |= 0x20;
        record1.flag |= 0x1; // is paired
        record2.flag |= 0x1;
        record1.flag |= 0x40; // first in pair
        record2.flag |= 0x80; // second in paired

        if (!hasFlagRC(record1) && hasFlagRC(record2) && (record1.beginPos + record1.tLen) < record2.beginPos)
        { // ---r1--->     <---r2---
            record1.tLen = record2.beginPos - record1.beginPos + length(record1.seq); // approximate fragment size
            record1.flag |= 0x2; // mapped in proper pair
            record2.flag |= 0x2;
        }
        else if (hasFlagRC(record1) && !hasFlagRC(record2) && (record2.beginPos + record2.tLen) < record1.beginPos)
        { // ---r2--->     <---r1---
            record1.tLen = record1.beginPos - record2.beginPos + length(record2.seq); // approximate fragment size
            record1.flag |= 0x2; // mapped in proper pair
            record2.flag |= 0x2;
        }

        #pragma omp critical
        write(pseudo_bamfile, record1, bamIOContext, Sam());
        #pragma omp critical
        write(pseudo_bamfile, record2, bamIOContext, Sam());
        //writeRecord(bamFileOut, record);
    }
}

struct ConsensusConfig
{
    bool verbose{false};
    bool fix_indels{false};
    bool only_proper_pairs{true};
    unsigned buffer{150};

    // qualitiy statistics
    double baseQ_mean{0}; // will be overriden once when reading in reads
    double baseQ_std{1};  // will be overriden once when reading in reads
    double mappQ_mean{0}; // will be overriden every round after mapping
    double mappQ_std{1};  // will be overriden every round after mapping
    double alpha{1}; // scaling factor such that mapping and base quality are comparable
    double mean_coverage{0};
    double min_coverage{4};

    /* [used in fill_profiles
     * Returns the value to be added to the profile. A profile consists of counts
     * for the letters A,C,T,G,- weighted by their quality. The profile is later
     * on used to determine the new reference sequence.*/
    double add_to_profile(double const & baseQ,
                        double const & mappQ,
                        bool is_in_propair_pair) const
    {
        if (is_in_propair_pair)
            return (baseQ + (mappQ / alpha)); // return double the score
        return ((baseQ + (mappQ / alpha)) / 2);
    }

    /* [used in build_consensus_from_profile]
     * Determines whether short read information is to be trusted
     * TRUE: if base with given score will be substituted in consensus sequence
     * FALSE: else*/
    bool score_passes_threshold(double const & score) const
    {
        /* return true if score is at least as good as half the mean coverage
         * of bases with half as good quality as the mean quality*/
        return score >= (max(min_coverage, mean_coverage/2) * (baseQ_mean + mappQ_mean/alpha) / 4);
    }
};

template <typename TStore>
void compute_mappQ_stats(ConsensusConfig & config,
                         TStore const & store)
{
    if (length(store.alignQualityStore) == 0)
        return;

    double avg{0.0};
    double std{0.0};

    // compute average
    for (auto const & qual : store.alignQualityStore)
        avg += qual.score;
    avg = avg / length(store.alignQualityStore);

    // compute standard deviation
    for (auto const & qual : store.alignQualityStore)
        std += std::pow((qual.score - avg), 2);
    std = std::sqrt(std / length(store.alignQualityStore));

    config.mappQ_mean = avg;
    config.mappQ_std = std;
}

void compute_baseQ_stats(ConsensusConfig & config,
                         String<Dna5QString> const & reads1,
                         String<Dna5QString> const & reads2)
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

template <typename TStore>
inline void fill_profiles(String<ProfileChar<Dna5, double> > & quali_profile,
                          String<ProfileChar<Dna5, double> > & cover_profile,
                          TStore const & store,
                          ConsensusConfig const & config)
{
    typedef typename Value<typename TStore::TAlignedReadStore>::Type   TAlignedRead;
    typedef typename Value<typename TStore::TAlignQualityStore>::Type  TAlignQuality;
    typedef Gaps<typename TStore::TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >  TReadGaps;
    typedef typename TAlignedRead::TId TAlignedReadId;

    typename TStore::TReadSeq readSeq;
    typename TAlignedRead::TPos begin;

    SEQAN_ASSERT_EQ(length(quali_profile), length(cover_profile));

    // construct readId to alignId map ...
    // which I can't find in the fragmentstore but need it for mate retrieval
    String<TAlignedReadId> readId2alignId;
    resize(readId2alignId, length(store.alignedReadStore));
    for (TAlignedReadId i = 0; i < length(store.alignedReadStore); ++i)
    {
        TAlignedRead const & ar = store.alignedReadStore[i];
        readId2alignId[ar.readId] = i;
    }

    for (unsigned i = 0; i < length(store.alignedReadStore); ++i)
    {
        TAlignedRead  const & ar = store.alignedReadStore[i];              // aligned read
        TAlignQuality const & aq = store.alignQualityStore[i];             // aligned quality

        if (ar.contigId != 0) // not aligned to first contig, which is the consensus sequence
            continue;

        bool is_proper_pair{false};
        // check if read is aligned in proper pair
        if (ar.pairMatchId != TAlignedRead::INVALID_ID)
        {
            TAlignedRead const & am = store.alignedReadStore[readId2alignId[ar.pairMatchId]]; // aligned mate

            SEQAN_ASSERT_EQ(store.readNameStore[ar.readId], store.readNameStore[am.readId]);

            if ((ar.endPos < ar.beginPos /*(-)*/ && am.endPos > am.beginPos /*(+)*/ && am.endPos < ar.beginPos /*am maps before ar*/) ||
                (ar.endPos > ar.beginPos /*(+)*/ && am.endPos < am.beginPos /*(-)*/ && ar.endPos < am.beginPos /*ar maps before am*/)) // TODO incoorparate correct insert size
                is_proper_pair = true;
        }

        if (config.only_proper_pairs && !is_proper_pair)
            continue;

        readSeq = store.readSeqStore[ar.readId];
        if (ar.endPos < ar.beginPos)
        {
            reverseComplement(readSeq);
            begin = ar.endPos;
        }
        else
        {
            begin = ar.beginPos;
        }

        TReadGaps readGaps(readSeq, ar.gaps);

        // go through the read and fill quali_profile of contig
        for (unsigned row = 0; (row < length(readGaps) && (begin + row) < length(quali_profile)); ++row)
        {
            unsigned ord_idx{5};
            double baseQ = aq.score; // initialize base quality with maping quality

            if (!isGap(readGaps, row))
            {
                unsigned sourcePos = toSourcePosition(readGaps, row);
                auto w = readSeq[sourcePos];
                ord_idx = ordValue(w);
                baseQ = getQualityValue(readSeq[sourcePos]);
            }

            quali_profile[begin + row].count[ord_idx] += config.add_to_profile(baseQ, aq.score, is_proper_pair);
            cover_profile[begin + row].count[ord_idx] += 1;
        }
    }
}

template <typename TContigGaps>
inline Dna5String consensus_from_profile(String<ProfileChar<Dna5, double> > const & profile,
                                         TContigGaps const & contigGaps,
                                         ConsensusConfig const & config) // ebuffer at beginning and end
{
    SEQAN_ASSERT(length(profile) == length(contigGaps));

    Dna5String cns;
    unsigned begin = toViewPosition(contigGaps, config.buffer);
    unsigned end   = toViewPosition(contigGaps, length(source(contigGaps)) - config.buffer);

    unsigned confirmed{0};
    unsigned substitutions{0};
    unsigned insertions{0};
    unsigned deletions{0};

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
                ++insertions;
            else if (source(contigGaps)[toSourcePosition(contigGaps, i)] == Dna5(getMaxIndex(profile[i])))
                ++confirmed;
            else
                ++substitutions;
        }
        else if (!config.fix_indels) // if idx < 5 but indels should not be fixed
        {
            appendValue(cns, source(contigGaps)[toSourcePosition(contigGaps, i)]);
        }
        else // idx < 5 and fix_indels is true
        {
            if (!isGap(contigGaps, i))
                ++deletions;
        }
    }

    // at last, append config.buffer at the end
    append(cns, suffix(source(contigGaps),  length(source(contigGaps)) - config.buffer));

    if (!config.fix_indels)
        SEQAN_ASSERT_EQ(length(source(contigGaps)), length(cns));

    if (config.verbose)
        cout << "CHANGES ARE (one base): "
             << confirmed << " confirmed, "
             << substitutions << " substitutions, "
             << insertions << " insertions and "
             << deletions << " deletions." << endl;

    return cns;
}

template <typename TStore>
inline void print_profile(String<ProfileChar<Dna5, double> > const & profile, TStore const & store)
{
    typedef typename Value<typename TStore::TContigStore>::Type TContig;
    typedef Gaps<typename TContig::TContigSeq, AnchorGaps<typename TContig::TGapAnchors>> TContigGaps;

    TContigGaps contigGaps(store.contigStore[0].seq, store.contigStore[0].gaps);
    SEQAN_ASSERT(length(contigGaps) == length(profile));

    unsigned confirmed{0};
    unsigned substitutions{0};
    unsigned insertions{0};
    unsigned deletions{0};

    for (unsigned j = 0; j < length(profile[0].count); j++)
    {
        for (unsigned i = 0; i < length(profile); ++i)
        {
            cout << profile[i].count[j] << " ";
        }
        cout << endl;
    }

    stringstream prof;
    stringstream change;
    for (unsigned i = 0; i < length(profile); ++i)
    {
        if (isGap(contigGaps, i))
        {
            if (getMaxIndex(profile[i]) < 5) // not gap
            {
                prof << Dna5(getMaxIndex(profile[i]));
                change << "X";
                ++insertions;
            }
            else
            {
                prof << "-";
                change << "|";
            }
        }
        else
        {
            if (getMaxIndex(profile[i]) < 5) // not gap
            {
                prof << Dna5(getMaxIndex(profile[i]));

                if (Dna5(getMaxIndex(profile[i])) == (source(contigGaps))[toSourcePosition(contigGaps, i)])
                {
                    ++confirmed;
                    change << "|";
                }
                else
                {
                    ++substitutions;
                    change << "X";
                }
            }
            else
            {
                prof << "-";
                change << "X";
                ++deletions;
            }
        }
    }
    cout << endl;
    cout << prof.str() << endl;
    cout << change.str() << endl;
    cout << contigGaps << endl;
    cout << "confirmed: " << confirmed << endl;
    cout << "substitutions: " << substitutions << endl;
    cout << "insertions: " << insertions << endl;
    cout << "deletions: " << deletions << endl;

    cout << endl;
}

Dna5String polish(StringSet<Dna5QString> const & reads1,
                  StringSet<String<char>> const & ids1,
                  StringSet<Dna5QString> const & reads2,
                  StringSet<String<char>> const & ids2,
                  Dna5String ref,
                  CharString id,
                  ConsensusConfig & config)
{
    typedef FragmentStore<>                       TFragmentStore;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type    TContig;
    typedef Gaps<typename TContig::TContigSeq, AnchorGaps<typename TContig::TGapAnchors>> TContigGaps;

    if (config.verbose)
        cout << "### Mapping..." << endl;
    paired_mapping(reads1, ids1, reads2, ids2, ref, id); // creates pseudo.sam for FragmentStore

    TFragmentStore store;
    // loadContigs(store, "consensus.fa");
    appendValue(store.contigNameStore, id);
    resize(store.contigStore, 1);
    TContig & contig = back(store.contigStore);
    contig.usage = 0;
    contig.fileId = 0;
    contig.fileBeginPos = 0;
    contig.seq = ref;
    SEQAN_ASSERT(store.contigStore[0].seq == ref);
    BamFileIn bamfile("pseudo.sam");
    readRecords(store, bamfile);

    if (config.verbose)
        cout << "### Fill profiles..." << endl;
    TContigGaps contigGaps(store.contigStore[0].seq, store.contigStore[0].gaps);
    String<ProfileChar<Dna5, double> > quali_profile;
    String<ProfileChar<Dna5, double> > cover_profile;
    resize(quali_profile, length(contigGaps));
    resize(cover_profile, length(contigGaps));
    compute_mappQ_stats(config, store);
    config.alpha = config.mappQ_mean/config.baseQ_mean;
    fill_profiles(quali_profile, cover_profile, store, config);

    double mean_coverage{0};
    for (auto const & p : cover_profile)
    {
        mean_coverage += totalCount(p);
    }
    config.mean_coverage = mean_coverage/length(cover_profile);

    if (config.verbose)
        cout << "### MappQ avg: " << config.mappQ_mean << " and std: " << config.mappQ_std
             << ". Cov: " << config.mean_coverage << endl;

    if (config.verbose)
        cout << "### Build new consensus... length before: " << length(ref) << endl;
    // Do not correct bases to the left and right of the consensus sequences otherwise
    // we extend the alignment with low coverage.
    // Also skip 50 bp of the beginning because reads will map poorly there
    ref = consensus_from_profile(quali_profile,
                                 contigGaps,
                                 config);

    // UNCOMMENT THIS FOR PRETTY DEBUGGING OUTPUT:
    // print_profile(quali_profile, store);
    // cout << "### Layout alignment for debug and print profile..." << endl;
    // AlignedReadLayout layout2;
    // layoutAlignment(layout2, store);
    // printAlignment(std::cout, layout2, store, 0, 0, 1550, 0, 46);

    return ref;
}

Dna5String polish_to_perfection(StringSet<Dna5QString> const & reads1,
                                StringSet<String<char>> const & ids1,
                                StringSet<Dna5QString> const & reads2,
                                StringSet<String<char>> const & ids2,
                                Dna5String ref,
                                CharString id,
                                ConsensusConfig & config)
{
    Dna5String old_ref;

    unsigned round{1};

    if (config.verbose)
        std::cout << "Number of short reads: " << (length(reads1) * 2)
                  << ". Base Quality mean: " << config.baseQ_mean
                  << " and std: " << config.baseQ_std << std::endl;

    while (ref != old_ref && round < 20)
    {
        old_ref = ref; // store prior result
        ref = polish(reads1, ids1, reads2, ids2, ref, id, config);
        ++round;
        // break;
    }

    // after all substitutions has been corrected, correct insertions/deletions a few times with only proper pairs
    config.fix_indels = true;
    for (unsigned i= 0; i < 10; ++i)
    {
        ref = polish(reads1, ids1, reads2, ids2, ref, id, config);
        ++round;
        // break;
    }

    config.only_proper_pairs = false;
    old_ref = ""; // reset
    while (ref != old_ref && round < 50)
    {
        old_ref = ref; // store prior result
        ref = polish(reads1, ids1, reads2, ids2, ref, id, config);
        ++round;
        // break;
    }

    return ref;
}
