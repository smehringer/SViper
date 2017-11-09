#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

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


BamAlignmentRecord map_single_read(Dna5String read,       // copy for gaps object
                                   CharString const & id,
                                   CharString qual,       // copy for reversing
                                   Dna5String ref)        // copy for gaps object ?
{
    typedef String<CigarElement<char, unsigned>> TCigarString;
    typedef Gaps<Dna5String, ArrayGaps> TGaps;

    int const MATCH = 3;
    int const MISMATCH = -2;
    int const GAP_OPEN = -3;
    int const GAP_EXT = -1;

    Dna5String readRC = read; // copy for reverse complementing in place
    Dna5String refRC  = ref;   // copy for assigning to gaps as ref

    reverseComplement(readRC);

    TGaps gapsRef(ref);
    TGaps gapsRead(read);
    TGaps gapsRefRC(refRC);
    TGaps gapsReadRC(readRC);

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
        record.qual = qual;
        record.mapQ = mapQ;
    }
    else // reverse strand
    {
        getCigarString(cigar, gapsRefRC, gapsReadRC);

        CharString qual_cp = qual;
        reverse(qual_cp);

        double mapP = (((double)scoreRC / (double)(length(readRC)*MATCH)));
        unsigned short mapQ = mapP * 100;

        record.beginPos = countLeadingGaps(gapsReadRC);
        record.flag = 16;
        record.seq = readRC;
        record.qual = qual_cp;
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

inline void paired_mapping(StringSet<Dna5String> const & reads,
                           StringSet<String<char>> const & ids,
                           StringSet<String<char>> const & quals,
                           StringSet<Dna5String> const & reads2,
                           StringSet<String<char>> const & ids2,
                           StringSet<String<char>> const & quals2,
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

    for (unsigned i = 0; i < length(reads); ++i)
    {
        if (ids[i] != ids2[i])
            throw std::ios_base::failure("paired end data are not correct");

        BamAlignmentRecord record1 = map_single_read(reads[i], ids[i], quals[i], ref);
        BamAlignmentRecord record2 = map_single_read(reads2[i], ids2[i], quals2[i], ref);

        record1.pNext = record2.beginPos;
        record2.pNext = record2.beginPos;
        record1.rNextId = 0; // ref is the same
        record2.rNextId = 0; // ref is the same

        // right now, the records only have a flag for reverse set
        if (hasFlagRC(record2))
            record1.flag |= 0x20; // mate is reversed
        if (hasFlagRC(record1))
            record2.flag |= 0x20;
        record1.flag |= 0x1; // is paired
        record2.flag |= 0x1;
        record1.flag |= 0x40; // first in pair
        record2.flag |= 0x80; // second in paired

        if ((!hasFlagRC(record1) && hasFlagRC(record2) && (record1.beginPos + record1.tLen) < record2.beginPos) ||
            (hasFlagRC(record1) && !hasFlagRC(record2) && (record2.beginPos + record2.tLen) < record1.beginPos)) // ignore insertsize for now
        {
            record1.tLen = abs(record2.beginPos - record1.beginPos) + length(reads[i]); // approximate fragment size
            record1.flag |= 0x2; // mapped in proper pair
            record2.flag |= 0x2;
        }

        write(pseudo_bamfile, record1, bamIOContext, Sam());
        write(pseudo_bamfile, record2, bamIOContext, Sam());
        //writeRecord(bamFileOut, record);
    }
}

struct ConsensusConfig
{
    bool fix_indels{false};
    bool only_proper_pairs{true};
    unsigned buffer{150};
};

template <typename TStore>
inline void fill_profiles(String<ProfileChar<Dna5, double> > & quality_profile,
                                String<ProfileChar<Dna5, double> > & coverage_profile,
                                TStore const & store,
                                ConsensusConfig const & config)
{
    typedef typename Value<typename TStore::TAlignedReadStore>::Type   TAlignedRead;
    typedef typename Value<typename TStore::TAlignQualityStore>::Type  TAlignQuality;
    typedef Gaps<typename TStore::TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >  TReadGaps;
    typedef typename TAlignedRead::TId TAlignedReadId;

    SEQAN_ASSERT_EQ(length(quality_profile), length(coverage_profile));

    typename TStore::TReadSeq readSeq;
    typename TAlignedRead::TPos begin;

    // construct readId to alignId map ... which I can't find in the fragmentstore but need it for mate retrieval
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
        TAlignedRead  const & am = store.alignedReadStore[readId2alignId[ar.pairMatchId]]; // aligned mate
        TAlignQuality const & aq = store.alignQualityStore[i];             // aligned quality

        SEQAN_ASSERT_EQ(store.readNameStore[ar.readId], store.readNameStore[am.readId]);

        if (ar.contigId != 0) // not aligned to first contig, which is the consensus sequence
            continue;

        if (config.only_proper_pairs)
            if (!(ar.endPos < ar.beginPos /* ar (-)*/ && am.endPos > am.beginPos /*am (+)*/ && am.endPos < ar.beginPos /*am maps before ar*/) &&
                !(ar.endPos > ar.beginPos /* ar (+)*/ && am.endPos < am.beginPos /*am (-)*/ && ar.endPos < am.beginPos /*ar maps before am*/)) // TODO incoorparate correct insert size
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

        // go through the read and fill quality_profile of contig
        for (unsigned row = 0; (row < length(readGaps) && (begin + row) < length(quality_profile)); ++row)
        {
            unsigned ord_idx{5};
            double qual{(double)aq.score/100};

            if (!isGap(readGaps, row))
            {
                unsigned sourcePos = toSourcePosition(readGaps, row);
                auto w = readSeq[sourcePos];
                ord_idx = ordValue(w);
                //qual = (qual + (getQualityValue(readSeq[sourcePos]) / 100) ) / 2; // TODO:: incoorperate base qual
            }

            quality_profile[begin + row].count[ord_idx] += qual;
            coverage_profile[begin + row].count[ord_idx] += 1;

            // one liner: (bad for debugging)
            // profile[begin + row].count[ordValue(readGaps[row])] += ((getQualityValue(readSeq[0]) + aq.score) / 200);
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
        if (!config.fix_indels) // if only substitutions are allowed
            if (isGap(contigGaps, i)) // do not insert bases
                continue;

        int idx = getMaxIndex(profile[i]);

        if (idx < 5)  // is not gap TODO replace by seqan alphabet size or something like that
        {
            appendValue(cns, Dna5(getMaxIndex(profile[i])));
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

Dna5String polish(StringSet<Dna5String> const & reads1,
                  StringSet<String<char>> const & ids1,
                  StringSet<String<char>> const & quals1,
                  StringSet<Dna5String> const & reads2,
                  StringSet<String<char>> const & ids2,
                  StringSet<String<char>> const & quals2,
                  Dna5String ref,
                  CharString id,
                  ConsensusConfig const & config)
{
    typedef FragmentStore<>                       TFragmentStore;
    typedef typename TFragmentStore::TContigStore TContigStore;
    typedef typename Value<TContigStore>::Type    TContig;
    typedef Gaps<typename TContig::TContigSeq, AnchorGaps<typename TContig::TGapAnchors>> TContigGaps;

    cout << "### Mapping..." << endl;
    paired_mapping(reads1, ids1, quals1, reads2, ids2, quals2, ref, id); // creates pseudo.sam for FragmentStore

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

    cout << "### Create Fragment Store..." << endl;
    BamFileIn bamfile("pseudo.sam");
    readRecords(store, bamfile);

    cout << "### Fill profiles..." << endl;
    TContigGaps contigGaps(store.contigStore[0].seq, store.contigStore[0].gaps);
    String<ProfileChar<Dna5, double> > quality_profile;
    String<ProfileChar<Dna5, double> > coverage_profile;
    resize(quality_profile, length(contigGaps));
    resize(coverage_profile, length(contigGaps));
    fill_profiles(quality_profile, coverage_profile, store, config);

    cout << "## stats" << endl;
    double mean_coverage{0};
    for (auto const & p : coverage_profile)
    {
        mean_coverage += totalCount(p);
    }
    mean_coverage = mean_coverage/length(coverage_profile);
    std::cout << mean_coverage << endl;


    cout << "### Build new consensus... length before: " << length(ref) << endl;
    // Do not correct bases to the left and right of the consensus sequences otherwise
    // we extend the alignment with low coverage.
    // Also skip 50 bp of the beginning because reads will map poorly there
    ref = consensus_from_profile(quality_profile,
                                 contigGaps,
                                 config);

    // cout << "### Layout alignment for debug and print profile..." << endl;
    // AlignedReadLayout layout2;
    // layoutAlignment(layout2, store);
    // printAlignment(std::cout, layout2, store, 0, 0, 1550, 0, 46);

    return ref;
}
