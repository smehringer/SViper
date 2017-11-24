#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random> // random_shuffle

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#include <mateHunter.h>

template <typename file_type>
bool open_file_success(file_type & file, const char * name)
{
    if (!seqan::open(file, name))
    {
        std::cerr << "ERROR: Could not open file " << name << std::endl;
        return 0;
    }
    return 1;
}

struct bamRecordNameLess
{
    bool operator()(BamAlignmentRecord lhs, BamAlignmentRecord rhs) const
    {
        return lhs.qName < rhs.qName;
    }
};

void view_bam(std::vector<seqan::BamAlignmentRecord> & records,
              seqan::BamFileIn & bam_file,
              seqan::BamIndex<seqan::Bai> const & bam_index,
              seqan::CharString ref_name,
              int start,
              int end,
              bool long_reads)
{
    int rID = 0;
    bool hasAlignments = false;
    seqan::BamAlignmentRecord record{};

    // seqan jummpToRegion jumps to a bin containing the start position
    // since long reads often start way before this, we don?t want to miss those
    int start_reading = start;
    if (long_reads)
        start_reading = std::min(0, start - 10000);

    // get reference contig id
    if (!seqan::getIdByName(rID, seqan::contigNamesCache(seqan::context(bam_file)), ref_name))
    {
        std::cerr << "[ERROR]: view_bam - Reference sequence named "
                  << ref_name << " is not present in bam file." << std::endl;
        return;
    }

    if (!seqan::jumpToRegion(bam_file, hasAlignments, rID, start_reading, start_reading+1, bam_index))
    {
        std::cerr << "[ERROR]: view_bam - Could not jump to " << start_reading << std::endl;
        return;
    }

    if (!hasAlignments)
    {
        std::cerr << "[ERROR]: view_bam - No alignments for reference region ["
                  << start_reading << "-" << end << "]." << std::endl;
        return;
    }

    while (!seqan::atEnd(bam_file))
    {
        seqan::readRecord(record, bam_file);

        if (record.beginPos >= end)
            break;

        // check if read is in region since binning of BAM file might be off
        if (record.beginPos + getAlignmentLengthInRef(record) < start)
            continue;

        if (!hasFlagQCNoPass(record) && !hasFlagDuplicate(record) && // passes QC
            record.mapQ >= 15)                                       // only fairly unique hits
            records.push_back(record);
    }
}

Dna5String append_ref_flanks(Dna5String const & seq,
                             FaiIndex const & faiIndex,
                             unsigned fai_idx,
                             int ref_length,
                             int start,
                             int end,
                             int length)
{
    SEQAN_ASSERT(start <= end);

    Dna5String leftFlank;
    Dna5String rightFlank;

    readRegion(leftFlank, faiIndex, fai_idx, max(0, start - length), max(0, start));
    readRegion(rightFlank, faiIndex, fai_idx, min(end, ref_length), min(end + length, ref_length));

    Dna5String new_seq = leftFlank;
    append(new_seq, seq);
    append(new_seq, rightFlank);

    return new_seq;
}

void records_to_read_pairs(StringSet<Dna5QString> & reads1,
                           StringSet<Dna5QString> & reads2,
                           vector<BamAlignmentRecord> & records,
                           seqan::BamFileIn & bam_file,
                           seqan::BamIndex<seqan::Bai> const & bam_index)
{
    // sort records by name
    std::sort(records.begin(), records.end(), bamRecordNameLess());
    Dna5QString dummy_seq{}; // read sequence for dummy record is an empty string

    CharString id1; // stores last seen id

    for (unsigned i = 0; i < length(records); ++i)
    {
        Dna5QString seq = records[i].seq;
        assignQualities(seq, records[i].qual);

        if (hasFlagRC(records[i]))
            reverseComplement(seq);

        if (length(reads1) == length(reads2)) // new read instead of potential mate
        {
            id1 = records[i].qName;
            appendValue(reads1, seq);
        }
        else if (records[i].qName == id1) // |ids1| = |ids2| + 1
        {
            appendValue(reads2, seq);
        }
        else // |ids1| = |ids2| + 1 but records[i] is not the former reads mate
        {
            // This case means that there is no mate for records[i -1] in records.
            // This only happens for two type of paired reads:
            // 1) The read is in a proper pair, but the mate is outside of our
            //    region of interest.
            // 2) The read is an abbarrant mapping pair with his mate unmapped
            //    or mapped to another chromosome.
            // We ONLY want to hunt for mates of pairs of the 2. case, since those
            // can be interesting, e.g. for novel insertions.

            if (!hasFlagAllProper(records[i - 1]))
            {
                // hunt mate for read in records[i - 1]
                BamAlignmentRecord mate = mateHunt(records[i - 1], bam_file, bam_index);

                if (mate.qName == records[i - 1].qName) // mateHunt was successfull
                {
                    Dna5QString mseq = mate.seq;
                    assignQualities(mseq, mate.qual);

                    if (hasFlagRC(mate))
                        reverseComplement(mseq);

                    appendValue(reads2, mseq);
                }
                else
                {
                    // append dummy mate
                    appendValue(reads2, dummy_seq);
                }
            }
            else
            {
                // append dummy mate
                appendValue(reads2, dummy_seq);
            }

            // append rec to reads1. his mate is still be out there
            id1 = records[i].qName;
            appendValue(reads1, seq);
        }
    }

    // mazbe the last record had no mate. Correct this
    if (length(reads1) != length(reads2))
    {
        // hunt mate for read in records[i - 1]
        BamAlignmentRecord mate = mateHunt(records[length(records) - 1], bam_file, bam_index);

        if (mate.qName == records[length(records) - 1].qName) // mateHunt was successfull
        {
            Dna5QString mseq = mate.seq;
            assignQualities(mseq, mate.qual);
            appendValue(reads2, mseq);
        }
        else
        {
            // append dummy mate
            appendValue(reads2, dummy_seq);
        }
    }
}

template <typename container_type>
void subsample(container_type & c1, container_type & c2, unsigned N)
{
    SEQAN_ASSERT_EQ(length(c1), length(c2));

    container_type c1_new;
    container_type c2_new;
    resize(c1_new, N);
    resize(c2_new, N);

    std::vector<unsigned int> indices(length(c1));
    std::iota(indices.begin(), indices.end(), 0);
    std::random_shuffle(indices.begin(), indices.end());

    // now fill c1/2_new with the first N indices randomly chosen
    for (unsigned idx = 0; idx < N; ++idx)
    {
        c1_new[idx] = c1[indices[idx]];
        c2_new[idx] = c2[indices[idx]];
    }

    c1 = c1_new;
    c2 = c2_new;
}
