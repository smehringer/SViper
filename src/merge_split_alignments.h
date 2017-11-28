#pragma once

#include <iostream>
#include <seqan/bam_io.h>

#include <helper_functions.h>
#include <basics.h>

using namespace seqan;

struct bamRecordQualityLess
{
    bool operator()(BamAlignmentRecord lhs, BamAlignmentRecord rhs) const
    {
        if (!hasFlagSupplementary(lhs) && !hasFlagSecondary(lhs)) // is primary
            return true;
        if (!hasFlagSupplementary(rhs) && !hasFlagSecondary(rhs)) // is primary
            return false;
        return lhs.qual >= rhs.qual;
    }
};

// This function computes the end position of the mapping
int compute_map_end_pos(unsigned map_begin_pos,
                        String<CigarElement<char, unsigned>> & cigar)
{
    int len{0}; // tracks read length so far
    for (auto ce : cigar) //
        if (ce.operation & 4) // D or M, not I nor S
            len += ce.count;
    return map_begin_pos + len;
}


unsigned compute_fragment_length(String<CigarElement<char, unsigned>> & cigar)
{
    unsigned len{0}; // tracks read length so far
    for (auto ce : cigar) //
        if (ce.operation != 'D')
            len += ce.count;
    return len;
}

unsigned get_read_begin_and_alter_cigar(BamAlignmentRecord & record)
{
    unsigned read_begin_pos{0};
    if ((record.cigar[0]).operation == 'H') // hard clipping
    {
        read_begin_pos += (record.cigar[0]).count;
        erase(record.cigar, 0);
    }
    if ((record.cigar[0]).operation == 'S') // soft clipping
    {
        read_begin_pos += (record.cigar[0]).count;
        erase(record.cigar, 0);
    }
    return read_begin_pos;
}

unsigned get_read_end_and_alter_cigar(BamAlignmentRecord & record)
{
    unsigned read_end_pos{length(record.seq) - 1};
    if ((record.cigar[length(record.cigar) - 1]).operation == 'H') // hard clipping
    {
        read_end_pos -= (record.cigar[length(record.cigar) - 1]).count;
        erase(record.cigar, length(record.cigar) - 1);
    }
    if ((record.cigar[length(record.cigar) - 1]).operation == 'S') // soft clipping
    {
        read_end_pos -= (record.cigar[length(record.cigar) - 1]).count;
        erase(record.cigar, length(record.cigar) - 1);
    }
    return read_end_pos;
}

void truncate_cigar_right(BamAlignmentRecord & record, int until)
{
    unsigned cigar_pos{0};
    int read_pos{0};
    int ref_pos{record.beginPos};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&until] (int /*ref*/, int read) {return read >= (until - 1);});

    if (cigar_pos == 0) // no advancement happened. only case: until=1
    {
        erase(record.cigar, cigar_pos + 1, length(record.cigar)); // erase the rest
        record.cigar[0].count = until;
        return;
    }

    erase(record.cigar, cigar_pos, length(record.cigar)); // erase the rest

    if (read_pos == (until -1))
    {
        appendValue(record.cigar, CigarElement<char, unsigned>('I', 1));
    }
    if (read_pos > (until -1))
    {
        (record.cigar[length(record.cigar) - 1]).count -= (read_pos - until);
    }

    if ((record.cigar[length(record.cigar) - 1]).operation != 'M' &&
        (record.cigar[length(record.cigar) - 1]).operation != 'I') // meaning == H,S or D
    {
        if ((record.cigar[length(record.cigar) - 1]).operation == 'D')
        {   // continuing deletion in record must be added to the laarge DEL event
            erase(record.cigar, length(record.cigar) - 1);
        }
        else // meaning == H or S
        {
            (record.cigar[length(record.cigar) - 1]).operation = 'I'; // add those bases
        }
    }
}

void transform_cigar_right_into_insertion(BamAlignmentRecord & record, int until)
{
    unsigned cigar_pos{0};
    int read_pos{0};
    int ref_pos{record.beginPos};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&until] (int ref, int /*read*/) {return ref >= until;});

//    if (cigar_pos == 0) // no advancement happened. only case: until=1
//    {
//        erase(record.cigar, cigar_pos + 1, length(record.cigar)); // erase the rest
//        record.cigar[0].count = until;
//        return;
//    }
    --cigar_pos; // move to cigar operation that over stepped
    unsigned insertion_size{0};

    if (ref_pos == until)
    {
        --cigar_pos;
    }
    if (ref_pos > until)
    {
        if ((record.cigar[cigar_pos]).operation == 'M')
            insertion_size += (ref_pos - until);
        (record.cigar[cigar_pos]).count -= (ref_pos - until);
    }

    if ((record.cigar[cigar_pos]).operation == 'I')
    {
        --cigar_pos;
    }

    // count all events from cigar pos to have one big insertion
    for (unsigned idx = cigar_pos + 1; idx < length(record.cigar); ++idx)
        if ((record.cigar[idx]).operation != 'D') // only add isnerted and matched bases
            insertion_size += (record.cigar[idx]).count;

    erase(record.cigar, cigar_pos + 1, length(record.cigar)); // erase the cigar operations

    appendValue(record.cigar, CigarElement<char, unsigned>('I', insertion_size));
}

void truncate_cigar_left(BamAlignmentRecord & record, int until)
{
    unsigned num_truncated_bases{0};

    unsigned cigar_pos{0};
    int ref_pos{record.beginPos};
    int curr_read_pos{0};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     curr_read_pos,
                     record.cigar,
                     [&until] (int /*ref_pos*/, int read_pos) {return read_pos >= (until + 1);});

    num_truncated_bases += ref_pos - record.beginPos;
    erase(record.cigar, 0, cigar_pos - 1); // erase the overlapping beginning, Note that erase(0,4) = erase [0,4)

    if (curr_read_pos == (until + 1))
    {
        erase(record.cigar, 0); // default erase(..., pos, pos+1)
    }
    else if (curr_read_pos > (until + 1))
    {
        if ((record.cigar[0]).operation == 'M')
            num_truncated_bases -= (curr_read_pos - until - 1);

        (record.cigar[0]).count = (curr_read_pos - until - 1);
    }

    if (length(record.cigar) > 0           &&
        (record.cigar[0]).operation != 'M' &&
        (record.cigar[0]).operation != 'I') // meaning == H,S or D
    {
        if ((record.cigar[0]).operation == 'D')
        {   // continuing deletion in record must be added to the laarge DEL event
            num_truncated_bases += (record.cigar[0]).count;
            erase(record.cigar, 0);
        }
        else // meaning == H or S
        {
            (record.cigar[0]).operation = 'I'; // add those bases
        }
    }

    record.beginPos = record.beginPos + num_truncated_bases;
}

void transform_cigar_left_into_insertion(BamAlignmentRecord & record, int until)
{
    unsigned insertion_size{0};

    unsigned cigar_pos{0};
    int ref_pos{record.beginPos};
    int curr_read_pos{0};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     curr_read_pos,
                     record.cigar,
                     [&until] (int ref_pos, int /*read_pos*/) {return ref_pos >= until;});

    --cigar_pos; // go the former operation that did the overstepping

    if (ref_pos == until)
    {
        ++cigar_pos;
    }
    else if (ref_pos > until)
    {
        if ((record.cigar[cigar_pos]).operation == 'M')
            insertion_size += (record.cigar[cigar_pos]).count - (ref_pos - until);
        (record.cigar[cigar_pos]).count = (ref_pos - until);
    }

    if ((record.cigar[cigar_pos]).operation == 'I') // meaning == H,S or D
    {
        ++cigar_pos;
    }

    // count all events from cigar pos to have one big insertion
    for (unsigned idx = 0; idx < cigar_pos; ++idx)
        if ((record.cigar[idx]).operation != 'D') // only add isnerted and matched bases
            insertion_size += (record.cigar[idx]).count;

    erase(record.cigar, 0, cigar_pos); // erase the overlapping beginning,
    insert(record.cigar, 0, CigarElement<char, unsigned>('I', insertion_size));
}

BamAlignmentRecord merge_record_group(vector<BamAlignmentRecord> & record_group)
{
    SEQAN_ASSERT(record_group.size() > 0);

    // we now have all reads with the same name, given that the file was sorted
    // sort alignment by quality but put primary alignment on top.
    std::sort(record_group.begin(), record_group.end(), bamRecordQualityLess());

    BamAlignmentRecord final_record = record_group[0];
    // now check for supplementary alignment that can be merged

#ifndef NDEBUG
        if (record_group.size() > 1)
        {
            std::cout << "Merging group of " << record_group.size() << " with name " << final_record.qName << "\t" << endl;
            for (auto const & rec : record_group)
                cerr << "  -> "<< rec.qName << " " << rec.flag << " " << rec.rID << " " << rec.beginPos << " " << rec.mapQ << endl;
        }
#endif

        if (final_record.qName == "final_chr1_4939457_92")
            std::cerr << "oh oh";

    for (unsigned i = 1; i < record_group.size(); ++i)
    {
        // check if supplementary is fit for merging
        if (hasFlagSecondary(record_group[i])) // do not merge with secondary alignments
            continue;
        if (!hasFlagSupplementary(record_group[i])) // only want to merge supplementary alignments
            continue;
        if (final_record.rID != record_group[i].rID) // reference must be the same
            continue;
        if (hasFlagRC(final_record) != hasFlagRC(record_group[i])) // orientation must be the same
            continue;

        BamAlignmentRecord prim_record{final_record}; // copy so we do not screw with the final one right wawy
        BamAlignmentRecord supp_record{record_group[i]};

        if (supp_record.beginPos < prim_record.beginPos)
        { // ---supp-- before --final--
#ifndef NDEBUG
                cerr << "  -> ---supp-- before --final--" << endl;
#endif
            if (compute_map_end_pos(final_record.beginPos, final_record.cigar) <= compute_map_end_pos(supp_record.beginPos, supp_record.cigar))
                continue; // do not merge supp alignments that are completely covered

            // Alter cigar string
            // 1) remove soft clipping at the left end from the primary (right) alignment
            unsigned final_read_begin = get_read_begin_and_alter_cigar(prim_record);

            // 2) Cut supp (left) alignment at the right end such that the read
            // sequence does not overlap.
            truncate_cigar_right(supp_record, final_read_begin);
            if (length(supp_record.cigar) == 0) // no interesting bases left in supp alignment
                continue;

            // 3) Check if a Deletion or a Insertion must be added to connect the two alignments.
            //    This depends on the reference positions.
            if (compute_map_end_pos(supp_record.beginPos, supp_record.cigar) < prim_record.beginPos) // DEL
            {
                // 4) concatenate cropped cigar string to one big one with a deletion inside
                int deletion_size = prim_record.beginPos - compute_map_end_pos(supp_record.beginPos, supp_record.cigar);
                appendValue(supp_record.cigar, CigarElement<char, unsigned>('D', deletion_size));
            }
            else if (compute_map_end_pos(supp_record.beginPos, supp_record.cigar) > prim_record.beginPos)// INS
            {
                // 4) Compute bases that must inserted for the alignment to fit
                transform_cigar_right_into_insertion(supp_record, prim_record.beginPos);
            }

            append(supp_record.cigar, prim_record.cigar);

            // replace final variables
            prim_record.beginPos = supp_record.beginPos;
            prim_record.cigar = supp_record.cigar;
        }
        else if (prim_record.beginPos < supp_record.beginPos)
        { // ---final-- before --supp--
#ifndef NDEBUG
                cerr << "  -> ---final-- before --supp--" << endl;
#endif
            if (compute_map_end_pos(final_record.beginPos, final_record.cigar) >= compute_map_end_pos(supp_record.beginPos, supp_record.cigar))
                continue; // do not merge supp alignments that are completely covered

            // Alter cigar string
            // 1) remove soft clipping at the right end from the primary (left) alignment
            unsigned final_read_end = get_read_end_and_alter_cigar(prim_record);

            // 2) Cut supp (right) alignment at the left end such that the read
            // sequence does not overlap.
            truncate_cigar_left(supp_record, final_read_end);
            if (length(supp_record.cigar) == 0) // no interesting bases left in supp alignment
                continue;

            // 3) Check if a Deletion or a Insertion must be added to connect the two alignments.
            //    This depends on the reference positions.
            if (compute_map_end_pos(prim_record.beginPos, prim_record.cigar) < supp_record.beginPos) // DEL
            {
                // 4) concatenate cropped cigar string to one big one with a deletion inside
                int deletion_size = supp_record.beginPos - compute_map_end_pos(prim_record.beginPos, prim_record.cigar);
                appendValue(prim_record.cigar, CigarElement<char, unsigned>('D', deletion_size));
            }
            else
            {
                // 4) Compute bases that must inserted for the alignment to fit
                transform_cigar_left_into_insertion(supp_record, compute_map_end_pos(prim_record.beginPos, prim_record.cigar));
            }

            append(prim_record.cigar, supp_record.cigar);
        }

        if (length(prim_record.seq) != compute_fragment_length(prim_record.cigar))
        { // screwing with cigar was not successfull. Do not replace final_record
            cerr << "[ERROR] CIGAR and sequence length don't match: "
                 << length(prim_record.seq) << "(seq length) != "
                 << compute_fragment_length(prim_record.cigar)
                 << "(cigar length). Do not use merged record." << endl;
#ifndef NDEBUG
            throw std::exception();
#endif
        }
        else // screwing with cigar was successfull. Replace final_record
        {
            final_record = prim_record;
        }
    }

    return final_record;
}

vector<BamAlignmentRecord>
merge_alignments(vector<BamAlignmentRecord> const & records)
{
    vector<BamAlignmentRecord> merged_records;
    vector<BamAlignmentRecord> record_group;

    for (auto const & rec : records)
    {
        if (!record_group.empty() &&
            (record_group[record_group.size() - 1]).qName != rec.qName)
        {
            BamAlignmentRecord merged_record = merge_record_group(record_group);
            merged_records.push_back(merged_record);
            record_group.clear();
            record_group.push_back(rec);
        }
        else
        {
            record_group.push_back(rec);
        }

#ifndef NDEBUG
        if (!record_group.empty() && // avoid segementation fault
            (record_group[record_group.size() - 1]).qName > rec.qName) // check if records are sorted
        std::cout << "[ERROR] Merging alignments - Records are not sorted by name. "
                  << " Merging will not be successfull in this case." << endl;
#endif

    }
    //process last group
    BamAlignmentRecord merged_record = merge_record_group(record_group);
    merged_records.push_back(merged_record);

    return merged_records;
}
