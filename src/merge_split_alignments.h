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

unsigned truncate_cigar_left(BamAlignmentRecord & record, int until)
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

    if ((record.cigar[0]).operation != 'M' &&
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

    return num_truncated_bases;
}

unsigned truncate_cigar_left_ins(BamAlignmentRecord & record, int until)
{
    unsigned cigar_pos{0};
    int ref_pos{record.beginPos};
    int read_pos{0};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&until] (int ref, int /*read*/) {return ref >= until;});

    erase(record.cigar, 0, cigar_pos - 1); // erase the overlapping beginning

    if (ref_pos == (until + 1))
    {
        erase(record.cigar, 0);
    }
    else if (ref_pos > (until + 1))
    {
        (record.cigar[0]).count = ref_pos - until - 1;
    }

    return (read_pos - (ref_pos - until - 1));
}

//void merge_two_records(BamAlignmentRecord & prim, BamAlignmentRecord & supp, )
//{
//
//}


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

        if (final_record.qName == "final_chr6_110052570_19150")
            cerr << "hiiii";

    for (unsigned i = 1; i < record_group.size(); ++i)
    {
        BamAlignmentRecord supp_record{record_group[i]};

        if (hasFlagSecondary(supp_record)) // do not merge with secondary alignments
            continue;
        if (!hasFlagSupplementary(supp_record)) // only want to merge supplementary alignments
            continue;
        if (final_record.rID != supp_record.rID) // reference must be the same
            continue;
        if (hasFlagRC(final_record) != hasFlagRC(supp_record)) // orientation must be the same
            continue;

        if (compute_map_end_pos(supp_record.beginPos, supp_record.cigar) < final_record.beginPos)
        { // ---supp-- DEL --final--
#ifndef NDEBUG
                cerr << "  -> ---supp-- DEL --final--" << endl;
#endif
            // Alter cigar string
            // 1) remove soft clipping at the left end from the right alignment
            unsigned final_read_begin = get_read_begin_and_alter_cigar(final_record);
            // 2) Cut left alignment at the right end such that the read
            // seequence does not overlap.
            truncate_cigar_right(supp_record, final_read_begin);
            // 3) append clipped bases, because there can't be clipping in the middle

            // 4) concatenate cropped cigar string to one big one with a deletion inside
            int deletion_size = final_record.beginPos - compute_map_end_pos(supp_record.beginPos, supp_record.cigar);
            appendValue(supp_record.cigar, CigarElement<char, unsigned>('D', deletion_size));
            append(supp_record.cigar, final_record.cigar);

            // replace final variables
            final_record.beginPos = supp_record.beginPos;
            final_record.cigar = supp_record.cigar;

            // update mapping info needed for evaluataion
            BamTagsDict fin_tagsDict(final_record.tags);
            BamTagsDict sup_tagsDict(supp_record.tags);
            int fin_id{-1};
            int sup_id{-1};
            findTagKey(fin_id, fin_tagsDict, "NM");
            findTagKey(sup_id, sup_tagsDict, "NM");
            int fin_nm{0};
            int sup_nm{0};
            if (fin_id != -1)
                extractTagValue(fin_nm, fin_tagsDict, fin_id);
            if (sup_id != -1)
                extractTagValue(sup_nm, sup_tagsDict, sup_id);
            setTagValue(fin_tagsDict, "NM", fin_nm + sup_nm + deletion_size);
        }
        else if (compute_map_end_pos(final_record.beginPos, final_record.cigar) < supp_record.beginPos)
        { // ---final-- DEL --supp--
#ifndef NDEBUG
                cerr << "  -> ---final-- DEL --supp--" << endl;
#endif
            // Alter cigar string
            // 1) remove soft clipping at the right end from the primary (left) alignment
            unsigned final_read_end = get_read_end_and_alter_cigar(final_record);

            // 2) Cut supp (right) alignment at the left end such that the read
            // sequence does not overlap.
            unsigned num_truncated_bases = truncate_cigar_left(supp_record, final_read_end);
            // 3) append clipped bases, because there can't be clipping in the middle

            // 4) concatenate cropped cigar string to one big one with a deletion inside
            int deletion_size = supp_record.beginPos + num_truncated_bases - compute_map_end_pos(final_record.beginPos, final_record.cigar);
            appendValue(final_record.cigar, CigarElement<char, unsigned>('D', deletion_size));
            append(final_record.cigar, supp_record.cigar);

            // update mapping info needed for evaluataion
            BamTagsDict fin_tagsDict(final_record.tags);
            BamTagsDict sup_tagsDict(supp_record.tags);
            int fin_id{-1};
            int sup_id{-1};
            findTagKey(fin_id, fin_tagsDict, "NM");
            findTagKey(sup_id, sup_tagsDict, "NM");
            int fin_nm{0};
            int sup_nm{0};
            if (fin_id != -1)
                extractTagValue(fin_nm, fin_tagsDict, fin_id);
            if (sup_id != -1)
                extractTagValue(sup_nm, sup_tagsDict, sup_id);
            setTagValue(fin_tagsDict, "NM", fin_nm + sup_nm + deletion_size);
        }
//        else if (compute_map_end_pos(final_record.beginPos, final_record.cigar) > supp_record.beginPos &&
//                 final_record.beginPos < supp_record.beginPos)
//        { // ---final-- INS --supp-
//
//            // Alter cigar string
//            // 1) remove soft clipping at the right end from the left alignment
//            unsigned final_read_end = get_read_end_and_alter_cigar(final_record);
//            // 2) Cut right alignment at the left end such that the read
//            // seequence does not overlap.
//            unsigned supp_read_start = truncate_cigar_left_ins(supp_record, compute_map_end_pos(final_record.beginPos, final_record.cigar));
//            // 3) add an insertion
//            appendValue(final_record.cigar, CigarElement<char, unsigned>('I', supp_read_start - final_read_end - 1));
//            // 4) concatenate cropped cigar string to one big one with a deletion inside
//            append(final_record.cigar, supp_record.cigar);
//
//            // update mapping info needed for evaluataion
//            BamTagsDict fin_tagsDict(final_record.tags);
//            BamTagsDict sup_tagsDict(supp_record.tags);
//            int fin_id{-1};
//            int sup_id{-1};
//            findTagKey(fin_id, fin_tagsDict, "NM");
//            findTagKey(sup_id, sup_tagsDict, "NM");
//            int fin_nm{0};
//            int sup_nm{0};
//            if (fin_id != -1)
//                extractTagValue(fin_nm, fin_tagsDict, fin_id);
//            if (sup_id != -1)
//                extractTagValue(sup_nm, sup_tagsDict, sup_id);
//            setTagValue(fin_tagsDict, "NM", fin_nm + sup_nm);
//        }

#ifndef NDEBUG
        if (length(final_record.seq) != compute_fragment_length(final_record.cigar))
        {
            cerr << "[ERROR] CIGAR and sequence length don't match: "
                 << length(final_record.seq) << "(seq length) != "
                 << compute_fragment_length(final_record.cigar)
                 << "(cigar length)." << endl;
            throw std::exception();
        }
#endif
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
