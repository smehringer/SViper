#pragma once

#include <iostream>
#include <tuple>
#include <regex>
#include <string>
#include <sstream>

#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

using namespace std;
using namespace seqan;

bool verbose{false};
double const DEV_POS = 0.9;
double const DEV_SIZE = 0.8;


enum SV_TYPE
{
    UNKOWN=0,
    DEL=1,
    INS=2,
    DUP=3,
    INV=4,
    TRA=5
};

struct Variant
{
    Variant(const string & line)
    {
        stringstream ss(line);

        // read variant line into member variables
        ss >> ref_chrom;
        ss >> ref_pos;
        if (ss.fail()) // e.g. when quality was '.' and cast to double failed
            throw std::iostream::failure("ERROR when reading vcf file. Reference position could not be read for line " + line);
        ss >> id;
        ss >> ref_seq;
        ss >> alt_seq;
        ss >> quality;
        if (ss.fail()) // e.g. when quality was '.' and cast to double failed
            ss.clear();
        ss >> filter;
        ss >> info;
        ss >> format;
        ss >> samples;

        // check if all information are given


        if (ref_chrom.empty())
            throw std::iostream::failure("ERROR when reading vcf file. ref_chrom was not provided.");
        if (ref_pos == -1)
             throw std::iostream::failure("ERROR when reading vcf file. ref_pos was not provided.");
        if (id.empty())
            throw std::iostream::failure("ERROR when reading vcf file. id was not provided.");
        if (ref_seq.empty())
            throw std::iostream::failure("ERROR when reading vcf file. ref_seq was not provided.");
        if (alt_seq.empty())
            throw std::iostream::failure("ERROR when reading vcf file. alt_seq was not provided.");
        if (filter.empty())
            throw std::iostream::failure("ERROR when reading vcf file. filter was not provided.");
        if (info.empty())
            throw std::iostream::failure("ERROR when reading vcf file. info was not provided.");
        if (format.empty())
            throw std::iostream::failure("ERROR when reading vcf file. format was not provided.");
        if (samples.empty())
            throw std::iostream::failure("ERROR when reading vcf file. samples was not provided.");

        // determine sv_type from alt_seq
        // this assumes the vcf format to be have the <TYPE> tags
        if (alt_seq == "<DEL>")
            sv_type = SV_TYPE::DEL;
        else if (alt_seq == "<INS>")
            sv_type = SV_TYPE::INS;
        else if (alt_seq == "<DUP>")
            sv_type = SV_TYPE::DUP;
        else if (alt_seq == "<INV>")
            sv_type = SV_TYPE::INV;
        else if (alt_seq == "<TRA>")
            sv_type = SV_TYPE::TRA;
        else
            sv_type = SV_TYPE::UNKOWN;

        // determine END position from END info tag
        auto n = info.find("END=");
        if (n !=  std::string::npos)
        {
            stringstream ss(info.substr(n+4, info.find(';', n) - n - 4));
            ss >> ref_pos_end;                 // store end temporarily
            if (ss.fail())
                throw std::iostream::failure("ERROR when reading vcf file. END value "+
                                             info.substr(n+4, info.find(';', n) - n - 4)+
                                             " of variant " + ref_chrom + ":" + to_string(ref_pos) +
                                             " could not be read.");
            if (ref_pos_end < ref_pos)
                throw std::iostream::failure("ERROR when reading vcf file. END value "+
                                             info.substr(n+4, info.find(';', n) - n - 4)+
                                             " of variant " + ref_chrom + ":" + to_string(ref_pos) +
                                             " is larger than start pos " + to_string(ref_pos));

        }
        else
        {
            throw std::iostream::failure(string("ERROR when reading vcf file.") +
                                         " END tag not found in info field of variant "+
                                         ref_chrom + ":" + to_string(ref_pos) + "." +
                                         " INFO: " + info + "."
                                         );
        }

        n = info.find("SVLEN=");
        if (n != std::string::npos &&
            info.substr(n+6, info.find(';', n) - n - 6) != "NA" &&
            info.substr(n+6, info.find(';', n) - n - 6) != ".")
        {
            stringstream ss(info.substr(n+6, info.find(';', n) - n - 6));
            ss >> sv_length;
            if (ss.fail())
                throw std::iostream::failure("ERROR when reading vcf file. SVLEN value "+
                                             info.substr(n+6, info.find(';', n) - n - 6)+
                                             " of variant " + ref_chrom + ":" + to_string(ref_pos) +
                                             " could not be read.");
        }
        else // if SVLEN is not in the info, use END (will not work for insertions...)
        {
            sv_length = ref_pos_end - ref_pos; // calculate actual length
        }

    }

    Variant() = default;
    Variant(const Variant&) = default;
    Variant(Variant&&) = default;
    Variant& operator=(const Variant&) = default;
    Variant& operator=(Variant&&) = default;

    SV_TYPE sv_type;
    int     sv_length{-1};
    string  ref_chrom;
    int     ref_pos{-1};
    int     ref_pos_end{-1};
    string  id;
    string  ref_seq;
    string  alt_seq;
    double  quality{0};
    string  filter;
    string  info;
    string  format;
    string  samples;

    void write(std::ostream & stream)
    {
        stream << ref_chrom
               << "\t" << ref_pos
               << "\t" << id
               << "\t" << ref_seq
               << "\t" << alt_seq
               << "\t" << quality
               << "\t" << filter
               << "\t" << info
               << "\t" << format
               << "\t" << samples
               << "\n";
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

struct customLess
{
    bool operator()(BamAlignmentRecord lhs, BamAlignmentRecord rhs) const
    {
        if (!hasFlagSupplementary(lhs) && !hasFlagSecondary(lhs))
            return false;
        return lhs.qual < rhs.qual;
    }
};

struct recordQualityLess
{
    bool operator()(BamAlignmentRecord lhs, BamAlignmentRecord rhs) const
    {
        return lhs.qual < rhs.qual;
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

template <typename lambda_type>
void advance_in_cigar(unsigned & cigar_pos,
                      int & ref_pos,
                      int & read_pos,
                      String<CigarElement<char, unsigned>> const & cigar,
                      lambda_type && stop_criterion)
{
    while (cigar_pos < length(cigar))
    {
        if (stop_criterion(ref_pos, read_pos))
            break;

        if ((cigar[cigar_pos]).operation == 'M')
        {
            read_pos += (cigar[cigar_pos]).count;
            ref_pos  += (cigar[cigar_pos]).count;
        }
        else if ((cigar[cigar_pos]).operation == 'I' ||
                 (cigar[cigar_pos]).operation == 'S' ||
                 (cigar[cigar_pos]).operation == 'H')
        {
            read_pos += (cigar[cigar_pos]).count;
        }
        else // D
        {
            ref_pos += (cigar[cigar_pos]).count;
        }
        ++cigar_pos;
    }
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

    erase(record.cigar, cigar_pos, length(record.cigar)); // erase the rest

    if (read_pos > (until -1))
    {
        (record.cigar[length(record.cigar) - 1]).count -= (read_pos - until);
    }
}

unsigned truncate_cigar_left(BamAlignmentRecord & record, int until)
{
    unsigned num_truncated_bases{0};

    unsigned cigar_pos{0};
    int ref_pos{record.beginPos};
    int read_pos{0};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&until] (int /*ref*/, int read) {return read >= (until + 1);});

    num_truncated_bases += ref_pos - record.beginPos;
    erase(record.cigar, 0, cigar_pos - 1); // erase the overlapping beginning

    if (read_pos == (until + 1))
    {
        erase(record.cigar, 0);
    }
    else if (read_pos > (until + 1))
    {
        if ((record.cigar[0]).operation == 'M')
            num_truncated_bases -= (read_pos - until - 1);
        (record.cigar[0]).count = (read_pos - until - 1); // why -1 ?
    }
    return num_truncated_bases;
}

BamAlignmentRecord process_record_group(vector<BamAlignmentRecord> & record_group)
{
    // cout << (record_group[0]).qName << "\t" << record_group.size() << endl;

    // we now have all reads with the same name, given that the file was sorted
    // sort alignment by quality but put primary alignment on top.
    std::sort(record_group.begin(), record_group.end(), customLess());
    BamAlignmentRecord final_record = record_group[0];

    if (verbose)
        cout << "Merging group of " << record_group.size() << " with name " << final_record.qName << "\t" << endl;

    // now check for supplementary alignment that can be merged
    for (unsigned i = 1; i < record_group.size(); ++i)
    {
        BamAlignmentRecord supp_record{record_group[i]};

        if (!hasFlagSupplementary(supp_record)) // only want to merge supplementary alignments (not secondary)
            continue;
        if (final_record.rID != supp_record.rID) // reference must be the same
            continue;
        if ((hasFlagRC(final_record) && !hasFlagRC(supp_record)) || // orientation must be the same
            (!hasFlagRC(final_record) && hasFlagRC(supp_record)))
            continue;

        if (compute_map_end_pos(supp_record.beginPos, supp_record.cigar) < final_record.beginPos)
        {
            // supp before final, non overlapping
            if (verbose)
                cout << "Merging supplementary alignment on the left to primary on the right." << endl;

            // Alter cigar string
            // 1) remove soft clipping at the left end from the right alignment
            unsigned final_read_begin = get_read_begin_and_alter_cigar(final_record);
            // 2) Cut left alignment at the right end such that the read
            // seequence does not overlap.
            truncate_cigar_right(supp_record, final_read_begin);
            // 3) append clipped bases, because there can't be clipping in the middle
            if ((supp_record.cigar[length(supp_record.cigar) - 1]).operation != 'M' &&
                (supp_record.cigar[length(supp_record.cigar) - 1]).operation != 'I') // meaning == H,S or D
            {
                (supp_record.cigar[length(supp_record.cigar) - 1]).operation = 'M'; // add those bases
            }
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
        {
            // final before supp
            if (verbose)
                cout << "Merging supplementary alignment on the right to primary on the left." << endl;

            // Alter cigar string
            // 1) remove soft clipping at the right end from the left alignment
            unsigned final_read_end = get_read_end_and_alter_cigar(final_record);
            // 2) Cut right alignment at the left end such that the read
            // seequence does not overlap.
            unsigned num_truncated_bases = truncate_cigar_left(supp_record, final_read_end);
            // 3) append clipped bases, because there can't be clipping in the middle
            if ((supp_record.cigar[0]).operation != 'M' &&
                (supp_record.cigar[0]).operation != 'I') // meaning == H,S or D
            {
                (supp_record.cigar[0]).operation = 'M'; // add those bases
            }
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
    }

    if (verbose)
        if (length(final_record.seq) != compute_fragment_length(final_record.cigar))
            cerr << "[ERROR] CIGAR and sequence length don't match: "
                 << length(final_record.seq) << "(seq length) != "
                 << compute_fragment_length(final_record.cigar)
                 << "(cigar length)." << endl;

    return final_record;
}

bool is_same_sv_type(char cigar, SV_TYPE type)
{
    if ((cigar == 'D') && (type == SV_TYPE::DEL))
        return true;
    else if ((cigar == 'I') && (type == SV_TYPE::INS))
        return true;
    return false;
}

bool is_supporting(BamAlignmentRecord const & record, Variant const & variant)
{
    bool is_supporting{false};

    // first advance to region of interest (- buffer)
    unsigned cigar_pos{0};
    int read_pos{0};
    int ref_pos{record.beginPos};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&variant] (int ref, int /*read*/) {return ref >= (variant.ref_pos - DEV_POS * variant.sv_length);});

    // now look for variant (DEL/INS) until region of interest (+ buffer) is surpassed
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&] (int ref, int /*read*/)
                     {
                        if (is_same_sv_type((record.cigar[cigar_pos]).operation, variant.sv_type))
                        {
                            if ((record.cigar[cigar_pos]).count >= variant.sv_length - (DEV_SIZE * variant.sv_length) &&
                                (record.cigar[cigar_pos]).count <= variant.sv_length + (DEV_SIZE * variant.sv_length))
                            {
                                is_supporting = true;
                                return true; // stop criterion met.
                            }
                        }
                        return ref > (variant.ref_pos + variant.sv_length + DEV_POS * variant.sv_length);
                     }
                     );
    return is_supporting;
}

bool refine_variant(BamAlignmentRecord const & record, Variant & variant)
{
    bool has_variant{false};

    // first advance to region of interest (- buffer)
    unsigned cigar_pos{0};
    int read_pos{0};
    int ref_pos{record.beginPos};
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&variant] (int ref, int /*read*/) {return ref >= (variant.ref_pos - DEV_POS * variant.sv_length);});

    // now look for variant (DEL/INS) until region of interest (+ buffer) is surpassed
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&] (int ref, int /*read*/)
                     {
                        if (is_same_sv_type((record.cigar[cigar_pos]).operation, variant.sv_type))
                        {
                            if ((record.cigar[cigar_pos]).count >= variant.sv_length - (DEV_SIZE * variant.sv_length) &&
                                (record.cigar[cigar_pos]).count <= variant.sv_length + (DEV_SIZE * variant.sv_length))
                            {

                                has_variant = true;
                                return true; // stop criterion met.
                            }
                        }
                        return ref > (variant.ref_pos + variant.sv_length + DEV_POS * variant.sv_length);
                     }
                     );

    if (has_variant) // then advance_in_cigar stopped at this variant
    {
        variant.ref_pos = ref_pos;
        variant.sv_length = (record.cigar[cigar_pos]).count;



        // refine info
        std::regex svlen_re("SVLEN=[0-9]*;");
        std::regex end_re("END=[0-9]*;");
        std::regex seq_re("SEQ=[A-Za-z]*;");

        if (std::regex_search(variant.info, svlen_re))
            variant.info = std::regex_replace(variant.info, svlen_re, std::string("SVLEN=" + std::to_string(variant.sv_length) + ";"));
        else
            variant.info.append(std::string(";SVLEN=" + std::to_string(variant.sv_length)));

        if (variant.sv_type == SV_TYPE::DEL)
        {
            if (std::regex_search(variant.info, end_re))
                variant.info = std::regex_replace(variant.info, end_re, std::string("END=" + std::to_string(record.beginPos + variant.sv_length) + ";"));
            else
                variant.info.append(std::string(";END=" + std::to_string(record.beginPos + variant.sv_length)));
        }

        if (variant.sv_type == SV_TYPE::INS)
        {
            stringstream ss;
            ss << "SEQ=" << infix(record.seq, read_pos,  read_pos + variant.sv_length);

            if (std::regex_search(variant.info, seq_re))
               variant.info = std::regex_replace(variant.info, seq_re, ss.str() + ";");
            else
                variant.info.append(string(";" + ss.str()));
        }

    }

    return has_variant;
}

tuple<int, int> get_read_region_boundaries(BamAlignmentRecord const & record,
                                           int ref_region_begin,
                                           int ref_region_end)
{
    int read_region_begin;
    int read_region_end;

    if (ref_region_begin >= ref_region_end)
    {
        cerr << "[GET READ REGION ERROR]"
             << " for record " << record.qName
             << " Start (" << ref_region_begin
             << ") >= End (" << ref_region_end << ")" << endl;
             return make_tuple(0, 0);
    }

    // advance to begin of region of interest
    unsigned cigar_pos{0};
    int read_pos{-1};                // -1 because ref position is 0 based
    int ref_pos{record.beginPos -1}; // -1 because ref position is 0 based

    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&ref_region_begin] (int ref, int /*read*/) {return ref >= ref_region_begin;});

    read_region_begin = (read_pos > 0 ) ? read_pos : 0; // in case mapping pos is inside ref span

    // advance to end of region of interest
    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&ref_region_end] (int ref, int /*read*/) {return ref >= ref_region_end;});

    // calculate region end pos in read
    --cigar_pos;
    if ((record.cigar[cigar_pos]).operation == 'M')
    {
        read_region_end = read_pos - (ref_pos - ref_region_end - 1);
    }
    else if ((record.cigar[cigar_pos]).operation == 'D')
    {
        read_region_end = read_pos + 1;
    }
    else // I
    {
        read_region_end = read_pos; // should never happen.. but just in case
    }

    return make_tuple(read_region_begin, read_region_end);
}

inline DnaString build_consensus(StringSet<DnaString> const & seqs,
                                 vector<double> const & mapQ) // mapping qualities
{
    Align<DnaString> align;
    resize(rows(align), length(seqs));
    for (unsigned i = 0; i < length(seqs); ++i)
        assignSource(row(align, i), seqs[i]);

    globalMsaAlignment(align, SimpleScore(5, -3, -1, -3));

    DnaString consensus;
    String<ProfileChar<Dna, double> > profile;
    resize(profile, 1);

    for (unsigned i = 0; i < length(row(align, 0)); ++i)
    {
        for (unsigned rowNo = 0; rowNo < length(seqs); ++rowNo)
            profile[0].count[ordValue(getValue(row(align, rowNo), i))] += (1.0 * mapQ[rowNo] / 100);

        int idx = getMaxIndex(profile[0]);

        if (idx < 4)  // is not gap
            appendValue(consensus, Dna(getMaxIndex(profile[0])));

        // clear profile
        //(Note: clear(profile) compiles but has not the correct effect
        for (auto & c : profile[0].count)
            c = 0;
    }

    return consensus;
}

double compute_identity_measure(BamAlignmentRecord const & record)
{
    BamTagsDict tagsDict(record.tags);
    int id{-1};
    findTagKey(id, tagsDict, "NM");

    if (id == -1)
        return -1;

    double nm{0};
    extractTagValue(nm, tagsDict, id);

    return (1 - (nm/(length(record.seq) - 10000))) * 100;
}

double compute_fuzzyness_measure(BamAlignmentRecord const & record)
{
    unsigned count{0};

    // go over the whole cigar and add up the number of insertions and deletions
    unsigned cigar_pos{0};
    int read_pos{-1};                // -1 because ref position is 0 based
    int ref_pos{record.beginPos -1}; // -1 because ref position is 0 based

    advance_in_cigar(cigar_pos,
                     ref_pos,
                     read_pos,
                     record.cigar,
                     [&] (int /**/, int /**/)
                     {
                        if ((record.cigar[cigar_pos]).operation != 'M')
                            ++count; // count insertions and deletions
                        return false; // never stop. Go over whole cigar
                     });

    return (1.0 - ((double)count / ((length(record.seq) - 10000) * 0.3 / 1.3) )) * 100.0;
}

void assign_quality(BamAlignmentRecord const & record,
                    Variant & variant,
                    bool true_variant)
{
    double identity = compute_identity_measure(record);
    double fuzzyness = compute_fuzzyness_measure(record);

    // subtract variant length from score so to not penalize the variant
    // that one actually expects (e.g. deletion of 1000bp's otherwise screws up
    // the edit distance)
    if (true_variant)
        identity += ((double)variant.sv_length / (length(record.seq) - 10000)) * 100;

    // write seperate scores into info field of the variant
    ostringstream ss;
    ss << ";IDENTITIY_SCORE=" << identity
       << ";FUZZYNESS_SCORE=" << fuzzyness;
    variant.info.append(ss.str());

    variant.quality = (identity + fuzzyness) / 2;
}

Variant evaluate_alignment(BamAlignmentRecord const & record,
                           map<string, Variant> const & variant_map)
{
    // in the pipeline the final.fa sequnces will have the chr position and if
    // of the variant they correspond to in their name plit by '_'.
    // This way they can be easily matched to the former variant.
    stringstream ss;
    ss << record.qName;

    // extract read information
    string name;
    string chrom;
    string pos;
    string id;

    getline(ss, name, '_');
    getline(ss, chrom, '_');
    getline(ss, pos, '_');
    getline(ss, id, '_');

    // now find the corresponding var
    auto it = variant_map.find(id);

    if (it == variant_map.end())
    {
        cerr << "Variant with id '" << id
             << "' could not be found in vcf file." << endl;
        return Variant(string(chrom + "\t" + pos + "\t" + id +
                              "\tN\tN\t.\tNOTFOUND\t*\t*\t*"));
    }

    Variant var = it->second;

    if (var.ref_chrom != chrom || to_string(var.ref_pos) != pos)
    {
        cerr << "Variant with id '" << id << "' does not match with the read"
             << " name information: "
             << var.ref_chrom << ":" << var.ref_pos << " "
             << "vs. "
             << chrom << ":" << pos << "." << endl;
    }

    if (!refine_variant(record, var))
    {
        var.filter = "FAIL2";
        assign_quality(record, var, false);
    }
    else
    {
        assign_quality(record, var, true);
    }

    return var;
}
