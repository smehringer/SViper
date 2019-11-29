#pragma once

#include <sstream>
#include <regex>

#include <sviper/basics.h>

namespace sviper
{
double const DEV_SIZE = 0.6;
double const DEV_POS = 0.5;

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
    // Constructor, destructor, and assignment
    // -------------------------------------------------------------------------
    Variant() = default;
    Variant(const Variant&) = default;
    Variant(Variant&&) = default;
    Variant& operator=(const Variant&) = default;
    Variant& operator=(Variant&&) = default;

    Variant(const std::string & line)
    {
        std::stringstream ss(line);
        std::string dummy{};

        // read variant line into member variables
        ss >> ref_chrom;
        ss >> ref_pos;
        if (ss.fail()) // e.g. when pos was '*' and cast to int failed
            throw_verbose_exception("VCF-ERROR. ref_pos could not be read for line " + line);
        ss >> id;
        ss >> ref_seq;
        ss >> alt_seq;
        ss >> quality;
        if (ss.fail()) // e.g. when quality was '.' and cast to double failed
        {
            ss.clear();
            //ss >> dummy;
        }
        ss >> filter;
        ss >> info;
        ss >> format;
        while (!ss.eof()) // read the rest
        {
            ss >> dummy;
            samples.push_back(dummy);
        }

        // check if all information are given
        if (ref_chrom.empty())
            throw_verbose_exception("VCF-ERROR. ref_chrom is empty.\nLine:" + line + "\n");
        if (ref_pos == -1)
             throw_verbose_exception("VCF-ERROR. ref_pos is empty.\nLine:" + line + "\n");
        if (id.empty())
            throw_verbose_exception("VCF-ERROR. id is empty.\nLine:" + line + "\n");
        if (ref_seq.empty())
            throw_verbose_exception("VCF-ERROR. ref_seq is empty.\nLine:" + line + "\n");
        if (alt_seq.empty())
            throw_verbose_exception("VCF-ERROR. alt_seq is empty.\nLine:" + line + "\n");
        if (filter.empty())
            throw_verbose_exception("VCF-ERROR. filter is empty.\nLine:" + line + "\n");
        if (info.empty())
            throw_verbose_exception("VCF-ERROR. info is empty.\nLine:" + line + "\n");
        if (format.empty())
            throw_verbose_exception("VCF-ERROR. format is empty.\nLine:" + line + "\n");
        if (samples.empty())
            throw_verbose_exception("VCF-ERROR. samples is empty.\nLine:" + line + "\n");

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

        if (sv_type == SV_TYPE::DEL || sv_type == SV_TYPE::INS) // no check for variants that are not processed anyway
        {
            // determine END position from END info tag
            auto const end_n = info.find("END=");
            auto const len_n = info.find("SVLEN=");

            if (end_n != std::string::npos)
            {
                std::stringstream ss(info.substr(end_n+4, info.find(';', end_n) - end_n - 4));
                ss >> ref_pos_end;                 // store end temporarily
                if (ss.fail())
                    throw_verbose_exception("VCF-ERROR. END value "+
                                            info.substr(end_n+4, info.find(';', end_n) - end_n - 4)+
                                            " of variant " + ref_chrom + ":" + std::to_string(ref_pos) +
                                            " could not be read.");
            }

            if (len_n != std::string::npos &&
                info.substr(len_n+6, info.find(';', len_n) - len_n - 6) != "NA" &&
                info.substr(len_n+6, info.find(';', len_n) - len_n - 6) != ".")
            {
                std::stringstream ss(info.substr(len_n+6, info.find(';', len_n) - len_n - 6));
                ss >> sv_length;
                if (ss.fail())
                    throw_verbose_exception("VCF-ERROR. SVLEN value "+
                                            info.substr(len_n+6, info.find(';', len_n) - len_n - 6)+
                                            " of variant " + ref_chrom + ":" + std::to_string(ref_pos) +
                                            " could not be read.");
                sv_length = std::abs(sv_length); // some tools report a negative length since bases were deleted

                if (end_n == std::string::npos) // no end tag but sv_len tag
                {
                    if (sv_type == SV_TYPE::DEL)
                        ref_pos_end = ref_pos + sv_length - 1;
                    else // IMPORTANT: this is correct for <INS> but undefined for all other types which are currently not supported
                        ref_pos_end = ref_pos;
                }
            }
            else // if SVLEN is not in the info, use END (will not work for insertions...)
            {
                if (end_n == std::string::npos)
                {
                    throw_verbose_exception("VCF-ERROR."
                                            "neither END nor SVLEN tag not found in info field.");
                }

                sv_length = ref_pos_end - ref_pos; // calculate actual length
            }
        }
    }

    // Member variables
    // -------------------------------------------------------------------------
    SV_TYPE sv_type{};
    int     sv_length{-1};
    std::string  ref_chrom{};
    int     ref_pos{-1};
    int     ref_pos_end{-1};
    std::string  id{};
    std::string  ref_seq{};
    std::string  alt_seq{};
    double  quality{0};
    std::string  filter{};
    std::string  info{};
    std::string  format{};
    std::vector<std::string>  samples{};

    // Public Member Function
    // -------------------------------------------------------------------------
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
               << "\t" << format;
        for (auto sample : samples)
            stream << "\t" << sample;
        stream << "\n";
    }

private:
    // Private functions
    // -------------------------------------------------------------------------
    // Make error messages more expressive by adding the full VCF record info
    void throw_verbose_exception(std::string const & what)
    {
        std::ostringstream os{};
        os << what << "Read: ";
        (*this).write(os);
        throw std::iostream::failure(os.str());
    }
};

inline bool is_same_sv_type(char cigar_operation, SV_TYPE type)
{
    if ((cigar_operation == 'D') && (type == SV_TYPE::DEL))
        return true;
    else if ((cigar_operation == 'I') && (type == SV_TYPE::INS))
        return true;
    return false;
}

inline bool record_supports_variant(seqan::BamAlignmentRecord const & record, Variant const & variant)
{
    bool is_supporting{false};

    // first advance to beginning of the region of interest.
    // The region of interest starts at the allowed pos deviation.
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

inline bool refine_variant(seqan::BamAlignmentRecord const & record, Variant & variant)
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

    if (has_variant) // if true, then advance_in_cigar stopped at this cigar_pos
    {
        variant.ref_pos = ref_pos;
        variant.sv_length = (record.cigar[cigar_pos]).count;

        if (variant.sv_type == SV_TYPE::INS)
            variant.ref_pos_end = ref_pos;
        else // DEL
            variant.ref_pos_end = ref_pos + variant.sv_length;

        // refine info
        std::regex svlen_neg_re("SVLEN=-[0-9]*");
        std::regex svlen_re("SVLEN=[0-9]*");
        std::regex end_re("END=[0-9]*");
        std::regex seq_re("SEQ=[A-Za-z]*");

        if (std::regex_search(variant.info, svlen_neg_re))
            variant.info = std::regex_replace(variant.info, svlen_neg_re, std::string("SVLEN=-" + std::to_string(variant.sv_length)));
        else if (std::regex_search(variant.info, svlen_re))
            variant.info = std::regex_replace(variant.info, svlen_re, std::string("SVLEN=" + std::to_string(variant.sv_length)));
        else
            variant.info.append(std::string(";SVLEN=" + std::to_string(variant.sv_length)));

        if (std::regex_search(variant.info, end_re))
            variant.info = std::regex_replace(variant.info, end_re, std::string("END=" + std::to_string(variant.ref_pos_end)));
        else
            variant.info.append(std::string(";END=" + std::to_string(variant.ref_pos_end)));

        if (variant.sv_type == SV_TYPE::INS)
        {
            std::stringstream ss;
            ss << "SEQ=" << seqan::infix(record.seq, read_pos, read_pos + variant.sv_length);

            if (std::regex_search(variant.info, seq_re))
               variant.info = std::regex_replace(variant.info, seq_re, ss.str());
            else
                variant.info.append(std::string(";" + ss.str()));
        }

    }

    return has_variant;
}
} // namespace sviper
