#pragma once

#include <sstream>
#include <regex>

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
        if (ss.fail()) // e.g. when pos was '*' and cast to int failed
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

bool is_same_sv_type(char cigar_operation, SV_TYPE type)
{
    if ((cigar_operation == 'D') && (type == SV_TYPE::DEL))
        return true;
    else if ((cigar_operation == 'I') && (type == SV_TYPE::INS))
        return true;
    return false;
}

bool record_supports_variant(BamAlignmentRecord const & record, Variant const & variant)
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

bool refine_variant(seqan::BamAlignmentRecord const & record, Variant & variant)
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
                variant.info = std::regex_replace(variant.info, end_re, std::string("END=" + std::to_string(variant.ref_pos + variant.sv_length) + ";"));
            else
                variant.info.append(std::string(";END=" + std::to_string(variant.ref_pos + variant.sv_length)));
        }

        if (variant.sv_type == SV_TYPE::INS)
        {
            stringstream ss;
            ss << "SEQ=" << seqan::infix(record.seq, read_pos, read_pos + variant.sv_length);

            if (std::regex_search(variant.info, seq_re))
               variant.info = std::regex_replace(variant.info, seq_re, ss.str() + ";");
            else
                variant.info.append(string(";" + ss.str()));
        }

    }

    return has_variant;
}
