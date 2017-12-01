#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>
#include <algorithm>

#include <basics.h>
#include <variant.h>

using namespace std;

// global variable of the allowed deviation until SV's are still considered to be
// the same
const double ALLOWED_POS_DEVIATION = 0.95; // percent of the sv_length
const double ALLOWED_LENGTH_DEVIATION = 0.8;
const double SCORE_DEVIATION = 0.05;
const double QUALITY_CUTOFF = 65.0;


struct Statistics
{
    vector<unsigned> total_numbers{0, 0, 0, 0, 0, 0};
    vector<double> sum_scores_after{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    vector<unsigned> after_could_not_be_scored{0, 0, 0, 0, 0, 0};
    unsigned variant_in_before_but_not_in_after{0};
};

// IMPORTANT: Those operators are unusual and can e.g. not be used for sorting
bool operator ==(const Variant& lhs, const Variant& rhs)
{
    int length_dev = max(80.0, (ALLOWED_LENGTH_DEVIATION * ((lhs.sv_length + rhs.sv_length)/2)));
    int pos_dev = max(100.0, (ALLOWED_POS_DEVIATION * ((lhs.sv_length + rhs.sv_length)/2)));

    return (lhs.ref_chrom == rhs.ref_chrom) &&
           (lhs.sv_type == rhs.sv_type) &&
           (abs(lhs.sv_length - rhs.sv_length) < length_dev) &&
           (abs(lhs.ref_pos - rhs.ref_pos) < pos_dev);
}

bool operator <(const Variant& lhs, const Variant& rhs)
{
    return lhs.ref_pos < rhs.ref_pos;
}

bool operator >(const Variant& lhs, const Variant& rhs)
{
    return lhs.ref_pos > rhs.ref_pos;
}

void add_pairing_stats(Statistics & stats,
                       Variant const & curr_var_after)
{
    SV_TYPE t{curr_var_after.sv_type}; // the type is used to index the vectors

    ++stats.total_numbers[t];
    stats.sum_scores_after[t] += curr_var_after.quality;

    if (curr_var_after.quality == 0.0)
        ++stats.after_could_not_be_scored[t];
}

unsigned sum(vector<unsigned> v)
{
    unsigned sum{0};
    for (auto n : v)
        sum+=n;
    return sum;
}

int main(int argc, char ** argv)
{
    unsigned return_code{0};
    if (argc == 1 || string(argv[1]) == "-h" || string(argv[1]) == "-help" || string(argv[1]) == "--help")
    {
        cout << endl
             << "ectractThemAll - A tiny helper to extract reads by ID from a SAM file." << endl
             << "----------------------------------------------------------------------" << endl
             << endl
             << "USAGE: SORTED_IDS_FILE  SORTED_SAM_FILE  OUT_FILE_NAME" << endl
             << endl
             << "SORTED_IDS_FILE\t An alphabetically sorted file with one id per line." << endl
             << "SORTED_SAM_FILE\t An alphabetically sorted sam file containing all the reads." << endl
             << "OUT_FILE_NAME\t Provide an file name for the ouput sam file." << endl
             << endl
             << "Note: The script will notify you in case not all ids were found." << endl
             << "      It is important that the files are alphabetically ordered otherwise" << endl
             << "      it wont work. In case you sorted the file with samtools sort and" << endl
             << "      there are errors, try sorting it with unix sort." << endl
             << "----------------------------------------------------------------------"
             << endl;
             return 0;
    }
    else if (argc < 2 || argc > 3)
    {
        cout << endl
             << "USAGE: SORTED_IDS_FILE  SORTED_SAM_FILE  OUT_FILE_NAME" << endl
             << endl
             << "Or type -h/--help for more information." << endl
             << endl;
             return 0;
    }

    bool golden{false}; // for readability on further if clauses
    ifstream vcf_file_after(argv[1]);  // must be sorted by position grouped by chromosome
    ifstream vcf_file_golden;
    ofstream log_file("comparing_vcf_files.log");

    if (argc == 3) // golden vcf was provided
    {
        golden = true;
        vcf_file_golden.open(argv[2]); // must be sorted by position grouped by chromosome
    }

    if (!vcf_file_after.is_open())
    {
        cerr << "[ERROR] Could not open file " << argv[1] << endl;
        return 1;
    }

    if (golden && !vcf_file_golden.is_open())
    {
        cerr << "[ERROR] Could not open file " << argv[2] << endl;
        vcf_file_after.close();
        return 1;
    }

    Variant curr_var_after;  // stores information of the current vcf line in file vcf_after
    Variant curr_var_golden; // stores information of the current vcf line in file vcf_golden
    Statistics stats;
    Statistics golden_stats; // only used of golden vcf is provided

    // Initialization:
    // -------------------------------------------------------------------------
    // Skip all headers and init variant with first line in file
    string dummy_line{"#"}; // used throughout processing to read into

    while (dummy_line[0] == '#')
        getline(vcf_file_after, dummy_line);
    curr_var_after = Variant(dummy_line);

    if (golden)
    {
        dummy_line = "#";
        while (dummy_line[0] == '#')
            getline(vcf_file_golden, dummy_line);
        curr_var_golden = Variant(dummy_line);
    }

    // Iteration:
    // -------------------------------------------------------------------------
    // Now go through all vcf files trying to pair the variants
    // stayed the same = total - impr - wors - notScored
    vector<vector<unsigned>> bp_deviation_after_to_golden{{}, {}, {}, {}, {}, {}}; // for golden
    vector<unsigned> variants_unique_to_golden{0, 0, 0, 0, 0, 0};

    while (true) // will be broken if one of the vcf file is at end
    {
        while (golden &&
               !(curr_var_golden == curr_var_after) && // != is not equal to < so both must be tested
               (curr_var_golden < curr_var_after))
        {
            log_file << "2"
                     << "\t" << "."
                     << "\t" << "."
                     << "\t" << "."
                     << "\t" << "."
                     << "\t" << "."
                     << "\t" << "."
                     << "\t" << "FN:" << curr_var_golden.sv_type << ":" << curr_var_golden.ref_chrom << ":" << curr_var_golden.ref_pos << ":" << curr_var_golden.sv_length
                     << endl;

            ++variants_unique_to_golden[curr_var_golden.sv_type];

            if (!getline(vcf_file_golden, dummy_line))
            {
                golden = false;
            }
            else
            {
                curr_var_golden = Variant(dummy_line);
            }

        }

        if (curr_var_after.filter == "PASS" && // variant passed polishing
            curr_var_after.quality >= QUALITY_CUTOFF)
        {
            string golden_info{"NO_INFO"};

            if (golden)
                golden_info = "FP"; // default, to be changed below

            if (golden && curr_var_golden == curr_var_after)
            {
                assert(curr_var_after.ref_chrom == curr_var_after.ref_chrom);
                assert(curr_var_after.sv_type == curr_var_after.sv_type);

                golden_info = "TP:" + to_string(curr_var_golden.ref_pos);
                add_pairing_stats(golden_stats, curr_var_after);
                (bp_deviation_after_to_golden[curr_var_golden.sv_type]).push_back(abs(curr_var_golden.ref_pos - curr_var_after.ref_pos));

                if (!getline(vcf_file_golden, dummy_line))
                    golden = false;
                else
                    curr_var_golden = Variant(dummy_line);
            }

            log_file << "0"
                     << "\t" << curr_var_after.sv_type
                     << "\t" << curr_var_after.ref_chrom
                     << "\t" << curr_var_after.ref_pos
                     << "\t" << curr_var_after.sv_length
                     << "\t" << curr_var_after.quality
                     << "\t" << golden_info
                     << endl;

            add_pairing_stats(stats, curr_var_after);
        }
        else // variant failed polishing
        {
            string golden_info{"NO_INFO"};

            if (golden)
                golden_info = "TN"; // default. to be changed in next if clause

            if (golden && curr_var_golden == curr_var_after)
            {
                golden_info = "FN:" + to_string(curr_var_golden.ref_pos);
                ++golden_stats.variant_in_before_but_not_in_after;

                if (!getline(vcf_file_golden, dummy_line))
                    golden = false;
                else
                    curr_var_golden = Variant(dummy_line);
            }

            log_file << "1"
                     << "\t" << curr_var_after.sv_type
                     << "\t" << curr_var_after.ref_chrom
                     << "\t" << curr_var_after.ref_pos
                     << "\t" << curr_var_after.sv_length
                     << "\t" << curr_var_after.quality
                     << "\t" << golden_info
                     << endl;

            ++stats.variant_in_before_but_not_in_after;
        }

        // update variant
        if (!getline(vcf_file_after, dummy_line))
            break;
        curr_var_after = Variant(dummy_line);
    }

    // Finalizing:
    // -------------------------------------------------------------------------

    // print stats
    cout.precision(4);
    cout << "------------------------------------------------------------------------------------------------------------------" << endl
         << "                                              General Statistics" << endl
         << "------------------------------------------------------------------------------------------------------------------" << endl
         << endl
         << "==================================================================================================================" << endl
         << "     what             |\tUNKOWN\t|\tDEL\t|\tINS\t|\tDUP\t|\tINV\t|\tTRA\t|" << endl
         << "------------------------------------------------------------------------------------------------------------------" << endl
         << "Total number          |\t";
    for (auto num : stats.total_numbers)
        cout << num << "\t|\t";
    cout << endl;
    cout << "Average Score         |\t";
    for (unsigned i = 0; i < stats.total_numbers.size(); ++i)
        if (stats.total_numbers[i] > 0)
            cout << ((stats.sum_scores_after[i])/(stats.total_numbers[i])) << "\t|\t";
        else
            cout << "0\t|\t";
    cout << endl;
    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Not scored (After)    |\t";
    for (auto num : stats.after_could_not_be_scored)
        cout << num << "\t|\t";
    cout << endl;

    cout << "==================================================================================================================" << endl;
    cout << "Number of Variants that were \"polished away\": " << stats.variant_in_before_but_not_in_after << endl
         << "------------------------------------------------------------------------------------------------------------------" << endl;


    cout << "------------------------------------------------------------------------------------------------------------------" << endl
         << "                                                 Golden Statistics" << endl
         << "------------------------------------------------------------------------------------------------------------------" << endl
         << endl
         << "==================================================================================================================" << endl
         << "     what             |\tUNKOWN\t|\tDEL\t|\tINS\t|\tDUP\t|\tINV\t|\tTRA\t|" << endl
         << "------------------------------------------------------------------------------------------------------------------" << endl
         << "True Positives        |\t";
    for (auto num : golden_stats.total_numbers)
        cout << num << "\t|\t";
    cout << endl;
    cout << "Average Score         |\t";
    for (unsigned i = 0; i < golden_stats.total_numbers.size(); ++i)
        if (golden_stats.total_numbers[i] > 0)
            cout << ((golden_stats.sum_scores_after[i])/(golden_stats.total_numbers[i])) << "\t|\t";
        else
            cout << "0\t|\t";
    cout << endl;
    cout << "Mean Bp Dev After     |\t";
    for (unsigned i = 0; i < golden_stats.total_numbers.size(); ++i)
    {
        if ((bp_deviation_after_to_golden[i]).size() > 0)
        {
            std::nth_element((bp_deviation_after_to_golden[i]).begin(), (bp_deviation_after_to_golden[i]).begin() + (bp_deviation_after_to_golden[i]).size()/2, (bp_deviation_after_to_golden[i]).end());
            cout << (bp_deviation_after_to_golden[i])[(bp_deviation_after_to_golden[i]).size()/2]
                 << "/" << ((sum(bp_deviation_after_to_golden[i])+1)/(golden_stats.total_numbers[i]+1)-1) << "\t|\t";
        }
        else
        {
            cout << "-\t|\t";
        }
    }
    cout << endl;

    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "Not scored (After)    |\t";
    for (auto num : golden_stats.after_could_not_be_scored)
        cout << num << "\t|\t";
    cout << endl;

    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "False Negatives      |\t";
    for (auto num : variants_unique_to_golden)
        cout << num << "\t|\t";
    cout << endl;

    cout << "==================================================================================================================" << endl;
    cout << "Number of Variants  that were \"polished away\": " << golden_stats.variant_in_before_but_not_in_after << endl
         << endl
         << "------------------------------------------------------------------------------------------------------------------" << endl;


	vcf_file_after.close();
	if (golden)
    	vcf_file_golden.close();

    return return_code;
}
