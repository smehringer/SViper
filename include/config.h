#pragma once

#include <seqan/align.h>

//!\brief Command line option struct
struct CmdOptions
{
    bool verbose{false};
    bool veryVerbose{false};
    bool output_polished_bam{false};
    unsigned threads{0};
    int flanking_region{400}; // size of flanking region for breakpoints
    int mean_coverage_of_short_reads{36}; // original coverage
    double mean_insert_size_of_short_reads{280.054}; // original coverage
    double stdev_insert_size_of_short_reads{145.162};
    int length_of_short_reads{150}; // original coverage
    std::string long_read_file_name;
    std::string short_read_file_name;
    std::string candidate_file_name;
    std::string output_prefix;
    std::string reference_file_name;
    std::string log_file_name;

    CmdOptions(std::string long_, std::string short_, std::string candidate_, std::string ref_) :
        threads(std::thread::hardware_concurrency()),
        long_read_file_name(long_),
        short_read_file_name(short_),
        candidate_file_name(candidate_),
        output_prefix(candidate_ + "_polished"),
        reference_file_name(ref_) {}

    CmdOptions() = default;
    CmdOptions(const CmdOptions&) = default;
    CmdOptions(CmdOptions&&) = default;
    CmdOptions& operator=(const CmdOptions&) = default;
    CmdOptions& operator=(CmdOptions&&) = default;
};

/*! A global struct containing all the important information.
 * This struct holds information needed throughout the polishing process and
 * also includes the most important functions (thresholds, computations etc.).
 * Most changes in the important steps of the polishing routine might be made
 * here.
 */
struct SViperConfig
{
    bool verbose{false};

    // options to be set through the command line options
    int flanking_region; // size of flanking region for breakpoints
    int mean_coverage_of_short_reads; // original coverage
    double mean_insert_size_of_short_reads;
    double stdev_insert_size_of_short_reads;
    int length_of_short_reads;

    unsigned ref_flank_length{500}; // length to flank to the consensus sequence with the reference

    // mapping penalties for the final alignment
    double const MM{90000};   // Match
    double const MX{-80000};  // Mismatch
    double const GE{-1};      // Gap extension
    double const GO{-100000}; // Gap open

    // polishing parameter
    bool fix_indels{false};

    // quality statistics
    double baseQ_mean{0}; // will be overridden once when reading in reads
    double baseQ_std{1};  // will be overridden once when reading in reads
    double mappQ_mean{0}; // will be overridden every round after mapping
    double mappQ_std{1};  // will be overridden every round after mapping
    double alpha{1}; // scaling factor such that mapping and base quality are comparable
    double mean_coverage{0}; // will be overridden every round after mapping
    double min_coverage{4};
    std::vector<unsigned> cov_profile;

    /* Polishing statistics (total base count).
     * Note: The total number can exceed the length easily since bases can be
     * changed multiple times during different rounds. */
    unsigned substituted_bases{0};
    unsigned inserted_bases{0};
    unsigned deleted_bases{0};
    unsigned rounds{0};

    SViperConfig(CmdOptions & options) :
        verbose(options.verbose),
        flanking_region(options.flanking_region),
        mean_coverage_of_short_reads(options.mean_coverage_of_short_reads),
        mean_insert_size_of_short_reads(options.mean_insert_size_of_short_reads),
        stdev_insert_size_of_short_reads(options.stdev_insert_size_of_short_reads),
        length_of_short_reads(options.length_of_short_reads)
    {}

    bool insert_size_in_range(int end_pos_forward, int begin_pos_reverse) const
    {
        double range_start = std::max(0.0, mean_insert_size_of_short_reads - 3 * stdev_insert_size_of_short_reads);
        double range_end   = std::min(mean_insert_size_of_short_reads + 3 * mean_insert_size_of_short_reads,
                                      mean_insert_size_of_short_reads + 3 * stdev_insert_size_of_short_reads);

        return (begin_pos_reverse - end_pos_forward > range_start) &&
               (begin_pos_reverse - end_pos_forward < range_end);
    }

    /* [used in fill_profiles
     * Returns the value to be added to the profile. A profile consists of counts
     * for the letters A,C,T,G,- weighted by their quality. The profile is later
     * on used to determine the new reference sequence.*/
    double add_to_profile(double const & baseQ,
                          double const & mappQ,
                          bool is_in_propair_pair) const
    {
        if (is_in_propair_pair)
            return (baseQ * mappQ); // return double the score
        return ((baseQ * mappQ) / 2);
    }

    /* [used in build_consensus_from_profile]
     * Determines whether short read information is to be trusted
     * TRUE: if base with given score will be substituted in consensus sequence
     * FALSE: else*/
    bool score_passes_threshold(double const & score, unsigned const pos) const
    {
        /* return true if score is at least as good as half the coverage at this pos
         * of bases with half as good quality as the mean quality. */
        return score >= (std::max(min_coverage, (double)cov_profile[pos]/2) * (baseQ_mean * mappQ_mean));
    }

    /* Allowed Length Deviation.
     * Returns the minimum and maximum size a variant might have to be considered
     * the same variant. E.g. a deletion of size 1000 is considered polished but
     * the same, when the length is between 800-1200. */
    bool sv_length_passes_threshold(unsigned length, unsigned golden_length)
    {
        return (length >= golden_length - (0.8 * golden_length) &&
                length <= golden_length + (0.8 * golden_length));
    }
};
