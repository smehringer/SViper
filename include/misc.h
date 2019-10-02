#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <memory>
#include <thread>
#include <limits>

#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/align.h>
#include <seqan/graph_msa.h>

struct file_info{
    std::ofstream                           log_file;
    CmdOptions                              options;
    std::vector<std::unique_ptr<BamFileIn>> long_read_file_handles;
    BamHeader                               long_read_header;   // The bam header object needed to fill bam context
    BamIndex<Bai>                           long_read_bai;
    std::vector<std::unique_ptr<BamFileIn>> short_read_file_handles;
    BamHeader                               short_read_header;  // The bam header object needed to fill bam context
    BamIndex<Bai>                           short_read_bai;
    std::vector<std::unique_ptr<FaiIndex>>  faidx_file_handles;
    std::vector<seqan::BamAlignmentRecord>  polished_reads;
};
