// variant_caller.hpp
#ifndef VARIANT_CALLER_H
#define VARIANT_CALLER_H

#include <string>
#include <map>
#include <variant>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>

class VariantCaller {
public:
    VariantCaller();
    ~VariantCaller();
    bool help;
    bool parse_arguments(int argc, char **argv);
    void print_usage();
    bool run();

private:
    // Options
    bool keep_read_extremities;
    std::string base_csv_file;
    std::string indel_csv_file;
    std::string prefix_out;
    int min_qual;
    bool is_r1_rev;
    bool is_r2_rev;
    std::string bam_input;
    std::string fasta_reference;
    float min_freq;
    std::string call_strand;
    int min_count;
    int min_alt_count;
    bool skip_secondary;
    bool skip_duplicate;
    int max_n_pileup;

    // Data structures
    faidx_t *fai;
    samFile *in;
    bam_hdr_t *header;

    void process_pileup();
    static int pileup_callback(uint32_t tid, hts_pos_t pos, int n, const bam_pileup1_t *plp, void *data);
    void update_counts(uint32_t tid, int pos, int n, const bam_pileup1_t *plp, std::map<int, std::pair<int, int>>& del_depth);

    // Counting structures
    struct BaseCounts {
        int depth_fw;
        int depth_rv;
        int A_fwd, A_rev;
        int T_fwd, T_rev;
        int C_fwd, C_rev;
        int G_fwd, G_rev;
        int N_fwd, N_rev;
    };
    struct IndelCounts {
        std::map<std::string, int> fwd_counts;
        std::map<std::string, int> rev_counts;
    };
    std::map<std::string, BaseCounts> base_counts_map;
    std::map<std::string, IndelCounts> indel_counts_map;

    // Helper methods
    bool is_read1(const bam1_t *b);
    bool is_read2(const bam1_t *b);
    bool is_duplicate(const bam1_t *b);
    bool is_secondary(const bam1_t *b);
    std::string tid_to_name(int32_t tid);

    // Output methods
    void write_base_csv();
    void write_indel_csv();
    void call_variants();
    void make_a_call(
        const std::string& region, const int& pos,
        int& total_count, int& total_count_fw, int& total_count_rv,
        int& alt_count, int& alt_count_fw, int& alt_count_rv,
        char& ref, const std::string& alt_base,
        std::ofstream& both_ofs, std::ofstream& fwd_ofs, std::ofstream& rev_ofs);

};

#endif // VARIANT_CALLER_H