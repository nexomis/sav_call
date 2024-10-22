// variant_caller.hpp
#ifndef VARIANT_CALLER_H
#define VARIANT_CALLER_H

#include <string>
#include <map>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>

class VariantCaller {
public:
    VariantCaller();
    ~VariantCaller();
    bool parse_arguments(int argc, char **argv);
    void print_usage();
    bool run();

private:
    // Options
    bool count_read_extremities;
    std::string base_csv_file;
    std::string indel_csv_file;
    std::string called_variant_file;
    int min_qual;
    bool is_r1_rev;
    bool is_r2_rev;
    std::string bam_input;
    std::string fasta_reference;
    float min_freq;
    std::string call_strand;
    int min_count;

    // Data structures
    faidx_t *fai;
    samFile *in;
    bam_hdr_t *header;

    void process_pileup();
    static int pileup_callback(uint32_t tid, hts_pos_t pos, int n, const bam_pileup1_t *plp, void *data);
    void update_counts(uint32_t tid, hts_pos_t pos, int n, const bam_pileup1_t *plp);

    // Counting structures
    struct BaseCounts {
        int depth;
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
    std::string tid_to_name(int32_t tid);

    // Output methods
    void write_base_csv();
    void write_indel_csv();
    void call_variants();

};

#endif // VARIANT_CALLER_H