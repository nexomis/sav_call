// variant_caller.cpp
#include "VariantCaller.hpp"
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cctype>
#include <cstring>

VariantCaller::VariantCaller()
    : help(false),
      keep_read_extremities(false),
      min_qual(0),
      is_r1_rev(false),
      is_r2_rev(true),
      min_freq_fw(0.02),
      min_freq_rv(0.02),
      min_freq_both(0.02),
      call_strand("both"),
      min_count_fw(20),
      min_count_rv(20),
      min_count_both(20),
      min_alt_count_fw(10),
      min_alt_count_rv(10),
      min_alt_count_both(10),
      skip_secondary(true),
      skip_duplicate(true),
      max_n_pileup(1000000),
      fai(nullptr),
      in(nullptr),
      header(nullptr)
{
}

VariantCaller::~VariantCaller() {
    if (header) bam_hdr_destroy(header);
    if (in) sam_close(in);
    if (fai) fai_destroy(fai);
}

void VariantCaller::print_usage() {
    std::cout << "Usage: sav_call [options] --bam <input.bam> --reference <reference.fasta>       " << std::endl
              << "                                                                                " << std::endl
              << "  Version 0.2.0                                                                 " << std::endl
              << "                                                                                " << std::endl
              << "Required arguments:                                                             " << std::endl
              << "  --prefix-out <string>                Prefix for output files with called      " << std::endl
              << "                                         variants. 3 files for forward, reverse " << std::endl
              << "                                         and both strands, respectively:        " << std::endl
              << "                                         - fwd: <prefix>.fwd.csv                " << std::endl
              << "                                         - rev: <prefix>.rev.csv                " << std::endl
              << "                                         - both: <prefix>.both.csv              " << std::endl
              << "  --bam <input.bam>                    Input BAM file                           " << std::endl
              << "  --reference <reference.fasta>        Reference FASTA used for mapping         " << std::endl
              << "                                                                                " << std::endl
              << "Optional arguments:                                                             " << std::endl
              << "  --base-csv <file>                    Base counts and depths for both forward  " << std::endl
              << "                                         and reverse strand with all pileup     " << std::endl
              << "                                         position to CSV file. Position without " << std::endl
              << "                                         reads are not reported. Gaps in reads  " << std::endl
              << "                                         are accounted in depths                " << std::endl
              << "  --indel-csv <file>                   Output indel counts to CSV file          " << std::endl
              << "  --keep-read-extremities              Keep read extremities (default: false)   " << std::endl
              << "  --keep-duplicate                     Keep duplicate reads (default: false)    " << std::endl
              << "  --keep-secondary                     Keep secondary mapping (default: false)  " << std::endl
              << "  --min-qual <int>                     Minimum base quality to count            " << std::endl
              << "                                         (default: 0)                           " << std::endl
              << "  --R1-strand <forward|reverse>        Set R1 strand (default: forward)         " << std::endl
              << "  --R2-strand <forward|reverse>        Set R2 strand (default: reverse)         " << std::endl
              << "  --bam <input.bam>                    Input BAM file (required)                " << std::endl
              << "  --reference <reference.fasta>        Reference FASTA file (required)          " << std::endl
              << "  --min-freq <float;float;float>       Minimum frequency to call a variant      " << std::endl
              << "                                         (default: 0.02;0.02;0.02)              " << std::endl
              << "                                         fwd;rev;both                           " << std::endl
              << "  --min-count <int;int;int>              Minimum counts to call a variant       " << std::endl
              << "                                         (default: 20;20;20)                    " << std::endl
              << "                                         fwd;rev;both                           " << std::endl
              << "  --min-alt-count <int;int;int>        Minimum alternative counts to call       " << std::endl
              << "                                         (default: 10;10;10)                    " << std::endl
              << "                                         fwd;rev;both                           " << std::endl
              << "  --max-n-pileup <int>                 Maximum reads in pileup                  " << std::endl
              << "                                         (default: 1000000)                     " << std::endl;
}

bool VariantCaller::parse_arguments(int argc, char **argv) {
    static struct option long_options[] = {
        {"help", no_argument, 0, 0},
        {"keep-read-extremities", no_argument, 0, 0},
        {"skip-duplicate", no_argument, 0, 0},
        {"skip-secondary", no_argument, 0, 0},
        {"base-csv", required_argument, 0, 0},
        {"indel-csv", required_argument, 0, 0},
        {"prefix-out", required_argument, 0, 0},
        {"min-qual", required_argument, 0, 0},
        {"R1-strand", required_argument, 0, 0},
        {"R2-strand", required_argument, 0, 0},
        {"bam", required_argument, 0, 0},
        {"reference", required_argument, 0, 0},
        {"min-freq", required_argument, 0, 0},
        {"call-strand", required_argument, 0, 0},
        {"min-count", required_argument, 0, 0},
        {"min-alt-count", required_argument, 0, 0},
        {"max-n-pileup", required_argument, 0, 0},
        {0, 0, 0, 0}
    };

    std::string R1_strand = "forward";
    std::string R2_strand = "reverse";

    int option_index = 0;
    opterr = 0;  // Suppress getopt error messages
    while (true) {
        int c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1)  // End of options
            break;
        if (c != 0) {
            std::cerr << "Unknown option: " << c << std::endl;
            return false;
        }
        std::string opt_name = long_options[option_index].name;
        if (opt_name == "keep-read-extremities") {
            keep_read_extremities = true;
        } else if (opt_name == "base-csv") {
            base_csv_file = optarg;
        } else if (opt_name == "indel-csv") {
            indel_csv_file = optarg;
        } else if (opt_name == "prefix-out") {
            prefix_out = optarg;
        } else if (opt_name == "min-qual") {
            min_qual = std::stoi(optarg);
        } else if (opt_name == "R1-strand") {
            R1_strand = optarg;
        } else if (opt_name == "R2-strand") {
            R2_strand = optarg;
        } else if (opt_name == "bam") {
            bam_input = optarg;
        } else if (opt_name == "reference") {
            fasta_reference = optarg;
        } else if (opt_name == "min-freq") {
            split_three(optarg, min_freq_fw, min_freq_rv, min_freq_both);
        } else if (opt_name == "min-count") {
            split_three(optarg, min_count_fw, min_count_rv, min_count_both);
        } else if (opt_name == "min-alt-count") {
            split_three(optarg, min_alt_count_fw, min_alt_count_rv, min_alt_count_both);
        } else if (opt_name == "max-n-pileup") {
            max_n_pileup = std::stoi(optarg);
        } else if (opt_name == "keep-secondary") {
            skip_secondary = false;
        } else if (opt_name == "help") {
            help = true;
        } else if (opt_name == "keep-duplicate") {
            skip_duplicate = false;
        } else {
            std::cerr << "Unknown option: " << opt_name << std::endl;
            return false;
        }
    }

    if ((bam_input.empty() || fasta_reference.empty()) && ! help) {
        std::cerr << "Error: BAM input and reference FASTA are required." << std::endl;
        return false;
    }

    if (R1_strand != "forward" && R1_strand != "reverse") {
        std::cerr << "Error: Invalid R1-strand option." << std::endl;
        return false;
    }

    if (R2_strand != "forward" && R2_strand != "reverse") {
        std::cerr << "Error: Invalid R2-strand option." << std::endl;
        return false;
    }

    if (call_strand != "forward" && call_strand != "reverse" && call_strand != "both") {
        std::cerr << "Error: Invalid call-strand option." << std::endl;
        return false;
    }

    if (R1_strand == "forward") {
        is_r1_rev = false;
    } else {
        is_r1_rev = true;
    }
    if (R2_strand == "forward") {
        is_r2_rev = false;
    } else {
        is_r2_rev = true;
    }

    return true;
}

bool VariantCaller::run() {
    // Open BAM file
    in = sam_open(bam_input.c_str(), "r");
    if (!in) {
        std::cerr << "Failed to open BAM file: " << bam_input << std::endl;
        return false;
    }

    // Read header
    header = sam_hdr_read(in);
    if (!header) {
        std::cerr << "Failed to read BAM header" << std::endl;
        return false;
    }

    // Load reference
    fai = fai_load(fasta_reference.c_str());
    if (!fai) {
        std::cerr << "Failed to load reference FASTA: " << fasta_reference << std::endl;
        return false;
    }

    // Process pileup
    process_pileup();

    // Write outputs
    if (!base_csv_file.empty())
        write_base_csv();
    if (!indel_csv_file.empty())
        write_indel_csv();
    call_variants();

    return true;
}

void VariantCaller::process_pileup() {
    bam_plp_t iter;
    bam_plp_auto_f fill_func = [](void *data, bam1_t *b) -> int {
        VariantCaller *vc = static_cast<VariantCaller*>(data);
        return sam_read1(vc->in, vc->header, b);
    };

    iter = bam_plp_init(fill_func, this);
    bam_plp_set_maxcnt(iter, max_n_pileup);

    int tid;
    int pos;
    const bam_pileup1_t *plp;
    std::map<int,std::pair<int,int>> del_depth; // to record depth of deletion
    int n_plp;

    while ((plp = bam_plp_auto(iter, &tid, &pos, &n_plp)) != nullptr) {
        update_counts(tid, pos, n_plp, plp, del_depth);
    }

    bam_plp_destroy(iter);
}

void VariantCaller::update_counts(uint32_t tid, int pos, int n_plp, const bam_pileup1_t *plp, std::map<int, std::pair<int, int> >& del_depth) {
    std::string chrom = tid_to_name(tid);
    std::string key = chrom + ":" + std::to_string(pos + 1);
    BaseCounts &counts = base_counts_map[key];

    counts.depth_fw = 0;
    counts.depth_rv = 0;

    if (del_depth.find(pos) != del_depth.end()) {
        counts.depth_fw += del_depth[pos].first;
        counts.depth_rv += del_depth[pos].second;
        del_depth.erase(pos);
    }

    for (int i = 0; i < n_plp; ++i) {
        const bam_pileup1_t *p = &plp[i];
        const bam1_t *b = p->b;

        // Skip deletions and reference skips
        if (p->is_del || p->is_refskip)
            continue;

        if (is_duplicate(b) && skip_duplicate)
            continue;

        if (is_secondary(b) && skip_secondary)
            continue;

        if (! keep_read_extremities){
            if (p->is_head || p->is_tail){
                continue;
            }
        }
            

        // Get base quality
        uint8_t *q = bam_get_qual(b);
        if (q[p->qpos] < min_qual)
            continue;

        // Get base
        uint8_t *s = bam_get_seq(b);
        int base = bam_seqi(s, p->qpos);
        char base_char = seq_nt16_str[base];

        // Determine read labels and strands
        bool is_forward_strand = ! bam_is_rev(b);
        
        if ((is_read1(b) && is_r1_rev) || (is_read2(b) && is_r2_rev)) {
            is_forward_strand = ! is_forward_strand;
        }

        if (is_forward_strand)
            counts.depth_fw++;
        else
            counts.depth_rv++;

        switch (base_char) {
            case 'A':
                if (is_forward_strand)
                    counts.A_fwd++;
                else
                    counts.A_rev++;
                break;
            case 'C':
                if (is_forward_strand)
                    counts.C_fwd++;
                else
                    counts.C_rev++;
                break;
            case 'G':
                if (is_forward_strand)
                    counts.G_fwd++;
                else
                    counts.G_rev++;
                break;
            case 'T':
                if (is_forward_strand)
                    counts.T_fwd++;
                else
                    counts.T_rev++;
                break;
            default:
                if (is_forward_strand)
                    counts.N_fwd++;
                else
                    counts.N_rev++;
                break;
        }

        // Process indels
        if (p->indel != 0) {
            // Indels
            int indel_len = p->indel;
            std::string seq = "";
            if (indel_len > 0) {
                // Insertion
                for (int j = 1; j <= indel_len; ++j) {
                    int indel_base = bam_seqi(s, p->qpos + j);
                    seq += "=ACMGRSVTWYHKDBN"[indel_base];
                }
            } else {
                // Deletion
                int ref_len;
                char *ref_seq = faidx_fetch_seq(fai, chrom.c_str(), pos + 1, pos - indel_len, &ref_len);
                if (ref_seq) {
                    seq = std::string(ref_seq, ref_seq + ref_len);
                    free(ref_seq);
                }
                for (int i = 0; i < - indel_len; i++) {
                    int new_pos = pos + i + 1;
                    if (del_depth.find(new_pos) == del_depth.end()) {
                        del_depth[new_pos] = std::make_pair(0,0);
                    }
                    if (is_forward_strand)
                        del_depth[new_pos].first++;
                    else
                        del_depth[new_pos].second++;
                }
            }

            if (indel_len > 0)
                seq = "+" + seq;
            else
                seq = "-" + seq;

            // Update indel counts
            IndelCounts &indel_counts = indel_counts_map[key];
            if (indel_counts.fwd_counts.find(seq) == indel_counts.fwd_counts.end()) {
                indel_counts.fwd_counts[seq] = 0;
                indel_counts.rev_counts[seq] = 0;
            }
            if (is_forward_strand)
                indel_counts.fwd_counts[seq]++;
            else
                indel_counts.rev_counts[seq]++;

        }
    }
}

bool VariantCaller::is_read1(const bam1_t *b) {
    return (b->core.flag & BAM_FREAD1) != 0;
}

bool VariantCaller::is_duplicate(const bam1_t *b) {
    return (b->core.flag & BAM_FDUP) != 0;
}

bool VariantCaller::is_secondary(const bam1_t *b) {
    return (b->core.flag & BAM_FSECONDARY) != 0;
}

bool VariantCaller::is_read2(const bam1_t *b) {
    return (b->core.flag & BAM_FREAD2) != 0;
}

std::string VariantCaller::tid_to_name(int32_t tid) {
    return std::string(header->target_name[tid]);
}

void VariantCaller::write_base_csv() {
    std::ofstream ofs(base_csv_file);
    if (!ofs) {
        std::cerr << "Error opening base CSV file for writing: " << base_csv_file << std::endl;
        return;
    }
    // Output header
    ofs << "region;pos;ref;depth;depth_fw;depth_rv;A;a;T;t;C;c;G;g;N;n\n";

    for (int tid = 0; tid < header->n_targets; ++tid) {
        const char* chrom_name = header->target_name[tid];
        int64_t chrom_length = header->target_len[tid];
        for (int64_t pos = 1; pos <= chrom_length; ++pos) {
            const std::string key = std::string(chrom_name) + ":" + std::to_string(pos);
            const BaseCounts &counts = base_counts_map[key];
            int ref_len;
            char *ref_base = faidx_fetch_seq(fai, chrom_name, pos - 1, pos - 1, &ref_len);
            char ref = (ref_base && ref_len == 1) ? ref_base[0] : 'N';
            if (ref_base) free(ref_base);

            ofs << chrom_name << ';' << pos << ';' << ref << ';' 
                << counts.depth_fw + counts.depth_rv << ';'
                << counts.depth_fw << ';' << counts.depth_rv << ';'
                << counts.A_fwd << ';' << counts.A_rev << ';'
                << counts.T_fwd << ';' << counts.T_rev << ';'
                << counts.C_fwd << ';' << counts.C_rev << ';'
                << counts.G_fwd << ';' << counts.G_rev << ';'
                << counts.N_fwd << ';' << counts.N_rev << '\n';
        }
    }

    ofs.close();
}

void VariantCaller::write_indel_csv() {
    std::ofstream ofs(indel_csv_file);
    if (!ofs) {
        std::cerr << "Error opening indel CSV file for writing: " << indel_csv_file << std::endl;
        return;
    }

    // Output header
    ofs << "region;pos;ref;alt;forward;reverse\n";

    for (const auto &entry : indel_counts_map) {
        const std::string &key = entry.first;
        const IndelCounts &counts = entry.second;
        // Split key into chrom and pos
        size_t delim_pos = key.find(':');
        std::string chrom = key.substr(0, delim_pos);
        int position = std::stoi(key.substr(delim_pos + 1));

        // Get reference base
        int ref_len;
        char *ref_base = faidx_fetch_seq(fai, chrom.c_str(), position - 1, position - 1, &ref_len);
        char ref = (ref_base && ref_len == 1) ? ref_base[0] : 'N';
        if (ref_base) free(ref_base);

        // Output indels
        for (const auto &indel_entry : counts.fwd_counts) {
            const std::string &alt = indel_entry.first;
            int fwd_count = indel_entry.second;
            int rev_count = counts.rev_counts.count(alt) ? counts.rev_counts.at(alt) : 0;

            ofs << chrom << ';' << position << ';' << ref << ';' << alt << ';'
                << fwd_count << ';' << rev_count << '\n';
        }
    }

    ofs.close();
}

void VariantCaller::call_variants() {
    std::string both_file = prefix_out + ".both.csv";
    std::string rev_file = prefix_out + ".rev.csv";
    std::string fwd_file = prefix_out + ".fwd.csv";

    std::ofstream both_ofs(both_file);
    std::ofstream rev_ofs(rev_file);
    std::ofstream fwd_ofs(fwd_file);

    if (!(both_ofs || rev_ofs || fwd_ofs)) {
        std::cerr << "Error opening output file for writing: " << std::endl;
        return;
    }

    // Output header
    both_ofs << "region;pos;ref;alt;depth;freq" << std::endl;
    rev_ofs << "region;pos;ref;alt;depth;freq" << std::endl;
    fwd_ofs << "region;pos;ref;alt;depth;freq" << std::endl;

    for (int tid = 0; tid < header->n_targets; ++tid) {
        const char* chrom_name = header->target_name[tid];
        int64_t chrom_length = header->target_len[tid];
        for (int64_t pos = 0; pos < chrom_length; ++pos) {
            const std::string key = std::string(chrom_name) + ":" + std::to_string(pos);
            const BaseCounts &counts = base_counts_map[key];

            // Get reference base
            int ref_len;
            char *ref_base = faidx_fetch_seq(fai, chrom_name, pos - 1, pos - 1, &ref_len);
            char ref = (ref_base && ref_len == 1) ? ref_base[0] : 'N';
            if (ref_base) free(ref_base);

            int total_count_fw = counts.depth_fw;
            int total_count_rv = counts.depth_rv;
            int total_count = total_count_fw + total_count_rv;

            // For each possible alt base
            char bases[] = {'A', 'C', 'G', 'T', 'N'};
            for (char alt_base : bases) {
                if (alt_base == std::toupper(ref))
                    continue;  // Skip reference base

                int alt_count_fw = 0, alt_count_rv = 0;

                switch (alt_base) {
                    case 'A':
                        alt_count_fw = counts.A_fwd;
                        alt_count_rv = counts.A_rev;
                        break;
                    case 'C':
                        alt_count_fw = counts.C_fwd;
                        alt_count_rv = counts.C_rev;
                        break;
                    case 'G':
                        alt_count_fw = counts.G_fwd;
                        alt_count_rv = counts.G_rev;
                        break;
                    case 'T':
                        alt_count_fw = counts.T_fwd;
                        alt_count_rv = counts.T_rev;
                        break;
                    case 'N':
                        alt_count_fw = counts.N_fwd;
                        alt_count_rv = counts.N_rev;
                        break;
                }

                int alt_count = alt_count_fw + alt_count_rv;

                std::string alt_str = "";
                alt_str += alt_base;

                make_a_call(chrom_name, pos,
                            total_count, total_count_fw, total_count_rv,
                            alt_count, alt_count_fw, alt_count_rv,
                            ref, alt_str,
                            both_ofs, fwd_ofs, rev_ofs); 

            }

            if (indel_counts_map.find(key) != indel_counts_map.end()) {
                const IndelCounts &indel_counts =  indel_counts_map[key];
                for (const auto &indel_entry : indel_counts.fwd_counts) {
                    const std::string &alt_base = indel_entry.first;
                    int alt_count_fw = indel_entry.second;
                    int alt_count_rv = indel_counts.rev_counts.count(alt_base) ? indel_counts.rev_counts.at(alt_base) : 0;
                    int alt_count = alt_count_fw + alt_count_rv;

                    make_a_call(chrom_name, pos,
                                total_count, total_count_fw, total_count_rv,
                                alt_count, alt_count_fw, alt_count_rv,
                                ref, alt_base,
                                both_ofs, fwd_ofs, rev_ofs);

                }
            }

        }
    }

    both_ofs.close();
    rev_ofs.close();
    fwd_ofs.close();
}

void VariantCaller::make_a_call(
    const std::string& region, const int& pos,
    int& total_count, int& total_count_fw, int& total_count_rv,
    int& alt_count, int& alt_count_fw, int& alt_count_rv,
    char& ref, const std::string& alt_base,
    std::ofstream& both_ofs, std::ofstream& fwd_ofs, std::ofstream& rev_ofs) {
    
    bool pass_count = (alt_count >= min_alt_count_both) && (total_count >= min_count_both);
    bool pass_count_fw = (alt_count_fw >= min_alt_count_fw) && (total_count_fw >= min_count_fw);
    bool pass_count_rv = (alt_count_rv >= min_alt_count_rv) && (total_count_rv >= min_count_rv);
    float freq = 0;
    float freq_fw = 0;
    float freq_rv = 0;
    if (total_count > 0)
        freq = static_cast<float>(alt_count) / total_count;
    if (total_count_fw > 0)
        freq_fw = static_cast<float>(alt_count_fw) / total_count_fw;
    if (total_count_rv > 0)
        freq_rv = static_cast<float>(alt_count_rv) / total_count_rv;
    bool pass_freq = (freq >= min_freq_both);
    bool pass_freq_rv = (freq_rv >= min_freq_fw);
    bool pass_freq_fw = (freq_fw >= min_freq_rv);

    if (pass_freq && pass_count){
        both_ofs << region << ';' << pos << ';' << ref << ';' << alt_base << ';'
        << total_count << ';' << freq << std::endl;
    }
    if (pass_freq_rv && pass_count_rv){
        rev_ofs << region << ';' << pos << ';' << ref << ';' << alt_base << ';'
        << total_count_rv << ';' << freq_rv << std::endl;
    }
    if (pass_freq_fw && pass_count_fw){
        fwd_ofs << region << ';' << pos << ';' << ref << ';' << alt_base << ';'
        << total_count_fw << ';' << freq_fw << std::endl;
    }
}

void VariantCaller::split_three(const std::string& input_str, 
                std::string& val1, 
                std::string& val2, 
                std::string& val3) 
{
    size_t pos1 = input_str.find(';');
    if (pos1 == std::string::npos) {
        throw std::invalid_argument("Input string does not contain enough delimiters.");
    }

    size_t pos2 = input_str.find(';', pos1 + 1);
    if (pos2 == std::string::npos) {
        throw std::invalid_argument("Input string does not contain enough delimiters.");
    }

    val1 = input_str.substr(0, pos1);
    val2 = input_str.substr(pos1 + 1, pos2 - pos1 - 1);
    val3 = input_str.substr(pos2 + 1);
}


void VariantCaller::split_three(const std::string& input_str, 
                               int& val1, 
                               int& val2, 
                               int& val3) 
{
    std::string s1, s2, s3;

    // Call the string-based split function
    try {
        split_three(input_str, s1, s2, s3);
    }
    catch (const std::invalid_argument& e) {
        throw; // Re-throw the exception for the caller to handle
    }

    // Convert strings to ints
    try {
        val1 = std::stoi(s1);
    }
    catch (const std::exception& e) {
        throw std::invalid_argument("Conversion to int failed for val1: " + s1);
    }

    try {
        val2 = std::stoi(s2);
    }
    catch (const std::exception& e) {
        throw std::invalid_argument("Conversion to int failed for val2: " + s2);
    }

    try {
        val3 = std::stoi(s3);
    }
    catch (const std::exception& e) {
        throw std::invalid_argument("Conversion to int failed for val3: " + s3);
    }
}

void VariantCaller::split_three(const std::string& input_str, 
                               float& val1, 
                               float& val2, 
                               float& val3) 
{
    std::string s1, s2, s3;

    // Call the string-based split function
    try {
        split_three(input_str, s1, s2, s3);
    }
    catch (const std::invalid_argument& e) {
        throw; // Re-throw the exception for the caller to handle
    }

    // Convert strings to floats
    try {
        val1 = std::stof(s1);
    }
    catch (const std::exception& e) {
        throw std::invalid_argument("Conversion to float failed for val1: " + s1);
    }

    try {
        val2 = std::stof(s2);
    }
    catch (const std::exception& e) {
        throw std::invalid_argument("Conversion to float failed for val2: " + s2);
    }

    try {
        val3 = std::stof(s3);
    }
    catch (const std::exception& e) {
        throw std::invalid_argument("Conversion to float failed for val3: " + s3);
    }
}
