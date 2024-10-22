// variant_caller.cpp
#include "VariantCaller.hpp"
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <cctype>
#include <cstring>

VariantCaller::VariantCaller()
    : count_read_extremities(false),
      min_qual(0),
      is_r1_rev(false),
      is_r2_rev(true),
      min_freq(0.0),
      call_strand("both"),
      min_count(0),
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
    std::cout << "Usage: variant_caller [options] --bam <input.bam> --reference <reference.fasta>\n"
              << "Options:\n"
              << "  --count-read-extremities        Count read extremities (default: false)\n"
              << "  --base-csv <file>               Output base counts to CSV file\n"
              << "  --indel-csv <file>              Output indel counts to CSV file\n"
              << "  --called-variant-csv <file>     Output called variants to CSV file\n"
              << "  --min-qual <int>                Minimum base quality to count\n"
              << "  --R1-strand <forward|reverse>   Set R1 strand (default: forward)\n"
              << "  --R2-strand <forward|reverse>   Set R2 strand (default: reverse)\n"
              << "  --bam <input.bam>               Input BAM file (required)\n"
              << "  --reference <reference.fasta>   Reference FASTA file (required)\n"
              << "  --min-freq <float>              Minimum frequency to call a variant\n"
              << "  --call-strand <forward|reverse|both>  Strand to apply thresholds (default: both)\n"
              << "  --min-count <int>               Minimum count to report in outputs\n";
}

bool VariantCaller::parse_arguments(int argc, char **argv) {
    static struct option long_options[] = {
        {"count-read-extremities", no_argument, 0, 0},
        {"base-csv", required_argument, 0, 0},
        {"indel-csv", required_argument, 0, 0},
        {"called-variant-csv", required_argument, 0, 0},
        {"min-qual", required_argument, 0, 0},
        {"R1-strand", required_argument, 0, 0},
        {"R2-strand", required_argument, 0, 0},
        {"bam", required_argument, 0, 0},
        {"reference", required_argument, 0, 0},
        {"min-freq", required_argument, 0, 0},
        {"call-strand", required_argument, 0, 0},
        {"min-count", required_argument, 0, 0},
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
        if (opt_name == "count-read-extremities") {
            count_read_extremities = true;
        } else if (opt_name == "base-csv") {
            base_csv_file = optarg;
        } else if (opt_name == "indel-csv") {
            indel_csv_file = optarg;
        } else if (opt_name == "called-variant-csv") {
            called_variant_file = optarg;
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
            min_freq = std::stof(optarg);
        } else if (opt_name == "call-strand") {
            call_strand = optarg;
        } else if (opt_name == "min-count") {
            min_count = std::stoi(optarg);
        } else {
            std::cerr << "Unknown option: " << opt_name << std::endl;
            return false;
        }
    }

    if (bam_input.empty() || fasta_reference.empty()) {
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
    if (!called_variant_file.empty())
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
    bam_plp_set_maxcnt(iter, 50000);  // Set a reasonable max depth to avoid excessive memory usage TODO set as option

    int tid;
    int pos;
    const bam_pileup1_t *plp;
    int n_plp;

    while ((plp = bam_plp_auto(iter, &tid, &pos, &n_plp)) != nullptr) {
        update_counts(tid, pos, n_plp, plp);
    }

    bam_plp_destroy(iter);
}

void VariantCaller::update_counts(uint32_t tid, hts_pos_t pos, int n_plp, const bam_pileup1_t *plp) {
    std::string chrom = tid_to_name(tid);
    std::string key = chrom + ":" + std::to_string(pos + 1);
    BaseCounts &counts = base_counts_map[key];
    counts.depth = n_plp;

    for (int i = 0; i < n_plp; ++i) {
        const bam_pileup1_t *p = &plp[i];
        const bam1_t *b = p->b;

        // Skip deletions and reference skips
        if (p->is_del || p->is_refskip){
            counts.depth--;
            continue;
        }

        if ((p->is_head || p->is_tail) && (!count_read_extremities)){
            counts.depth--;
            continue;
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
            case 'N':
                if (is_forward_strand)
                    counts.N_fwd++;
                else
                    counts.N_rev++;
                break;
            default:
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
            }

            // For reverse strand, reverse complement
            if (!is_forward_strand) {
                std::string revcomp_seq = "";
                for (auto it = seq.rbegin(); it != seq.rend(); ++it) {
                    char c = *it;
                   switch (std::toupper(c)) {
                        case 'A': revcomp_seq += 'T'; break;
                        case 'C': revcomp_seq += 'G'; break;
                        case 'G': revcomp_seq += 'C'; break;
                        case 'T': revcomp_seq += 'A'; break;
                        default: revcomp_seq += 'N';
                    }
                }
                seq = revcomp_seq;
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
    ofs << "region;pos;ref;depth;A;a;T;t;C;c;G;g;N;n\n";

    for (int tid = 0; tid < header->n_targets; ++tid) {
        const char* chrom_name = header->target_name[tid];
        int64_t chrom_length = header->target_len[tid];
        for (int64_t pos = 0; pos < chrom_length; ++pos) {
            const std::string key = std::string(chrom_name) + ":" + std::to_string(pos);
            const BaseCounts &counts = base_counts_map[key];
            int ref_len;
            char *ref_base = faidx_fetch_seq(fai, chrom_name, pos - 1, pos - 1, &ref_len);
            char ref = (ref_base && ref_len == 1) ? ref_base[0] : 'N';
            if (ref_base) free(ref_base);

            ofs << chrom_name << ';' << pos << ';' << ref << ';' << counts.depth << ';'
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

            if (fwd_count + rev_count >= min_count) {
                ofs << chrom << ';' << position << ';' << ref << ';' << alt << ';'
                    << fwd_count << ';' << rev_count << '\n';
            }
        }
    }

    ofs.close();
}

void VariantCaller::call_variants() {
    std::ofstream ofs(called_variant_file);
    if (!ofs) {
        std::cerr << "Error opening called variant CSV file for writing: " << called_variant_file << std::endl;
        return;
    }

    // Output header
    ofs << "region;pos;ref;alt;depth;freq;depth_fw;freq_fw;depth_rv;freq_rv\n";

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

            int total_count_fw = counts.A_fwd + counts.C_fwd + counts.G_fwd + 
                               counts.T_fwd + counts.N_fwd ;

            int total_count_rv = counts.A_rev + counts.C_rev + counts.G_rev + 
                               counts.T_rev + counts.N_rev ;

            int total_count = total_count_fw + total_count_rv;

            if (call_strand == "forward") {
                if (total_count_fw < min_count)
                    continue;
            } else if (call_strand == "reverse") {
                if (total_count_rv < min_count)
                    continue;
            } else {
                if (total_count < min_count)
                    continue;
            }

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
                if (alt_count == 0)
                    continue;

                float freq = 0;
                float freq_fw = 0;
                float freq_rv = 0;
                
                if (total_count > 0)
                    freq = static_cast<float>(alt_count) / total_count;
                if (total_count_fw > 0)
                    freq_fw = static_cast<float>(alt_count_fw) / total_count_fw;
                if (total_count_rv > 0)
                    freq_rv = static_cast<float>(alt_count_rv) / total_count_rv;

                bool pass_freq = (freq >= min_freq);

                // Apply call_strand filter
                if (call_strand == "forward")
                    pass_freq = (freq_fw >= min_freq);
                else if (call_strand == "reverse")
                    pass_freq = (freq_rv >= min_freq);

                if (!pass_freq)
                    continue;

                ofs << chrom_name << ';' << pos << ';' << ref << ';' << alt_base << ';'
                    << total_count << ';' << freq << ';'
                    << total_count_fw << ';' << freq_fw << ';'
                    << total_count_rv << ';' << freq_rv << '\n';
            }
        }
    }

    ofs.close();
}

