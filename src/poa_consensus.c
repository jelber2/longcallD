#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "call_var_main.h"
#include "align.h"
#include "poa_consensus.h"

// abPOA used an adaptive band (wb = -1); poa-consensus's default (50 bp floor
// plus adaptive widening) is the equivalent.
#define LONGCALLD_POA_BAND_WIDTH 50

static size_t *alloc_lens(int n, const int *lens) {
    size_t *out = (size_t*)malloc((size_t)n * sizeof(size_t));
    if (out == NULL) return NULL;
    for (int i = 0; i < n; ++i) out[i] = (size_t)lens[i];
    return out;
}

// Copy one POA consensus + its read cluster into plain C allocations.  Returns
// 0 on success.
static int copy_consensus(const poa_consensus_output_t *out, int ci,
                          const int *read_ids,
                          int *cons_len, uint8_t **cons_seq,
                          int *clu_n_seqs, int **clu_read_ids) {
    size_t clen = out->cons_lens[ci];
    uint8_t *seq = (uint8_t*)malloc(clen * sizeof(uint8_t));
    if (clen > 0 && seq == NULL) return -1;
    if (clen > 0) memcpy(seq, out->cons_seqs[ci], clen);
    *cons_len = (int)clen;
    *cons_seq = seq;

    if (clu_n_seqs != NULL && clu_read_ids != NULL) {
        size_t n = out->cluster_n_reads[ci];
        int *ids = (int*)malloc((n > 0 ? n : 1) * sizeof(int));
        if (ids == NULL) { free(seq); return -1; }
        for (size_t j = 0; j < n; ++j) {
            size_t ri = out->cluster_read_indices[ci][j];
            ids[j] = read_ids != NULL ? read_ids[ri] : (int)ri;
        }
        *clu_n_seqs = (int)n;
        *clu_read_ids = ids;
    }
    return 0;
}

int poa_consensus_aln_cons(const call_var_opt_t *opt,
                           int n_reads, int *read_ids,
                           uint8_t **read_seqs, int *read_lens,
                           int max_n_cons,
                           int *cons_lens, uint8_t **cons_seqs,
                           int *clu_n_seqs, int **clu_read_ids) {
    if (n_reads <= 0) return 0;
    if (max_n_cons < 1) max_n_cons = 1;

    size_t *lens = alloc_lens(n_reads, read_lens);
    if (lens == NULL) return 0;

    poa_consensus_output_t out;
    memset(&out, 0, sizeof(out));
    int ret = poa_consensus_run((size_t)n_reads,
                                (const uint8_t *const *)read_seqs, lens,
                                1 /* multi-allele */, 1 /* min_reads */,
                                opt->match, -opt->mismatch,
                                -opt->gap_open1, -opt->gap_ext1,
                                LONGCALLD_POA_BAND_WIDTH, opt->min_af,
                                (size_t)max_n_cons, &out);
    free(lens);

    if (ret != 0) {
        if (out.error != NULL)
            fprintf(stderr, "poa-consensus: %s\n", out.error);
        poa_consensus_output_free(&out);
        return 0;
    }

    int n_cons = (int)out.n_cons;
    if (n_cons > max_n_cons) n_cons = max_n_cons;
    for (int i = 0; i < n_cons; ++i) {
        if (copy_consensus(&out, i, read_ids,
                           cons_lens + i, cons_seqs + i,
                           clu_n_seqs != NULL ? clu_n_seqs + i : NULL,
                           clu_read_ids != NULL ? clu_read_ids + i : NULL) != 0) {
            n_cons = i;
            break;
        }
    }
    poa_consensus_output_free(&out);
    return n_cons;
}

int poa_consensus_partial_aln_cons(const call_var_opt_t *opt,
                                   int sampling_reads,
                                   int n_reads, int *read_ids,
                                   uint8_t **read_seqs, int *read_lens,
                                   int *read_full_cover,
                                   int *cons_len, uint8_t **cons_seq,
                                   int *clu_n_seqs, int **clu_read_ids) {
    *cons_len = 0; *cons_seq = NULL;
    if (clu_n_seqs != NULL) *clu_n_seqs = 0;
    if (clu_read_ids != NULL) *clu_read_ids = NULL;
    if (n_reads <= 0) return 0;

    uint8_t **trim_seqs = (uint8_t**)malloc((size_t)n_reads * sizeof(uint8_t*));
    int *trim_lens = (int*)malloc((size_t)n_reads * sizeof(int));
    int *trim_ids = (int*)malloc((size_t)n_reads * sizeof(int));
    if (trim_seqs == NULL || trim_lens == NULL || trim_ids == NULL) {
        free(trim_seqs); free(trim_lens); free(trim_ids);
        return 0;
    }

    int n_valid = 0;
    // `collect_partial_aln_beg_end` needs a fully spanning target read.  If
    // read 0 only covers one side, feed the reads untrimmed instead of
    // aborting (poa-consensus is semi-global and copes with the overhangs).
    int can_trim = LONGCALLD_NOISY_IS_BOTH_COVER(read_full_cover[0]) != 0;
    for (int i = 0; i < n_reads; ++i) {
        const uint8_t *seq = read_seqs[i];
        int len = read_lens[i];
        if (len <= 0) continue;
        if (i != 0 && can_trim) {
            int ref_beg, ref_end, read_beg, read_end;
            // Same pre-processing as the old abPOA subgraph path: keep only
            // the part of each read that aligns to read 0.
            if (collect_partial_aln_beg_end(opt, sampling_reads,
                                            read_seqs[0], read_lens[0], read_full_cover[0],
                                            read_seqs[i], read_lens[i], read_full_cover[i],
                                            &ref_beg, &ref_end, &read_beg, &read_end) == 0) {
                continue;
            }
            int seq_beg_cut = read_beg - 1;
            int seq_end_cut = read_lens[i] - read_end;
            seq = read_seqs[i] + seq_beg_cut;
            len = read_lens[i] - seq_beg_cut - seq_end_cut;
            if (len <= 0) continue;
        }
        trim_seqs[n_valid] = (uint8_t*)seq;
        trim_lens[n_valid] = len;
        trim_ids[n_valid] = read_ids != NULL ? read_ids[i] : i;
        ++n_valid;
    }

    if (n_valid == 0) {
        free(trim_seqs); free(trim_lens); free(trim_ids);
        return 0;
    }

    size_t *lens = alloc_lens(n_valid, trim_lens);
    poa_consensus_output_t out;
    memset(&out, 0, sizeof(out));
    int ret = 0;
    if (lens != NULL) {
        ret = poa_consensus_run((size_t)n_valid,
                                (const uint8_t *const *)trim_seqs, lens,
                                0 /* one consensus */, 1 /* min_reads */,
                                opt->match, -opt->mismatch,
                                -opt->gap_open1, -opt->gap_ext1,
                                LONGCALLD_POA_BAND_WIDTH, opt->min_af,
                                1, &out);
    }
    free(lens);
    free(trim_seqs); free(trim_lens); free(trim_ids);

    if (ret != 0 || out.n_cons == 0) {
        if (out.error != NULL)
            fprintf(stderr, "poa-consensus: %s\n", out.error);
        poa_consensus_output_free(&out);
        return 0;
    }

    if (copy_consensus(&out, 0, NULL, cons_len, cons_seq, NULL, NULL) != 0) {
        poa_consensus_output_free(&out);
        return 0;
    }
    // A haplotype is one allele: every read in this haplotype belongs to the
    // same cluster, including reads that were left out of the POA build.
    if (clu_n_seqs != NULL && clu_read_ids != NULL) {
        int *ids = (int*)malloc((size_t)n_reads * sizeof(int));
        if (ids == NULL) {
            free(*cons_seq); *cons_seq = NULL; *cons_len = 0;
            poa_consensus_output_free(&out);
            return 0;
        }
        for (int i = 0; i < n_reads; ++i) ids[i] = read_ids != NULL ? read_ids[i] : i;
        *clu_n_seqs = n_reads;
        *clu_read_ids = ids;
    }
    poa_consensus_output_free(&out);
    return 1;
}
