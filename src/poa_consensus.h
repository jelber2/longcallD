#ifndef LONGCALLD_POA_CONSENSUS_H
#define LONGCALLD_POA_CONSENSUS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct call_var_opt_t;

// Mirrors PoaConsensusOutput in poa-consensus/src/ffi.rs.
typedef struct poa_consensus_output_t {
    size_t n_cons;
    size_t *cons_lens;
    uint8_t **cons_seqs;
    size_t *cluster_n_reads;
    size_t **cluster_read_indices;
    char *error;
} poa_consensus_output_t;

int poa_consensus_run(size_t n_reads,
                      const uint8_t *const *read_seqs,
                      const size_t *read_lens,
                      int multi,
                      size_t min_reads,
                      int match_score,
                      int mismatch_score,
                      int gap_open,
                      int gap_extend,
                      size_t band_width,
                      double min_allele_freq,
                      size_t max_n_cons,
                      poa_consensus_output_t *out);
void poa_consensus_output_free(poa_consensus_output_t *out);

// Build up to max_n_cons consensus sequences from all reads.  On success
// returns n_cons and fills cons_lens/cons_seqs plus per-cluster read ids.
// read_ids are the caller's read identifiers; returned clu_read_ids are already
// translated through read_ids.
int poa_consensus_aln_cons(const struct call_var_opt_t *opt,
                           int n_reads, int *read_ids,
                           uint8_t **read_seqs, int *read_lens,
                           int max_n_cons,
                           int *cons_lens, uint8_t **cons_seqs,
                           int *clu_n_seqs, int **clu_read_ids);

// Haplotype-specific consensus.  Reads are first trimmed to the part that
// aligns to read_seqs[0] (the same partial-alignment pre-processing abPOA
// needed), then a single consensus is built.  Returns 1 and fills cons_len/
// cons_seq and a single cluster containing all n_reads, or 0 on failure.
int poa_consensus_partial_aln_cons(const struct call_var_opt_t *opt,
                                   int sampling_reads,
                                   int n_reads, int *read_ids,
                                   uint8_t **read_seqs, int *read_lens,
                                   int *read_full_cover,
                                   int *cons_len, uint8_t **cons_seq,
                                   int *clu_n_seqs, int **clu_read_ids);

#ifdef __cplusplus
}
#endif

#endif
