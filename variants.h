/**
 * variants.h
 *
 * Header file defining data structures and function prototypes for variant calling
 * in the rpoB rifampicin resistance scanner.
 *
 * Identifies codon substitutions (SNPs), insertions, deletions, frameshifts and
 * codons with ambiguous base calls, mapping them to codon coordinates, amino acid
 * substitutions, and flagging mutations within the Rifampicin Resistance-Determining
 * Region (RRDR).
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */
#ifndef VARIANTS_H
#define VARIANTS_H

#include "mutdb.h"

// Variant types
#define VAR_SNP        0   // One or more bases changed inside a single codon
#define VAR_INSERTION  1   // Patient has extra bases not in the reference
#define VAR_DELETION   2   // Reference bases missing from the patient
#define VAR_AMBIGUOUS  3   // Codon contains N or another IUPAC ambiguity code

// M. tuberculosis rpoB RRDR boundaries (codons 426-452)
#define RRDR_START_CODON 426
#define RRDR_END_CODON   452

// Longest indel sequence we store for display (longer ones are truncated)
#define MAX_INDEL_SEQ 15

typedef struct {
    int codon_num;           // 1-based codon coordinate in reference
    char ref_codon[4];       // Reference codon (e.g., "TCG")
    char alt_codon[4];       // Patient codon (e.g., "TTG")
    char ref_aa;             // Reference amino acid
    char alt_aa;             // Mutated amino acid ('?' for ambiguous codons and indels)
    int nuc_pos;             // 1-based reference position: first changed/deleted base,
                             // or the base just before an insertion
    int indel_len;           // Number of inserted or deleted bases (indels only)
    char indel_seq[MAX_INDEL_SEQ + 1]; // Inserted or deleted bases (indels only)
    int is_frameshift;       // 1 if an indel's length is not a multiple of 3
    int is_synonymous;       // 1 if silent change, 0 otherwise
    int in_rrdr;             // 1 if any affected codon lies within codons 426-452
    int type;                // VAR_SNP, VAR_INSERTION, VAR_DELETION or VAR_AMBIGUOUS
    MutEntry *db_hit;        // Cached pointer to mutation database entry (or NULL)
} Variant;

typedef struct {
    Variant *variants;
    int count;
    int capacity;
    int cov_start_pos;       // 1-based start of patient coverage on reference (0 = none)
    int cov_end_pos;         // 1-based end of patient coverage on reference (0 = none)
    int unaligned_bases;     // Patient bases outside the aligned region (flanks, ignored)
} VariantList;

/**
 * Walk an alignment and list every difference between patient and reference,
 * in reference coordinates. Looks each substitution up in the mutation database.
 * Returns NULL on error. The caller must free the result with free_variant_list().
 */
VariantList *call_variants(const char *ref_aligned, const char *qry_aligned,
                           int align_len, MutDB *db);

/**
 * Free a VariantList and its array of variants.
 */
void free_variant_list(VariantList *vl);

#endif // VARIANTS_H
