/**
 * report.h — Report generation
 *
 * Takes the variant list and mutation database lookup results
 * and prints a structured clinical-style report.
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#ifndef REPORT_H
#define REPORT_H

#include "variants.h"
#include "mutdb.h"

// Number of reference bases in the RRDR (codons 426-452 = 81 bp)
#define RRDR_LENGTH_BP ((RRDR_END_CODON - RRDR_START_CODON + 1) * 3)

// Everything the report needs to know besides the variants
typedef struct
{
    const char *sample_id;     // Name/ID of the patient sample
    const char *ref_id;        // Name/ID of the reference sequence
    double identity_pct;       // % identity over aligned, unambiguous bases
    int aligned_bases;         // Patient bases aligned to reference bases
    int patient_length;        // Total patient bases
    int reverse_complemented;  // 1 if the patient was read on the opposite strand
    int alignment_ok;          // 0 if the sample does not align well to rpoB
    int rrdr_resolved;         // RRDR bases (0-81) read as A, C, G or T
} ReportInfo;

/**
 * Print the full analysis report and the final verdict.
 */
void print_report(const ReportInfo *info, const VariantList *vl);

/**
 * Calculate alignment identity percentage.
 * Counts matching columns / columns where the patient has an unambiguous
 * base (A, C, G, T) aligned to a reference base, × 100.
 */
double calc_identity(const char *ref_aligned, const char *qry_aligned, int align_len);

/**
 * Count alignment columns where a patient base sits on a reference base.
 */
int count_aligned_bases(const char *ref_aligned, const char *qry_aligned, int align_len);

/**
 * Count how many of the 81 RRDR reference bases (codons 426-452) are
 * covered by an unambiguous patient base (A, C, G or T). Gaps, N and
 * other IUPAC codes do not count. 81 means the RRDR is fully resolved.
 */
int count_rrdr_resolved(const char *ref_aligned, const char *qry_aligned, int align_len);

#endif // REPORT_H
