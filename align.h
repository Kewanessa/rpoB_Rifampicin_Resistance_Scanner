/**
 * align.h
 *
 * Header file defining the interface and data structures for pairwise sequence
 * alignment using the Needleman-Wunsch algorithm (semi-global variant)
 * for the rpoB rifampicin resistance scanner.
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#ifndef ALIGN_H
#define ALIGN_H

// Scoring constants for sequence alignment
#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_SCORE -2

/**
 * Struct representing the result of an alignment between two sequences.
 *
 * ref_aligned: Aligned reference sequence (may include gap characters '-')
 * qry_aligned: Aligned query sequence (may include gap characters '-')
 * length:      Length of the aligned strings (both strings have identical length)
 * score:       Alignment score computed by Needleman-Wunsch dynamic programming
 */
typedef struct
{
    char *ref_aligned;
    char *qry_aligned;
    int length;
    int score;
} AlignResult;

/**
 * Performs semi-global Needleman-Wunsch alignment between a reference sequence
 * and a query sequence.
 *
 * Semi-global: leading and trailing gaps in either sequence are free, so a
 * short amplicon can sit anywhere along the gene and flanking DNA beyond the
 * gene ends costs nothing. Internal gaps cost GAP_SCORE per base. Indels are
 * kept contiguous and placed at their left-most equivalent position.
 * Refuses (returns NULL) if the DP grid would be unreasonably large.
 *
 * @param ref      Reference sequence string
 * @param ref_len  Length of the reference sequence
 * @param qry      Query sequence string
 * @param qry_len  Length of the query sequence
 * @return         Pointer to dynamically allocated AlignResult, or NULL on error.
 *                 Caller is responsible for freeing the result using free_align_result.
 */
AlignResult *align_sequences(const char *ref, int ref_len, const char *qry, int qry_len);

/**
 * Frees memory dynamically allocated for an AlignResult and its aligned strings.
 *
 * @param result   Pointer to AlignResult to free (safe to call with NULL).
 */
void free_align_result(AlignResult *result);

#endif // ALIGN_H
