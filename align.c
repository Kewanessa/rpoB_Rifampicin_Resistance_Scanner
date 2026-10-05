/**
 * align.c — Needleman-Wunsch sequence alignment
 *
 * THE BIG IDEA:
 * You have two DNA strings (reference and patient). They might differ
 * in length, have insertions or deletions, or start at different points.
 * You can't just compare them position-by-position.
 *
 * The Needleman-Wunsch algorithm finds the BEST way to line them up
 * by building a 2D grid (a "dynamic programming matrix").
 *
 * HOW THE GRID WORKS (tiny example):
 *
 *   Suppose ref = "ACGT" (rows) and qry = "CGT" (columns)
 *
 *          -     C     G     T
 *    -  [  0] [  0] [  0] [  0]    ← first row: zeros (free end gaps)
 *    A  [  0] [ -1] [ -1] [ -1]
 *    C  [  0] [  1] [ -1] [ -2]
 *    G  [  0] [ -1] [  2] [  0]
 *    T  [  0] [ -1] [  0] [  3]    ← best end cell (score 3)
 *
 *   Each cell = best of three options:
 *     DIAGONAL: align ref[i] with qry[j] (match +1 or mismatch -1)
 *     UP:       gap in the query (deletion, -2)
 *     LEFT:     gap in the reference (insertion, -2)
 *
 *   Best alignment:   ref  A C G T
 *                     qry  - C G T     (the patient read starts at the C)
 *
 *   We use SEMI-GLOBAL alignment: the first row and column are zeros
 *   (free end gaps), so a short patient sequence won't be penalized
 *   for not covering the whole reference.
 *
 * MEMORY NOTE:
 *   For 3.5kb × 3.5kb sequences, the grid has ~12 million cells.
 *   That MUST go on the heap (malloc), not the stack.
 *   We store it as one flat array, indexed as [i * (qry_len+1) + j].
 *   Indexes are computed with size_t so they cannot overflow, and
 *   MAX_CELLS stops us from trying to build an absurdly large grid.
 *
 * TRACEBACK:
 *   After filling the grid, we follow arrows backward from the best
 *   end cell to reconstruct the alignment. A cell can often be reached
 *   equally well from more than one direction, so we store EVERY optimal
 *   direction as a bit flag (DIAG, UP, LEFT) in a separate 1-byte array.
 *
 *   When walking back we keep extending a gap that is already open
 *   whenever that is still optimal. With a linear gap penalty, a 3 bp
 *   deletion scores the same whether it is one block (---) or split
 *   (-G--), so without this rule an in-frame codon deletion could be
 *   reported as two frameshifts. Otherwise ties are broken
 *   DIAGONAL > UP > LEFT, which also places an indel inside a repeat
 *   at its left-most equivalent position (like `bcftools norm`).
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "align.h"

// Direction bit flags for traceback (several can be set in one cell)
#define DIAG 1
#define UP   2
#define LEFT 4

// Largest grid we are willing to allocate: 50 million cells (~250 MB).
// Enough for the 3,519 bp reference against a ~14 kb patient sequence.
#define MAX_CELLS 50000000UL

/**
 * Return the maximum of three integers.
 */
static int max3(int a, int b, int c)
{
    int max = a;
    if (b > max) max = b;
    if (c > max) max = c;
    return max;
}

AlignResult *align_sequences(const char *ref, int ref_len, const char *qry, int qry_len)
{
    if (ref == NULL || qry == NULL || ref_len <= 0 || qry_len <= 0)
    {
        fprintf(stderr, "Error: cannot align empty sequences\n");
        return NULL;
    }

    // Matrix dimensions: (ref_len+1) rows × (qry_len+1) columns
    size_t rows = (size_t) ref_len + 1;
    size_t cols = (size_t) qry_len + 1;

    // --- Refuse grids that would use too much memory ---
    size_t total_cells = rows * cols;
    if (total_cells > MAX_CELLS)
    {
        fprintf(stderr, "Error: sequences too long to align (%d bp vs %d bp).\n"
                "Expected an rpoB gene or amplicon (at most ~%lu bp against this reference).\n",
                qry_len, ref_len, (unsigned long) (MAX_CELLS / rows - 1));
        return NULL;
    }

    // --- Allocate the score matrix (one flat block on the heap) ---
    int *score = malloc(total_cells * sizeof(int));
    if (score == NULL)
    {
        fprintf(stderr, "Error: cannot allocate score matrix (%lu cells)\n",
                (unsigned long) total_cells);
        return NULL;
    }

    // Traceback direction matrix (1 byte per cell to save memory)
    unsigned char *trace = malloc(total_cells * sizeof(unsigned char));
    if (trace == NULL)
    {
        fprintf(stderr, "Error: cannot allocate traceback matrix\n");
        free(score);
        return NULL;
    }

    // --- Initialize first row and column ---
    // Semi-global: first row = 0 (free end gaps for query start)
    // First column = 0 (free end gaps for reference start)
    for (size_t j = 0; j < cols; j++)
    {
        score[j] = 0;  // First row: no penalty for starting mid-query
        trace[j] = LEFT;
    }
    for (size_t i = 0; i < rows; i++)
    {
        score[i * cols] = 0;  // First column: no penalty for starting mid-reference
        trace[i * cols] = UP;
    }
    trace[0] = 0; // Origin cell

    // --- Fill the matrix ---
    for (size_t i = 1; i < rows; i++)
    {
        for (size_t j = 1; j < cols; j++)
        {
            // Diagonal: match or mismatch
            int match = (ref[i - 1] == qry[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE;
            int diag_score = score[(i - 1) * cols + (j - 1)] + match;

            // Up: gap in query (deletion from reference)
            int up_score = score[(i - 1) * cols + j] + GAP_SCORE;

            // Left: gap in reference (insertion in query)
            int left_score = score[i * cols + (j - 1)] + GAP_SCORE;

            int best = max3(diag_score, up_score, left_score);
            score[i * cols + j] = best;

            // Remember every direction that achieves the best score
            unsigned char moves = 0;
            if (diag_score == best) moves |= DIAG;
            if (up_score == best) moves |= UP;
            if (left_score == best) moves |= LEFT;
            trace[i * cols + j] = moves;
        }
    }

    // --- Find the best endpoint (semi-global) ---
    // Check last row (reference fully consumed) and last column (query fully consumed)
    size_t best_i = rows - 1;
    size_t best_j = cols - 1;
    int best_score = score[best_i * cols + best_j];

    // Check last row (best alignment ending anywhere in query)
    for (size_t j = 0; j < cols; j++)
    {
        if (score[(rows - 1) * cols + j] > best_score)
        {
            best_score = score[(rows - 1) * cols + j];
            best_i = rows - 1;
            best_j = j;
        }
    }
    // Check last column (best alignment ending anywhere in reference)
    for (size_t i = 0; i < rows; i++)
    {
        if (score[i * cols + (cols - 1)] > best_score)
        {
            best_score = score[i * cols + (cols - 1)];
            best_i = i;
            best_j = cols - 1;
        }
    }

    // --- Traceback ---
    // Build the alignment in reverse, then flip it
    // Maximum possible alignment length = ref_len + qry_len
    size_t max_align = (size_t) ref_len + (size_t) qry_len;
    char *ref_aln = malloc(max_align + 1);
    char *qry_aln = malloc(max_align + 1);
    if (ref_aln == NULL || qry_aln == NULL)
    {
        fprintf(stderr, "Error: cannot allocate alignment buffers\n");
        free(score);
        free(trace);
        free(ref_aln);
        free(qry_aln);
        return NULL;
    }

    size_t aln_len = 0;
    // Add trailing gaps if we didn't end at the corner
    // (semi-global: free end gaps at the end too).
    // Traceback accumulates right-to-left before the entire alignment is reversed.
    // Therefore, trailing bases must be added from right-to-left (highest index down
    // to best endpoint) so that reversing the alignment places them in forward order.
    for (size_t k = rows - 1; k > best_i; k--)
    {
        ref_aln[aln_len] = ref[k - 1];
        qry_aln[aln_len] = '-';
        aln_len++;
    }
    for (size_t k = cols - 1; k > best_j; k--)
    {
        ref_aln[aln_len] = '-';
        qry_aln[aln_len] = qry[k - 1];
        aln_len++;
    }

    // Traceback starting from best endpoint
    size_t i = best_i;
    size_t j = best_j;
    unsigned char prev = DIAG; // Direction taken in the previous step

    // Follow the traceback arrows back to the origin
    while (i > 0 || j > 0)
    {
        unsigned char moves;
        if (i == 0)
        {
            moves = LEFT; // Free leading gap in the reference
        }
        else if (j == 0)
        {
            moves = UP;   // Free leading gap in the query
        }
        else
        {
            moves = trace[i * cols + j];
        }

        // Pick one optimal direction (see TRACEBACK in the header comment)
        unsigned char dir;
        if (prev != DIAG && (moves & prev))
        {
            dir = prev;   // Keep extending the gap that is already open
        }
        else if (moves & DIAG)
        {
            dir = DIAG;
        }
        else if (moves & UP)
        {
            dir = UP;
        }
        else
        {
            dir = LEFT;
        }

        if (dir == DIAG)
        {
            ref_aln[aln_len] = ref[i - 1];
            qry_aln[aln_len] = qry[j - 1];
            i--;
            j--;
        }
        else if (dir == UP)
        {
            ref_aln[aln_len] = ref[i - 1];
            qry_aln[aln_len] = '-';
            i--;
        }
        else // LEFT
        {
            ref_aln[aln_len] = '-';
            qry_aln[aln_len] = qry[j - 1];
            j--;
        }
        aln_len++;
        prev = dir;
    }

    // --- Reverse the alignment (we built it backward) ---
    for (size_t k = 0; k < aln_len / 2; k++)
    {
        char tmp;
        tmp = ref_aln[k]; ref_aln[k] = ref_aln[aln_len - 1 - k]; ref_aln[aln_len - 1 - k] = tmp;
        tmp = qry_aln[k]; qry_aln[k] = qry_aln[aln_len - 1 - k]; qry_aln[aln_len - 1 - k] = tmp;
    }

    // Null-terminate
    ref_aln[aln_len] = '\0';
    qry_aln[aln_len] = '\0';

    // Clean up the matrices
    free(score);
    free(trace);

    // --- Package the result ---
    AlignResult *result = malloc(sizeof(AlignResult));
    if (result == NULL)
    {
        fprintf(stderr, "Error: out of memory\n");
        free(ref_aln);
        free(qry_aln);
        return NULL;
    }

    result->ref_aligned = ref_aln;
    result->qry_aligned = qry_aln;
    result->length = (int) aln_len;
    result->score = best_score;

    return result;
}

void free_align_result(AlignResult *result)
{
    if (result != NULL)
    {
        free(result->ref_aligned);
        free(result->qry_aligned);
        free(result);
    }
}
