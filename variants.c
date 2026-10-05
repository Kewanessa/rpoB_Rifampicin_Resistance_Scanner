/**
 * variants.c — Variant calling from aligned sequences
 *
 * WHAT THIS DOES:
 * After alignment, we have two strings of equal length:
 *   ref_aligned:  A C G T - T C G A A ...
 *   qry_aligned:  A C G T T T T G A A ...
 *                         ^   ^
 *                    insertion  SNP
 *
 * We walk through these column by column, tracking our position
 * in the reference (ignoring gaps in the reference = insertions).
 * We build codons from the reference-coordinate bases, translate
 * them, and report any amino acid differences.
 *
 * KEY CONCEPTS:
 * - "Reference position" = where we are in the original reference,
 *   skipping over gaps in the reference.
 * - "Codon number" = (ref_position / 3) + 1  (1-based)
 * - "Aligned window" = from the first to the last column where a
 *   patient base sits on a reference base. Reference outside it was
 *   not sequenced, and patient bases outside it hang over the ends of
 *   the reference (primers, flanking DNA), so neither is a variant.
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "variants.h"
#include "translate.h"

#define INITIAL_VAR_CAPACITY 64

static int add_variant(VariantList *vl, const Variant *v)
{
    if (vl->count >= vl->capacity)
    {
        int new_cap = vl->capacity * 2;
        Variant *new_arr = realloc(vl->variants, new_cap * sizeof(Variant));
        if (new_arr == NULL)
        {
            fprintf(stderr, "Error: out of memory growing variant list\n");
            return -1;
        }
        vl->variants = new_arr;
        vl->capacity = new_cap;
    }
    vl->variants[vl->count++] = *v;
    return 0;
}

/**
 * qsort comparator: order variants by reference position, then by type.
 */
static int compare_variants(const void *a, const void *b)
{
    const Variant *va = a;
    const Variant *vb = b;
    if (va->nuc_pos != vb->nuc_pos)
    {
        return (va->nuc_pos > vb->nuc_pos) - (va->nuc_pos < vb->nuc_pos);
    }
    return (va->type > vb->type) - (va->type < vb->type);
}

/**
 * Return 1 if c is an unambiguous DNA base (A, C, G or T).
 */
static int is_acgt(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

/**
 * Return 1 if any codon from first_codon to last_codon lies in the RRDR.
 */
static int overlaps_rrdr(int first_codon, int last_codon)
{
    return first_codon <= RRDR_END_CODON && last_codon >= RRDR_START_CODON;
}

/**
 * Record an insertion or deletion of len bases.
 *
 * ref_index: 0-based reference index of the first deleted base, or of
 *            the reference base just before an insertion
 * bases:     the deleted (reference) or inserted (patient) bases
 */
static int add_indel(VariantList *vl, int type, int ref_index, int len,
                     const char *bases, const char *ref_bases, int ref_len)
{
    Variant v;
    memset(&v, 0, sizeof(Variant));
    v.type = type;
    v.nuc_pos = ref_index + 1;
    v.indel_len = len;
    v.is_frameshift = (len % 3 != 0);

    // Keep the bases for the report, shortened with "..." if very long
    if (len <= MAX_INDEL_SEQ)
    {
        memcpy(v.indel_seq, bases, len);
        v.indel_seq[len] = '\0';
    }
    else
    {
        memcpy(v.indel_seq, bases, MAX_INDEL_SEQ - 3);
        strcpy(v.indel_seq + MAX_INDEL_SEQ - 3, "...");
    }

    // Codons touched: a deletion spans its own bases; an insertion lies
    // between the base at ref_index and the next one
    int first_codon = ref_index / 3 + 1;
    int last_codon = (type == VAR_DELETION) ? (ref_index + len - 1) / 3 + 1
                                            : (ref_index + 1) / 3 + 1;
    v.codon_num = (type == VAR_DELETION) ? first_codon : last_codon;
    v.in_rrdr = overlaps_rrdr(first_codon, last_codon);

    // Reference codon at the reported position (for information only)
    int codon_start = (v.codon_num - 1) * 3;
    if (codon_start + 3 <= ref_len)
    {
        memcpy(v.ref_codon, &ref_bases[codon_start], 3);
        v.ref_codon[3] = '\0';
        v.ref_aa = translate_codon(v.ref_codon);
    }
    else
    {
        v.ref_aa = '?';
    }
    v.alt_aa = '?';
    v.db_hit = NULL;

    return add_variant(vl, &v);
}

/**
 * Compare every fully sequenced codon and record substitutions
 * (synonymous or not) and codons with ambiguous base calls.
 */
static int add_codon_changes(VariantList *vl, const char *ref_bases, const char *qry_bases,
                             int ref_len, int cov_start, int cov_end, MutDB *db)
{
    int num_codons = ref_len / 3;
    for (int c = 0; c < num_codons; c++)
    {
        int codon_start = c * 3;
        int codon_num = c + 1;

        // Skip codons outside the sequenced coverage window
        if (codon_start < cov_start || (codon_start + 2) > cov_end)
        {
            continue;
        }

        // If codon contains gaps, the deletion step already reported it
        if (qry_bases[codon_start] == '-' ||
            qry_bases[codon_start + 1] == '-' ||
            qry_bases[codon_start + 2] == '-')
        {
            continue;
        }

        char ref_codon[4] = {ref_bases[codon_start], ref_bases[codon_start + 1],
                             ref_bases[codon_start + 2], '\0'};
        char qry_codon[4] = {qry_bases[codon_start], qry_bases[codon_start + 1],
                             qry_bases[codon_start + 2], '\0'};

        if (strcmp(ref_codon, qry_codon) == 0)
        {
            continue;
        }

        Variant v;
        memset(&v, 0, sizeof(Variant));
        v.codon_num = codon_num;
        memcpy(v.ref_codon, ref_codon, 4);
        memcpy(v.alt_codon, qry_codon, 4);
        v.nuc_pos = codon_start + 1;
        v.in_rrdr = overlaps_rrdr(codon_num, codon_num);
        v.db_hit = NULL;

        // Bacterial start codon: position 1 GTG/TTG translates to Methionine
        v.ref_aa = translate_codon(ref_codon);
        if (codon_num == 1 && is_start_codon(ref_codon))
        {
            v.ref_aa = 'M';
        }

        if (!is_acgt(qry_codon[0]) || !is_acgt(qry_codon[1]) || !is_acgt(qry_codon[2]))
        {
            // N or a mixed base call (e.g. Y = C/T): the amino acid is unknown
            v.type = VAR_AMBIGUOUS;
            v.alt_aa = '?';
        }
        else
        {
            v.type = VAR_SNP;
            v.alt_aa = translate_codon(qry_codon);
            if (codon_num == 1 && is_start_codon(qry_codon))
            {
                v.alt_aa = 'M';
            }
            v.is_synonymous = (v.ref_aa == v.alt_aa && v.ref_aa != '?');

            // Cache database hit during variant creation. The entry must
            // describe the same reference amino acid as this sample.
            if (!v.is_synonymous)
            {
                MutEntry *hit = lookup_mutation(db, codon_num, v.alt_aa);
                if (hit != NULL && hit->ref_aa == v.ref_aa)
                {
                    v.db_hit = hit;
                }
            }
        }

        if (add_variant(vl, &v) != 0)
        {
            return -1;
        }
    }
    return 0;
}

VariantList *call_variants(const char *ref_aligned, const char *qry_aligned,
                           int align_len, MutDB *db)
{
    if (!ref_aligned || !qry_aligned || align_len <= 0)
    {
        return NULL;
    }

    VariantList *vl = malloc(sizeof(VariantList));
    if (!vl) return NULL;

    vl->variants = malloc(INITIAL_VAR_CAPACITY * sizeof(Variant));
    if (!vl->variants)
    {
        free(vl);
        return NULL;
    }
    vl->count = 0;
    vl->capacity = INITIAL_VAR_CAPACITY;
    vl->cov_start_pos = 0;
    vl->cov_end_pos = 0;
    vl->unaligned_bases = 0;

    // Step 1: Find the aligned window (first and last columns where a
    // patient base sits on a reference base)
    int first_col = -1;
    int last_col = -1;
    for (int i = 0; i < align_len; i++)
    {
        if (ref_aligned[i] != '-' && qry_aligned[i] != '-')
        {
            if (first_col == -1) first_col = i;
            last_col = i;
        }
    }

    // Patient bases outside the window are flanking sequence, not variants
    for (int i = 0; i < align_len; i++)
    {
        if (qry_aligned[i] != '-' && (first_col == -1 || i < first_col || i > last_col))
        {
            vl->unaligned_bases++;
        }
    }

    if (first_col == -1)
    {
        return vl; // Nothing aligned: no coverage and no variants
    }

    // Step 2: Project alignment onto reference coordinates
    int total_ref_len = 0;
    for (int i = 0; i < align_len; i++)
    {
        if (ref_aligned[i] != '-') total_ref_len++;
    }

    char *ref_bases = malloc(total_ref_len + 1);
    char *qry_bases = malloc(total_ref_len + 1);
    if (!ref_bases || !qry_bases)
    {
        free(ref_bases);
        free(qry_bases);
        free_variant_list(vl);
        return NULL;
    }

    int ref_idx = 0;
    int cov_start_ref = -1;
    int cov_end_ref = -1;

    for (int i = 0; i < align_len; i++)
    {
        if (ref_aligned[i] != '-')
        {
            ref_bases[ref_idx] = ref_aligned[i];

            // Only capture patient sequence within the aligned window
            if (i >= first_col && i <= last_col)
            {
                qry_bases[ref_idx] = qry_aligned[i]; // A base, or '-' if deleted
                if (cov_start_ref == -1) cov_start_ref = ref_idx;
                cov_end_ref = ref_idx;
            }
            else
            {
                qry_bases[ref_idx] = ' '; // Explicitly marks unsequenced regions
            }
            ref_idx++;
        }
    }
    ref_bases[total_ref_len] = '\0';
    qry_bases[total_ref_len] = '\0';

    vl->cov_start_pos = cov_start_ref + 1; // 1-based
    vl->cov_end_pos = cov_end_ref + 1;

    int ok = 1;

    // Step 3: Deletions = runs of '-' in the patient inside the window.
    // The window starts and ends on patient bases, so every run is internal.
    int pos = cov_start_ref;
    while (ok && pos <= cov_end_ref)
    {
        if (qry_bases[pos] != '-')
        {
            pos++;
            continue;
        }
        int start = pos;
        while (qry_bases[pos] == '-')
        {
            pos++;
        }
        ok = add_indel(vl, VAR_DELETION, start, pos - start,
                       &ref_bases[start], ref_bases, total_ref_len) == 0;
    }

    // Step 4: Insertions = runs of reference gaps inside the window
    int last_ref = cov_start_ref - 1; // Index of the last reference base passed
    int col = first_col;
    while (ok && col <= last_col)
    {
        if (ref_aligned[col] != '-')
        {
            last_ref++;
            col++;
            continue;
        }
        int start = col;
        while (ref_aligned[col] == '-') // Stops at last_col at the latest
        {
            col++;
        }
        ok = add_indel(vl, VAR_INSERTION, last_ref, col - start,
                       &qry_aligned[start], ref_bases, total_ref_len) == 0;
    }

    // Step 5: Codon substitutions and ambiguous codons
    if (ok)
    {
        ok = add_codon_changes(vl, ref_bases, qry_bases, total_ref_len,
                               cov_start_ref, cov_end_ref, db) == 0;
    }

    free(ref_bases);
    free(qry_bases);

    if (!ok)
    {
        free_variant_list(vl);
        return NULL;
    }

    // Report variants in reference order
    qsort(vl->variants, vl->count, sizeof(Variant), compare_variants);
    return vl;
}

void free_variant_list(VariantList *vl)
{
    if (vl)
    {
        free(vl->variants);
        free(vl);
    }
}
