/**
 * report.c — Report generation and classification
 *
 * THE REPORT STRUCTURE:
 * 1. Header: sample ID, reference, strand, alignment identity %,
 *    aligned bases, sequenced span, RRDR coverage
 * 2. Variants sorted into categories (each variant lands in exactly one):
 *    - Category 1: Known resistance-associated: a WHO group 1-2 entry in
 *      the DB, or any other amino acid change / in-frame indel in the RRDR
 *      (WHO 2023 "additional grading rule", see follows_rrdr_rule)
 *    - Category 2: Uncharacterized RRDR variant (nonsense or frameshift)
 *    - Category 3: Other non-synonymous variants and indels outside the RRDR,
 *      plus RRDR variants the DB grades "not associated" (WHO group 4-5)
 *    - Category 4: Synonymous (DNA changed but protein didn't)
 *    - Ambiguous codons: N or mixed base calls, amino acid unknown
 * 3. Overall verdict
 * 4. Disclaimer
 *
 * THE VERDICT LOGIC (the first rule that applies wins):
 * - Poor alignment → "Inconclusive, sample does not align well to rpoB"
 * - Any Cat 1 → "Resistance-associated mutation detected"
 * - Else any Cat 2 → "Uncharacterized RRDR variant, phenotypic testing advised"
 * - Else if any RRDR base is missing, a gap or ambiguous → "Inconclusive"
 * - Else → "No known resistance mutation in RRDR"
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */
#include <stdio.h>
#include <string.h>
#include "report.h"

// Report categories (see above)
#define CAT_KNOWN      1
#define CAT_RRDR       2
#define CAT_OTHER      3
#define CAT_SYNONYMOUS 4
#define CAT_AMBIGUOUS  5

static const char *confidence_label(int tier)
{
    switch (tier)
    {
        case CONF_ASSOC_W_R:     return "Associated with resistance";
        case CONF_ASSOC_W_R_INT: return "Associated with resistance (interim)";
        case CONF_UNCERTAIN:     return "Uncertain significance";
        case CONF_NOT_ASSOC_INT: return "Not associated with resistance (interim)";
        case CONF_NOT_ASSOC:     return "Not associated with resistance";
        default:                 return "Unknown";
    }
}

/**
 * Return 1 if c is an unambiguous DNA base (A, C, G or T).
 */
static int is_base(char c)
{
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

/**
 * Return 1 if the variant has a WHO group 1-2 ("associated with
 * resistance") entry in the mutation database.
 */
static int has_resistance_entry(const Variant *v)
{
    return v->type == VAR_SNP && v->db_hit != NULL && v->db_hit->confidence <= CONF_ASSOC_W_R_INT;
}

/**
 * Return 1 if the database grades the variant "not associated with
 * resistance" (WHO group 4-5).
 */
static int is_not_associated(const Variant *v)
{
    return v->db_hit != NULL && v->db_hit->confidence >= CONF_NOT_ASSOC_INT;
}

/**
 * WHO 2023 "additional grading rule" for rifampicin: any non-silent
 * RRDR mutation is assumed to confer resistance (group 2, interim)
 * unless there is evidence to the contrary. We apply it to amino acid
 * substitutions and in-frame indels. Nonsense and frameshift changes
 * stay "uncharacterized": rpoB is essential, so they usually point to a
 * sequencing error rather than a resistant strain.
 */
static int follows_rrdr_rule(const Variant *v)
{
    if (!v->in_rrdr)
    {
        return 0;
    }
    // Evidence to the contrary: the database grades it "not associated"
    if (is_not_associated(v))
    {
        return 0;
    }
    if (v->type == VAR_SNP)
    {
        return !v->is_synonymous && v->alt_aa != '*' && v->alt_aa != '?';
    }
    if (v->type == VAR_INSERTION || v->type == VAR_DELETION)
    {
        return !v->is_frameshift;
    }
    return 0;
}

/**
 * Decide which report category a variant belongs to.
 * Only WHO groups 1-2 count as known resistance. An RRDR variant the
 * database grades "not associated" (groups 4-5) is not uncharacterized,
 * so it is listed with the other non-synonymous variants.
 */
static int category(const Variant *v)
{
    if (v->type == VAR_AMBIGUOUS)
    {
        return CAT_AMBIGUOUS;
    }
    if (v->type == VAR_SNP && v->is_synonymous)
    {
        return CAT_SYNONYMOUS;
    }
    if (has_resistance_entry(v) || follows_rrdr_rule(v))
    {
        return CAT_KNOWN;
    }
    return (v->in_rrdr && !is_not_associated(v)) ? CAT_RRDR : CAT_OTHER;
}

/**
 * Print the exact DNA change of a codon substitution in c. notation
 * (positions within the rpoB gene), e.g. "c.1349C>T" for one base or
 * "c.1303_1304delinsTT" when several bases of the codon changed.
 */
static void print_dna_change(const Variant *v)
{
    int first = -1;
    int last = -1;
    for (int k = 0; k < 3; k++)
    {
        if (v->ref_codon[k] != v->alt_codon[k])
        {
            if (first == -1) first = k;
            last = k;
        }
    }
    if (first == -1)
    {
        return;
    }
    if (first == last)
    {
        printf("c.%d%c>%c", v->nuc_pos + first, v->ref_codon[first], v->alt_codon[first]);
    }
    else
    {
        printf("c.%d_%ddelins%.*s", v->nuc_pos + first, v->nuc_pos + last,
               last - first + 1, &v->alt_codon[first]);
    }
}

/**
 * Print the start of a line describing a codon substitution, e.g.
 * "  * S450L (TCG -> TTG, c.1349C>T)"
 */
static void print_substitution(const Variant *v)
{
    printf("  * %c%d%c (%s -> %s, ", v->ref_aa, v->codon_num, v->alt_aa,
           v->ref_codon, v->alt_codon);
    print_dna_change(v);
    printf(")");
}

/**
 * Print the start of a line describing an insertion or deletion, e.g.
 * "  * In-frame deletion at codon 435: c.1303_1305delGAC (3 bp)"
 * Positions are 1-based within the rpoB gene (c. = coding DNA).
 */
static void print_indel(const Variant *v)
{
    const char *frame = v->is_frameshift ? "Frameshift" : "In-frame";
    if (v->type == VAR_DELETION)
    {
        printf("  * %s deletion at codon %d: c.%d", frame, v->codon_num, v->nuc_pos);
        if (v->indel_len > 1)
        {
            printf("_%d", v->nuc_pos + v->indel_len - 1);
        }
        printf("del%s (%d bp)", v->indel_seq, v->indel_len);
    }
    else
    {
        printf("  * %s insertion at codon %d: c.%d_%dins%s (%d bp)",
               frame, v->codon_num, v->nuc_pos, v->nuc_pos + 1, v->indel_seq, v->indel_len);
    }
}

double calc_identity(const char *ref_aligned, const char *qry_aligned, int align_len)
{
    int matches = 0;
    int compared = 0;

    for (int i = 0; i < align_len; i++)
    {
        // Only compare reference bases against unambiguous patient bases
        if (ref_aligned[i] == '-' || !is_base(qry_aligned[i])) continue;
        compared++;
        if (ref_aligned[i] == qry_aligned[i]) matches++;
    }

    return compared > 0 ? (100.0 * matches) / compared : 0.0;
}

int count_aligned_bases(const char *ref_aligned, const char *qry_aligned, int align_len)
{
    int aligned = 0;
    for (int i = 0; i < align_len; i++)
    {
        if (ref_aligned[i] != '-' && qry_aligned[i] != '-') aligned++;
    }
    return aligned;
}

int count_rrdr_resolved(const char *ref_aligned, const char *qry_aligned, int align_len)
{
    int rrdr_start = (RRDR_START_CODON - 1) * 3; // 1275 (0-based)
    int rrdr_end = (RRDR_END_CODON * 3) - 1;     // 1355 (0-based)

    // Aligned window: first and last columns with a patient base on a reference base
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

    int ref_pos = 0;
    int resolved = 0;
    for (int i = 0; i < align_len; i++)
    {
        if (ref_aligned[i] != '-')
        {
            // A clear base, or a deletion inside the sequenced region,
            // counts as resolved; end gaps, N and IUPAC codes do not
            int deleted = (qry_aligned[i] == '-' && i > first_col && i < last_col);
            if (ref_pos >= rrdr_start && ref_pos <= rrdr_end &&
                (is_base(qry_aligned[i]) || deleted))
            {
                resolved++;
            }
            ref_pos++;
        }
    }
    return resolved;
}

void print_report(const ReportInfo *info, const VariantList *vl)
{
    printf("\n========================================================\n");
    printf("  rpoB Rifampicin Resistance Scan Report\n");
    printf("========================================================\n");
    printf("  Sample:         %s\n", info->sample_id);
    printf("  Reference:      %s\n", info->ref_id);
    printf("  Strand:         %s\n",
           info->reverse_complemented ? "reverse complement (auto-detected)" : "forward");
    printf("  Identity:       %.2f%%\n", info->identity_pct);
    printf("  Aligned:        %d of %d patient bases\n", info->aligned_bases, info->patient_length);
    if (vl->cov_start_pos > 0)
    {
        printf("  Sequenced Span: bp %d - %d\n", vl->cov_start_pos, vl->cov_end_pos);
    }
    else
    {
        printf("  Sequenced Span: none\n");
    }
    if (info->rrdr_resolved == RRDR_LENGTH_BP)
    {
        printf("  RRDR Coverage:  Complete (%d/%d bp)\n", info->rrdr_resolved, RRDR_LENGTH_BP);
    }
    else
    {
        printf("  RRDR Coverage:  INCOMPLETE (%d/%d bp resolved)\n",
               info->rrdr_resolved, RRDR_LENGTH_BP);
    }
    if (!info->alignment_ok)
    {
        printf("  WARNING:        low identity or few aligned bases - is this M. tuberculosis rpoB?\n");
    }
    printf("========================================================\n\n");

    // Count of variants in each category (index = category number)
    int counts[CAT_AMBIGUOUS + 1] = {0};
    int rule_count = 0;   // Known resistance by the WHO RRDR rule only
    int indel_count = 0;
    int truncating = 0;   // Frameshifts and premature stop codons
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        counts[category(v)]++;
        if (category(v) == CAT_KNOWN && !has_resistance_entry(v))
        {
            rule_count++;
        }
        if (v->type == VAR_INSERTION || v->type == VAR_DELETION)
        {
            indel_count++;
            truncating += v->is_frameshift;
        }
        if (v->type == VAR_SNP && v->alt_aa == '*' && v->ref_aa != '*')
        {
            truncating++;
        }
    }

    // --- Section 1: Known Resistance Mutations ---
    printf("--- KNOWN RESISTANCE-ASSOCIATED MUTATIONS ---\n");
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) != CAT_KNOWN)
        {
            continue;
        }
        if (has_resistance_entry(v))
        {
            print_substitution(v);
            printf(" - %s%s\n", confidence_label(v->db_hit->confidence),
                   v->in_rrdr ? "" : " (outside RRDR)");
            if (v->db_hit->source[0] != '\0')
            {
                printf("    Source: %s\n", v->db_hit->source);
            }
        }
        else
        {
            if (v->type == VAR_SNP)
            {
                print_substitution(v);
            }
            else
            {
                print_indel(v);
            }
            printf(" - %s\n", confidence_label(CONF_ASSOC_W_R_INT));
            printf("    Source: WHO 2023 additional grading rule (non-silent RRDR mutation)\n");
        }
    }
    if (counts[CAT_KNOWN] == 0) printf("  (none)\n");

    // --- Section 2: Uncharacterized RRDR Variants ---
    printf("\n--- UNCHARACTERIZED RRDR VARIANTS ---\n");
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) != CAT_RRDR)
        {
            continue;
        }
        if (v->type == VAR_SNP)
        {
            print_substitution(v);
            printf(" - significance unknown");
            if (v->db_hit != NULL)
            {
                printf(" (database: %s)", confidence_label(v->db_hit->confidence));
            }
            printf("\n");
        }
        else
        {
            print_indel(v);
            printf(" - significance unknown%s\n",
                   v->is_frameshift ? ", reading frame disrupted" : "");
        }
    }
    if (counts[CAT_RRDR] == 0) printf("  (none)\n");

    // --- Section 3: Other Non-Synonymous Variants ---
    printf("\n--- OTHER NON-SYNONYMOUS VARIANTS ---\n");
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) != CAT_OTHER)
        {
            continue;
        }
        if (v->type == VAR_SNP)
        {
            print_substitution(v);
            if (v->db_hit != NULL)
            {
                printf(" (database: %s)", confidence_label(v->db_hit->confidence));
            }
            printf("%s\n", v->in_rrdr ? " (in RRDR)" : "");
        }
        else
        {
            print_indel(v);
            printf("%s\n", v->is_frameshift ? " - reading frame disrupted" : "");
        }
    }
    if (counts[CAT_OTHER] == 0) printf("  (none)\n");

    // --- Section 4: Synonymous Mutations ---
    printf("\n--- SYNONYMOUS CHANGES ---\n");
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) == CAT_SYNONYMOUS)
        {
            printf("  * ");
            print_dna_change(v);
            printf(" (codon %d %s -> %s, aa=%c unchanged)%s\n",
                   v->codon_num, v->ref_codon, v->alt_codon, v->ref_aa,
                   v->in_rrdr ? " (in RRDR)" : "");
        }
    }
    if (counts[CAT_SYNONYMOUS] == 0) printf("  (none)\n");

    // --- Section 5: Ambiguous Codons ---
    printf("\n--- AMBIGUOUS CODONS (N or mixed base calls) ---\n");
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) == CAT_AMBIGUOUS)
        {
            printf("  * codon %d (%s -> %s, ", v->codon_num, v->ref_codon, v->alt_codon);
            print_dna_change(v);
            printf(") - amino acid cannot be determined%s\n", v->in_rrdr ? " (in RRDR)" : "");
        }
    }
    if (counts[CAT_AMBIGUOUS] == 0) printf("  (none)\n");

    // --- Summary & Clinical Verdict ---
    printf("\n--- SUMMARY ---\n");
    printf("  Variants found:    %d total\n", vl->count);
    printf("    Known resistance: %d", counts[CAT_KNOWN]);
    if (rule_count > 0)
    {
        printf(" (%d by the WHO RRDR rule)", rule_count);
    }
    printf("\n");
    printf("    Unknown RRDR:     %d\n", counts[CAT_RRDR]);
    printf("    Other variants:   %d\n", counts[CAT_OTHER]);
    printf("    Synonymous:       %d\n", counts[CAT_SYNONYMOUS]);
    printf("    Ambiguous codons: %d\n", counts[CAT_AMBIGUOUS]);
    printf("    (Indels among these: %d)\n", indel_count);
    if (vl->unaligned_bases > 0)
    {
        printf("  Ignored %d patient base(s) outside the aligned region (flanks)\n",
               vl->unaligned_bases);
    }
    if (truncating > 0)
    {
        printf("  QC WARNING: %d frameshift/premature stop variant(s). rpoB is essential,\n"
               "  so these usually mean a sequencing or assembly error.\n", truncating);
    }

    printf("\n========================================================\n");
    printf("  VERDICT: ");
    if (!info->alignment_ok)
    {
        printf("INCONCLUSIVE - sample does not align well to the rpoB reference\n");
    }
    else if (counts[CAT_KNOWN] > 0)
    {
        printf("RIFAMPICIN RESISTANCE-ASSOCIATED MUTATION DETECTED\n");
    }
    else if (counts[CAT_RRDR] > 0)
    {
        printf("UNCHARACTERIZED RRDR VARIANT - phenotypic testing advised\n");
    }
    else if (info->rrdr_resolved < RRDR_LENGTH_BP)
    {
        printf("INCONCLUSIVE - RRDR not fully covered or contains ambiguous bases\n");
    }
    else
    {
        printf("No known resistance mutation detected in RRDR\n");
    }
    printf("========================================================\n");

    printf("\n  DISCLAIMER: Research and educational use only.\n");
    printf("  Not validated for standalone clinical diagnostics.\n");
    printf("  A negative result does not rule out rifampicin resistance: mutations\n");
    printf("  outside the sequenced region, minority populations hidden in a\n");
    printf("  consensus sequence and non-rpoB mechanisms are not detected.\n\n");
}
