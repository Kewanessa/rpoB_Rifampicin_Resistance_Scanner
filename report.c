/**
 * report.c — Report generation and classification
 *
 * THE REPORT STRUCTURE (each part has a heading and a rule, see ui.h):
 * 1. Sample: sample ID, reference, strand
 *    Alignment & coverage: identity %, aligned bases, sequenced span,
 *    RRDR coverage
 * 2. Variants sorted into categories (each variant lands in exactly one).
 *    Each category is a section whose heading shows how many it holds:
 *    - Category 1: Known resistance-associated: a WHO group 1-2 entry in
 *      the DB, or any other amino acid change / in-frame indel in the RRDR
 *      (WHO 2023 "additional grading rule", see follows_rrdr_rule)
 *    - Category 2: Uncharacterized RRDR variant (nonsense or frameshift)
 *    - Category 3: Other non-synonymous variants and indels outside the RRDR,
 *      plus RRDR variants the DB grades "not associated" (WHO group 4-5)
 *    - Category 4: Synonymous (DNA changed but protein didn't)
 *    - Ambiguous codons: N or mixed base calls, amino acid unknown
 * 3. Summary counts and QC warnings
 * 4. Overall verdict, framed and coloured (red = resistance, yellow =
 *    needs a closer look, green = no known resistance)
 * 5. Disclaimer
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
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include "report.h"
#include "ui.h"

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

// Return 1 if c is an unambiguous DNA base (A, C, G or T).
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
 * Print a codon substitution (without a newline), e.g.
 * "S450L (TCG -> TTG, c.1349C>T)"
 */
static void print_substitution(const Variant *v)
{
    printf("%c%d%c (%s -> %s, ", v->ref_aa, v->codon_num, v->alt_aa,
           v->ref_codon, v->alt_codon);
    print_dna_change(v);
    printf(")");
}

/**
 * Print an insertion or deletion (without a newline), e.g.
 * "In-frame deletion at codon 435: c.1303_1305delGAC (3 bp)"
 * Positions are 1-based within the rpoB gene (c. = coding DNA).
 */
static void print_indel(const Variant *v)
{
    const char *frame = v->is_frameshift ? "Frameshift" : "In-frame";
    if (v->type == VAR_DELETION)
    {
        printf("%s deletion at codon %d: c.%d", frame, v->codon_num, v->nuc_pos);
        if (v->indel_len > 1)
        {
            printf("_%d", v->nuc_pos + v->indel_len - 1);
        }
        printf("del%s (%d bp)", v->indel_seq, v->indel_len);
    }
    else
    {
        printf("%s insertion at codon %d: c.%d_%dins%s (%d bp)",
               frame, v->codon_num, v->nuc_pos, v->nuc_pos + 1, v->indel_seq, v->indel_len);
    }
}

// Colour of a category: red = resistance, yellow = needs a closer look
static const char *category_color(int cat)
{
    switch (cat)
    {
        case CAT_KNOWN:     return ui.red;
        case CAT_RRDR:      return ui.yellow;
        case CAT_AMBIGUOUS: return ui.yellow;
        default:            return "";
    }
}

// Print the bullet that starts a variant line, in its category's colour
static void print_bullet(int cat)
{
    printf("    %s•%s ", category_color(cat), ui.reset);
}

/**
 * Start a variant entry that takes several lines: a blank line between
 * it and the entry above, so each one stands apart, then its bullet.
 * 'shown' counts the entries printed so far in this section.
 */
static void start_entry(int cat, int *shown)
{
    if (*shown > 0)
    {
        printf("\n");
    }
    (*shown)++;
    print_bullet(cat);
}

/**
 * Print a line under a variant (its interpretation, source...),
 * indented below the variant name and wrapped to the page width.
 * 'style' is a ui colour, or "" for plain text.
 */
static void print_detail(const char *style, const char *format, ...) UI_PRINTF(2, 3);

static void print_detail(const char *style, const char *format, ...)
{
    char text[256];
    va_list args;
    va_start(args, format);
    vsnprintf(text, sizeof(text), format, args);
    va_end(args);

    printf("      %s", style);
    ui_wrap(text, 6);
    printf("%s", ui.reset);
}

// Shown in a section that has no variants
static void print_none(void)
{
    printf("    %sNone found%s\n", ui.dim, ui.reset);
}

/**
 * Print a highlighted warning after a blank line, e.g.
 * "QC WARNING: 1 frameshift/premature stop variant(s)..."
 * Long text wraps and lines up after the label.
 */
static void print_warning(const char *label, const char *format, ...) UI_PRINTF(2, 3);

static void print_warning(const char *label, const char *format, ...)
{
    char text[256];
    va_list args;
    va_start(args, format);
    vsnprintf(text, sizeof(text), format, args);
    va_end(args);

    printf("\n    %s%s%s%s ", ui.bold, ui.yellow, label, ui.reset);
    ui_wrap(text, 4 + (int) strlen(label) + 1);
}

/**
 * Print one row of the summary counts, e.g. "Known resistance: 1".
 * A non-zero count stands out in 'color'; zero is dimmed. 'note' is
 * printed after the count unless it is "".
 */
static void print_count(const char *label, int count, const char *color, const char *note)
{
    if (count > 0)
    {
        printf("      %-18s%s%s%d%s", label, ui.bold, color, count, ui.reset);
    }
    else
    {
        printf("      %-18s%s%d%s", label, ui.dim, count, ui.reset);
    }
    if (note[0] != '\0')
    {
        printf("  %s", note);
    }
    printf("\n");
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
    // Title, set apart from the progress log above it
    printf("\n\n");
    ui_rule("═", "");
    printf("  %srpoB Rifampicin Resistance Scan Report%s\n", ui.bold, ui.reset);
    ui_rule("═", "");

    ui_heading("SAMPLE", -1, "");
    ui_field("Sample:", "%s", info->sample_id);
    ui_field("Reference:", "%s", info->ref_id);
    ui_field("Strand:", "%s",
             info->reverse_complemented ? "reverse complement (auto-detected)" : "forward");

    ui_heading("ALIGNMENT & COVERAGE", -1, "");
    ui_field("Identity:", "%.2f%%", info->identity_pct);
    ui_field("Aligned:", "%d of %d patient bases", info->aligned_bases, info->patient_length);
    if (vl->cov_start_pos > 0)
    {
        ui_field("Sequenced Span:", "bp %d - %d", vl->cov_start_pos, vl->cov_end_pos);
    }
    else
    {
        ui_field("Sequenced Span:", "none");
    }
    if (info->rrdr_resolved == RRDR_LENGTH_BP)
    {
        ui_field("RRDR Coverage:", "%sComplete%s (%d/%d bp)", ui.green, ui.reset,
                 info->rrdr_resolved, RRDR_LENGTH_BP);
    }
    else
    {
        ui_field("RRDR Coverage:", "%sINCOMPLETE%s (%d/%d bp resolved)", ui.yellow, ui.reset,
                 info->rrdr_resolved, RRDR_LENGTH_BP);
    }
    if (!info->alignment_ok)
    {
        print_warning("WARNING:",
                      "low identity or few aligned bases -\nis this M. tuberculosis rpoB?");
    }

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

    // Entries printed so far in the current section
    int shown = 0;

    // Section 1: Known Resistance Mutations
    ui_heading("KNOWN RESISTANCE-ASSOCIATED MUTATIONS", counts[CAT_KNOWN], ui.red);
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) != CAT_KNOWN)
        {
            continue;
        }
        start_entry(CAT_KNOWN, &shown);
        if (has_resistance_entry(v))
        {
            print_substitution(v);
            printf("\n");
            print_detail("", "%s%s", confidence_label(v->db_hit->confidence),
                         v->in_rrdr ? "" : " (outside RRDR)");
            if (v->db_hit->source[0] != '\0')
            {
                print_detail(ui.dim, "Source: %s", v->db_hit->source);
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
            printf("\n");
            print_detail("", "%s", confidence_label(CONF_ASSOC_W_R_INT));
            print_detail(ui.dim,
                         "Source: WHO 2023 additional grading rule (non-silent RRDR mutation)");
        }
    }
    if (counts[CAT_KNOWN] == 0) print_none();

    // Section 2: Uncharacterized RRDR Variants
    ui_heading("UNCHARACTERIZED RRDR VARIANTS", counts[CAT_RRDR], ui.yellow);
    shown = 0;
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) != CAT_RRDR)
        {
            continue;
        }
        start_entry(CAT_RRDR, &shown);
        if (v->type == VAR_SNP)
        {
            print_substitution(v);
            printf("\n");
            if (v->db_hit != NULL)
            {
                print_detail("", "Significance unknown (database: %s)",
                             confidence_label(v->db_hit->confidence));
            }
            else
            {
                print_detail("", "Significance unknown");
            }
        }
        else
        {
            print_indel(v);
            printf("\n");
            print_detail("", "Significance unknown%s",
                         v->is_frameshift ? ", reading frame disrupted" : "");
        }
    }
    if (counts[CAT_RRDR] == 0) print_none();

    // Section 3: Other Non-Synonymous Variants (mostly one line each and
    // there can be many, so no blank lines between them)
    ui_heading("OTHER NON-SYNONYMOUS VARIANTS", counts[CAT_OTHER], "");
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) != CAT_OTHER)
        {
            continue;
        }
        print_bullet(CAT_OTHER);
        if (v->type == VAR_SNP)
        {
            print_substitution(v);
            printf("%s\n", v->in_rrdr ? " (in RRDR)" : "");
            if (v->db_hit != NULL)
            {
                print_detail("", "Database: %s", confidence_label(v->db_hit->confidence));
            }
        }
        else
        {
            print_indel(v);
            printf("\n");
            if (v->is_frameshift)
            {
                print_detail("", "Reading frame disrupted");
            }
        }
    }
    if (counts[CAT_OTHER] == 0) print_none();

    // Section 4: Synonymous Mutations (one line each, so no blank lines between)
    ui_heading("SYNONYMOUS CHANGES", counts[CAT_SYNONYMOUS], "");
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) == CAT_SYNONYMOUS)
        {
            print_bullet(CAT_SYNONYMOUS);
            print_dna_change(v);
            printf(" (codon %d %s -> %s, aa=%c unchanged)%s\n",
                   v->codon_num, v->ref_codon, v->alt_codon, v->ref_aa,
                   v->in_rrdr ? " (in RRDR)" : "");
        }
    }
    if (counts[CAT_SYNONYMOUS] == 0) print_none();

    // Section 5: Ambiguous Codons (one line each; a sample with many Ns has many)
    ui_heading("AMBIGUOUS CODONS (N or mixed base calls, amino acid unknown)",
               counts[CAT_AMBIGUOUS], ui.yellow);
    for (int i = 0; i < vl->count; i++)
    {
        const Variant *v = &vl->variants[i];
        if (category(v) == CAT_AMBIGUOUS)
        {
            print_bullet(CAT_AMBIGUOUS);
            printf("codon %d (%s -> %s, ", v->codon_num, v->ref_codon, v->alt_codon);
            print_dna_change(v);
            printf(")%s\n", v->in_rrdr ? " (in RRDR)" : "");
        }
    }
    if (counts[CAT_AMBIGUOUS] == 0) print_none();

    // Summary: the total, then the count in each section
    ui_heading("SUMMARY", -1, "");
    printf("    %-19s%s%d total%s (indels: %d)\n", "Variants found:",
           ui.bold, vl->count, ui.reset, indel_count);
    char rule_note[48] = "";
    if (rule_count > 0)
    {
        snprintf(rule_note, sizeof(rule_note), "(%d by the WHO RRDR rule)", rule_count);
    }
    print_count("Known resistance:", counts[CAT_KNOWN], ui.red, rule_note);
    print_count("Unknown RRDR:", counts[CAT_RRDR], ui.yellow, "");
    print_count("Other variants:", counts[CAT_OTHER], "", "");
    print_count("Synonymous:", counts[CAT_SYNONYMOUS], "", "");
    print_count("Ambiguous codons:", counts[CAT_AMBIGUOUS], ui.yellow, "");
    if (vl->unaligned_bases > 0)
    {
        printf("\n    Ignored %d patient base(s) outside the aligned region (flanks)\n",
               vl->unaligned_bases);
    }
    if (truncating > 0)
    {
        print_warning("QC WARNING:", "%d frameshift/premature stop variant(s).\nrpoB is "
                      "essential, so these usually mean a sequencing or assembly error.",
                      truncating);
    }

    // Verdict: the first rule that applies wins (see the top of this file)
    const char *verdict;
    const char *color;
    if (!info->alignment_ok)
    {
        verdict = "INCONCLUSIVE - sample does not align well to the rpoB reference";
        color = ui.yellow;
    }
    else if (counts[CAT_KNOWN] > 0)
    {
        verdict = "RIFAMPICIN RESISTANCE-ASSOCIATED MUTATION DETECTED";
        color = ui.red;
    }
    else if (counts[CAT_RRDR] > 0)
    {
        verdict = "UNCHARACTERIZED RRDR VARIANT - phenotypic testing advised";
        color = ui.yellow;
    }
    else if (info->rrdr_resolved < RRDR_LENGTH_BP)
    {
        verdict = "INCONCLUSIVE - RRDR not fully covered or contains ambiguous bases";
        color = ui.yellow;
    }
    else
    {
        verdict = "No known resistance mutation detected in RRDR";
        color = ui.green;
    }

    // The verdict line itself stays "  VERDICT: ..." so scripts can find it
    printf("\n\n");
    ui_rule("━", color);
    printf("\n  %sVERDICT:%s %s%s%s%s\n\n", ui.bold, ui.reset, ui.bold, color, verdict, ui.reset);
    ui_rule("━", color);

    ui_heading("DISCLAIMER", -1, "");
    printf("    ");
    ui_wrap("Research and educational use only. Not validated for standalone "
            "clinical diagnostics. A negative result does not rule out rifampicin "
            "resistance: mutations outside the sequenced region, minority "
            "populations hidden in a consensus sequence and non-rpoB mechanisms "
            "are not detected.", 4);
    printf("\n");
}
