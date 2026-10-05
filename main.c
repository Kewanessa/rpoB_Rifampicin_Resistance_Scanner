/**
 * main.c — rpoB Rifampicin Resistance Scanner
 *
 * CS50 Final Project
 *
 * Usage: ./amr_scan <patient.fasta> <reference.fasta> <mutations.tsv>
 *
 * 1. Parse the patient FASTA → get patient DNA string
 * 2. Parse the reference FASTA → get reference DNA string
 * 3. Translate the reference to protein → check it is a complete gene
 * 4. Load the mutation database from TSV
 * 5. Validate the database against the reference protein
 * 6. Align patient to reference (Needleman-Wunsch), retrying with the
 *    reverse complement if the patient was read from the other strand
 * 7. Call variants (find every codon that differs)
 * 8. Look up each variant in the database
 * 9. Print the report with verdict
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta.h"
#include "align.h"
#include "translate.h"
#include "variants.h"
#include "mutdb.h"
#include "report.h"

// An alignment is trusted only if the patient matches the reference
// closely (M. tuberculosis rpoB is >99% identical to H37Rv) over at
// least MIN_ALIGNED_BASES bases (or half of a shorter sequence).
#define MIN_IDENTITY_PCT  90.0
#define MIN_ALIGNED_BASES 100

/**
 * Translate an entire DNA sequence to protein.
 * Returns a malloc'd string of amino acids (caller must free).
 * Handles the bacterial start codon (GTG/TTG → M).
 */
static char *translate_sequence(const char *dna, int dna_len)
{
    int num_codons = dna_len / 3;
    char *protein = malloc(num_codons + 1);
    if (protein == NULL)
    {
        return NULL;
    }

    for (int i = 0; i < num_codons; i++)
    {
        char codon[4] = {dna[i * 3], dna[i * 3 + 1], dna[i * 3 + 2], '\0'};
        protein[i] = translate_codon(codon);

        // Handle bacterial start codon at position 1
        if (i == 0 && is_start_codon(codon))
        {
            protein[i] = 'M';
        }
    }
    protein[num_codons] = '\0';
    return protein;
}

/**
 * Check that the reference is a complete protein-coding gene: a whole
 * number of codons, a start codon, a stop codon at the end and no
 * stop or unreadable codon in between. Codon numbers (and so the RRDR,
 * codons 426-452) are only meaningful if base 1 is the start codon.
 */
static int is_complete_gene(const char *dna, int dna_len, const char *protein)
{
    int protein_len = strlen(protein);
    if (dna_len % 3 != 0 || protein_len < 2)
    {
        return 0;
    }
    char first[4] = {dna[0], dna[1], dna[2], '\0'};
    if (!is_start_codon(first) || protein[protein_len - 1] != '*')
    {
        return 0;
    }
    for (int i = 0; i < protein_len - 1; i++)
    {
        if (protein[i] == '*' || protein[i] == '?')
        {
            return 0;
        }
    }
    return 1;
}

/**
 * Decide whether an alignment can be trusted (see MIN_IDENTITY_PCT).
 */
static int alignment_is_reliable(double identity, int aligned_bases, int patient_len)
{
    int needed = MIN_ALIGNED_BASES;
    if (patient_len / 2 < needed)
    {
        needed = patient_len / 2;
    }
    return aligned_bases > 0 && aligned_bases >= needed && identity >= MIN_IDENTITY_PCT;
}

int main(int argc, char *argv[])
{
    // ============================================================
    // Step 0: Check command-line arguments
    // ============================================================
    if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
    {
        printf("rpoB Rifampicin Resistance Scanner\n");
        printf("CS50 Final Project — kewanessa\n\n");
        printf("Usage: %s <patient.fasta> <reference.fasta> <mutations.tsv>\n", argv[0]);
        printf("       %s -h | --help\n", argv[0]);
        printf("       %s -v | --version\n\n", argv[0]);
        printf("Arguments:\n");
        printf("  patient.fasta    DNA sequence from patient isolate (full gene or RRDR fragment,\n");
        printf("                   either strand)\n");
        printf("  reference.fasta  H37Rv rpoB reference (Rv0667, NC_000962.3:759807-763325)\n");
        printf("  mutations.tsv    Mutation database (tab-separated, WHO 2023 catalogue)\n");
        return 0;
    }

    if (argc == 2 && (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0))
    {
        printf("rpoB Rifampicin Resistance Scanner version 1.1.0\n");
        printf("Author: kewanessa (CS50 Final Project)\n");
        return 0;
    }

    if (argc != 4)
    {
        fprintf(stderr, "Usage: %s <patient.fasta> <reference.fasta> <mutations.tsv>\n", argv[0]);
        fprintf(stderr, "\nArguments:\n");
        fprintf(stderr, "  patient.fasta    DNA sequence from patient isolate\n");
        fprintf(stderr, "  reference.fasta  H37Rv rpoB reference (Rv0667)\n");
        fprintf(stderr, "  mutations.tsv    Mutation database (tab-separated)\n");
        fprintf(stderr, "\nUse '%s --help' for additional information.\n", argv[0]);
        return 1;
    }

    const char *patient_file = argv[1];
    const char *reference_file = argv[2];
    const char *mutations_file = argv[3];

    // ============================================================
    // Step 1: Parse the FASTA files
    // ============================================================
    printf("Reading patient sequence from '%s'...\n", patient_file);
    FastaRecord *patient = read_fasta(patient_file);
    if (patient == NULL)
    {
        return 1;
    }
    printf("  → ID: %s\n  → Length: %d bp\n", patient->id, patient->length);

    printf("Reading reference sequence from '%s'...\n", reference_file);
    FastaRecord *reference = read_fasta(reference_file);
    if (reference == NULL)
    {
        free_fasta_record(patient);
        return 1;
    }
    printf("  → ID: %s\n  → Length: %d bp\n", reference->id, reference->length);

    // ============================================================
    // Step 2: Translate the reference to protein (for DB validation)
    // ============================================================
    printf("Translating reference to protein...\n");
    char *ref_protein = translate_sequence(reference->sequence, reference->length);
    if (ref_protein == NULL)
    {
        fprintf(stderr, "Error: failed to translate reference\n");
        free_fasta_record(patient);
        free_fasta_record(reference);
        return 1;
    }

    // Codon numbering only works if the reference is the whole gene
    if (!is_complete_gene(reference->sequence, reference->length, ref_protein))
    {
        fprintf(stderr, "Error: '%s' is not a complete protein-coding gene.\n"
                "Use the H37Rv rpoB coding sequence (3,519 bp, start to stop codon) so that\n"
                "codon numbers follow M. tuberculosis numbering.\n", reference_file);
        free(ref_protein);
        free_fasta_record(patient);
        free_fasta_record(reference);
        return 1;
    }

    // Remove the stop codon from the protein
    int protein_len = strlen(ref_protein) - 1;
    ref_protein[protein_len] = '\0';

    printf("  → Protein length: %d aa\n", protein_len);

    // ============================================================
    // Step 3: Load and validate the mutation database
    // ============================================================
    printf("Loading mutation database from '%s'...\n", mutations_file);
    MutDB *db = load_mutdb(mutations_file);
    if (db == NULL)
    {
        free(ref_protein);
        free_fasta_record(patient);
        free_fasta_record(reference);
        return 1;
    }

    // Self-check: every database entry must match the reference. A mismatch
    // means the wrong numbering system (e.g. E. coli) or the wrong reference,
    // and the results would be wrong, so we stop.
    int db_errors = validate_mutdb(db, ref_protein, protein_len);
    if (db_errors != 0)
    {
        fprintf(stderr, "Error: %d database entries do not match the reference protein.\n"
                "Check that the database uses M. tuberculosis (not E. coli) numbering\n"
                "and that the reference is the H37Rv rpoB gene.\n", db_errors);
        free(ref_protein);
        free_mutdb(db);
        free_fasta_record(patient);
        free_fasta_record(reference);
        return 1;
    }

    // ============================================================
    // Step 4: Align the sequences
    // ============================================================
    printf("Aligning sequences (%d bp vs %d bp)...\n",
           patient->length, reference->length);
    printf("  (this may take a moment for full-length genes)\n");

    AlignResult *aln = align_sequences(reference->sequence, reference->length,
                                       patient->sequence, patient->length);
    if (aln == NULL)
    {
        fprintf(stderr, "Error: alignment failed\n");
        free(ref_protein);
        free_mutdb(db);
        free_fasta_record(patient);
        free_fasta_record(reference);
        return 1;
    }

    double identity = calc_identity(aln->ref_aligned, aln->qry_aligned, aln->length);
    int aligned = count_aligned_bases(aln->ref_aligned, aln->qry_aligned, aln->length);
    int reversed = 0;

    // A read from the opposite DNA strand aligns badly: try its reverse complement
    if (!alignment_is_reliable(identity, aligned, patient->length))
    {
        char *rc = reverse_complement(patient->sequence, patient->length);
        AlignResult *rc_aln = NULL;
        if (rc != NULL)
        {
            rc_aln = align_sequences(reference->sequence, reference->length,
                                     rc, patient->length);
            free(rc);
        }
        if (rc_aln != NULL)
        {
            double rc_identity = calc_identity(rc_aln->ref_aligned, rc_aln->qry_aligned,
                                               rc_aln->length);
            int rc_aligned = count_aligned_bases(rc_aln->ref_aligned, rc_aln->qry_aligned,
                                                 rc_aln->length);
            if (alignment_is_reliable(rc_identity, rc_aligned, patient->length))
            {
                printf("  → Sample matches the opposite strand: using its reverse complement\n");
                free_align_result(aln);
                aln = rc_aln;
                identity = rc_identity;
                aligned = rc_aligned;
                reversed = 1;
            }
            else
            {
                free_align_result(rc_aln);
            }
        }
    }

    printf("  → Alignment length: %d columns, score: %d\n", aln->length, aln->score);

    // ============================================================
    // Step 5: Call variants
    // ============================================================
    printf("Calling variants...\n");
    VariantList *variants = call_variants(aln->ref_aligned, aln->qry_aligned, aln->length, db);
    if (variants == NULL)
    {
        fprintf(stderr, "Error: variant calling failed\n");
        free(ref_protein);
        free_align_result(aln);
        free_mutdb(db);
        free_fasta_record(patient);
        free_fasta_record(reference);
        return 1;
    }

    // ============================================================
    // Step 6: Print the report
    // ============================================================
    ReportInfo info;
    info.sample_id = patient->id;
    info.ref_id = reference->id;
    info.identity_pct = identity;
    info.aligned_bases = aligned;
    info.patient_length = patient->length;
    info.reverse_complemented = reversed;
    info.alignment_ok = alignment_is_reliable(identity, aligned, patient->length);
    info.rrdr_resolved = count_rrdr_resolved(aln->ref_aligned, aln->qry_aligned, aln->length);
    print_report(&info, variants);

    // ============================================================
    // Step 7: Clean up all memory
    // ============================================================
    free(ref_protein);
    free_variant_list(variants);
    free_align_result(aln);
    free_mutdb(db);
    free_fasta_record(patient);
    free_fasta_record(reference);

    return 0;
}
