/**
 * mutdb.c — Mutation database with hash table
 *
 * HOW THE HASH TABLE WORKS (CS50 week 5 concept):
 *
 * We want to look up mutations fast. Instead of searching through
 * a list every time, we use a hash table:
 *
 * 1. Take the mutation key (e.g., codon=450, alt_aa='L')
 * 2. Run it through a hash function → gives a number 0..210
 * 3. That number is the "bucket" in our array
 * 4. Store the entry in that bucket's linked list
 *
 * Lookup is O(1) on average instead of O(n).
 *
 *   table[0]  → NULL
 *   table[1]  → [S450L] → NULL
 *   table[2]  → [H445Y] → [H445D] → NULL   ← two entries hashed to same bucket
 *   ...
 *   table[210] → NULL
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mutdb.h"

#define MAX_LINE 1024

// Columns: gene, codon, ref_aa, alt_aa, drug, confidence, source (optional)
#define MIN_FIELDS 6
#define MAX_FIELDS 7

/**
 * Hash function: combines codon number and alternate amino acid
 * into a hash table index.
 */
static unsigned int hash_mutation(int codon, char alt_aa)
{
    // Simple hash: codon * 31 + ascii value of alt_aa
    // (unsigned arithmetic, so a huge codon number cannot overflow)
    unsigned int h = (unsigned int) codon * 31u + (unsigned char) alt_aa;
    return h % HASH_TABLE_SIZE;
}

/**
 * Return 1 if c is a one-letter amino acid code or '*' (stop codon).
 */
static int is_amino_acid(char c)
{
    return c != '\0' && strchr("ACDEFGHIKLMNPQRSTVWY*", c) != NULL;
}

/**
 * Split a line into tab-separated fields, keeping empty fields.
 * Each tab is replaced by '\0' and fields[] points into the line.
 * Returns the number of fields (max_fields + 1 means "too many").
 */
static int split_tabs(char *line, char *fields[], int max_fields)
{
    int count = 0;
    char *start = line;
    while (1)
    {
        if (count == max_fields)
        {
            return max_fields + 1;
        }
        fields[count++] = start;
        char *tab = strchr(start, '\t');
        if (tab == NULL)
        {
            return count;
        }
        *tab = '\0';
        start = tab + 1;
    }
}

/**
 * Read a whole field as a whole number. Returns 1 on success, 0 if the
 * field is empty, has extra characters (e.g. "1.5") or is out of range.
 */
static int parse_int(const char *text, int *value)
{
    char *end;
    errno = 0;
    long number = strtol(text, &end, 10);
    if (end == text || *end != '\0' || errno == ERANGE || number < INT_MIN || number > INT_MAX)
    {
        return 0;
    }
    *value = (int) number;
    return 1;
}

MutDB *load_mutdb(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: cannot open mutation database '%s'\n", filename);
        return NULL;
    }

    // Allocate and initialize the database
    MutDB *db = malloc(sizeof(MutDB));
    if (db == NULL)
    {
        fprintf(stderr, "Error: out of memory\n");
        fclose(fp);
        return NULL;
    }
    for (int i = 0; i < HASH_TABLE_SIZE; i++)
    {
        db->table[i] = NULL;
    }
    db->count = 0;

    char line[MAX_LINE];
    int line_num = 0;
    int other_genes = 0;

    while (fgets(line, MAX_LINE, fp) != NULL)
    {
        line_num++;

        // A line longer than the buffer: ignore the rest of it
        int len = strlen(line);
        if (len > 0 && line[len - 1] != '\n' && !feof(fp))
        {
            int ch;
            while ((ch = fgetc(fp)) != '\n' && ch != EOF)
            {
                // Discard
            }
        }

        // Strip trailing newline
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
        {
            line[--len] = '\0';
        }

        // Skip empty lines and comment lines
        if (len == 0 || line[0] == '#')
        {
            continue;
        }

        // Skip header line (starts with "gene")
        if (strncmp(line, "gene", 4) == 0)
        {
            continue;
        }

        // Parse: gene<tab>codon<tab>ref_aa<tab>alt_aa<tab>drug<tab>confidence[<tab>source]
        char *fields[MAX_FIELDS];
        int num_fields = split_tabs(line, fields, MAX_FIELDS);
        int codon;
        int confidence;
        if (num_fields < MIN_FIELDS || num_fields > MAX_FIELDS ||
            !parse_int(fields[1], &codon) || !parse_int(fields[5], &confidence))
        {
            fprintf(stderr, "Warning: skipping malformed line %d in '%s'\n", line_num, filename);
            continue;
        }
        const char *gene = fields[0];
        const char *ref_aa_str = fields[2];
        const char *alt_aa_str = fields[3];
        const char *drug = fields[4];
        const char *source = (num_fields == MAX_FIELDS) ? fields[6] : "";

        // This scanner only looks at rpoB
        if (strcmp(gene, "rpoB") != 0)
        {
            other_genes++;
            continue;
        }

        // Check the values make sense (one-letter amino acids, WHO group 1-5)
        char ref_aa = toupper((unsigned char) ref_aa_str[0]);
        char alt_aa = toupper((unsigned char) alt_aa_str[0]);
        if (codon < 1 || strlen(ref_aa_str) != 1 || strlen(alt_aa_str) != 1 ||
            !is_amino_acid(ref_aa) || !is_amino_acid(alt_aa) || ref_aa == alt_aa ||
            confidence < CONF_ASSOC_W_R || confidence > CONF_NOT_ASSOC)
        {
            fprintf(stderr, "Warning: skipping invalid entry on line %d in '%s'\n", line_num, filename);
            continue;
        }

        // The same mutation listed twice: keep the first one
        if (lookup_mutation(db, codon, alt_aa) != NULL)
        {
            fprintf(stderr, "Warning: duplicate entry %c%d%c on line %d in '%s' ignored\n",
                    ref_aa, codon, alt_aa, line_num, filename);
            continue;
        }

        // Create a new entry (snprintf truncates long text safely)
        MutEntry *entry = malloc(sizeof(MutEntry));
        if (entry == NULL)
        {
            fprintf(stderr, "Error: out of memory\n");
            fclose(fp);
            free_mutdb(db);
            return NULL;
        }

        snprintf(entry->gene, sizeof(entry->gene), "%s", gene);
        entry->codon = codon;
        entry->ref_aa = ref_aa;
        entry->alt_aa = alt_aa;
        snprintf(entry->drug, sizeof(entry->drug), "%s", drug);
        entry->confidence = confidence;
        snprintf(entry->source, sizeof(entry->source), "%s", source);

        // Insert into hash table (prepend to the linked list at this bucket)
        unsigned int bucket = hash_mutation(codon, entry->alt_aa);
        entry->next = db->table[bucket];
        db->table[bucket] = entry;
        db->count++;
    }

    fclose(fp);

    if (other_genes > 0)
    {
        printf("Note: skipped %d entries for genes other than rpoB\n", other_genes);
    }

    // An empty database would make every sample look susceptible
    if (db->count == 0)
    {
        fprintf(stderr, "Error: no valid rpoB mutations found in '%s'\n", filename);
        free_mutdb(db);
        return NULL;
    }

    printf("Loaded %d mutation(s) from database\n", db->count);
    return db;
}

MutEntry *lookup_mutation(MutDB *db, int codon, char alt_aa)
{
    if (db == NULL)
    {
        return NULL;
    }

    unsigned int bucket = hash_mutation(codon, alt_aa);
    MutEntry *entry = db->table[bucket];

    // Walk the linked list in this bucket
    while (entry != NULL)
    {
        if (entry->codon == codon && entry->alt_aa == alt_aa)
        {
            return entry;
        }
        entry = entry->next;
    }

    return NULL; // Not found
}

int validate_mutdb(MutDB *db, const char *ref_protein, int protein_len)
{
    if (db == NULL || ref_protein == NULL)
    {
        return -1;
    }

    int mismatches = 0;

    // Walk every entry in every bucket
    for (int i = 0; i < HASH_TABLE_SIZE; i++)
    {
        MutEntry *entry = db->table[i];
        while (entry != NULL)
        {
            // Check: is the database's ref_aa correct for this codon?
            int pos = entry->codon - 1; // Convert to 0-based
            if (pos < 0 || pos >= protein_len)
            {
                fprintf(stderr, "WARNING: database entry %c%d%c — codon %d is out of range "
                        "(protein has %d residues)\n",
                        entry->ref_aa, entry->codon, entry->alt_aa,
                        entry->codon, protein_len);
                mismatches++;
            }
            else if (ref_protein[pos] != entry->ref_aa)
            {
                fprintf(stderr, "WARNING: database says position %d should be '%c', "
                        "but reference protein has '%c'. "
                        "Possible numbering mismatch (MTB vs E.coli)!\n",
                        entry->codon, entry->ref_aa, ref_protein[pos]);
                mismatches++;
            }
            entry = entry->next;
        }
    }

    if (mismatches == 0)
    {
        printf("Database validation: all entries match the reference protein ✓\n");
    }
    else
    {
        fprintf(stderr, "Database validation: %d mismatch(es) found — check numbering!\n", mismatches);
    }

    return mismatches;
}

void free_mutdb(MutDB *db)
{
    if (db == NULL)
    {
        return;
    }

    // Free every linked list in every bucket
    for (int i = 0; i < HASH_TABLE_SIZE; i++)
    {
        MutEntry *entry = db->table[i];
        while (entry != NULL)
        {
            MutEntry *next = entry->next;
            free(entry);
            entry = next;
        }
    }
    free(db);
}
