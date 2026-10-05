/**
 * mutdb.h — Mutation database loader and lookup
 *
 * Loads a TSV (tab-separated values) file of known resistance mutations
 * and stores them in a hash table for fast lookup.
 *
 * The hash table is the CS50 week-5 concept: an array of linked lists,
 * keyed by a hash of the mutation (codon number + new amino acid,
 * e.g. 450 + 'L' for S450L).
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#ifndef MUTDB_H
#define MUTDB_H

// Hash table size (prime number for better distribution)
#define HASH_TABLE_SIZE 211

// Confidence tiers = the five WHO 2023 catalogue grading groups
#define CONF_ASSOC_W_R      1  // "Associated with resistance" — high confidence
#define CONF_ASSOC_W_R_INT  2  // "Associated with resistance – interim" — medium
#define CONF_UNCERTAIN      3  // "Uncertain significance"
#define CONF_NOT_ASSOC_INT  4  // "Not associated with resistance – interim"
#define CONF_NOT_ASSOC      5  // "Not associated with resistance"

// A single mutation entry from the database
typedef struct MutEntry
{
    char gene[16];          // Gene name (e.g., "rpoB")
    int codon;              // Codon number (e.g., 450)
    char ref_aa;            // Reference amino acid (e.g., 'S')
    char alt_aa;            // Mutant amino acid (e.g., 'L')
    char drug[32];          // Drug name (e.g., "rifampicin")
    int confidence;         // Confidence tier (WHO group 1-5)
    char source[64];        // Source (e.g., "WHO 2023")
    struct MutEntry *next;  // Next entry in the linked list (for hash collisions)
} MutEntry;

// The mutation database: a hash table of MutEntry linked lists
typedef struct
{
    MutEntry *table[HASH_TABLE_SIZE];
    int count; // Total number of entries loaded
} MutDB;

/**
 * Load the mutation database from a TSV file.
 * Expected columns (tab-separated):
 *   gene  codon  ref_aa  alt_aa  drug  confidence  [source]
 *
 * Only rpoB rows are used. 'ref_aa' and 'alt_aa' are one-letter amino
 * acid codes ('*' = stop) and 'confidence' is the WHO group 1-5.
 * Invalid or duplicate rows are skipped with a warning.
 *
 * Returns a pointer to a MutDB, or NULL on failure (including a file
 * with no valid entries). The caller must free it with free_mutdb().
 */
MutDB *load_mutdb(const char *filename);

/**
 * Look up a mutation in the database.
 * Takes a codon number and an alternate amino acid.
 * Returns the matching MutEntry, or NULL if not found.
 */
MutEntry *lookup_mutation(MutDB *db, int codon, char alt_aa);

/**
 * Validate the database against a reference protein sequence.
 * For each entry, check that the reference amino acid at that codon
 * matches what the database expects. This catches numbering mismatches
 * (e.g., M. tuberculosis vs. E. coli numbering).
 *
 * ref_protein: the translated reference protein (1172 amino acids)
 * protein_len: length of the protein
 *
 * Returns the number of mismatches found. Prints warnings to stderr.
 */
int validate_mutdb(MutDB *db, const char *ref_protein, int protein_len);

/**
 * Free all memory used by the mutation database.
 */
void free_mutdb(MutDB *db);

#endif // MUTDB_H
