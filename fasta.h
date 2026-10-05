/**
 * fasta.h — FASTA file parser
 *
 * Reads a FASTA file and returns the first sequence record.
 * Handles: Windows (\r\n) and old Mac (\r) line endings, lines of any
 * length, blank lines, lowercase bases, a UTF-8 byte order mark,
 * missing headers, empty files, and IUPAC ambiguity codes.
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#ifndef FASTA_H
#define FASTA_H

// A single FASTA record: an ID string and a DNA/protein sequence
typedef struct
{
    char *id;       // The header line (without the '>')
    char *sequence; // The sequence, all uppercase, no whitespace
    int length;     // Length of the sequence
} FastaRecord;

/**
 * Parse the first record from a FASTA file.
 * Returns a pointer to a FastaRecord on success, or NULL on failure.
 * Prints an error message to stderr on failure, and a note if the file
 * holds more than one record (only the first is used).
 * The caller must free the result with free_fasta_record().
 */
FastaRecord *read_fasta(const char *filename);

/**
 * Free a FastaRecord and all its internal memory.
 */
void free_fasta_record(FastaRecord *record);

/**
 * Check if a character is a valid DNA base (A, C, G, T)
 * or a valid IUPAC ambiguity code (N, R, Y, S, W, K, M, B, D, H, V).
 * Input should be uppercase.
 * Returns 1 if valid, 0 if not.
 */
int is_valid_base(char c);

/**
 * Return the reverse complement of a DNA sequence (IUPAC codes are
 * complemented too, e.g. R <-> Y). Used when a patient sequence was
 * read from the opposite DNA strand.
 * Returns a malloc'd string (caller must free), or NULL if out of memory.
 */
char *reverse_complement(const char *seq, int len);

#endif // FASTA_H
