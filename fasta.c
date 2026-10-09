/**
 * fasta.c — FASTA file parser implementation
 *
 * How it works:
 * 1. Opens the file and reads it one character at a time with fgetc(),
 *    so lines of any length are handled
 * 2. Lines starting with '>' are headers — we grab the ID
 * 3. All other non-blank lines are sequence data — we append them
 *    to a growing buffer using realloc (doubling strategy)
 * 4. We uppercase everything and validate each character
 * 5. Returns a FastaRecord struct with id, sequence, and length
 *
 * The "growing buffer" pattern:
 *   We don't know how long the sequence is ahead of time.
 *   So we start with capacity=1024 bytes, and every time we'd
 *   overflow, we double the capacity with realloc(). This is
 *   the same idea as CS50's "resizable array" from week 5.
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasta.h"

// Initial buffer capacity (will double as needed)
#define INITIAL_CAPACITY 1024

// Longest header we keep (anything after this is ignored)
#define MAX_ID 1024

int is_valid_base(char c)
{
    // Standard bases
    if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
    {
        return 1;
    }
    // IUPAC ambiguity codes (N = any, R = purine, Y = pyrimidine, etc.)
    if (c == 'N' || c == 'R' || c == 'Y' || c == 'S' || c == 'W' ||
        c == 'K' || c == 'M' || c == 'B' || c == 'D' || c == 'H' || c == 'V')
    {
        return 1;
    }
    return 0;
}

/**
 * Complement of one base, including IUPAC ambiguity codes.
 */
static char complement_base(char c)
{
    switch (c)
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'R': return 'Y'; // A/G  <-> C/T
        case 'Y': return 'R';
        case 'K': return 'M'; // G/T  <-> A/C
        case 'M': return 'K';
        case 'B': return 'V'; // not A <-> not T
        case 'V': return 'B';
        case 'D': return 'H'; // not C <-> not G
        case 'H': return 'D';
        default:  return c;   // S, W and N are their own complements
    }
}

char *reverse_complement(const char *seq, int len)
{
    char *rc = malloc(len + 1);
    if (rc == NULL)
    {
        return NULL;
    }
    for (int i = 0; i < len; i++)
    {
        rc[len - 1 - i] = complement_base(seq[i]);
    }
    rc[len] = '\0';
    return rc;
}

/**
 * Free everything read_fasta has allocated so far and close the file.
 */
static void abandon(FILE *fp, char *seq_buf, FastaRecord *record)
{
    fclose(fp);
    free(seq_buf);
    free(record);
}

FastaRecord *read_fasta(const char *filename)
{
    // --- Step 1: Open the file ---
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: cannot open file '%s'\n", filename);
        return NULL;
    }

    // --- Step 2: Allocate the record ---
    FastaRecord *record = malloc(sizeof(FastaRecord));
    if (record == NULL)
    {
        fprintf(stderr, "Error: out of memory\n");
        fclose(fp);
        return NULL;
    }
    record->id = NULL;
    record->sequence = NULL;
    record->length = 0;

    // --- Step 3: Allocate the growing sequence buffer ---
    int capacity = INITIAL_CAPACITY;
    char *seq_buf = malloc(capacity);
    if (seq_buf == NULL)
    {
        fprintf(stderr, "Error: out of memory\n");
        abandon(fp, NULL, record);
        return NULL;
    }
    int seq_len = 0;

    // --- Step 4: Skip a UTF-8 byte order mark (EF BB BF), if any ---
    int c = fgetc(fp);
    if (c == 0xEF)
    {
        int c2 = fgetc(fp);
        int c3 = fgetc(fp);
        if (c2 != 0xBB || c3 != 0xBF)
        {
            fprintf(stderr, "Error: invalid byte 0xEF on line 1, column 1 of '%s'\n", filename);
            abandon(fp, seq_buf, record);
            return NULL;
        }
    }
    else if (c != EOF)
    {
        ungetc(c, fp);
    }

    // --- Step 5: Read the file one character at a time ---
    char id[MAX_ID];
    int id_len = 0;
    int found_header = 0;  // Seen this record's '>' line
    int in_header = 0;     // Currently inside a '>' line
    int line_start = 1;    // Next character starts a new line
    int more_records = 0;  // File holds another record after this one
    int line = 1;
    int column = 0;

    while ((c = fgetc(fp)) != EOF)
    {
        // Treat \r\n (Windows) and a lone \r (old Mac) like \n
        if (c == '\r')
        {
            int next = fgetc(fp);
            if (next != '\n' && next != EOF)
            {
                ungetc(next, fp);
            }
            c = '\n';
        }
        if (c == '\n')
        {
            in_header = 0;
            line_start = 1;
            line++;
            column = 0;
            continue;
        }
        column++;

        // Header line?
        if (line_start && c == '>')
        {
            if (found_header || seq_len > 0)
            {
                // A second record starts here — we only read the first one
                more_records = 1;
                break;
            }
            found_header = 1;
            in_header = 1;
            line_start = 0;
            continue;
        }
        line_start = 0;

        // Store the header (everything after the '>'), up to MAX_ID - 1 characters
        if (in_header)
        {
            if (id_len < MAX_ID - 1)
            {
                id[id_len++] = (char) c;
            }
            continue;
        }

        // Skip whitespace within sequence lines
        if (c == ' ' || c == '\t')
        {
            continue;
        }

        // Validate the character (c is 0-255 from fgetc, so toupper is safe)
        char base = (char) toupper(c);
        if (!is_valid_base(base))
        {
            if (isprint(c))
            {
                fprintf(stderr, "Error: invalid character '%c' on line %d, column %d of '%s'\n",
                        c, line, column, filename);
            }
            else
            {
                fprintf(stderr, "Error: invalid byte 0x%02X on line %d, column %d of '%s'\n",
                        c, line, column, filename);
            }
            abandon(fp, seq_buf, record);
            return NULL;
        }

        // --- Grow the buffer if needed (doubling strategy) ---
        // We need room for this char + the null terminator
        if (seq_len + 1 >= capacity)
        {
            // Doubling past INT_MAX would overflow (a file of over 1 GB)
            if (capacity > INT_MAX / 2)
            {
                fprintf(stderr, "Error: sequence in '%s' is too long\n", filename);
                abandon(fp, seq_buf, record);
                return NULL;
            }
            capacity *= 2;
            char *new_buf = realloc(seq_buf, capacity);
            if (new_buf == NULL)
            {
                fprintf(stderr, "Error: out of memory while growing buffer\n");
                abandon(fp, seq_buf, record);
                return NULL;
            }
            seq_buf = new_buf;
        }

        seq_buf[seq_len++] = base;
    }

    if (ferror(fp))
    {
        fprintf(stderr, "Error: could not read '%s'\n", filename);
        abandon(fp, seq_buf, record);
        return NULL;
    }
    fclose(fp);

    // --- Step 6: Validate what we got ---
    if (seq_len == 0)
    {
        fprintf(stderr, "Error: no sequence data found in '%s'\n", filename);
        free(seq_buf);
        free(record);
        return NULL;
    }
    if (more_records)
    {
        fprintf(stderr, "Note: '%s' holds more than one sequence; only the first was analysed\n",
                filename);
    }

    // Trim trailing spaces from the header
    while (id_len > 0 && isspace((unsigned char) id[id_len - 1]))
    {
        id_len--;
    }

    // If there was no header (or an empty one), use the file's base name
    const char *name = id;
    if (id_len == 0)
    {
        const char *slash = strrchr(filename, '/');
        name = (slash != NULL) ? slash + 1 : filename;
        id_len = strlen(name);
    }
    record->id = malloc(id_len + 1);
    if (record->id == NULL)
    {
        fprintf(stderr, "Error: out of memory\n");
        free(seq_buf);
        free(record);
        return NULL;
    }
    memcpy(record->id, name, id_len);
    record->id[id_len] = '\0';

    // Null-terminate the sequence
    seq_buf[seq_len] = '\0';

    // --- Step 7: Package the result ---
    record->sequence = seq_buf;
    record->length = seq_len;

    return record;
}

void free_fasta_record(FastaRecord *record)
{
    if (record != NULL)
    {
        free(record->id);
        free(record->sequence);
        free(record);
    }
}
