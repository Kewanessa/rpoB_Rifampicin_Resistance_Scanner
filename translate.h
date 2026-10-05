/**
 * translate.h
 *
 * Header file defining interfaces for codon translation and start codon
 * detection for an rpoB scanner.
 *
 * Supports translating 3-character nucleotide codons into single-letter
 * amino acid abbreviations and identifying bacterial start codons (ATG, GTG, TTG).
 *
 * kewanessa
 * CS50 Final Project
 */

#ifndef TRANSLATE_H
#define TRANSLATE_H

/**
 * Translates a 3-character codon string into a single-letter amino acid code.
 *
 * Translates standard genetic code triplets (e.g., "ATG" -> 'M', "GAG" -> 'E').
 * Returns '*' for standard stop codons (TAA, TAG, TGA).
 *
 * @param codon Pointer to a null-terminated 3-character string representing a codon.
 * @return Single-letter uppercase amino acid code, '*' for stop codons,
 *         or '?' for unknown, invalid, or ambiguous codons.
 */
char translate_codon(const char *codon);

/**
 * Determines whether a codon is a valid start codon.
 *
 * Recognizes bacterial start codons: ATG, GTG, or TTG.
 *
 * @param codon Pointer to a null-terminated 3-character string representing a codon.
 * @return 1 if codon is ATG, GTG, or TTG; 0 otherwise.
 */
int is_start_codon(const char *codon);

#endif // TRANSLATE_H
