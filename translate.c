/**
 * translate.c — Codon translation table
 *
 * HOW TRANSLATION WORKS:
 * DNA is read in groups of 3 bases called "codons." Each codon
 * maps to one amino acid (or a stop signal). There are 64 possible
 * codons (4 bases × 4 bases × 4 bases = 64).
 *
 * This file contains the standard genetic code as a lookup table.
 * We convert each base to a number (A=0, C=1, G=2, T=3) and use
 * those to index into a flat array of 64 amino acids.
 *
 * Example: ATG → A=0, T=3, G=2 → index = 0*16 + 3*4 + 2 = 14 → 'M' (Met)
 *
 * kewanessa
 * CS50 Final Project
 * AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
 */

#include <ctype.h>
#include <string.h>

#include "translate.h"

/**
 * The standard genetic code, laid out in a 64-entry table.
 * Index = (first_base * 16) + (second_base * 4) + third_base
 * where A=0, C=1, G=2, T=3.
 *
 * Amino acid letters:
 *   A=Ala, C=Cys, D=Asp, E=Glu, F=Phe, G=Gly, H=His, I=Ile,
 *   K=Lys, L=Leu, M=Met, N=Asn, P=Pro, Q=Gln, R=Arg, S=Ser,
 *   T=Thr, V=Val, W=Trp, Y=Tyr, *=Stop
 */

//                   Second base:  A          C          G          T
//                   Third base: ACGT       ACGT       ACGT       ACGT
static const char CODON_TABLE[64] = {
    // First base = A (0)
    'K', 'N', 'K', 'N',   // AA_: AAA=K, AAC=N, AAG=K, AAT=N
    'T', 'T', 'T', 'T',   // AC_: ACA=T, ACC=T, ACG=T, ACT=T
    'R', 'S', 'R', 'S',   // AG_: AGA=R, AGC=S, AGG=R, AGT=S
    'I', 'I', 'M', 'I',   // AT_: ATA=I, ATC=I, ATG=M, ATT=I

    // First base = C (1)
    'Q', 'H', 'Q', 'H',   // CA_: CAA=Q, CAC=H, CAG=Q, CAT=H
    'P', 'P', 'P', 'P',   // CC_: CCA=P, CCC=P, CCG=P, CCT=P
    'R', 'R', 'R', 'R',   // CG_: CGA=R, CGC=R, CGG=R, CGT=R
    'L', 'L', 'L', 'L',   // CT_: CTA=L, CTC=L, CTG=L, CTT=L

    // First base = G (2)
    'E', 'D', 'E', 'D',   // GA_: GAA=E, GAC=D, GAG=E, GAT=D
    'A', 'A', 'A', 'A',   // GC_: GCA=A, GCC=A, GCG=A, GCT=A
    'G', 'G', 'G', 'G',   // GG_: GGA=G, GGC=G, GGG=G, GGT=G
    'V', 'V', 'V', 'V',   // GT_: GTA=V, GTC=V, GTG=V, GTT=V

    // First base = T (3)
    '*', 'Y', '*', 'Y',   // TA_: TAA=*, TAC=Y, TAG=*, TAT=Y
    'S', 'S', 'S', 'S',   // TC_: TCA=S, TCC=S, TCG=S, TCT=S
    '*', 'C', 'W', 'C',   // TG_: TGA=*, TGC=C, TGG=W, TGT=C
    'L', 'F', 'L', 'F',   // TT_: TTA=L, TTC=F, TTG=L, TTT=F
};

/**
 * Convert a nucleotide base to a numeric index.
 * A=0, C=1, G=2, T=3, anything else=-1
 */
static int base_to_index(char c)
{
    switch (toupper((unsigned char) c))
    {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1; // N, gap, or invalid
    }
}

char translate_codon(const char *codon)
{
    if (codon == NULL || strlen(codon) < 3)
    {
        return '?';
    }

    int b1 = base_to_index(codon[0]);
    int b2 = base_to_index(codon[1]);
    int b3 = base_to_index(codon[2]);

    // If any base is ambiguous (N, gap, etc.), we can't translate
    if (b1 < 0 || b2 < 0 || b3 < 0)
    {
        return '?';
    }

    // Index into the 64-entry table
    int index = b1 * 16 + b2 * 4 + b3;
    return CODON_TABLE[index];
}

int is_start_codon(const char *codon)
{
    if (codon == NULL || strlen(codon) < 3)
    {
        return 0;
    }

    // Bacteria use ATG, GTG, and TTG as start codons
    // (All are read as Met at the start of translation)
    char c[4];
    for (int i = 0; i < 3; i++)
    {
        c[i] = toupper((unsigned char) codon[i]);
    }
    c[3] = '\0';

    return (strcmp(c, "ATG") == 0 ||
            strcmp(c, "GTG") == 0 ||
            strcmp(c, "TTG") == 0);
}
