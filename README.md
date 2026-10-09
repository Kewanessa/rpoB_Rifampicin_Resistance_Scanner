# rpoB Rifampicin Resistance Scanner

#### Video Demo: <URL HERE>

#### Description:

By **kewanessa**, CS50x Final Project.

`amr_scan` is a command-line tool written in C. It reads the DNA sequence of the *rpoB* gene from a *Mycobacterium tuberculosis* sample and reports whether it carries mutations known to cause resistance to **rifampicin**, one of the two most important tuberculosis drugs. Rifampicin resistance is the main marker of multidrug-resistant TB, so detecting it quickly decides which treatment a patient gets.

About 95% of rifampicin resistance comes from mutations in an 81 bp stretch of *rpoB* called the **Rifampicin Resistance-Determining Region (RRDR)**, codons 426 to 452. The most common one worldwide is **S450L**. `amr_scan` aligns the sample to the reference gene, finds every substitution, insertion and deletion, looks each one up in a table built from the **WHO 2023 mutation catalogue**, and prints a report that ends in a verdict. It accepts a full-length gene or a short PCR amplicon, read from either DNA strand.

> For research and education only. This is not a diagnostic device.

## Usage

```bash
make             # build ./amr_scan
make test        # run the 33-test suite
make asan        # run the tests under AddressSanitizer and UndefinedBehaviorSanitizer
make valgrind    # check for memory leaks

./amr_scan <patient.fasta> <reference.fasta> <mutations.tsv>
./amr_scan tests/test_S450L.fasta ncbi_data/ncbi_dataset/data/gene.fna mutations.tsv
```

Part of the report for that sample:

```
  Strand:         forward
  Identity:       99.97%
  RRDR Coverage:  Complete (81/81 bp)

--- KNOWN RESISTANCE-ASSOCIATED MUTATIONS ---
  * S450L (TCG -> TTG, c.1349C>T) - Associated with resistance
    Source: WHO 2023

  VERDICT: RIFAMPICIN RESISTANCE-ASSOCIATED MUTATION DETECTED
```

The program exits with code 0 when it prints a report and 1 on any error (missing file, invalid FASTA, bad database, wrong reference).

## Files and design choices

**`main.c`** runs the pipeline: read both FASTA files, translate the reference, load and check the database, align, call variants, print the report, and free every allocation, on error paths too. It refuses a reference that is not a complete gene (a start codon, and a stop codon only at the end), because the RRDR is defined by codon numbers and extra flanking DNA would silently shift every one of them. If the sample aligns poorly, it tries again with the reverse complement, since a read can come from either DNA strand.

**`fasta.c` / `fasta.h`** read the first record of a FASTA file one character at a time into a buffer that doubles with `realloc` when it is full (the resizable array from week 5), so lines can be any length. They accept lowercase, Windows and old Mac line endings, a UTF-8 byte order mark and IUPAC ambiguity codes (N, R, Y…), and reject anything else with its exact line and column.

**`align.c` / `align.h`** implement **Needleman-Wunsch** alignment (+1 match, −1 mismatch, −2 per gap base). I made it *semi-global*: leading and trailing gaps are free, so a 200 bp amplicon is not punished for the 3,300 reference bases it does not cover. A 3,519 × 3,519 grid has 12.4 million cells, so the matrices are flat arrays on the heap (`score[i * cols + j]`), with a size limit for absurd inputs. One subtle bug shaped the traceback: with a linear gap penalty, a 3 bp deletion scores the same as a 1 bp plus a 2 bp deletion, so an in-frame codon deletion was reported as two frameshifts. The traceback now stores *every* optimal direction as bit flags and keeps extending an open gap whenever that is still optimal.

**`translate.c` / `translate.h`** translate codons with a 64-entry lookup table indexed by `b1 * 16 + b2 * 4 + b3` (A=0, C=1, G=2, T=3) instead of a long `switch`. The H37Rv *rpoB* gene starts with `TTG`, which bacteria read as methionine, so the start codon is special-cased to avoid a false "L1M" call.

**`variants.c` / `variants.h`** walk the alignment in reference coordinates and list synonymous and non-synonymous substitutions, in-frame and frameshift insertions and deletions, and codons with ambiguous bases. Only the stretch the sample actually covers is analysed, so unsequenced reference and primer overhangs are never reported as variants.

**`mutdb.c` / `mutdb.h`** load `mutations.tsv` into a **hash table** with chaining (211 buckets, keyed by codon number and new amino acid) for O(1) average lookups. Malformed, duplicate and non-*rpoB* rows are skipped with a warning. Every entry is then checked against the reference protein. Older papers number *rpoB* by the *E. coli* gene, 81 codons higher (S531L instead of S450L), and a database in that numbering would silently miss every mutation, so a mismatch stops the program.

**`report.c` / `report.h`** put each variant in exactly one category and print the verdict. The first rule that applies wins:

| Condition | Verdict |
|---|---|
| Sample does not align well to *rpoB* | INCONCLUSIVE |
| A WHO group 1–2 mutation, or any amino acid change or in-frame indel in the RRDR | RESISTANCE-ASSOCIATED MUTATION DETECTED |
| A nonsense or frameshift change in the RRDR | UNCHARACTERIZED RRDR VARIANT |
| Part of the RRDR missing or ambiguous | INCONCLUSIVE |
| Otherwise | No known resistance mutation detected in RRDR |

The second row follows WHO's 2023 rule that any non-silent RRDR mutation confers resistance unless there is evidence to the contrary, such as a database entry graded "not associated with resistance". Nonsense and frameshift changes stay "uncharacterized": *rpoB* is essential, so they usually mean a sequencing error, and the report adds a QC warning. A mixed base call such as Y (C or T) at codon 450 could hide S450L, so ambiguous RRDR bases make the result inconclusive rather than negative.

**`mutations.tsv`** lists the 103 *rpoB* missense mutations that the WHO 2023 catalogue grades group 1 (23) or group 2 (80), in *M. tuberculosis* numbering. **`tools/extract_who_rpoB.py`** rebuilds it from the official WHO Excel file with only the Python standard library, so every row can be traced to its source.

**`ncbi_data/`** holds the H37Rv *rpoB* reference gene and protein as downloaded from NCBI, with checksums.

**`Makefile`** builds with strict flags (`-Wall -Wextra -Werror -std=c11`) and provides the test, sanitizer and Valgrind targets.

**`tests/`** holds `run_tests.sh` and the FASTA and database files it uses. Each of the 33 tests checks the exit code, the exact verdict and key report lines. They cover known mutations (S450L, H445Y, D435V, D435F, a double mutant, I491F outside the RRDR), in-frame and frameshift indels, amplicons, a reverse-strand read, lowercase input with Windows line endings, flanking DNA, synonymous and nonsense changes, Ns and mixed bases in the RRDR, unrelated DNA, and error cases such as an empty or malformed file, an *E. coli*-numbered database and a reference that is not a whole gene. All tests pass with the normal build and under the sanitizers, and Valgrind reports no leaks.

## Data sources

- **Reference**: NCBI Datasets, *M. tuberculosis* H37Rv *rpoB* (Gene ID 888164, locus Rv0667, `NC_000962.3:759807-763325`, 3,519 bp; protein `NP_215181.1`, 1,172 aa).
- **Mutations**: WHO 2023 catalogue, file `WHO-UCN-TB-2023.7-eng.xlsx` from <https://github.com/GTB-tbsequencing/mutation-catalogue-2023> (licence ODC-By 1.0). Filter: rifampicin, *rpoB*, missense, group 1 or 2. Rebuild with `python3 tools/extract_who_rpoB.py WHO-UCN-TB-2023.7-eng.xlsx > mutations.tsv`.

## Limitations

- It reads one consensus sequence, not raw FASTQ reads, so it cannot measure low-frequency (heteroresistant) populations.
- It covers only *rpoB* and rifampicin, not other TB drugs (*katG*, *inhA*, *gyrA*…).
- In-frame indels are graded by the RRDR rule rather than looked up in the catalogue, and are reported at their left-most position (WHO names use the right-most).
- The aligner uses linear gap penalties; affine penalties would model long indels better.

## AI assistance

In line with CS50's academic honesty policy, I disclose that Anthropic's Claude (through Claude Code) helped during the final review of this project. It reviewed the code for bugs and helped fix the ones it found (indels split into frameshifts by the traceback, ambiguous RRDR bases, flanking DNA reported as insertions, an out-of-bounds read, database validation and parsing, and "not associated" database grades in the RRDR); wrote `tools/extract_who_rpoB.py` to rebuild `mutations.tsv` from the WHO catalogue; extended the test suite; and helped revise and shorten this README. Source files touched by this work say so in their header comment.

## License and references

The code is released under the MIT License (see [`LICENSE`](LICENSE)). The WHO catalogue data is used under ODC-By 1.0.

1. World Health Organization. *Catalogue of mutations in Mycobacterium tuberculosis complex and their association with drug resistance*, 2nd ed. Geneva: WHO; 2023. ISBN 978-92-4-008241-0.
2. Needleman SB, Wunsch CD. A general method applicable to the search for similarities in the amino acid sequence of two proteins. *J Mol Biol*. 1970;48(3):443-453.
3. Cole ST, Brosch R, Parkhill J, et al. Deciphering the biology of *Mycobacterium tuberculosis* from the complete genome sequence. *Nature*. 1998;393(6685):537-544.
