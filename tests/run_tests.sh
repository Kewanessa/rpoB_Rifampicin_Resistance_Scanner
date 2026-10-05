#!/bin/bash
# run_tests.sh — Test suite for rpoB Rifampicin Resistance Scanner
#
# kewanessa
# CS50 Final Project
# AI assistance: reviewed and revised with Anthropic Claude (Claude Code); see README.
#
# Run from anywhere: bash tests/run_tests.sh
# Or use: make test
# To test another build: PROG=./amr_scan_asan bash tests/run_tests.sh
#
# Every test checks the exit code, the exact VERDICT line and, where it
# matters, specific lines of the report (variant names, positions...).
# Any AddressSanitizer / UndefinedBehaviorSanitizer message is a failure.

# Run from the project root so the relative paths below always work
cd "$(dirname "$0")/.." || exit 1

PROG="${PROG:-./amr_scan}"
REF="ncbi data/ncbi_dataset/data/gene.fna"
DB="mutations.tsv"
TESTS_DIR="tests"
PASS=0
FAIL=0

# The five possible verdicts
RESISTANT="RIFAMPICIN RESISTANCE-ASSOCIATED MUTATION DETECTED"
UNCHARACTERIZED="UNCHARACTERIZED RRDR VARIANT - phenotypic testing advised"
NEGATIVE="No known resistance mutation detected in RRDR"
INCOMPLETE="INCONCLUSIVE - RRDR not fully covered or contains ambiguous bases"
POOR_ALIGNMENT="INCONCLUSIVE - sample does not align well to the rpoB reference"

# Colors for output (only when printing to a terminal)
if [ -t 1 ]; then
    GREEN=$'\033[0;32m'
    RED=$'\033[0;31m'
    NC=$'\033[0m' # No Color
else
    GREEN=""
    RED=""
    NC=""
fi

if [ ! -x "$PROG" ]; then
    echo "Error: '$PROG' not found. Run 'make' first."
    exit 1
fi

pass() {
    printf '%sPASS%s\n' "$GREEN" "$NC"
    PASS=$((PASS + 1))
}

# fail <reason lines...>
fail() {
    printf '%sFAIL%s\n' "$RED" "$NC"
    printf '    %s\n' "$@"
    FAIL=$((FAIL + 1))
}

# Did a sanitizer (make asan) report a memory or undefined-behaviour error?
sanitizer_error() {
    grep -qE "AddressSanitizer|LeakSanitizer|runtime error:" <<< "$1"
}

# run_test <description> <patient file> <expected verdict> [text...]
# Each extra argument is text that must appear in the report,
# or, if it starts with '!', text that must NOT appear.
run_test() {
    local description="$1"
    local input="$2"
    local expected="$3"
    shift 3

    printf '  Test: %s... ' "$description"

    local output exit_code verdict text
    output=$("$PROG" "$input" "$REF" "$DB" 2>&1)
    exit_code=$?
    verdict=$(sed -n 's/^  VERDICT: //p' <<< "$output")

    if [ "$exit_code" -ne 0 ]; then
        fail "expected exit code 0, got $exit_code" "$(tail -3 <<< "$output")"
        return
    fi
    if sanitizer_error "$output"; then
        fail "sanitizer reported an error"
        return
    fi
    if [ "$verdict" != "$expected" ]; then
        fail "expected verdict: $expected" "got verdict:      $verdict"
        return
    fi
    for text in "$@"; do
        if [ "${text:0:1}" = "!" ]; then
            if grep -qF -- "${text:1}" <<< "$output"; then
                fail "report should not contain: ${text:1}"
                return
            fi
        elif ! grep -qF -- "$text" <<< "$output"; then
            fail "report should contain: $text"
            return
        fi
    done
    pass
}

# expect_error <description> <expected message> <program arguments...>
# The program must exit with code 1 and print the expected message.
expect_error() {
    local description="$1"
    local message="$2"
    shift 2

    printf '  Test: %s... ' "$description"

    local output exit_code
    output=$("$PROG" "$@" 2>&1)
    exit_code=$?

    if [ "$exit_code" -ne 1 ]; then
        fail "expected exit code 1, got $exit_code"
        return
    fi
    if sanitizer_error "$output"; then
        fail "sanitizer reported an error"
        return
    fi
    if ! grep -qF -- "$message" <<< "$output"; then
        fail "expected message: $message"
        return
    fi
    pass
}

# expect_output <description> <expected text> <program arguments...>
# The program must exit with code 0 and print the expected text.
expect_output() {
    local description="$1"
    local message="$2"
    shift 2

    printf '  Test: %s... ' "$description"

    local output exit_code
    output=$("$PROG" "$@" 2>&1)
    exit_code=$?

    if [ "$exit_code" -ne 0 ]; then
        fail "expected exit code 0, got $exit_code"
        return
    fi
    if ! grep -qF -- "$message" <<< "$output"; then
        fail "expected output: $message"
        return
    fi
    pass
}

echo ""
echo "================================================"
echo "  rpoB Scanner Test Suite"
echo "================================================"

echo ""
echo "Known resistance mutations"

run_test "Wild type (no mutations)" \
    "$TESTS_DIR/test_wildtype.fasta" "$NEGATIVE" \
    "Variants found:    0 total" \
    "RRDR Coverage:  Complete (81/81 bp)"

run_test "S450L mutation (most common)" \
    "$TESTS_DIR/test_S450L.fasta" "$RESISTANT" \
    "* S450L (TCG -> TTG, c.1349C>T) - Associated with resistance" \
    "Known resistance: 1"

run_test "H445Y mutation" \
    "$TESTS_DIR/test_H445Y.fasta" "$RESISTANT" \
    "* H445Y (CAC -> TAC, c.1333C>T) - Associated with resistance"

run_test "D435V mutation" \
    "$TESTS_DIR/test_D435V.fasta" "$RESISTANT" \
    "* D435V (GAC -> GTC, c.1304A>T) - Associated with resistance"

run_test "D435F (two bases changed in one codon)" \
    "$TESTS_DIR/test_D435F.fasta" "$RESISTANT" \
    "* D435F (GAC -> TTC, c.1303_1304delinsTT) - Associated with resistance"

run_test "Double mutation (S450L + H445Y)" \
    "$TESTS_DIR/test_double.fasta" "$RESISTANT" \
    "* H445Y (CAC -> TAC, c.1333C>T)" \
    "* S450L (TCG -> TTG, c.1349C>T)" \
    "Known resistance: 2"

run_test "I491F outside the RRDR (WHO borderline)" \
    "$TESTS_DIR/test_I491F.fasta" "$RESISTANT" \
    "* I491F (ATC -> TTC, c.1471A>T) - Associated with resistance (outside RRDR)" \
    "Source: WHO 2023; borderline"

run_test "In-frame deletion in RRDR (p.Asp435del)" \
    "$TESTS_DIR/test_deletion.fasta" "$RESISTANT" \
    "* In-frame deletion at codon 435: c.1303_1305delGAC (3 bp)" \
    "Source: WHO 2023 additional grading rule" \
    "!Frameshift"

run_test "In-frame insertion in RRDR" \
    "$TESTS_DIR/test_insertion.fasta" "$RESISTANT" \
    "* In-frame insertion at codon 435: c.1302_1303insGCC (3 bp)" \
    "!Frameshift"

run_test "RRDR amplicon carrying S450L" \
    "$TESTS_DIR/test_fragment_S450L.fasta" "$RESISTANT" \
    "* S450L (TCG -> TTG, c.1349C>T)" \
    "Sequenced Span: bp 1201 - 1400"

run_test "Reverse-strand read (S450L)" \
    "$TESTS_DIR/test_revcomp_S450L.fasta" "$RESISTANT" \
    "Strand:         reverse complement (auto-detected)" \
    "* S450L (TCG -> TTG, c.1349C>T)"

run_test "Lowercase, Windows line endings, no header" \
    "$TESTS_DIR/test_lowercase_crlf.fasta" "$RESISTANT" \
    "Sample:         test_lowercase_crlf.fasta" \
    "* S450L (TCG -> TTG, c.1349C>T)"

echo ""
echo "No resistance, uncharacterized and inconclusive results"

run_test "Synonymous change (no AA change)" \
    "$TESTS_DIR/test_synonymous.fasta" "$NEGATIVE" \
    "* c.600G>T (codon 200 CCG -> CCT, aa=P unchanged)" \
    "Synonymous:       1"

run_test "Synonymous mutation inside RRDR" \
    "$TESTS_DIR/test_rrdr_silent.fasta" "$NEGATIVE" \
    "* c.1299C>T (codon 433 TTC -> TTT, aa=F unchanged) (in RRDR)"

run_test "RRDR fragment input" \
    "$TESTS_DIR/test_rrdr_fragment.fasta" "$NEGATIVE" \
    "Sequenced Span: bp 1201 - 1400" \
    "RRDR Coverage:  Complete (81/81 bp)"

run_test "Flanking DNA around the gene is ignored" \
    "$TESTS_DIR/test_flanks.fasta" "$NEGATIVE" \
    "Variants found:    0 total" \
    "Ignored 50 patient base(s)" \
    "!insertion"

run_test "Nonsense mutation in RRDR (S450*)" \
    "$TESTS_DIR/test_nonsense.fasta" "$UNCHARACTERIZED" \
    "* S450* (TCG -> TAG, c.1349C>A) - significance unknown" \
    "QC WARNING"

run_test "Frameshift deletion in RRDR" \
    "$TESTS_DIR/test_frameshift.fasta" "$UNCHARACTERIZED" \
    "* Frameshift deletion at codon 447: c.1341delC (1 bp)"

run_test "Ns in RRDR (inconclusive)" \
    "$TESTS_DIR/test_rrdr_ns.fasta" "$INCOMPLETE" \
    "* codon 440 (CTG -> NNN" \
    "RRDR Coverage:  INCOMPLETE (78/81 bp resolved)"

run_test "Mixed base call at S450 (C/T = Y)" \
    "$TESTS_DIR/test_hetero_S450.fasta" "$INCOMPLETE" \
    "* codon 450 (TCG -> TYG, c.1349C>Y) - amino acid cannot be determined (in RRDR)"

run_test "Fragment that misses the RRDR" \
    "$TESTS_DIR/test_no_rrdr.fasta" "$INCOMPLETE" \
    "RRDR Coverage:  INCOMPLETE (0/81 bp resolved)"

run_test "Unrelated DNA" \
    "$TESTS_DIR/test_unrelated.fasta" "$POOR_ALIGNMENT" \
    "WARNING:"

run_test "Sequence made only of Ns" \
    "$TESTS_DIR/test_all_n.fasta" "$POOR_ALIGNMENT" \
    "Sequenced Span: none"

echo ""
echo "Error handling"

expect_error "Empty file" "no sequence data found" \
    "$TESTS_DIR/test_empty.fasta" "$REF" "$DB"

expect_error "Malformed file (invalid characters)" "invalid character '1' on line 2, column 5" \
    "$TESTS_DIR/test_malformed.fasta" "$REF" "$DB"

expect_error "Missing patient file" "cannot open file" \
    "$TESTS_DIR/no_such_file.fasta" "$REF" "$DB"

expect_error "Wrong number of arguments" "Usage:" \
    "$TESTS_DIR/test_S450L.fasta"

expect_error "Database with no valid entries" "no valid rpoB mutations found" \
    "$TESTS_DIR/test_S450L.fasta" "$REF" "$TESTS_DIR/db_empty.tsv"

expect_error "Database in E. coli numbering (S531L)" "do not match the reference protein" \
    "$TESTS_DIR/test_S450L.fasta" "$REF" "$TESTS_DIR/db_ecoli_numbering.tsv"

expect_error "Reference that is not a whole gene" "is not a complete protein-coding gene" \
    "$TESTS_DIR/test_S450L.fasta" "$TESTS_DIR/test_rrdr_fragment.fasta" "$DB"

expect_output "Help message" "Usage:" --help

expect_output "Version message" "version" --version

echo ""
echo "================================================"
echo "  Results: ${GREEN}$PASS passed${NC}, ${RED}$FAIL failed${NC}"
echo "================================================"
echo ""

if [ "$FAIL" -eq 0 ]; then
    exit 0
else
    exit 1
fi
