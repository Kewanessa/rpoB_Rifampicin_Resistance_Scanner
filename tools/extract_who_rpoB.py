#!/usr/bin/env python3
"""
extract_who_rpoB.py - Rebuild mutations.tsv from the official WHO catalogue

Reads the WHO "Catalogue of mutations in Mycobacterium tuberculosis complex
and their association with drug resistance, 2nd edition" (2023) Excel file
and prints every rpoB rifampicin missense mutation graded
"1) Assoc w R" or "2) Assoc w R - Interim" in the TSV format that
amr_scan reads.

Get the Excel file (WHO-UCN-TB-2023.7-eng.xlsx) from:
  https://github.com/GTB-tbsequencing/mutation-catalogue-2023
  (folder "Final Result Files")

Usage:
  python3 tools/extract_who_rpoB.py WHO-UCN-TB-2023.7-eng.xlsx > mutations.tsv

Only the Python standard library is used (an .xlsx file is a zip of XML).

kewanessa
CS50 Final Project
AI assistance: written with Anthropic Claude (Claude Code); see README.
"""

import re
import sys
import zipfile
import xml.etree.ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
SHEET_NAME = "Catalogue_master_file"

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# WHO "Additional grading criteria applied" -> short note for the source column
CRITERIA_NOTES = {
    "": "",
    "Borderline": "borderline",
    "RRDR": "RRDR rule",
    "Evidence from ALL dataset only": "ALL dataset only",
}


def column_index(cell_ref):
    """Convert an Excel cell reference such as 'AB12' to a 0-based column."""
    letters = re.match(r"[A-Z]+", cell_ref).group(0)
    index = 0
    for letter in letters:
        index = index * 26 + (ord(letter) - ord("A") + 1)
    return index - 1


def sheet_path(book, name):
    """Find the XML file inside the .xlsx that holds the named sheet."""
    workbook = ET.fromstring(book.read("xl/workbook.xml"))
    rels = ET.fromstring(book.read("xl/_rels/workbook.xml.rels"))
    targets = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
    for sheet in workbook.iter(NS + "sheet"):
        if sheet.attrib.get("name") == name:
            target = targets[sheet.attrib[REL_NS + "id"]].lstrip("/")
            return target if target.startswith("xl/") else "xl/" + target
    sys.exit(f"Error: sheet '{name}' not found")


def read_rows(xlsx_path):
    """Yield each row of the catalogue sheet as a list of strings."""
    with zipfile.ZipFile(xlsx_path) as book:
        shared = []
        with book.open("xl/sharedStrings.xml") as f:
            for _, el in ET.iterparse(f):
                if el.tag == NS + "si":
                    shared.append("".join(t.text or "" for t in el.iter(NS + "t")))
                    el.clear()

        with book.open(sheet_path(book, SHEET_NAME)) as f:
            for _, el in ET.iterparse(f):
                if el.tag != NS + "row":
                    continue
                cells = {}
                for cell in el.findall(NS + "c"):
                    value = cell.find(NS + "v")
                    if cell.attrib.get("t") == "s" and value is not None:
                        cells[column_index(cell.attrib["r"])] = shared[int(value.text)]
                    elif value is not None:
                        cells[column_index(cell.attrib["r"])] = value.text
                el.clear()
                width = max(cells) + 1 if cells else 0
                yield [cells.get(i, "") for i in range(width)]


def main():
    if len(sys.argv) != 2:
        sys.exit("Usage: python3 tools/extract_who_rpoB.py WHO-UCN-TB-2023.7-eng.xlsx")

    header = None
    entries = []
    for row in read_rows(sys.argv[1]):
        if header is None:
            # The real header is the row whose first cell is 'drug' (row 3)
            if row and row[0] == "drug":
                header = {name: i for i, name in enumerate(row)}
            continue

        def get(name):
            i = header[name]
            return row[i] if i < len(row) else ""

        if get("drug") != "Rifampicin" or get("gene") != "rpoB":
            continue
        if get("effect") != "missense_variant":
            continue
        grading = get("FINAL CONFIDENCE GRADING")
        if not grading.startswith(("1)", "2)")):
            continue

        match = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", get("mutation"))
        if match is None:
            sys.exit(f"Error: unexpected mutation format '{get('mutation')}'")
        ref_aa = THREE_TO_ONE[match.group(1)]
        codon = int(match.group(2))
        alt_aa = THREE_TO_ONE[match.group(3)]

        criterion = get("Additional grading criteria applied")
        if criterion not in CRITERIA_NOTES:
            sys.exit(f"Error: unknown grading criterion '{criterion}'")
        note = CRITERIA_NOTES[criterion]
        source = "WHO 2023" + ("; " + note if note else "")

        entries.append((codon, int(grading[0]), alt_aa, ref_aa, source))

    if header is None or not entries:
        sys.exit("Error: no rpoB rifampicin entries found - is this the WHO catalogue?")

    entries.sort()
    tier1 = sum(1 for e in entries if e[1] == 1)

    print("# rpoB rifampicin resistance mutation database")
    print("# kewanessa")
    print("# CS50 Final Project")
    print("# M. tuberculosis numbering (H37Rv, Rv0667)")
    print("# Columns: gene\tcodon\tref_aa\talt_aa\tdrug\tconfidence\tsource")
    print("# Confidence = WHO 2023 grading group: 1=Associated with resistance,")
    print("#   2=Associated with resistance - interim, 3=Uncertain significance,")
    print("#   4=Not associated with resistance - interim, 5=Not associated with resistance")
    print("# Source: WHO. Catalogue of mutations in Mycobacterium tuberculosis complex and")
    print("#   their association with drug resistance, 2nd ed. Geneva: WHO; 2023.")
    print("#   ISBN 9789240082410. Data licence: ODC-By 1.0.")
    print("# Extracted from WHO-UCN-TB-2023.7-eng.xlsx (sheet Catalogue_master_file),")
    print("#   https://github.com/GTB-tbsequencing/mutation-catalogue-2023")
    print("#   by tools/extract_who_rpoB.py. Filter: drug=Rifampicin, gene=rpoB,")
    print("#   effect=missense_variant, final grading group 1 or 2.")
    print(f"#   Result: {len(entries)} mutations ({tier1} in group 1, {len(entries) - tier1} in group 2).")
    print("#   WHO-graded in-frame indels are not single-codon substitutions and are")
    print("#   not listed; amr_scan reports every RRDR indel separately.")
    print("# Notes in the source column: borderline = WHO 'borderline' mutation;")
    print("#   RRDR rule = graded by WHO's RRDR additional grading criterion.")
    print("# IMPORTANT: This uses M. tuberculosis numbering, NOT E. coli numbering.")
    print("# E. coli S531 = M. tuberculosis S450, H526 = H445, D516 = D435, etc.")
    print("# Ref AAs verified against H37Rv protein NP_215181.1 (1172 aa)")
    print("gene\tcodon\tref_aa\talt_aa\tdrug\tconfidence\tsource")
    for codon, tier, alt_aa, ref_aa, source in entries:
        print(f"rpoB\t{codon}\t{ref_aa}\t{alt_aa}\trifampicin\t{tier}\t{source}")


if __name__ == "__main__":
    main()
