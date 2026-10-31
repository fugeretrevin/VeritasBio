#!/usr/bin/env python3
"""Rule out the D-8 class of defect across the WHOLE frozen corpus, not just CR9114.

THE DEFECT THIS SCRIPT HUNTS. CombinGym's clean CR9114 workbook applies an undocumented
filter -- `NOT(site 1 AND site 16) AND h1_mean is present` -- that drops 16,384 genotypes,
an entire quadrant of the 2^16 design space, regardless of affinity. It is mechanically
depth-correlated (a depth-k genotype carries two specified mutations with probability
k(k-1)/240), so it removed the fully matured antibody and moved a pre-registered covariate
by 6.0 % relative. See docs/cr9114-reprep.md.

That was found in the one landscape anyone happened to examine. The other nine clean
workbooks had never been diffed against their raw counterparts. Every one of them HAS a raw
counterpart at the same pinned commit `c5dccd51`, so ruling the defect in or out is cheap.

WHAT IS COMPUTED, per landscape:
  1. raw row count, clean row count, and the SET DIFFERENCE over canonical genotypes
     (both directions -- a clean row absent from raw would be just as alarming);
  2. for every missing row, WHY: absent measurement, a stated floor, a quality filter, or
     an undocumented structural exclusion;
  3. whether the attrition is DEPTH-CORRELATED, against the EXACT design-space expectation
     for that library -- C(n,k) for a binary library, and for the mixed alphabets here the
     coefficient of x^k in prod_i (1 + (a_i - 1) x), computed properly;
  4. where rows are missing, a SWEEP over all single sites and all site pairs, reporting
     false positives and false negatives, so a structural exclusion is either named to zero
     error or ruled out with numbers rather than with a story;
  5. whether the landscape's pre-registered covariate moves, by recomputing it over the
     raw universe with identical code, plus the depth-blind reweighting bound as a check.

PARSER. All genotypes go through scripts/analysis/genotype_canon.py, which expands mutation
strings to full library length filling wild type at omitted positions and ASSERTS equal
length. That module's validation gate (TrpB 168/2,242 = 7.4933 %, GB1 72/2,167 = 3.3226 %)
runs first here and this script refuses to continue if it fails.

Reproduce:
    python scripts/analysis/corpus_raw_fetch.py fetch
    python scripts/analysis/corpus_attrition.py
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))
sys.path.insert(0, str(REPO_ROOT / "scripts" / "analysis"))

from freeze_cross_protein_corpus import read_xlsx  # noqa: E402
from genotype_canon import (  # noqa: E402
    alphabet_design_space, canonicalise, depth_expectation, depth_of,
    reweighting_bound, shell_stats, validate,
)

CORPUS = REPO_ROOT / "downloads" / "combingym"

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}


def three_letter_to_tokens(mutant: str) -> str:
    """`p.Thr7Ser` -> `T7S`; `p.Thr7Ser;p.Gly3Glu` -> `T7S:G3E`. CreiLOV uses this form."""
    s = (mutant or "").strip()
    if not s or s.lower() == "wt":
        return "WT"
    out = []
    for part in s.replace(",", ";").split(";"):
        part = part.strip()
        if not part:
            continue
        if part.startswith("p."):
            part = part[2:]
        wt3, rest = part[:3], part[3:]
        digits = "".join(c for c in rest if c.isdigit())
        mut3 = rest[len(digits):]
        if wt3 not in THREE_TO_ONE or mut3 not in THREE_TO_ONE or not digits:
            raise ValueError(f"cannot convert three-letter mutation {mutant!r} (part {part!r})")
        out.append(f"{THREE_TO_ONE[wt3]}{digits}{THREE_TO_ONE[mut3]}")
    return ":".join(out)


# ==================================================================================
# Table readers
# ==================================================================================
def _num(cell: Any) -> Optional[float]:
    if cell is None:
        return None
    s = str(cell).strip()
    if s == "" or s.upper() == "NA" or s.upper() == "NAN":
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v if v == v else None  # reject NaN


def read_csv_rows(path: Path, delimiter: str = ",") -> List[Dict[str, str]]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def xlsx_rows(path: Path, sheet: Optional[str]) -> Tuple[List[str], List[List[Optional[str]]]]:
    _name, header, rows = read_xlsx(path, sheet)
    return header, rows


def col(header: Sequence[str], name: str) -> int:
    try:
        return list(header).index(name)
    except ValueError as exc:
        raise KeyError(f"{name!r} not in header {list(header)}") from exc


def cellat(row: Sequence[Any], i: int) -> Any:
    return row[i] if 0 <= i < len(row) else None


# ==================================================================================
# A Table is the audit's common currency: canonical key -> value (or None if absent)
# ==================================================================================
class Table:
    def __init__(self, label: str, entries: Dict[str, Optional[float]],
                 extra: Optional[Dict[str, Dict[str, Any]]] = None,
                 rows_in_sheet: Optional[int] = None,
                 blank_key_rows: int = 0, duplicate_keys: int = 0):
        self.label = label
        self.entries = entries
        self.extra = extra or {}
        self.rows_in_sheet = rows_in_sheet if rows_in_sheet is not None else len(entries)
        self.blank_key_rows = blank_key_rows
        self.duplicate_keys = duplicate_keys

    @property
    def keys(self):
        return set(self.entries)

    @property
    def measured(self):
        return {k for k, v in self.entries.items() if v is not None}

    def __len__(self) -> int:
        return len(self.entries)


def build_table(label: str, pairs: Sequence[Tuple[str, Optional[float]]],
                rows_in_sheet: int, blank_key_rows: int,
                extra: Optional[Dict[str, Dict[str, Any]]] = None) -> Table:
    entries: Dict[str, Optional[float]] = {}
    dup = 0
    for key, value in pairs:
        if key in entries:
            dup += 1
            continue
        entries[key] = value
    return Table(label, entries, extra, rows_in_sheet, blank_key_rows, dup)


# ==================================================================================
# Landscape adapters
# ==================================================================================
# Each adapter returns (clean Table, raw Table, library) where library is
# (positions-or-site-labels, wild-type canonical key).
Library = Tuple[List[str], str]


def _sequence_library(seqs: Sequence[str], wt_row: str) -> Tuple[List[int], str]:
    """Varying positions and parent residues, derived from the clean sequences themselves."""
    width = len(wt_row)
    varying = [i for i in range(width) if len({s[i] for s in seqs}) > 1]
    return [i + 1 for i in varying], "".join(wt_row[i] for i in varying)


def adapt_bnab(protein: str, raw_rel: str, clean_rel: str, sheet: str,
               value_col: str, sites: int, wt_seq: str) -> Tuple[Table, Table, Library]:
    """CR9114 / CR6261. Raw keys on an n-bit site indicator; clean keys on the sequence.

    The raw `genotype` column IS the canonical per-site key for a binary library, which is
    cleaner than parsing residue sequences -- `som_mut` is its popcount, asserted here.
    """
    raw_rows = read_csv_rows(CORPUS / raw_rel)
    raw_pairs: List[Tuple[str, Optional[float]]] = []
    raw_extra: Dict[str, Dict[str, Any]] = {}
    for r in raw_rows:
        g = (r["genotype"] or "").strip()
        assert len(g) == sites, f"{protein}: raw genotype {g!r} is not {sites} bits"
        bits = "".join((r[f"pos{i}"] or "").strip() for i in range(1, sites + 1))
        assert bits == g, f"{protein}: pos1..pos{sites} {bits!r} != genotype {g!r}"
        assert int(r["som_mut"]) == g.count("1"), f"{protein}: som_mut != popcount for {g!r}"
        raw_pairs.append((g, _num(r.get(value_col))))
        raw_extra[g] = r
    raw = build_table(f"{protein} raw", raw_pairs, len(raw_rows), 0, raw_extra)

    header, rows = xlsx_rows(CORPUS / clean_rel, sheet)
    gi, vi = col(header, "genotype"), col(header, value_col)
    ni = col(header, "n_mut")
    seqs = [r[gi] for r in rows if cellat(r, gi)]
    positions, parent = _sequence_library(seqs, wt_seq)
    assert len(positions) == sites, (
        f"{protein}: clean sequences vary at {len(positions)} positions, expected {sites}")
    clean_pairs: List[Tuple[str, Optional[float]]] = []
    blank = 0
    for r in rows:
        s = cellat(r, gi)
        if not s:
            blank += 1
            continue
        key = "".join("1" if s[p - 1] != wt_seq[p - 1] else "0" for p in positions)
        assert len(key) == sites
        nm = _num(cellat(r, ni))
        assert nm is None or int(nm) == key.count("1"), f"{protein}: n_mut != depth"
        clean_pairs.append((key, _num(cellat(r, vi))))
    clean = build_table(f"{protein} clean", clean_pairs, len(rows), blank)
    return clean, raw, ([str(p) for p in positions], "0" * sites)


def adapt_creilov() -> Tuple[Table, Table, Library]:
    header, rows = xlsx_rows(CORPUS / "Data/DMS/Clean/CreiLOV_clean.xlsx", "Sheet1")
    mi, gi, vi = col(header, "mutant"), col(header, "genotype"), col(header, "mean")
    ni = col(header, "n_mut")
    seqs = [r[gi] for r in rows if cellat(r, gi)]
    wt_row = next(r[gi] for r in rows if _num(cellat(r, ni)) == 0)
    positions, parent = _sequence_library(seqs, wt_row)
    # The mutant strings are three-letter and 1-based into the same string; cross-check.
    clean_pairs: List[Tuple[str, Optional[float]]] = []
    blank = 0
    for r in rows:
        s = cellat(r, gi)
        if not s:
            blank += 1
            continue
        key = "".join(s[p - 1] for p in positions)
        from_mut = canonicalise(three_letter_to_tokens(cellat(r, mi) or "WT"),
                                positions, parent)
        assert key == from_mut, (
            f"CreiLOV: sequence key {key!r} disagrees with mutant string "
            f"{cellat(r, mi)!r} -> {from_mut!r}")
        assert depth_of(key, parent) == int(_num(cellat(r, ni))), "CreiLOV: n_mut != depth"
        clean_pairs.append((key, _num(cellat(r, vi))))
    clean = build_table("CreiLOV clean", clean_pairs, len(rows), blank)

    rheader, rrows = xlsx_rows(CORPUS / "Data/DMS/Raw/CreiLOV_raw.xlsx", "Sheet1")
    rvi = col(rheader, "mean")
    raw_pairs = []
    rblank = 0
    for r in rrows:
        m = cellat(r, 0)
        if not m:
            rblank += 1
            continue
        raw_pairs.append((canonicalise(three_letter_to_tokens(m), positions, parent),
                          _num(cellat(r, rvi))))
    raw = build_table("CreiLOV raw", raw_pairs, len(rrows), rblank)
    return clean, raw, ([str(p) for p in positions], parent)


def adapt_eqfp611(value_col: str) -> Tuple[Table, Table, Library]:
    header, rows = xlsx_rows(CORPUS / "Data/DMS/Clean/eqFP611_clean.xlsx", "genodata")
    mi, vi, ni = col(header, "mutant"), col(header, value_col), col(header, "n_mut")
    wt_key = next(cellat(r, mi) for r in rows if _num(cellat(r, ni)) == 0)
    width = len(wt_key)
    clean_pairs, blank = [], 0
    for r in rows:
        k = cellat(r, mi)
        if not k:
            blank += 1
            continue
        assert len(k) == width, f"eqFP611: mutant {k!r} is not {width} characters"
        assert depth_of(k, wt_key) == int(_num(cellat(r, ni))), "eqFP611: n_mut != depth"
        clean_pairs.append((k, _num(cellat(r, vi))))
    clean = build_table(f"eqFP611 {value_col} clean", clean_pairs, len(rows), blank)

    # Raw has a two-row header: ['binary','amino acid',...] then
    # ['genotype','sequence','input','red','blue',None,'red','blue',None,'combined'].
    rheader, rrows = xlsx_rows(CORPUS / "Data/DMS/Raw/eqFP611_raw.xlsx", "genodata")
    inner = rrows[0]
    assert [c for c in inner if c] == ["genotype", "sequence", "input", "red", "blue",
                                       "red", "blue", "combined"], f"unexpected {inner}"
    # Columns 3,4 are counts; the brightness block is 6 (red), 7 (blue), 9 (combined).
    bright = {"red": 6, "blue": 7, "combined": 9}[value_col]
    raw_pairs, rblank = [], 0
    for r in rrows[1:]:
        k = (cellat(r, 1) or "").strip()
        if not k:
            rblank += 1
            continue
        assert len(k) == width
        raw_pairs.append((k, _num(cellat(r, bright))))
    raw = build_table(f"eqFP611 {value_col} raw", raw_pairs, len(rrows) - 1, rblank)
    return clean, raw, ([f"s{i+1}" for i in range(width)], wt_key)


def adapt_gb1() -> Tuple[Table, Table, Library]:
    header, rows = xlsx_rows(CORPUS / "Data/DMS/Clean/GB1_clean.xlsx", "GB1_clean")
    gi, vi, ni = col(header, "genotype"), col(header, "fitness"), col(header, "n_mut")
    seqs = [r[gi] for r in rows if cellat(r, gi)]
    wt_row = next(r[gi] for r in rows if _num(cellat(r, ni)) == 0)
    positions, parent = _sequence_library(seqs, wt_row)
    clean_pairs, blank = [], 0
    for r in rows:
        s = cellat(r, gi)
        if not s:
            blank += 1
            continue
        key = "".join(s[p - 1] for p in positions)
        assert len(key) == len(positions)
        clean_pairs.append((key, _num(cellat(r, vi))))
    clean = build_table("GB1 clean", clean_pairs, len(rows), blank)

    rheader, rrows = xlsx_rows(CORPUS / "Data/DMS/Raw/GB1_raw.xlsx", "Supplementary Dataset 1")
    rki, rvi = col(rheader, "Variants"), col(rheader, "Fitness")
    raw_pairs, rblank = [], 0
    for r in rrows:
        k = (cellat(r, rki) or "").strip()
        if not k:
            rblank += 1
            continue
        assert len(k) == len(positions), f"GB1: raw variant {k!r} is not 4 characters"
        raw_pairs.append((k, _num(cellat(r, rvi))))
    raw = build_table("GB1 raw", raw_pairs, len(rrows), rblank)
    return clean, raw, ([str(p) for p in positions], parent)


def adapt_hiv1() -> Tuple[Table, Table, Library]:
    header, rows = xlsx_rows(CORPUS / "Data/DMS/Clean/HIV-1_clean.xlsx", "HIV-1")
    gi, vi, ni = col(header, "genotype"), col(header, "fitness"), col(header, "n_mut")
    seqs = [r[gi] for r in rows if cellat(r, gi)]
    wt_row = (CORPUS / "Data/WT Sequence/HIV-1.fasta").read_text(encoding="utf-8")
    wt_row = "".join(l.strip() for l in wt_row.splitlines() if not l.startswith(">"))
    positions, parent = _sequence_library(seqs, wt_row)
    clean_pairs, blank = [], 0
    for r in rows:
        s = cellat(r, gi)
        if not s:
            blank += 1
            continue
        key = "".join(s[p - 1] for p in positions)
        assert depth_of(key, parent) == int(_num(cellat(r, ni))), "HIV-1: n_mut != depth"
        clean_pairs.append((key, _num(cellat(r, vi))))
    clean = build_table("HIV-1 clean", clean_pairs, len(rows), blank)

    rrows = read_csv_rows(CORPUS / "Data/DMS/Raw/HIV-1_raw.tsv", delimiter="\t")
    raw_pairs, rblank = [], 0
    for r in rrows:
        m = (r.get("mutation") or "").strip()
        if not m:
            rblank += 1
            continue
        # Mutation strings omit wild-type positions -- canonicalise, never zip.
        raw_pairs.append((canonicalise(m, positions, parent),
                          _num(r.get("relative_fitness"))))
    raw = build_table("HIV-1 raw", raw_pairs, len(rrows), rblank)
    return clean, raw, ([str(p) for p in positions], parent)


def adapt_rhla(value_col: str) -> Tuple[Table, Table, Library]:
    header, rows = xlsx_rows(CORPUS / "Data/DMS/Clean/RhlA_clean.xlsx", "Sheet1")
    gi, vi, ni = col(header, "genotype"), col(header, value_col), col(header, "num_mut")
    wt_key = next(cellat(r, gi) for r in rows if _num(cellat(r, ni)) == 0)
    width = len(wt_key)
    clean_pairs, blank = [], 0
    for r in rows:
        k = cellat(r, gi)
        if not k:
            blank += 1
            continue
        assert len(k) == width, f"RhlA: genotype {k!r} is not {width} characters"
        assert depth_of(k, wt_key) == int(_num(cellat(r, ni))), "RhlA: num_mut != depth"
        clean_pairs.append((k, _num(cellat(r, vi))))
    clean = build_table(f"RhlA {value_col} clean", clean_pairs, len(rows), blank)

    rheader, rrows = xlsx_rows(CORPUS / "Data/DMS/Raw/RhlA_raw.xlsx", "Sheet1")
    rki, rvi = col(rheader, "Sequence"), col(rheader, value_col)
    raw_pairs, rblank = [], 0
    for r in rrows:
        k = cellat(r, rki)
        if not k:
            rblank += 1
            continue
        assert len(k) == width
        raw_pairs.append((k, _num(cellat(r, rvi))))
    raw = build_table(f"RhlA {value_col} raw", raw_pairs, len(rrows), rblank)
    return clean, raw, ([f"s{i+1}" for i in range(width)], wt_key)


SACAS9_POS = [887, 888, 889, 985, 986, 988, 989, 991]
SPCAS9_POS = [661, 695, 848, 923, 924, 926, 1003, 1060]


def adapt_sacas9() -> Tuple[Table, Table, Library]:
    header, rows = xlsx_rows(CORPUS / "Data/DMS/Clean/SaCas9_clean.xlsx", "Sheet1")
    gi, vi, ni = col(header, "genotype"), col(header, "Mean"), col(header, "n_mut")
    wt_row = next(cellat(r, gi) for r in rows if _num(cellat(r, ni)) == 0)
    parent = "".join(wt_row[p - 1] for p in SACAS9_POS)
    clean_pairs, blank = [], 0
    for r in rows:
        s = cellat(r, gi)
        if not s:
            blank += 1
            continue
        key = "".join(s[p - 1] for p in SACAS9_POS)
        assert depth_of(key, parent) == int(_num(cellat(r, ni))), "SaCas9: n_mut != depth"
        clean_pairs.append((key, _num(cellat(r, vi))))
    clean = build_table("SaCas9 clean", clean_pairs, len(rows), blank)

    _rh, rrows = xlsx_rows(CORPUS / "Data/DMS/Raw/SaCas9_raw.xlsx", "Sheet1")
    hdr_i = next(i for i, r in enumerate(rrows) if r and str(cellat(r, 0)).strip() == "No.")
    inner = rrows[hdr_i]
    assert str(cellat(inner, 1)).strip() == "-".join(str(p) for p in SACAS9_POS), (
        f"SaCas9: raw position header {cellat(inner, 1)!r} does not name {SACAS9_POS}")
    # The clean `Mean` is over sg1..sg3 only; the raw table carries six sgRNAs. Use the
    # same three so the comparison is of row membership, not of phenotype definition.
    sg_cols = [2, 3, 4]
    raw_pairs, rblank = [], 0
    for r in rrows[hdr_i + 1:]:
        if not r or cellat(r, 0) in (None, ""):
            continue
        key = (cellat(r, 1) or "").replace("-", "")
        if not key:
            rblank += 1
            continue
        assert len(key) == len(SACAS9_POS)
        vals = [_num(cellat(r, j)) for j in sg_cols]
        value = sum(vals) / len(vals) if all(v is not None for v in vals) else None
        raw_pairs.append((key, value))
    raw = build_table("SaCas9 raw", raw_pairs, len(rrows) - hdr_i - 1, rblank)
    return clean, raw, ([str(p) for p in SACAS9_POS], parent)


def adapt_spcas9() -> Tuple[Table, Table, Library]:
    header, rows = xlsx_rows(CORPUS / "Data/DMS/Clean/SpCas9_clean.xlsx", "SpCas9")
    gi, vi, ni = col(header, "genotype"), col(header, "mean"), col(header, "n_mut")
    wt_row = next(cellat(r, gi) for r in rows if _num(cellat(r, ni)) == 0)
    parent = "".join(wt_row[p - 1] for p in SPCAS9_POS)
    clean_pairs, blank = [], 0
    for r in rows:
        s = cellat(r, gi)
        if not s:
            blank += 1  # the sheet is padded with 5,151 wholly empty rows
            continue
        key = "".join(s[p - 1] for p in SPCAS9_POS)
        assert depth_of(key, parent) == int(_num(cellat(r, ni))), "SpCas9: n_mut != depth"
        clean_pairs.append((key, _num(cellat(r, vi))))
    clean = build_table("SpCas9 clean", clean_pairs, len(rows), blank)

    _rh, rrows = xlsx_rows(CORPUS / "Data/DMS/Raw/SpCas9_raw.xlsx", "Enrichment Score")
    hdr_i = next(i for i, r in enumerate(rrows)
                 if r and str(cellat(r, 0)).strip() == "Cas9 variant #")
    inner = rrows[hdr_i]
    assert [str(cellat(inner, j)) for j in range(1, 9)] == [str(p) for p in SPCAS9_POS], (
        f"SpCas9: raw header does not name {SPCAS9_POS}")
    on5 = col(inner, "RFPsg5 ON")
    on8 = col(inner, "RFPsg8 ON")
    raw_pairs, rblank = [], 0
    raw_extra: Dict[str, Dict[str, Any]] = {}
    for r in rrows[hdr_i + 1:]:
        if not r or cellat(r, 0) in (None, ""):
            continue
        key = "".join(str(cellat(r, j) or "?") for j in range(1, 9))
        if "?" in key:
            rblank += 1
            continue
        a, b = _num(cellat(r, on5)), _num(cellat(r, on8))
        value = (a + b) / 2 if a is not None and b is not None else None
        raw_pairs.append((key, value))
        raw_extra[key] = {"RFPsg5 ON": cellat(r, on5), "RFPsg8 ON": cellat(r, on8)}
    raw = build_table("SpCas9 raw", raw_pairs, len(rrows) - hdr_i - 1, rblank, raw_extra)
    return clean, raw, ([str(p) for p in SPCAS9_POS], parent)


CR9114_WT = ("QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVT"
             "ITADKSTSTAYMELSSLRSEDTAVYYCARHGNYYYYYGMDVWGQGTTVTVSS")
CR6261_WT = ("QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVT"
             "ITADKSTSTAYMELSSLRSEDTAVYYCAKHMGYQLRETMDVWGQGTTVTVSS")

# key -> (adapter, sites, primary depth cut, pre-registered enrichment or None, level)
LANDSCAPES: Dict[str, Dict[str, Any]] = {
    "creilov": dict(adapt=adapt_creilov, sites=15, prereg=1.875, level="final_test"),
    "eqfp611_blue": dict(adapt=lambda: adapt_eqfp611("blue"), sites=13, prereg=1.854,
                         level="final_test"),
    "eqfp611_red": dict(adapt=lambda: adapt_eqfp611("red"), sites=13, prereg=0.634,
                        level="final_test"),
    "eqfp611_combined": dict(adapt=lambda: adapt_eqfp611("combined"), sites=13,
                             prereg=1.830, level="final_test"),
    "cr9114_h1": dict(
        adapt=lambda: adapt_bnab("CR9114", "Data/DMS/Raw/bnAbs_CR9114_raw.csv",
                                 "Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx",
                                 "CR91114-h1", "h1_mean", 16, CR9114_WT),
        sites=16, prereg=0.584, level="final_test"),
    "cr6261_h1": dict(
        adapt=lambda: adapt_bnab("CR6261-H1", "Data/DMS/Raw/bnAbs_CR6261_raw.csv",
                                 "Data/DMS/Clean/bnAbs_CR6261_H1_clean.xlsx",
                                 "CR6261-h1", "h1_mean", 11, CR6261_WT),
        sites=11, prereg=0.210, level="final_test"),
    "cr6261_h9": dict(
        adapt=lambda: adapt_bnab("CR6261-H9", "Data/DMS/Raw/bnAbs_CR6261_raw.csv",
                                 "Data/DMS/Clean/bnAbs_CR6261_H9_clean.xlsx",
                                 "CR6261", "h9_mean", 11, CR6261_WT),
        sites=11, prereg=0.000, level="final_test"),
    "hiv1_protease": dict(adapt=adapt_hiv1, sites=11, prereg=None, level="reserve"),
    "rhla_rls": dict(adapt=lambda: adapt_rhla("RLs"), sites=11, prereg=None,
                     level="reserve"),
    "rhla_pct_rha": dict(adapt=lambda: adapt_rhla("%Rha-(C8-C10)"), sites=11, prereg=None,
                         level="reserve"),
    "rhla_rha": dict(adapt=lambda: adapt_rhla("Rha-(C8-C10)"), sites=11, prereg=None,
                     level="reserve"),
    "sacas9": dict(adapt=adapt_sacas9, sites=8, prereg=None, level="reserve"),
    "spcas9": dict(adapt=adapt_spcas9, sites=8, prereg=None, level="reserve"),
    "gb1_combingym": dict(adapt=adapt_gb1, sites=4, prereg=None, level="development"),
}


# ==================================================================================
# The sweep: name a structural exclusion to zero error, or rule it out with numbers
# ==================================================================================
def score_composite(rule_sites: Sequence[int], dropped_measured: Sequence[str],
                    dropped_unmeasured: Sequence[str], kept: Sequence[str],
                    kept_unmeasured_count: int, parent: str) -> Dict[str, Any]:
    """Score the COMPOSITE rule, which is the form CR9114's filter actually takes:

        NOT(mutated at every site in rule_sites)  AND  the measurement is present

    Scored as a predictor of EXCLUSION over the whole raw table. This matters because the
    structural half alone cannot explain the unquantified rows outside the quadrant (on
    CR9114, 311 of them), so scoring only the structural half understates the rule. The
    composite is what must reach 0 FP and 0 FN before a filter is called named.
    """
    def hits_rule(g: str) -> bool:
        return all(g[i] != parent[i] for i in rule_sites)

    # False negative: a dropped genotype the composite rule fails to exclude. A dropped
    # row is excluded by the composite if it hits the structural rule OR is unmeasured.
    fn = sum(1 for g in dropped_measured if not hits_rule(g))
    fn += sum(1 for g in dropped_unmeasured if False)  # unmeasured are always excluded
    # False positive: a kept genotype the composite rule would have excluded -- either it
    # hits the structural rule, or it was kept despite having no measurement.
    fp = sum(1 for g in kept if hits_rule(g)) + kept_unmeasured_count
    return {"rule": "NOT(mutated at sites " + ", ".join(str(i) for i in rule_sites)
                    + ") AND the measurement is present",
            "fp": fp, "fn": fn, "symdiff": fp + fn, "exact": fp + fn == 0}


def sweep_sites(dropped: Sequence[str], kept: Sequence[str], parent: str,
                site_labels: Sequence[str]) -> List[Dict[str, Any]]:
    """Score every rule `mutated at site i` and `mutated at BOTH site i and site j`.

    This is the CR9114 sweep generalised to mixed alphabets: "mutated" means any residue
    other than the parent's at that site, which for a binary library is exactly bit = 1.
    A rule is scored as a predictor of EXCLUSION:
        false positive = a KEPT genotype the rule would have excluded
        false negative = a DROPPED genotype the rule fails to exclude
    Symmetric difference FP + FN == 0 means the rule reproduces the exclusion exactly.
    """
    n = len(parent)
    drop_mut = [[g[i] != parent[i] for i in range(n)] for g in dropped]
    keep_mut = [[g[i] != parent[i] for i in range(n)] for g in kept]
    out: List[Dict[str, Any]] = []
    for i in range(n):
        fn = sum(1 for m in drop_mut if not m[i])
        fp = sum(1 for m in keep_mut if m[i])
        out.append({"rule": f"mutated at site {site_labels[i]}", "fp": fp, "fn": fn,
                    "symdiff": fp + fn})
    for i in range(n):
        for j in range(i + 1, n):
            fn = sum(1 for m in drop_mut if not (m[i] and m[j]))
            fp = sum(1 for m in keep_mut if m[i] and m[j])
            out.append({"rule": f"mutated at BOTH sites {site_labels[i]} and "
                                f"{site_labels[j]}", "fp": fp, "fn": fn,
                        "symdiff": fp + fn})
    out.sort(key=lambda d: d["symdiff"])
    return out


# ==================================================================================
# Per-landscape audit
# ==================================================================================
def audit(key: str, spec: Dict[str, Any]) -> Dict[str, Any]:
    clean, raw, (site_labels, parent) = spec["adapt"]()
    sites = spec["sites"]
    cut = sites // 2
    out: Dict[str, Any] = {
        "landscape": key, "level": spec["level"], "sites": sites, "depth_cut": cut,
        "prereg_enrichment": spec["prereg"],
        "clean_rows_in_sheet": clean.rows_in_sheet,
        "clean_blank_key_rows": clean.blank_key_rows,
        "clean_keys": len(clean), "clean_duplicate_keys": clean.duplicate_keys,
        "raw_rows_in_sheet": raw.rows_in_sheet, "raw_blank_key_rows": raw.blank_key_rows,
        "raw_keys": len(raw), "raw_duplicate_keys": raw.duplicate_keys,
        "clean_measured": len(clean.measured), "raw_measured": len(raw.measured),
    }

    # ---- design space, exactly ----
    all_keys = sorted(raw.keys | clean.keys)
    space, alpha = alphabet_design_space(all_keys)
    expected = depth_expectation(alpha)
    out["alphabet_sizes"] = alpha
    out["design_space"] = space
    out["binary_library"] = all(a == 2 for a in alpha)
    out["raw_saturation"] = len(raw) / space if space else None

    # ---- the set difference, both directions ----
    in_raw_not_clean = sorted(raw.keys - clean.keys)
    in_clean_not_raw = sorted(clean.keys - raw.keys)
    out["in_raw_not_clean"] = len(in_raw_not_clean)
    out["in_clean_not_raw"] = len(in_clean_not_raw)

    # ---- WHY are they missing? ----
    unmeasured = sum(1 for k in in_raw_not_clean if raw.entries[k] is None)
    measured_dropped = [k for k in in_raw_not_clean if raw.entries[k] is not None]
    out["missing_with_no_measurement"] = unmeasured
    out["missing_WITH_a_measurement"] = len(measured_dropped)
    # The decisive quantity. A missing row with NO measurement could never have entered the
    # eligible universe under any convention, so it cannot move a covariate however
    # depth-correlated its absence happens to look. Only a RECOVERABLE row -- missing from
    # clean yet carrying a real measurement in raw, which is what the CR9114 quadrant was --
    # can move anything.
    out["recoverable_rows"] = len(measured_dropped)
    out["covariate_can_move"] = len(measured_dropped) > 0
    # The named rule, scored: "present in clean iff the measurement is present in raw".
    rule_fp = sum(1 for k in (raw.keys & clean.keys) if raw.entries[k] is None)
    rule_fn = len(measured_dropped)
    out["measurement_presence_rule"] = {
        "statement": "kept in clean iff the raw measurement is present",
        "false_positives_kept_but_unmeasured": rule_fp,
        "false_negatives_dropped_but_measured": rule_fn,
        "symdiff": rule_fp + rule_fn,
        "exact": rule_fp + rule_fn == 0,
    }

    # ---- depth-correlation of the clean/raw attrition ----
    def hist(keys) -> Dict[int, int]:
        h: Dict[int, int] = {}
        for k in keys:
            d = depth_of(k, parent)
            h[d] = h.get(d, 0) + 1
        return h

    raw_h, clean_h = hist(raw.keys), hist(clean.keys)
    retention: List[Tuple[int, int, int, int, Optional[float]]] = []
    for d in range(len(expected)):
        e = expected[d] if d < len(expected) else 0
        rr, cc = raw_h.get(d, 0), clean_h.get(d, 0)
        retention.append((d, e, rr, cc, (cc / rr) if rr else None))
    out["retention"] = retention
    big = [r for r in retention if r[2] >= 10 and r[4] is not None]
    if big:
        lo, hi = min(r[4] for r in big), max(r[4] for r in big)
        out["clean_vs_raw_retention_range"] = [lo, hi]
        out["clean_vs_raw_spread_pp"] = (hi - lo) * 100
        out["depth_correlated"] = (hi - lo) > 0.15
    else:
        out["clean_vs_raw_retention_range"] = None
        out["depth_correlated"] = False
    # Separately: is the RAW table itself depth-biased against the design space?
    rawbig = [(d, e, rr) for d, e, rr, _c, _r in retention if e >= 10]
    if rawbig:
        rr_ret = [rr / e for _d, e, rr in rawbig]
        out["raw_vs_space_retention_range"] = [min(rr_ret), max(rr_ret)]
        out["raw_vs_space_spread_pp"] = (max(rr_ret) - min(rr_ret)) * 100

    # ---- the sweep, only where rows are actually missing ----
    out["sweep"] = None
    out["composite"] = None
    if in_raw_not_clean:
        kept_shared = sorted(clean.keys & raw.keys)
        results = sweep_sites(in_raw_not_clean, kept_shared, parent, site_labels)
        out["sweep"] = {"n_rules": len(results), "best": results[:4],
                        "any_exact": any(r["symdiff"] == 0 for r in results)}
        # Now score the COMPOSITE of each candidate's site set with measurement-presence.
        # `rule_fp` above is the count of kept-but-unmeasured rows, which the composite
        # rule would wrongly exclude.
        dropped_measured = measured_dropped
        dropped_unmeasured = [k for k in in_raw_not_clean if raw.entries[k] is None]
        composites: List[Dict[str, Any]] = []
        n = len(parent)
        candidates: List[List[int]] = [[i] for i in range(n)]
        candidates += [[i, j] for i in range(n) for j in range(i + 1, n)]
        for sites_idx in candidates:
            c = score_composite(sites_idx, dropped_measured, dropped_unmeasured,
                                kept_shared, rule_fp, parent)
            c["sites"] = [site_labels[i] for i in sites_idx]
            composites.append(c)
        composites.sort(key=lambda d: d["symdiff"])
        out["composite"] = {
            "n_rules": len(composites), "best": composites[:4],
            "any_exact": any(c["exact"] for c in composites),
            "exact_rules": [c["sites"] for c in composites if c["exact"]],
        }

    # ---- is the attrition depth-correlated in a way that CAN bias anything? ----
    # Only recoverable rows (missing yet measured) can bias a covariate. Where there are
    # none, a depth-correlated membership difference is real but vacuous, because every
    # row it concerns was never eligible under any convention.
    out["depth_correlation_matters"] = bool(
        out.get("depth_correlated") and out["recoverable_rows"] > 0)

    # ---- censoring floor, per landscape (the dataset-inventory §3 correction) ----
    def floor_of(t: Table):
        vals = [v for v in t.entries.values() if v is not None]
        if not vals:
            return None
        lo = min(vals)
        return {"min": lo, "mass_at_min": sum(1 for v in vals if v == lo),
                "n": len(vals)}
    out["clean_floor"] = floor_of(clean)
    out["raw_floor"] = floor_of(raw)
    cf, rf = out["clean_floor"], out["raw_floor"]
    out["floor_identical_raw_vs_clean"] = bool(
        cf and rf and abs(cf["min"] - rf["min"]) < 1e-12)

    # ---- fidelity: do the shared rows carry the same value? ----
    shared = [k for k in (raw.keys & clean.keys)
              if raw.entries[k] is not None and clean.entries[k] is not None]
    worst = 0.0
    worst_key = None
    for k in shared:
        d = abs(raw.entries[k] - clean.entries[k])
        if d > worst:
            worst, worst_key = d, k
    out["shared_measured_rows"] = len(shared)
    out["max_abs_value_deviation"] = worst
    out["max_deviation_key"] = worst_key

    # ---- does the covariate move? ----
    def covariate(t: Table) -> Optional[Dict[str, Any]]:
        items = [(k, v) for k, v in t.entries.items() if v is not None]
        if not items:
            return None
        ids = [k for k, _v in items]
        vals = [v for _k, v in items]
        depths = [depth_of(k, parent) for k in ids]
        return shell_stats(vals, depths, sites, ids=ids, cut=cut)

    cstat, rstat = covariate(clean), covariate(raw)
    out["clean_covariate"] = cstat
    out["raw_covariate"] = rstat
    if cstat and rstat and cstat.get("enrichment") is not None:
        out["enrichment_clean"] = cstat["enrichment"]
        out["enrichment_raw_universe"] = rstat["enrichment"]
        if rstat["enrichment"] is not None:
            out["enrichment_delta"] = rstat["enrichment"] - cstat["enrichment"]
            out["enrichment_rel_delta"] = (
                (rstat["enrichment"] - cstat["enrichment"]) / cstat["enrichment"]
                if cstat["enrichment"] else None)
        out["reweighting_bound"] = reweighting_bound(
            cstat["per_depth"], expected, cut, cstat["base_rate"])
    return out


# ==================================================================================
# Reporting
# ==================================================================================
def fmt_pct(x: Optional[float]) -> str:
    return "-" if x is None else f"{x:.2%}"


def report(results: Dict[str, Dict[str, Any]]) -> None:
    for key, r in results.items():
        print()
        print("=" * 100)
        print(f"{key}   ({r['level']}, {r['sites']} sites, primary cut d<={r['depth_cut']})")
        print("=" * 100)
        a = r["alphabet_sizes"]
        kind = "binary" if r["binary_library"] else f"mixed alphabets {a}"
        print(f"  design space {r['design_space']:,}  ({kind})")
        print(f"  RAW   rows in sheet {r['raw_rows_in_sheet']:,}   keys {r['raw_keys']:,}"
              f"   with a measurement {r['raw_measured']:,}"
              f"   saturation of design space {fmt_pct(r['raw_saturation'])}")
        print(f"  CLEAN rows in sheet {r['clean_rows_in_sheet']:,}"
              f"   blank-key rows {r['clean_blank_key_rows']:,}"
              f"   keys {r['clean_keys']:,}   with a measurement {r['clean_measured']:,}")
        print(f"  SET DIFFERENCE   raw - clean = {r['in_raw_not_clean']:,}"
              f"   clean - raw = {r['in_clean_not_raw']:,}")
        if r["in_raw_not_clean"]:
            print(f"    of the missing: {r['missing_with_no_measurement']:,} have NO raw "
                  f"measurement; {r['missing_WITH_a_measurement']:,} DO carry one")
            print(f"    RECOVERABLE rows (missing yet measured -- the only kind that can "
                  f"move a covariate): {r['recoverable_rows']:,}")
            m = r["measurement_presence_rule"]
            verdict = "EXACT -- 0 FP, 0 FN" if m["exact"] else (
                f"NOT exact: {m['false_positives_kept_but_unmeasured']} FP, "
                f"{m['false_negatives_dropped_but_measured']} FN")
            print(f"    rule 'kept iff the raw measurement is present': {verdict}")
        print(f"  value fidelity on {r['shared_measured_rows']:,} shared rows: "
              f"max |raw - clean| = {r['max_abs_value_deviation']:.3g}")
        cf, rf = r["clean_floor"], r["raw_floor"]
        if cf and rf:
            print(f"  FLOOR   raw min {rf['min']:.6g} ({rf['mass_at_min']:,} rows)"
                  f"   clean min {cf['min']:.6g} ({cf['mass_at_min']:,} rows)"
                  f"   identical: {r['floor_identical_raw_vs_clean']}")
        print(f"  {'depth':<7}{'expected':>12}{'raw':>10}{'clean':>10}{'clean/raw':>12}"
              f"{'raw/expected':>14}")
        for d, e, rr, cc, ret in r["retention"]:
            if e == 0 and rr == 0 and cc == 0:
                continue
            print(f"  {d:<7}{e:>12,}{rr:>10,}{cc:>10,}"
                  f"{(f'{ret:.1%}' if ret is not None else '-'):>12}"
                  f"{(f'{rr/e:.1%}' if e else '-'):>14}")
        rng = r["clean_vs_raw_retention_range"]
        if rng:
            print(f"  clean/raw retention over depths with >=10 raw rows: "
                  f"{rng[0]:.1%} to {rng[1]:.1%}  (spread {r['clean_vs_raw_spread_pp']:.1f} pp)")
            print(f"  => {'DEPTH-CORRELATED' if r['depth_correlated'] else 'depth-blind'}")
        if r.get("raw_vs_space_retention_range"):
            lo, hi = r["raw_vs_space_retention_range"]
            print(f"  raw/design-space retention (library completeness, NOT a cleaning "
                  f"filter): {lo:.1%} to {hi:.1%}  (spread {r['raw_vs_space_spread_pp']:.1f} pp)")
        if r["sweep"]:
            print(f"  SWEEP over {r['sweep']['n_rules']} single-site and site-pair "
                  f"exclusion rules -- best four:")
            for s in r["sweep"]["best"]:
                print(f"    {s['rule']:<52} FP {s['fp']:>6}  FN {s['fn']:>6}  "
                      f"symdiff {s['symdiff']:>6}")
            print(f"    any structural rule exact on its own (0 FP, 0 FN)? "
                  f"{'YES' if r['sweep']['any_exact'] else 'NO'}")
        if r.get("composite"):
            c = r["composite"]
            print(f"  COMPOSITE sweep over {c['n_rules']} rules of the CR9114 form "
                  f"`NOT(mutated at sites ...) AND the measurement is present`:")
            for s in c["best"]:
                print(f"    sites {str(s['sites']):<28} FP {s['fp']:>6}  FN {s['fn']:>6}  "
                      f"symdiff {s['symdiff']:>6}")
            if c["any_exact"]:
                print(f"    => EXACT: {c['exact_rules']}  -- D-8-CLASS STRUCTURAL "
                      f"EXCLUSION, NAMED TO ZERO ERROR")
            else:
                print("    => no composite rule reaches 0 FP / 0 FN: NO structural "
                      "exclusion of the D-8 class")
        if r["in_raw_not_clean"] and r["recoverable_rows"] == 0:
            print("  NOTE: 0 recoverable rows, so the covariate CANNOT move whatever the "
                  "depth profile of the membership difference looks like.")
        if r.get("enrichment_clean") is not None:
            cs, rs = r["clean_covariate"], r["raw_covariate"]
            print(f"  COVARIATE   pre-registered {r['prereg_enrichment']}")
            print(f"    from CLEAN  universe {cs['universe']:,}  tier {cs['tier_size']:,}"
                  f"  shell {cs['shallow_hits']}/{cs['shallow_n']:,}"
                  f" = {cs['shallow_rate']:.4%}  enrichment {cs['enrichment']:.4f}x")
            print(f"    from RAW    universe {rs['universe']:,}  tier {rs['tier_size']:,}"
                  f"  shell {rs['shallow_hits']}/{rs['shallow_n']:,}"
                  f" = {rs['shallow_rate']:.4%}  enrichment {rs['enrichment']:.4f}x")
            d = r.get("enrichment_delta")
            if d is not None:
                rel = r.get("enrichment_rel_delta")
                rels = f"{rel:+.2%} relative" if rel is not None else "relative n/a (clean is 0)"
                print(f"    MOVES BY {d:+.4f} ({rels})"
                      f"   -> {'MOVES' if abs(d) > 5e-4 else 'DOES NOT MOVE'}")
            b = r.get("reweighting_bound")
            if b is not None:
                print(f"    depth-blind reweighting bound (the fallback when no raw table "
                      f"exists): {b:.4f}x")


def summary(results: Dict[str, Dict[str, Any]]) -> None:
    print()
    print("=" * 112)
    print("SUMMARY")
    print("=" * 112)
    head = (f"{'landscape':18}{'level':12}{'raw':>9}{'clean':>9}{'diff':>7}"
            f"  {'cause':<24}{'depth-corr':>14}{'covariate':>16}")
    print(head)
    print("-" * len(head))
    for key, r in results.items():
        diff = r["in_raw_not_clean"]
        if diff == 0:
            cause = "none - identical sets"
        elif r["recoverable_rows"] == 0:
            cause = "absent measurement"
        elif r.get("composite") and r["composite"]["any_exact"]:
            cause = "D-8 STRUCTURAL EXCL."
        else:
            cause = "unexplained"
        if diff == 0:
            dc = "-"
        elif not r["depth_correlated"]:
            dc = "no"
        elif r["recoverable_rows"] == 0:
            dc = "yes(vacuous)"
        else:
            dc = "YES"
        d = r.get("enrichment_delta")
        if r["prereg_enrichment"] is None:
            cov = "not pre-reg."
        elif d is None:
            cov = "?"
        elif abs(d) > 5e-4:
            cov = f"MOVES {d:+.4f}"
        else:
            cov = "does not move"
        print(f"{key:18}{r['level']:12}{r['raw_keys']:>9,}{r['clean_keys']:>9,}"
              f"{diff:>7,}  {cause:<24}{dc:>14}{cov:>16}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--landscape", action="append", default=None)
    parser.add_argument("--json-out", type=Path, default=None)
    args = parser.parse_args(argv)

    print("=" * 100)
    print("GATE 1: the parser reproduces the pre-registered TrpB and GB1 anchors")
    print("=" * 100)
    if not validate():
        print("\nPARSER GATE FAILED -- refusing to audit anything.", file=sys.stderr)
        return 1
    print("  both anchors reproduce exactly; proceeding.")

    keys = args.landscape or list(LANDSCAPES)
    results: Dict[str, Dict[str, Any]] = {}
    for key in keys:
        if key not in LANDSCAPES:
            print(f"unknown landscape {key!r}", file=sys.stderr)
            return 2
        print(f"... {key}", file=sys.stderr)
        results[key] = audit(key, LANDSCAPES[key])

    report(results)
    summary(results)
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(results, indent=2, default=str) + "\n",
                                 encoding="utf-8")
        print(f"\nwrote {args.json_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
