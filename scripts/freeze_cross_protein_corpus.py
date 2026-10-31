#!/usr/bin/env python3
"""
Freeze, verify and characterise the external corpus for the cross-protein evaluation.

This script exists so that `docs/preregistration-cross-protein.md` can be *checked*
rather than believed. It does four things and nothing else:

    fetch   Download the pinned CombinGym files into a local corpus directory.
    record  Write a lock file holding, for every pinned file, its git blob SHA-1
            (the authority: it comes from the pinned commit's tree) and the SHA-256
            of the bytes on disk.
    verify  Recompute both hashes and fail loudly on any mismatch, missing file or
            unexpected size. Exit status 0 means the corpus on disk is bit-identical
            to the pinned commit.
    stats   Compute the label-only, model-free landscape statistics the protocol
            requires to be published BEFORE any arm is run: eligible universe size,
            censoring mass, the achievable top-tier and diversity-adjusted ceilings,
            the depth-stratified hit-rate profile and the shallow-shell enrichment.

Why two hashes. The git blob SHA-1 is reproducible from the commit tree without
trusting this script's author: GitHub's tree API for the pinned commit lists exactly
these values, and `git hash-object` recomputes them. The SHA-256 is the ordinary
content digest, independent of git's framing. A file that matches both is the file
the protocol was frozen against.

Why no third-party dependencies. A freeze that needs pandas and openpyxl to be
verifiable is a freeze that silently changes when those packages change. Everything
here is the standard library, including a deliberately minimal .xlsx reader.

Usage
-----
    python scripts/freeze_cross_protein_corpus.py fetch
    python scripts/freeze_cross_protein_corpus.py record
    python scripts/freeze_cross_protein_corpus.py verify
    python scripts/freeze_cross_protein_corpus.py stats

Default corpus directory is `downloads/combingym/` (gitignored). Pass
--corpus-dir to point elsewhere.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import urllib.request
import zipfile
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple
from xml.etree import ElementTree

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CORPUS_DIR = REPO_ROOT / "downloads" / "combingym"
DEFAULT_LOCK = REPO_ROOT / "scripts" / "data" / "cross_protein_corpus.lock.json"
DEFAULT_STATS = REPO_ROOT / "scripts" / "data" / "cross_protein_corpus_stats.json"

# ==================================================================================
# The pin
# ==================================================================================
# CombinGym has no tags and no releases, so the pin is a commit SHA. It was the tip
# of `main` when this protocol was frozen. The tree SHA is recorded as well so that
# the file list below can be re-derived without trusting the list itself:
#   https://api.github.com/repos/sitonglab/CombinGym/git/trees/<tree_sha>?recursive=1
CORPUS = {
    "name": "CombinGym",
    "repo": "https://github.com/sitonglab/CombinGym",
    "commit": "c5dccd514a9acfe1089a5f1555faea98afd78019",
    "commit_date_utc": "2026-09-16T07:33:53Z",
    "tree_sha": "911f19c7ba82e112056f154d3e8b2073f70dc7d1",
    "raw_base": (
        "https://raw.githubusercontent.com/sitonglab/CombinGym/"
        "c5dccd514a9acfe1089a5f1555faea98afd78019/"
    ),
    "collection_size": "14 landscapes over 9 proteins",
    "frozen_on": "2026-10-31",
}

# (repo path, git blob SHA-1 from the pinned commit tree, size in bytes).
# Only the files the protocol actually reads are pinned: the curated DMS tables, the
# wild-type sequences they are scored against, and the dataset summary. MSAs,
# structures and baseline code are not part of this evaluation and are not pinned.
PINNED_FILES: Tuple[Tuple[str, str, int], ...] = (
    ("Data/Data_summary.xlsx", "45d6067d821f338df36045f56c4959fc3cffc86d", 13055),
    ("Data/DMS/Clean/bnAbs_CR6261_H1_clean.xlsx",
     "90864f0970ab50b6b9ee0623f1acecd8aac81d3b", 164703),
    ("Data/DMS/Clean/bnAbs_CR6261_H9_clean.xlsx",
     "97f2a084427e67b4d017654bb7ad158f55a2c9cc", 149858),
    ("Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx",
     "d9895e1f7e6a25fba4c6e2aced9843f7ad4a7000", 4403699),
    ("Data/DMS/Clean/CreiLOV_clean.xlsx",
     "9ad0d8ff3cf7c742280c84befcd33780a575b31d", 26263997),
    ("Data/DMS/Clean/eqFP611_clean.xlsx",
     "0d7423fe554453675fc2af10ed0320ffdfc21904", 616607),
    ("Data/DMS/Clean/GB1_clean.xlsx",
     "c4b9a7dd4eedceb10c952ebc897c69362f215a4d", 4350703),
    ("Data/DMS/Clean/HIV-1_clean.xlsx",
     "ad96e7c6b000bb9f1033d91fe297fb41ced62c5c", 47401),
    ("Data/DMS/Clean/RhlA_clean.xlsx",
     "0cb393531ab6fff5fa8f037221d50cf543429dcd", 287350),
    ("Data/DMS/Clean/SaCas9_clean.xlsx",
     "76b854e6ec16c09c99efd2e153f6f2a10154f280", 126881),
    ("Data/DMS/Clean/SpCas9_clean.xlsx",
     "00ae9e52824b835ddbd09c569630d813df7414b0", 113908),
    ("Data/WT Sequence/bnAbs_CR6261.fasta",
     "dbd6181366b401d8f9a0dff7218925b23a24d7dc", 141),
    ("Data/WT Sequence/bnAbs_CR9114.fasta",
     "3da4c71ad344f0025bf5a97df49a7fdbcb813bbe", 141),
    ("Data/WT Sequence/CreiLOV.fasta",
     "28203777abb3dd6949c81b07cc5ad9ef187c63f6", 134),
    ("Data/WT Sequence/eqFP611(mTagBFP2).fasta",
     "b73e838216a07d21aa5572a7090a3e9d8826f101", 249),
    ("Data/WT Sequence/GB1.fasta",
     "7f856d7aeaf1769580d01548ed3298f8dd129750", 66),
    ("Data/WT Sequence/HIV-1.fasta",
     "b7724dc1b1db587d26ccca9a7cd1f08b840cfbaa", 120),
    ("Data/WT Sequence/RhlA.fasta",
     "50e0093113f21fbbf192869fe18506d3a68635ba", 307),
    ("Data/WT Sequence/SaCas9(SaCas9-KKH).fasta",
     "8fc81dbe4cff2d4d048d25c120ae61945b2132ff", 1072),
    ("Data/WT Sequence/SpCas9.fasta",
     "90f650f6286b6697f353629fa3f048354a3d85f6", 1383),
)

# ==================================================================================
# Landscape registry
# ==================================================================================
# One entry per (protein, phenotype) landscape, which is the unit CombinGym's own
# benchmark tables use. `sites` is the number of library positions, taken from
# Data_summary.xlsx; `level` is this protocol's two-level validation assignment.
#
# `variant_col` must be a column whose strings are equal-length and differ only at
# library positions, so that plain Hamming distance over the string equals Hamming
# distance over the library. Full genotypes satisfy this (every non-library position
# is identical across rows); the shorter `mutant` combo strings satisfy it directly.
LANDSCAPES: Dict[str, Dict[str, Any]] = {
    "eqfp611_blue": {
        "protein": "eqFP611 / mTagBFP2", "family": "GFP-like fluorescent protein",
        "phenotype": "fluorescence_blue", "sites": 13, "level": "final_test",
        "file": "Data/DMS/Clean/eqFP611_clean.xlsx", "sheet": "genodata",
        "variant_col": "mutant", "depth_col": "n_mut", "value_col": "blue",
    },
    "eqfp611_red": {
        "protein": "eqFP611 / mTagBFP2", "family": "GFP-like fluorescent protein",
        "phenotype": "fluorescence_red", "sites": 13, "level": "final_test",
        "file": "Data/DMS/Clean/eqFP611_clean.xlsx", "sheet": "genodata",
        "variant_col": "mutant", "depth_col": "n_mut", "value_col": "red",
    },
    "eqfp611_combined": {
        "protein": "eqFP611 / mTagBFP2", "family": "GFP-like fluorescent protein",
        "phenotype": "fluorescence_combined", "sites": 13, "level": "final_test",
        "file": "Data/DMS/Clean/eqFP611_clean.xlsx", "sheet": "genodata",
        "variant_col": "mutant", "depth_col": "n_mut", "value_col": "combined",
    },
    "creilov": {
        "protein": "CreiLOV", "family": "LOV / flavin-binding photosensor",
        "phenotype": "fluorescence", "sites": 15, "level": "final_test",
        "file": "Data/DMS/Clean/CreiLOV_clean.xlsx", "sheet": "Sheet1",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "mean",
    },
    "cr9114_h1": {
        "protein": "CR9114", "family": "immunoglobulin heavy-chain variable domain",
        "phenotype": "binding_H1", "sites": 16, "level": "final_test",
        "file": "Data/DMS/Clean/bnAbs_CR9114_H1_clean.xlsx", "sheet": "CR91114-h1",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "h1_mean",
    },
    "cr6261_h1": {
        "protein": "CR6261", "family": "immunoglobulin heavy-chain variable domain",
        "phenotype": "binding_H1", "sites": 11, "level": "final_test",
        "file": "Data/DMS/Clean/bnAbs_CR6261_H1_clean.xlsx", "sheet": "CR6261-h1",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "h1_mean",
    },
    "cr6261_h9": {
        "protein": "CR6261", "family": "immunoglobulin heavy-chain variable domain",
        "phenotype": "binding_H9", "sites": 11, "level": "final_test",
        "file": "Data/DMS/Clean/bnAbs_CR6261_H9_clean.xlsx", "sheet": "CR6261",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "h9_mean",
    },
    "hiv1_protease": {
        "protein": "HIV-1 protease", "family": "retroviral aspartyl protease",
        "phenotype": "activity", "sites": 11, "level": "reserve",
        "file": "Data/DMS/Clean/HIV-1_clean.xlsx", "sheet": "HIV-1",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "fitness",
    },
    "rhla_rls": {
        "protein": "RhlA", "family": "rhamnosyltransferase",
        "phenotype": "activity_RLs", "sites": 11, "level": "reserve",
        "file": "Data/DMS/Clean/RhlA_clean.xlsx", "sheet": "Sheet1",
        "variant_col": "genotype", "depth_col": "num_mut", "value_col": "RLs",
    },
    "rhla_pct_rha": {
        "protein": "RhlA", "family": "rhamnosyltransferase",
        "phenotype": "activity_%Rha-(C8-C10)", "sites": 11, "level": "reserve",
        "file": "Data/DMS/Clean/RhlA_clean.xlsx", "sheet": "Sheet1",
        "variant_col": "genotype", "depth_col": "num_mut", "value_col": "%Rha-(C8-C10)",
    },
    "rhla_rha": {
        "protein": "RhlA", "family": "rhamnosyltransferase",
        "phenotype": "activity_Rha-(C8-C10)", "sites": 11, "level": "reserve",
        "file": "Data/DMS/Clean/RhlA_clean.xlsx", "sheet": "Sheet1",
        "variant_col": "genotype", "depth_col": "num_mut", "value_col": "Rha-(C8-C10)",
    },
    "sacas9": {
        "protein": "SaCas9 (KKH)", "family": "Cas9 endonuclease",
        "phenotype": "activity", "sites": 8, "level": "reserve",
        "file": "Data/DMS/Clean/SaCas9_clean.xlsx", "sheet": "Sheet1",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "Mean",
    },
    "spcas9": {
        "protein": "SpCas9", "family": "Cas9 endonuclease",
        "phenotype": "activity", "sites": 8, "level": "reserve",
        "file": "Data/DMS/Clean/SpCas9_clean.xlsx", "sheet": "SpCas9",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "mean",
    },
    "gb1_combingym": {
        "protein": "GB1", "family": "immunoglobulin-binding protein G domain B1",
        "phenotype": "binding", "sites": 4, "level": "development",
        "file": "Data/DMS/Clean/GB1_clean.xlsx", "sheet": "GB1_clean",
        "variant_col": "genotype", "depth_col": "n_mut", "value_col": "fitness",
    },
}

# Top-tier candidates, tried in order. The protocol deepens the tier when the 1%
# tier holds too few variants for the endpoint to resolve.
TIER_LADDER: Tuple[float, ...] = (0.01, 0.02, 0.05, 0.10)
MIN_TIER_MEMBERS = 30
# The absolute depth<=2 shell is reported only when it holds at least this many
# variants; below it a 1%-order hit rate cannot be estimated usefully.
MIN_SHELL_FOR_ABSOLUTE_CUT = 500
BUDGETS: Tuple[int, ...] = (25, 50, 75, 100, 125)
DIVERSITY_DISTANCES: Tuple[int, ...] = (1, 2, 3)


# ==================================================================================
# Hashing
# ==================================================================================
def sha256_file(path: Path) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def git_blob_sha1(path: Path) -> str:
    """Reproduces `git hash-object <path>`: SHA-1 over b'blob <len>\\0' + content."""
    size = path.stat().st_size
    digest = hashlib.sha1()
    digest.update(b"blob %d\0" % size)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ==================================================================================
# Minimal .xlsx reader (stdlib only)
# ==================================================================================
_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _col_index(ref: str) -> int:
    letters = re.match(r"([A-Z]+)", ref or "")
    if not letters:
        return 0
    index = 0
    for char in letters.group(1):
        index = index * 26 + (ord(char) - 64)
    return index - 1


def _shared_strings(archive: zipfile.ZipFile) -> List[str]:
    try:
        raw = archive.read("xl/sharedStrings.xml")
    except KeyError:
        return []
    out: List[str] = []
    for si in ElementTree.fromstring(raw).iter(f"{_NS}si"):
        out.append("".join(t.text or "" for t in si.iter(f"{_NS}t")))
    return out


def _sheet_path(archive: zipfile.ZipFile, sheet: Optional[str]) -> Tuple[str, str]:
    workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    rels = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {rel.get("Id"): rel.get("Target") for rel in rels}
    entries: List[Tuple[str, str]] = []
    for node in workbook.iter(f"{_NS}sheet"):
        target = targets.get(node.get(f"{_REL_NS}id"), "")
        target = target[1:] if target.startswith("/") else "xl/" + target.lstrip("./")
        entries.append((node.get("name") or "", target))
    if not entries:
        raise ValueError("workbook declares no sheets")
    if sheet is None:
        return entries[0]
    for name, target in entries:
        if name == sheet:
            return name, target
    raise ValueError(f"sheet {sheet!r} not found; available: {[n for n, _ in entries]}")


def read_xlsx(path: Path, sheet: Optional[str] = None
              ) -> Tuple[str, List[str], List[List[Optional[str]]]]:
    """Returns (sheet name, header row, data rows). Cells come back as strings or None.

    Only what this script needs: shared and inline strings, numbers, booleans. Dates
    are returned as their raw serial number, which is fine because no landscape
    column holds a date.
    """
    with zipfile.ZipFile(path) as archive:
        strings = _shared_strings(archive)
        name, target = _sheet_path(archive, sheet)
        rows: List[List[Optional[str]]] = []
        with archive.open(target) as handle:
            for _event, element in ElementTree.iterparse(handle, events=("end",)):
                if element.tag != f"{_NS}row":
                    continue
                cells: List[Optional[str]] = []
                for cell in element.iter(f"{_NS}c"):
                    index = _col_index(cell.get("r", ""))
                    while len(cells) < index:
                        cells.append(None)
                    kind = cell.get("t")
                    if kind == "inlineStr":
                        node = cell.find(f"{_NS}is")
                        value = "".join(t.text or "" for t in node.iter(f"{_NS}t")) if node is not None else None
                    else:
                        node = cell.find(f"{_NS}v")
                        value = node.text if node is not None else None
                        if kind == "s" and value is not None:
                            position = int(value)
                            value = strings[position] if position < len(strings) else None
                    cells.append(value)
                rows.append(cells)
                element.clear()
    if not rows:
        return name, [], []
    header = [(c or "").strip() for c in rows[0]]
    return name, header, rows[1:]


# ==================================================================================
# Statistics
# ==================================================================================
def wilson_interval(successes: int, trials: int, z: float = 1.959964) -> Tuple[float, float]:
    """Wilson score interval. Used because hit counts here are often single digits,
    where the normal approximation is not usable."""
    if trials == 0:
        return (float("nan"), float("nan"))
    phat = successes / trials
    denominator = 1.0 + z * z / trials
    centre = (phat + z * z / (2 * trials)) / denominator
    half = (z / denominator) * math.sqrt(phat * (1 - phat) / trials + z * z / (4 * trials * trials))
    return (max(0.0, centre - half), min(1.0, centre + half))


def hamming(left: str, right: str) -> int:
    return sum(1 for a, b in zip(left, right) if a != b)


def greedy_distinct(ordered: Sequence[str], min_distance: int) -> int:
    """Diversity-adjusted count: walk hits in fitness-descending order and keep one
    only when it is at least `min_distance` from everything already kept. This is
    the definition the protocol uses, applied to the full hit set to give a ceiling."""
    kept: List[str] = []
    for variant in ordered:
        if all(hamming(variant, other) >= min_distance for other in kept):
            kept.append(variant)
    return len(kept)


def landscape_stats(key: str, spec: Dict[str, Any], corpus_dir: Path) -> Dict[str, Any]:
    path = corpus_dir / spec["file"]
    out: Dict[str, Any] = {
        "key": key, "protein": spec["protein"], "family": spec["family"],
        "phenotype": spec["phenotype"], "sites": spec["sites"],
        "level": spec["level"], "file": spec["file"],
    }
    if not path.exists():
        out["status"] = "file_absent"
        return out
    if not spec.get("sheet") and not spec.get("value_col"):
        out["status"] = "registry_incomplete"
        out["note"] = "sheet/column names not yet resolved; see the protocol's open items"
        return out

    sheet_name, header, rows = read_xlsx(path, spec["sheet"])
    out["sheet"] = sheet_name
    try:
        variant_at = header.index(spec["variant_col"])
        value_at = header.index(spec["value_col"])
        depth_at = header.index(spec["depth_col"])
    except ValueError as exc:
        out["status"] = "column_missing"
        out["note"] = f"{exc}; header={header}"
        return out

    out["rows_in_sheet"] = len(rows)
    variants: List[str] = []
    values: List[float] = []
    depths: List[int] = []
    dropped_missing = 0
    for row in rows:
        def cell(i: int) -> Optional[str]:
            return row[i] if i < len(row) else None
        raw_variant, raw_value, raw_depth = cell(variant_at), cell(value_at), cell(depth_at)
        if raw_variant is None or raw_value is None or raw_depth is None or raw_value == "":
            dropped_missing += 1
            continue
        try:
            value = float(raw_value)
            depth = int(float(raw_depth))
        except ValueError:
            dropped_missing += 1
            continue
        if not math.isfinite(value):
            dropped_missing += 1
            continue
        variants.append(raw_variant.strip())
        values.append(value)
        depths.append(depth)

    n = len(values)
    out["dropped_missing_or_unparseable"] = dropped_missing
    out["eligible_universe"] = n
    if n == 0:
        out["status"] = "empty"
        return out

    widths = sorted({len(v) for v in variants})
    out["variant_string_widths"] = widths
    out["variant_width_uniform"] = len(widths) == 1
    out["duplicate_variant_ids"] = n - len(set(variants))

    lo, hi = min(values), max(values)
    at_min = sum(1 for v in values if v == lo)
    at_max = sum(1 for v in values if v == hi)
    out["value_min"] = lo
    out["value_max"] = hi
    out["distinct_values"] = len(set(values))
    out["mass_at_min"] = at_min
    out["frac_at_min"] = at_min / n
    out["mass_at_max"] = at_max
    out["frac_at_max"] = at_max / n
    out["censored_from_below"] = at_min / n >= 0.10
    out["censored_from_above"] = at_max / n >= 0.01

    # Tier selection: deepen until the tier holds enough members to resolve.
    order = sorted(range(n), key=lambda i: (-values[i], variants[i]))
    tier_used = None
    tier_size = 0
    for tier in TIER_LADDER:
        tier_size = math.ceil(tier * n)
        if tier_size >= MIN_TIER_MEMBERS:
            tier_used = tier
            break
    if tier_used is None:
        tier_used = TIER_LADDER[-1]
        tier_size = math.ceil(tier_used * n)
    out["tier_fraction"] = tier_used
    out["tier_deepened_from_1pct"] = tier_used != 0.01
    out["tier_size_nominal"] = tier_size

    boundary_value = values[order[tier_size - 1]]
    hits_exclusive = [i for i in order if values[i] > boundary_value]
    hits_inclusive = [i for i in order if values[i] >= boundary_value]
    out["tier_boundary_value"] = boundary_value
    out["hits_strictly_above_boundary"] = len(hits_exclusive)
    out["hits_including_boundary_ties"] = len(hits_inclusive)
    out["boundary_tie_mass"] = len(hits_inclusive) - len(hits_exclusive)
    # The nominal tier is used as the hit set; ties at the boundary are reported so a
    # reader can see when the cut is arbitrary.
    hits = order[:tier_size]
    hit_set = set(hits)
    out["hit_count_ceiling"] = len(hits)

    ordered_hit_ids = [variants[i] for i in hits]
    if out["variant_width_uniform"]:
        for d in DIVERSITY_DISTANCES:
            out[f"diversity_adjusted_ceiling_d{d}"] = greedy_distinct(ordered_hit_ids, d)
    else:
        out["diversity_adjusted_ceiling_note"] = (
            "variant strings are not of uniform width; Hamming distance is undefined "
            "as computed here and the diversity ceiling was not calculated"
        )

    out["achievable_at_budget"] = {
        str(b): {
            "raw_hits": min(b, out["hit_count_ceiling"]),
            "diversity_adjusted_d2": min(b, out.get("diversity_adjusted_ceiling_d2", 0)),
            "fraction_of_universe": b / n,
        }
        for b in BUDGETS
    }

    hit_depths = [depths[i] for i in hits]
    out["mean_depth_of_hit"] = sum(hit_depths) / len(hit_depths)
    out["max_depth_observed"] = max(depths)

    profile: Dict[str, Any] = {}
    for depth in range(0, max(depths) + 1):
        members = [i for i, d in enumerate(depths) if d == depth]
        if not members:
            continue
        hits_here = sum(1 for i in members if i in hit_set)
        low, high = wilson_interval(hits_here, len(members))
        profile[str(depth)] = {
            "variants": len(members), "hits": hits_here,
            "hit_rate": hits_here / len(members),
            "wilson95": [low, high],
        }
    out["depth_profile"] = profile

    base_rate = out["hit_count_ceiling"] / n

    def shell(cut: int) -> Dict[str, Any]:
        members = [i for i, d in enumerate(depths) if d <= cut]
        hits_here = sum(1 for i in members if i in hit_set)
        low, high = wilson_interval(hits_here, len(members)) if members else (float("nan"), float("nan"))
        return {
            "depth_cut": cut,
            "shell_variants": len(members),
            "shell_hits": hits_here,
            "shell_hit_rate": (hits_here / len(members)) if members else None,
            "enrichment_vs_base_rate": ((hits_here / len(members)) / base_rate)
                                       if members and base_rate > 0 else None,
            "wilson95_hit_rate": [low, high],
            "enrichment_wilson95": [low / base_rate, high / base_rate] if base_rate > 0 else None,
            "resolvable": len(members) >= MIN_SHELL_FOR_ABSOLUTE_CUT,
        }

    primary_cut = spec["sites"] // 2
    out["shallow_shell_primary"] = shell(primary_cut)
    out["shallow_shell_primary"]["definition"] = "depth <= floor(sites/2)"
    out["shallow_shell_absolute_d2"] = shell(2)
    out["shallow_shell_absolute_d2"]["definition"] = "depth <= 2 (absolute)"
    out["shallow_shell_quarter"] = shell(max(1, round(spec["sites"] / 4)))
    out["shallow_shell_quarter"]["definition"] = "depth <= round(sites/4)"
    out["base_rate"] = base_rate
    out["status"] = "ok"
    return out


# ==================================================================================
# Commands
# ==================================================================================
def cmd_fetch(args: argparse.Namespace) -> int:
    corpus_dir: Path = args.corpus_dir
    failures = 0
    for rel, blob, size in PINNED_FILES:
        target = corpus_dir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and git_blob_sha1(target) == blob:
            print(f"have  {rel}")
            continue
        url = CORPUS["raw_base"] + urllib.request.quote(rel)
        print(f"get   {rel}  ({size} bytes)")
        try:
            with urllib.request.urlopen(url, timeout=args.timeout) as response:
                payload = response.read()
            target.write_bytes(payload)
        except Exception as exc:  # noqa: BLE001 - any failure is reported, not raised
            print(f"FAIL  {rel}: {type(exc).__name__}: {exc}")
            failures += 1
            continue
        actual = git_blob_sha1(target)
        if actual != blob:
            print(f"FAIL  {rel}: blob sha1 {actual} != pinned {blob}")
            failures += 1
    print(f"\n{len(PINNED_FILES) - failures}/{len(PINNED_FILES)} files present and matching the pin.")
    return 1 if failures else 0


def cmd_record(args: argparse.Namespace) -> int:
    corpus_dir: Path = args.corpus_dir
    lock_path: Path = args.lock
    if lock_path.exists() and not args.force:
        print(f"{lock_path} already exists. A freeze is meant to be written once; "
              f"pass --force only if you intend to change the pin.", file=sys.stderr)
        return 2

    entries: List[Dict[str, Any]] = []
    missing: List[str] = []
    mismatched: List[str] = []
    for rel, blob, size in PINNED_FILES:
        path = corpus_dir / rel
        if not path.exists():
            missing.append(rel)
            entries.append({"path": rel, "git_blob_sha1": blob, "size_bytes": size,
                            "sha256": None, "sha256_status": "not_recorded_file_absent"})
            continue
        actual_blob = git_blob_sha1(path)
        digest, actual_size = sha256_file(path)
        if actual_blob != blob or actual_size != size:
            mismatched.append(rel)
            entries.append({"path": rel, "git_blob_sha1": blob, "size_bytes": size,
                            "sha256": None,
                            "sha256_status": "refused_local_copy_disagrees_with_pin",
                            "observed_git_blob_sha1": actual_blob,
                            "observed_size_bytes": actual_size})
            continue
        entries.append({"path": rel, "git_blob_sha1": blob, "size_bytes": size,
                        "sha256": digest, "sha256_status": "recorded"})

    document = {
        "schema": "veritasbio/cross-protein-corpus-lock/1",
        "corpus": CORPUS,
        "notes": [
            "git_blob_sha1 is authoritative: it is the value GitHub's tree API reports "
            "for the pinned commit, and `git hash-object` reproduces it.",
            "sha256 is recorded only for files whose local copy already matched the pin.",
            "Files listed with a non-'recorded' sha256_status were not verifiable at "
            "record time; `verify` still checks their blob SHA-1 against the pin.",
        ],
        "files": entries,
    }
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {lock_path}")
    print(f"  recorded sha256 for {sum(1 for e in entries if e['sha256'])}/{len(entries)} files")
    if missing:
        print(f"  absent locally ({len(missing)}): {', '.join(missing)}")
    if mismatched:
        print(f"  DISAGREE with the pin ({len(mismatched)}): {', '.join(mismatched)}")
    return 1 if mismatched else 0


def cmd_verify(args: argparse.Namespace) -> int:
    corpus_dir: Path = args.corpus_dir
    lock_path: Path = args.lock
    if not lock_path.exists():
        print(f"no lock file at {lock_path}; run `record` first", file=sys.stderr)
        return 2
    document = json.loads(lock_path.read_text(encoding="utf-8"))
    if document.get("corpus", {}).get("commit") != CORPUS["commit"]:
        print(f"PIN MISMATCH: lock pins {document.get('corpus', {}).get('commit')!r}, "
              f"this script pins {CORPUS['commit']!r}", file=sys.stderr)
        return 2

    failures: List[str] = []
    skipped: List[str] = []
    checked = 0
    for entry in document["files"]:
        rel = entry["path"]
        path = corpus_dir / rel
        if not path.exists():
            failures.append(f"{rel}: absent")
            continue
        actual_blob = git_blob_sha1(path)
        if actual_blob != entry["git_blob_sha1"]:
            failures.append(f"{rel}: blob sha1 {actual_blob} != {entry['git_blob_sha1']}")
            continue
        digest, size = sha256_file(path)
        if size != entry["size_bytes"]:
            failures.append(f"{rel}: size {size} != {entry['size_bytes']}")
            continue
        if entry.get("sha256"):
            if digest != entry["sha256"]:
                failures.append(f"{rel}: sha256 {digest} != {entry['sha256']}")
                continue
            checked += 1
        else:
            skipped.append(f"{rel} (blob sha1 ok; no sha256 in lock)")
            checked += 1
        print(f"ok    {rel}")
    print()
    if skipped:
        print(f"{len(skipped)} file(s) verified by blob SHA-1 only:")
        for item in skipped:
            print(f"  - {item}")
    if failures:
        print(f"FAILED {len(failures)} of {len(document['files'])} file(s):", file=sys.stderr)
        for item in failures:
            print(f"  - {item}", file=sys.stderr)
        return 1
    print(f"CORPUS VERIFIED: {checked}/{len(document['files'])} files match "
          f"{CORPUS['name']} @ {CORPUS['commit'][:12]}")
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    corpus_dir: Path = args.corpus_dir
    keys: Iterable[str] = args.landscape or LANDSCAPES.keys()
    results: Dict[str, Any] = {}
    for key in keys:
        if key not in LANDSCAPES:
            print(f"unknown landscape {key!r}", file=sys.stderr)
            return 2
        print(f"... {key}", file=sys.stderr)
        results[key] = landscape_stats(key, LANDSCAPES[key], corpus_dir)

    document = {
        "schema": "veritasbio/cross-protein-landscape-stats/1",
        "corpus": {k: CORPUS[k] for k in ("name", "repo", "commit", "commit_date_utc")},
        "conventions": {
            "direction": "higher is better for every landscape listed",
            "tier_ladder": list(TIER_LADDER),
            "min_tier_members": MIN_TIER_MEMBERS,
            "budgets": list(BUDGETS),
            "diversity_distances": list(DIVERSITY_DISTANCES),
            "shallow_shell_primary": "depth <= floor(sites/2)",
            "min_shell_for_absolute_cut": MIN_SHELL_FOR_ABSOLUTE_CUT,
        },
        "landscapes": results,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {args.out}")

    header = f"{'landscape':22} {'N':>8} {'tier':>6} {'hits':>6} {'d>=2':>6} {'shallow':>9} {'enrich':>7} {'floor%':>7}"
    print(header)
    print("-" * len(header))
    for key, row in results.items():
        if row.get("status") != "ok":
            print(f"{key:22} {row.get('status', '?')}")
            continue
        shell = row["shallow_shell_primary"]
        print(f"{key:22} {row['eligible_universe']:8d} {row['tier_fraction']:6.2%} "
              f"{row['hit_count_ceiling']:6d} {row.get('diversity_adjusted_ceiling_d2', -1):6d} "
              f"{shell['shell_hit_rate']:9.4%} {shell['enrichment_vs_base_rate']:7.2f} "
              f"{row['frac_at_min']:7.2%}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--corpus-dir", type=Path, default=DEFAULT_CORPUS_DIR,
                        help="Local mirror of the pinned repo paths (default: %(default)s).")
    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch = subparsers.add_parser("fetch", help="Download the pinned files.")
    fetch.add_argument("--timeout", type=int, default=600)
    fetch.set_defaults(func=cmd_fetch)

    record = subparsers.add_parser("record", help="Write the hash lock file.")
    record.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    record.add_argument("--force", action="store_true",
                        help="Overwrite an existing lock file.")
    record.set_defaults(func=cmd_record)

    verify = subparsers.add_parser("verify", help="Check the corpus against the lock file.")
    verify.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    verify.set_defaults(func=cmd_verify)

    stats = subparsers.add_parser(
        "stats", help="Compute the pre-run landscape statistics (ceilings, depth profile, "
                      "shallow-shell enrichment).")
    stats.add_argument("--out", type=Path, default=DEFAULT_STATS)
    stats.add_argument("--landscape", action="append", default=None,
                       help="Restrict to one landscape key; repeatable.")
    stats.set_defaults(func=cmd_stats)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
