#!/usr/bin/env python3
"""Canonicalise combinatorial-library genotypes to full library length, and validate it.

THE RULE THIS MODULE EXISTS TO ENFORCE. ProteinGym-style mutation strings OMIT wild-type
positions, so a genotype over a 16-site library may present as three tokens. Zipping a
token list against a position list, or truncating either to the other's length, silently
re-assigns mutations to the wrong sites. That bug inverted a conclusion earlier in this
project. Here every genotype is expanded to a string of exactly `len(positions)`
characters with the wild-type residue filled in at every omitted position, and equal
length is ASSERTED, never assumed.

VALIDATION GATE. `python scripts/analysis/genotype_canon.py` reproduces the two
pre-registered reference figures from the committed TrpB and GB1 files:

    TrpB shallow shell  168/2,242 = 7.4933 %
    GB1  shallow shell   72/2,167 = 3.3226 %

under the protocol's depth-0 convention: the wild type is RETAINED in the eligible
universe and the tier, and EXCLUDED from both the shallow and the deep shell. If either
figure does not reproduce, the parser is wrong and nothing computed with it is
trustworthy; the script exits non-zero.
"""
from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = REPO_ROOT / "scripts" / "data"

TOKEN_RE = re.compile(r"^([A-Za-z\*])(\d+)([A-Za-z\*])$")


class ParseError(ValueError):
    """Raised when a mutation string cannot be canonicalised without guessing."""


def split_tokens(mutant: str) -> List[str]:
    """Split a mutation string on ':' (ProteinGym) or '_' (CombinGym raw HIV-1).

    'WT' / 'wt' / '' denote the parent and yield no tokens.
    """
    s = (mutant or "").strip()
    if not s or s.upper() == "WT":
        return []
    sep = ":" if ":" in s else ("_" if "_" in s else None)
    parts = [s] if sep is None else s.split(sep)
    return [p.strip() for p in parts if p.strip()]


def parse_tokens(mutant: str) -> List[Tuple[str, int, str]]:
    """Parse into (wt_aa, 1-based position, mut_aa). Raises on anything unparseable."""
    out: List[Tuple[str, int, str]] = []
    for tok in split_tokens(mutant):
        m = TOKEN_RE.match(tok)
        if not m:
            raise ParseError(f"unparseable mutation token {tok!r} in {mutant!r}")
        out.append((m.group(1).upper(), int(m.group(2)), m.group(3).upper()))
    return out


def canonicalise(mutant: str, positions: Sequence[int], wt_residues: str) -> str:
    """Expand a mutation string to one character per library position, WT at omissions.

    `positions` are the library positions in a FIXED order; `wt_residues[i]` is the
    parent residue at `positions[i]`. The returned string has length len(positions) --
    asserted, not assumed -- so Hamming distance over it is Hamming distance over the
    library.
    """
    if len(wt_residues) != len(positions):
        raise ParseError(
            f"library definition is inconsistent: {len(positions)} positions but "
            f"{len(wt_residues)} wild-type residues")
    index = {pos: i for i, pos in enumerate(positions)}
    if len(index) != len(positions):
        raise ParseError(f"duplicate library positions: {list(positions)}")
    chars = list(wt_residues)
    for wt_aa, pos, mut_aa in parse_tokens(mutant):
        slot = index.get(pos)
        if slot is None:
            raise ParseError(
                f"position {pos} in {mutant!r} is not a library position "
                f"{sorted(positions)}")
        if wt_aa != wt_residues[slot]:
            raise ParseError(
                f"{mutant!r} asserts wild type {wt_aa!r} at position {pos}, but the "
                f"library's parent residue there is {wt_residues[slot]!r}")
        chars[slot] = mut_aa
    out = "".join(chars)
    # The assertion the module exists for.
    assert len(out) == len(positions), (
        f"canonical genotype {out!r} has length {len(out)}, expected {len(positions)}")
    return out


def depth_of(canonical: str, wt_residues: str) -> int:
    """Hamming distance to the parent over library positions. Lengths must be equal."""
    if len(canonical) != len(wt_residues):
        raise ParseError(
            f"length mismatch: genotype {len(canonical)} vs wild type "
            f"{len(wt_residues)} -- never zip or truncate these")
    return sum(1 for a, b in zip(canonical, wt_residues) if a != b)


def infer_library(mutants: Sequence[str]) -> Tuple[List[int], str]:
    """Derive (positions, wt_residues) from a column of mutation strings.

    Every token that names a position must agree on the parent residue there; a
    disagreement means the strings are not a single combinatorial library and is fatal
    rather than silently resolved.
    """
    parent: Dict[int, str] = {}
    for mutant in mutants:
        for wt_aa, pos, _mut_aa in parse_tokens(mutant):
            if pos in parent and parent[pos] != wt_aa:
                raise ParseError(
                    f"position {pos} is given parent residue {parent[pos]!r} by one row "
                    f"and {wt_aa!r} by {mutant!r}")
            parent.setdefault(pos, wt_aa)
    positions = sorted(parent)
    return positions, "".join(parent[p] for p in positions)


def alphabet_design_space(canonicals: Sequence[str]) -> Tuple[int, List[int]]:
    """Product of observed per-site alphabet sizes, and those sizes.

    For a binary library of n sites this returns 2**n and the exact depth-k expectation
    is C(n, k); for a mixed-alphabet library the design space is the product and the
    depth-k expectation must be computed by the per-site generating product, which
    `depth_expectation` does.
    """
    if not canonicals:
        return 0, []
    width = len(canonicals[0])
    sizes: List[int] = []
    for i in range(width):
        sizes.append(len({c[i] for c in canonicals}))
    space = 1
    for s in sizes:
        space *= s
    return space, sizes


def depth_expectation(alphabet_sizes: Sequence[int]) -> List[int]:
    """Exact count of genotypes at each depth for a mixed-alphabet product library.

    Each site contributes (1 + (a_i - 1) x) to the generating polynomial in x, where
    a_i is that site's alphabet size: one way to be wild type, a_i - 1 ways to be
    mutated. The coefficient of x^k is the number of depth-k genotypes. For all
    a_i = 2 this reduces to C(n, k), the binary case.
    """
    poly = [1]
    for a in alphabet_sizes:
        nxt = [0] * (len(poly) + 1)
        for k, coeff in enumerate(poly):
            nxt[k] += coeff
            nxt[k + 1] += coeff * (a - 1)
        poly = nxt
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


# ==================================================================================
# Shell statistics, exactly as freeze_cross_protein_corpus.py defines them
# ==================================================================================
TIER_LADDER: Tuple[float, ...] = (0.01, 0.02, 0.05, 0.10)
MIN_TIER_MEMBERS = 30


def shell_stats(values: Sequence[float], depths: Sequence[int], sites: int,
                ids: Optional[Sequence[str]] = None,
                cut: Optional[int] = None) -> Dict[str, object]:
    """Tier, hit set, shallow/deep shells and enrichment.

    Mirrors the freeze script: the tier is the top ceil(frac * N) rows by value
    descending with ties broken by the variant identifier, deepening along TIER_LADDER
    until it holds MIN_TIER_MEMBERS. Depth 0 is retained in the universe and the tier
    and excluded from both shells.
    """
    n = len(values)
    if n == 0:
        return {"universe": 0}
    keys = list(ids) if ids is not None else [str(i) for i in range(n)]
    order = sorted(range(n), key=lambda i: (-values[i], keys[i]))
    tier_used, tier_size = TIER_LADDER[-1], math.ceil(TIER_LADDER[-1] * n)
    for frac in TIER_LADDER:
        size = math.ceil(frac * n)
        if size >= MIN_TIER_MEMBERS:
            tier_used, tier_size = frac, size
            break
    hits = set(order[:tier_size])
    base_rate = tier_size / n
    cut = sites // 2 if cut is None else cut
    shallow = [i for i, d in enumerate(depths) if 1 <= d <= cut]
    deep = [i for i, d in enumerate(depths) if d > cut]
    shallow_hits = sum(1 for i in shallow if i in hits)
    deep_hits = sum(1 for i in deep if i in hits)
    shallow_rate = (shallow_hits / len(shallow)) if shallow else None
    return {
        "universe": n,
        "tier_fraction": tier_used,
        "tier_size": tier_size,
        "tier_boundary_value": values[order[tier_size - 1]],
        "base_rate": base_rate,
        "depth_cut": cut,
        "shallow_n": len(shallow),
        "shallow_hits": shallow_hits,
        "shallow_rate": shallow_rate,
        "enrichment": (shallow_rate / base_rate) if shallow_rate is not None else None,
        "deep_n": len(deep),
        "deep_hits": deep_hits,
        "deep_rate": (deep_hits / len(deep)) if deep else None,
        "depth0_rows": sum(1 for d in depths if d == 0),
        "max_depth": max(depths) if depths else None,
        "per_depth": _per_depth(depths, hits),
    }


def _per_depth(depths: Sequence[int], hits) -> Dict[int, Tuple[int, int]]:
    out: Dict[int, Tuple[int, int]] = {}
    for i, d in enumerate(depths):
        n, h = out.get(d, (0, 0))
        out[d] = (n + 1, h + (1 if i in hits else 0))
    return dict(sorted(out.items()))


def reweighting_bound(per_depth: Dict[int, Tuple[int, int]],
                      expected: Sequence[int], cut: int,
                      base_rate: Optional[float] = None) -> Optional[float]:
    """Depth-blind reweighting bound on shallow enrichment, for when no raw table exists.

    Hold each depth's OBSERVED hit rate fixed and reweight to the exact design-space
    population at that depth. BOTH the shallow shell and the BASE RATE are reweighted:
    the base rate is the hit rate over the whole reweighted universe, because a table
    that under-represents deep shells also understates the overall hit rate whenever the
    hit rate varies with depth. Dividing a reweighted shell rate by the table's own
    unreweighted base rate is the mistake that makes the bound look optimistic.

    On CR9114 this predicts 0.5622x against a true 0.5845x -- right direction and 3.8 %
    relative -- which is the figure docs/cr9114-reprep.md §6 records, and is how this
    implementation is checked. `base_rate` is accepted and ignored for call compatibility.
    """
    del base_rate  # the whole point is that the bound reweights the base rate itself
    shell_num = shell_den = 0.0
    all_num = all_den = 0.0
    for d, (n_d, h_d) in per_depth.items():
        if n_d == 0 or d >= len(expected):
            continue
        rate = h_d / n_d
        all_num += expected[d] * rate
        all_den += expected[d]
        if 1 <= d <= cut:
            shell_num += expected[d] * rate
            shell_den += expected[d]
    if shell_den == 0 or all_den == 0 or all_num == 0:
        return None
    reweighted_base = all_num / all_den
    if reweighted_base <= 0:
        return None
    return (shell_num / shell_den) / reweighted_base


# ==================================================================================
# The validation gate
# ==================================================================================
def _load_prepared(path: Path) -> Tuple[List[str], List[float]]:
    mutants: List[str] = []
    values: List[float] = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            raw = row.get("DMS_score")
            if raw is None or raw == "":
                continue
            mutants.append((row.get("mutant") or "").strip())
            values.append(float(raw))
    return mutants, values


def validate(verbose: bool = True) -> bool:
    """Reproduce the TrpB and GB1 anchors. Returns True only if BOTH match exactly."""
    anchors = [
        ("TrpB", "TRPB_THEMA_Johnston_2024.csv", 4, 159_999, 168, 2_242, 0.074933),
        ("GB1", "GB1_STRSG_Wu_2016.csv", 4, 149_361, 72, 2_167, 0.033226),
    ]
    ok = True
    for name, fname, sites, exp_n, exp_hits, exp_shell, exp_rate in anchors:
        path = DATA_DIR / fname
        if not path.exists():
            print(f"  {name}: FILE ABSENT -- checked {path}")
            ok = False
            continue
        mutants, values = _load_prepared(path)
        positions, wt = infer_library(mutants)
        canon = [canonicalise(m, positions, wt) for m in mutants]
        depths = [depth_of(c, wt) for c in canon]
        stats = shell_stats(values, depths, sites, ids=canon)
        got_n = stats["universe"]
        got_hits = stats["shallow_hits"]
        got_shell = stats["shallow_n"]
        got_rate = stats["shallow_rate"]
        good = (got_n == exp_n and got_hits == exp_hits and got_shell == exp_shell
                and abs(got_rate - exp_rate) < 5e-7)
        ok = ok and good
        if verbose:
            print(f"  {name:5} sites={sites} positions={positions} wt={wt!r}")
            print(f"        universe {got_n:,} (expected {exp_n:,})   "
                  f"tier {stats['tier_size']:,}")
            print(f"        shallow  {got_hits}/{got_shell} = {got_rate:.4%}   "
                  f"(expected {exp_hits}/{exp_shell} = {exp_rate:.4%})   "
                  f"{'PASS' if good else 'FAIL'}")
            print(f"        deep     {stats['deep_hits']}/{stats['deep_n']} = "
                  f"{stats['deep_rate']:.4%}   enrichment {stats['enrichment']:.4f}x")
    return ok


def main(argv: Optional[Sequence[str]] = None) -> int:
    print("=" * 96)
    print("PARSER VALIDATION GATE -- reproduce the pre-registered TrpB and GB1 anchors")
    print("=" * 96)
    print("  convention: depth 0 retained in the universe and the tier, excluded from "
          "both shells")
    if not validate():
        print("\nGATE FAILED. The parser is wrong; nothing computed with it is "
              "trustworthy.", file=sys.stderr)
        return 1
    print("\nGATE PASSED. Both anchors reproduce exactly.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
