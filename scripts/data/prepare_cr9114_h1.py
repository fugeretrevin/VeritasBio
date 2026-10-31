"""Re-prepare CR9114_H1_Phillips_2021.csv from Phillips et al.'s deposited source data.

WHY THIS EXISTS
---------------
The previously committed file held 48,840 of the 65,536 CR9114 genotypes. The attrition
was not random: it was the exact exclusion of one quadrant of the design space (every
genotype carrying BOTH mutation 1, F29S, and mutation 16, Y106S), which is mechanically
depth-correlated because the fraction of depth-k genotypes carrying two specified
mutations is C(14,k-2)/C(16,k) = k(k-1)/240. Every per-depth statistic computed on that
file was therefore shaped by the filter rather than by biology. See `docs/cr9114-reprep.md`.

THE AUTHORITATIVE SOURCE
------------------------
Phillips A.M. et al. (2021) eLife 10:e71393, Figure 1 - source data 1,
"CR9114 library -logK_D to H1, H3, and influenza B":

    https://cdn.elifesciences.org/articles/71393/elife-71393-fig1-data1-v3.csv
    12,258,680 bytes
    sha256 ee31f8fc26fce2730d2ae772e7a095aa334b4c6c5d8e7fbfe4b4a5a1c9a794fe

It holds all 2^16 = 65,536 genotypes, one row each, with `h1_mean` blank on the 442 that
were not quantified. CombinGym's `Data/DMS/Raw/bnAbs_CR9114_raw.csv` at the pinned commit
is BYTE-IDENTICAL to this file, so the corpus pin is preserved: this is not a new upstream,
it is the raw table CombinGym already ships, read instead of its derived clean workbook.

CONVENTIONS, matched to the committed CR6261 files and the corpus lock
---------------------------------------------------------------------
  columns           mutant,DMS_score,DMS_score_bin,stage
  mutated_sequence  dropped (wild types live in library_wt_sequences.json)
  mutant            ProteinGym notation over the 16 library positions, wild-type
                    positions OMITTED, tokens in ascending position order
  wild-type row     dropped (the parent is not a library member)
  DMS_score         h1_mean verbatim, repr(float) so it round-trips
  DMS_score_bin     1 when DMS_score >= the germline value, else 0
  stage             'train' on every row
  row order         DMS_score descending, ties broken by `mutant` ascending

PARSER VALIDATION
-----------------
`validate` reproduces the two pre-registered reference figures from the committed TrpB and
GB1 files before anything else runs: TrpB shallow 168/2,242 = 7.4933% and GB1 shallow
72/2,167 = 3.3226%, under the convention that depth 0 is excluded from both shells but
retained in the eligible universe. Genotypes are canonicalised to full library length with
wild type filled in at omitted positions, and equal length is ASSERTED - never zipped,
never truncated.

Usage
-----
    .venv\\Scripts\\python.exe -X utf8 -I scripts\\data\\prepare_cr9114_h1.py validate
    .venv\\Scripts\\python.exe -X utf8 -I scripts\\data\\prepare_cr9114_h1.py prepare --source <csv> [--out <csv>]
    .venv\\Scripts\\python.exe -X utf8 -I scripts\\data\\prepare_cr9114_h1.py stats [--file <csv>]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

REPO = Path(__file__).resolve().parent.parent.parent
DATA = REPO / "scripts" / "data"
TARGET = DATA / "CR9114_H1_Phillips_2021.csv"

SOURCE_SHA256 = "ee31f8fc26fce2730d2ae772e7a095aa334b4c6c5d8e7fbfe4b4a5a1c9a794fe"
SOURCE_URL = ("https://cdn.elifesciences.org/articles/71393/"
              "elife-71393-fig1-data1-v3.csv")

# The 16 library positions, 1-based over the 121-residue CR9114 heavy-chain variable
# domain, in the order the source data's `genotype` bitstring uses (left to right).
# Recovered from CombinGym's clean workbook by taking the sequence indices that vary,
# and cross-checked against the wild-type FASTA: every parent residue below equals the
# FASTA residue at that position.
CR9114_POSITIONS: Tuple[int, ...] = (
    29, 30, 31, 52, 57, 58, 59, 71, 74, 75, 76, 77, 84, 87, 95, 106)
CR9114_WT = "FSSITANTKSTSSRYY"
CR9114_ALT = "SNNSSTASIFSNNTFS"
CR9114_GERMLINE_SCORE = 8.424622015254

# Reference landscapes the parser must reproduce before it may be trusted.
VALIDATION = {
    "TrpB": {"file": "TRPB_THEMA_Johnston_2024.csv", "sites": 4,
             "shell_n": 2242, "shell_hits": 168, "shell_rate_pct": 7.4933},
    "GB1": {"file": "GB1_STRSG_Wu_2016.csv", "sites": 4,
            "shell_n": 2167, "shell_hits": 72, "shell_rate_pct": 3.3226},
}


# ======================================================================================
# Mutation-string handling. The one rule that matters: strings OMIT wild-type positions.
# ======================================================================================
def tokens(mutant: str) -> Optional[List[Tuple[int, str, str]]]:
    """'F29S:Y106S' -> [(29,'F','S'), (106,'Y','S')]. 'WT' / '' -> []."""
    s = (mutant or "").strip()
    if not s or s.upper() == "WT":
        return []
    out: List[Tuple[int, str, str]] = []
    for t in s.split(":"):
        if not t:
            continue
        digits = "".join(c for c in t if c.isdigit())
        if not digits or len(t) < 3:
            return None
        out.append((int(digits), t[0], t[-1]))
    return out


def canonical(mutant: str, positions: Sequence[int], wt: Dict[int, str]
              ) -> Optional[Tuple[str, ...]]:
    """Residue at EVERY library position, wild type filled in where the string omits one.

    A 16-site genotype may present as three tokens. Returning a tuple of fixed length
    len(positions) is what makes Hamming distance well defined; comparing raw token
    lists and zipping them silently truncates and understates distance. This inverted a
    conclusion earlier in this project, so the length is asserted below, not assumed.
    """
    toks = tokens(mutant)
    if toks is None:
        return None
    seen: Dict[int, str] = {}
    for pos, _w, mut in toks:
        if pos not in wt:
            return None                       # touches a non-library position
        seen[pos] = mut
    g = tuple(seen.get(p, wt[p]) for p in positions)
    assert len(g) == len(positions), (
        f"canonical genotype length {len(g)} != library size {len(positions)}")
    return g


def hamming(a: Sequence[str], b: Sequence[str]) -> int:
    if len(a) != len(b):
        raise ValueError(f"hamming on unequal lengths: {len(a)} vs {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


def depth_of(mutant: str) -> Optional[int]:
    toks = tokens(mutant)
    return None if toks is None else len(toks)


def greedy_diverse(hits: Sequence[Tuple[float, Tuple[str, ...]]], d: int) -> int:
    """Hits in fitness-descending order; keep one only when >= d from all kept."""
    kept: List[Tuple[str, ...]] = []
    for _, g in hits:
        if all(hamming(g, k) >= d for k in kept):
            kept.append(g)
    return len(kept)


# ======================================================================================
# Library inference, for the validation landscapes
# ======================================================================================
def infer_library(path: Path) -> Tuple[Tuple[int, ...], Dict[int, str]]:
    """The combinatorial positions and their wild-type residues, read from the data.

    The library sites are the heavily varied ones; stray single-position DMS rows that
    touch anything else are not library members. Same 25%-of-max rule as
    scripts/analysis/diversity.py so the two agree by construction.
    """
    counts: Counter = Counter()
    wt: Dict[int, str] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            toks = tokens(str(r.get("mutant", "")))
            if not toks:
                continue
            for pos, w, _m in toks:
                counts[pos] += 1
                wt.setdefault(pos, w)
    if not counts:
        return (), {}
    cut = max(counts.values()) * 0.25
    pos = tuple(sorted(p for p, n in counts.items() if n >= cut))
    return pos, {p: wt[p] for p in pos}


# ======================================================================================
# Statistics
# ======================================================================================
def landscape_stats(path: Path, sites: int,
                    positions: Optional[Tuple[int, ...]] = None,
                    wt: Optional[Dict[int, str]] = None,
                    distances: Sequence[int] = (2, 3)) -> Dict[str, object]:
    """Top-1% tier, shell enrichment at two cuts, and diversity-adjusted ceilings.

    Conventions, fixed to match scripts/freeze_cross_protein_corpus.py:
      - the eligible universe is every row with a parseable score;
      - the tier is the top ceil(0.01 * n) rows by score descending;
      - depth 0 is RETAINED in the universe and the tier, and EXCLUDED from both shells.
    """
    if positions is None or wt is None:
        positions, wt = infer_library(path)

    vals: List[float] = []
    depths: List[int] = []
    genos: List[Tuple[str, ...]] = []
    outside = 0
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                v = float(r["DMS_score"])
            except (TypeError, ValueError, KeyError):
                continue
            if not math.isfinite(v):
                continue
            m = str(r.get("mutant", ""))
            d = depth_of(m)
            g = canonical(m, positions, wt)
            if d is None or g is None:
                outside += 1
                continue
            vals.append(v)
            depths.append(d)
            genos.append(g)

    n = len(vals)
    order = sorted(range(n), key=lambda i: (-vals[i], genos[i]))
    tier = math.ceil(0.01 * n)
    hits = order[:tier]
    hit_set = set(hits)
    threshold = vals[hits[-1]]
    base = tier / n

    def shell(cut: int) -> Dict[str, object]:
        mem = [i for i, d in enumerate(depths) if 1 <= d <= cut]
        deep = [i for i, d in enumerate(depths) if d > cut]
        h = sum(1 for i in mem if i in hit_set)
        dh = sum(1 for i in deep if i in hit_set)
        rate = (h / len(mem)) if mem else None
        return {"cut": cut, "n": len(mem), "hits": h, "rate": rate,
                "enrichment": (rate / base) if rate is not None and base else None,
                "deep_n": len(deep), "deep_hits": dh,
                "deep_rate": (dh / len(deep)) if deep else None}

    ranked = sorted(((vals[i], genos[i]) for i in hits), key=lambda p: (-p[0], p[1]))
    out: Dict[str, object] = {
        "file": path.name, "n": n, "rows_outside_library": outside,
        "sites": sites, "positions": positions,
        "tier_size": tier, "threshold": threshold, "base_rate": base,
        "min": min(vals), "max": max(vals),
        "at_floor_7": sum(1 for v in vals if v == 7.0),
        "depth0_rows": sum(1 for d in depths if d == 0),
        "mean_hit_depth": sum(depths[i] for i in hits) / tier,
        "max_depth": max(depths),
        "primary": shell(sites // 2),
        "fixed_d2": shell(2),
        "depth_profile": {
            d: (sum(1 for i, dd in enumerate(depths) if dd == d),
                sum(1 for i, dd in enumerate(depths) if dd == d and i in hit_set))
            for d in range(0, max(depths) + 1)},
    }
    for dd in distances:
        out[f"ceiling_d{dd}"] = greedy_diverse(ranked, dd)
    return out


# ======================================================================================
# Commands
# ======================================================================================
def cmd_validate(args: argparse.Namespace) -> int:
    print("=" * 96)
    print("PARSER VALIDATION - reproduce the two pre-registered reference figures")
    print("=" * 96)
    print("  Convention: depth 0 retained in the eligible universe and the tier,")
    print("  excluded from both the shallow and the deep shell.\n")
    ok = True
    for name, spec in VALIDATION.items():
        path = DATA / str(spec["file"])
        if not path.exists():
            print(f"  ! {name}: {path} absent")
            ok = False
            continue
        s = landscape_stats(path, int(spec["sites"]))
        sh = s["primary"]                                      # type: ignore[index]
        got_n, got_h = sh["n"], sh["hits"]                     # type: ignore[index]
        got_pct = sh["rate"] * 100                             # type: ignore[operator]
        hit = (got_n == spec["shell_n"] and got_h == spec["shell_hits"]
               and abs(got_pct - float(spec["shell_rate_pct"])) < 5e-5)
        ok &= hit
        print(f"  {name:<6} universe {s['n']:>8,}  tier {s['tier_size']:>6,}  "
              f"shell {got_h}/{got_n} = {got_pct:.4f}%   "
              f"expected {spec['shell_hits']}/{spec['shell_n']} = "
              f"{spec['shell_rate_pct']}%   {'PASS' if hit else 'FAIL'}")
        print(f"         enrichment {sh['enrichment']:.4f}x   "        # type: ignore[index]
              f"deep rate {sh['deep_rate'] * 100:.4f}%   "             # type: ignore[index]
              f"library {s['positions']}")
    print()
    print("  => " + ("PARSER VALIDATED - downstream numbers may be trusted"
                     if ok else "PARSER WRONG - do not trust anything downstream"))
    return 0 if ok else 1


def cmd_prepare(args: argparse.Namespace) -> int:
    src: Path = args.source
    if not src.exists():
        print(f"source not found: {src}", file=sys.stderr)
        print(f"download it from {SOURCE_URL}", file=sys.stderr)
        return 2
    digest = hashlib.sha256(src.read_bytes()).hexdigest()
    print(f"source {src}")
    print(f"  {src.stat().st_size:,} bytes   sha256 {digest}")
    if digest != SOURCE_SHA256:
        print(f"  ! sha256 does not match the recorded source "
              f"({SOURCE_SHA256})", file=sys.stderr)
        if not args.allow_unpinned:
            return 2
    else:
        print("  sha256 matches the recorded eLife source data")

    wt_at = dict(zip(CR9114_POSITIONS, CR9114_WT))
    alt_at = dict(zip(CR9114_POSITIONS, CR9114_ALT))

    with open(src, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        rows = list(rdr)
    print(f"  {len(rows):,} source rows")
    assert len(rows) == 65_536, f"expected 65,536 source rows, got {len(rows):,}"

    prepared: List[Tuple[float, str]] = []
    blank = 0
    wt_rows = 0
    for r in rows:
        bits = r["genotype"].strip()
        assert len(bits) == 16, f"genotype {bits!r} is not 16 bits"
        raw = (r.get("h1_mean") or "").strip()
        if not raw:
            blank += 1
            continue
        v = float(raw)
        assert math.isfinite(v)
        # som_mut is the source's own mutation count; it must equal the bit count.
        assert int(r["som_mut"]) == bits.count("1"), (
            f"som_mut {r['som_mut']} != bit count for {bits}")
        if bits.count("1") == 0:
            wt_rows += 1
            continue                           # the parent is not a library member
        mutant = ":".join(
            f"{wt_at[p]}{p}{alt_at[p]}"
            for p, b in zip(CR9114_POSITIONS, bits) if b == "1")
        prepared.append((v, mutant))

    print(f"  {blank:,} rows with no h1_mean (not quantified) dropped")
    print(f"  {wt_rows} wild-type row dropped (repo convention)")
    print(f"  {len(prepared):,} rows prepared")

    prepared.sort(key=lambda p: (-p[0], p[1]))
    out: Path = args.out
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["mutant", "DMS_score", "DMS_score_bin", "stage"])
        for v, mutant in prepared:
            w.writerow([mutant, repr(v),
                        1 if v >= CR9114_GERMLINE_SCORE else 0, "train"])
    print(f"\nwrote {out}")
    print(f"  {out.stat().st_size:,} bytes   "
          f"sha256 {hashlib.sha256(out.read_bytes()).hexdigest()}")
    return 0


def cmd_stats(args: argparse.Namespace) -> int:
    wt_at = dict(zip(CR9114_POSITIONS, CR9114_WT))
    for path in args.file:
        if not path.exists():
            print(f"  ! {path} absent")
            continue
        s = landscape_stats(path, 16, CR9114_POSITIONS, wt_at)
        print("=" * 96)
        print(f"{path}")
        print("=" * 96)
        print(f"  variants (eligible universe)   {s['n']:,}")
        print(f"  rows outside the 16-site library {s['rows_outside_library']:,}")
        print(f"  score range                    [{s['min']:.6g}, {s['max']:.6g}]"
              f"   at the 7.0 floor: {s['at_floor_7']:,}")
        print(f"  top-1% tier size (hits)        {s['tier_size']:,}")
        print(f"  top-1% threshold               {s['threshold']:.14g}")
        print(f"  base rate                      {s['base_rate'] * 100:.4f}%")
        print(f"  mean depth of a top-1% hit     {s['mean_hit_depth']:.4f}")
        print(f"  max depth observed             {s['max_depth']}")
        print(f"  diversity-adjusted ceiling     d>=2: {s['ceiling_d2']}   "
              f"d>=3: {s['ceiling_d3']}")
        for key, label in (("primary", "PRIMARY  depth <= floor(16/2) = 8"),
                           ("fixed_d2", "FIXED    depth <= 2")):
            sh = s[key]                                        # type: ignore[index]
            print(f"\n  {label}")
            print(f"    shallow shell  {sh['hits']}/{sh['n']:,} = "          # type: ignore[index]
                  f"{sh['rate'] * 100:.4f}%" if sh['rate'] is not None     # type: ignore[index]
                  else f"    shallow shell  {sh['hits']}/{sh['n']:,}")      # type: ignore[index]
            if sh["enrichment"] is not None:                   # type: ignore[index]
                print(f"    enrichment     {sh['enrichment']:.4f}x")        # type: ignore[index]
            print(f"    deep shell     {sh['deep_hits']}/{sh['deep_n']:,} = "  # type: ignore[index]
                  f"{sh['deep_rate'] * 100:.4f}%")             # type: ignore[index]
        print("\n  depth profile  (hits / variants)")
        for d, (nv, nh) in sorted(s["depth_profile"].items()):  # type: ignore[union-attr]
            if nv:
                print(f"    {d:>3}  {nh:>5,} / {nv:>7,}  = {nh / nv * 100:6.3f}%")
        print()
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    v = sub.add_parser("validate", help="Reproduce the TrpB and GB1 reference figures.")
    v.set_defaults(func=cmd_validate)

    pr = sub.add_parser("prepare", help="Write the CSV from the eLife source data.")
    pr.add_argument("--source", type=Path, required=True,
                    help="elife-71393-fig1-data1-v3.csv")
    pr.add_argument("--out", type=Path, default=TARGET)
    pr.add_argument("--allow-unpinned", action="store_true",
                    help="Proceed even if the source sha256 does not match.")
    pr.set_defaults(func=cmd_prepare)

    st = sub.add_parser("stats", help="Report the pre-registered statistics.")
    st.add_argument("--file", type=Path, action="append", default=None)
    st.set_defaults(func=cmd_stats)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if getattr(args, "file", None) is None and args.cmd == "stats":
        args.file = [TARGET]
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
