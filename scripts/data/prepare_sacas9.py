"""Prepare the reserved KKH-SaCas9 landscape from the pinned CombinGym corpus.

WHY THIS EXISTS
---------------
This project holds no non-development landscape among its measured results: rule 2.4
has fired on HER2 and, on 2026-11-01, on CreiLOV. A new landscape must be reserved
before any generality claim is made again, and `docs/corpus-audit.md` §3 designates
SaCas9 `reserve` on the strength of a rare property: 1,296 of 1,296 rows, zero
dropped, **combinatorially complete** over a mixed 2.2.3.3.3.2.2.3 design. Every
genotype in the design space is present, at every depth, in both the clean and the raw
table, so there is no room for a filter to hide.

`docs/prereg-reserved-sacas9.md` is the sealing document. This script is its
instrument, and nothing here runs a method, an arm or a campaign.

THE LABEL, AND WHY IT IS ONE NAMED sgRNA AND NOT CombinGym's `Mean`
--------------------------------------------------------------------
`docs/corpus-audit.md` §5.2 records that CombinGym's clean workbook reports a `Mean`
column which **silently averages three of the six sgRNAs the raw table measures**
(sg1, sg2, sg3; the raw table also carries sg5, sg6, sg7). Two measurements made here
settle what the label should be, and both are assay structure rather than any method's
behaviour:

1. **The six readouts are close to rank-independent.** Spearman rho over the 1,295
   mutants is 0.126-0.151 among sg1/sg2/sg3 and |rho| <= 0.079 between the sg1-3 group
   and the sg5-7 group; within sg5/sg6/sg7 it is |rho| <= 0.018. "SaCas9 activity" is
   therefore not one latent phenotype - it is six nearly independent landscapes over one
   genotype set - and a mean over them is a composite no assay measured, correlated only
   0.32-0.52 with each of its own constituents.

2. **On sg1, sg2 and sg3 the parent is already at the top, so there is nothing to
   optimise.** The KKH parent sits at the 99.77th, 100.00th and 99.61st percentile of
   the 1,295 mutants; exactly 3, 0 and 5 mutants improve on it, and on sg2 the parent
   exceeds *every* mutant in the library. Every member of a top-5% tier defined on those
   readouts is WORSE than the starting enzyme. A benchmark that rewards finding variants
   worse than wild type is not a benchmark, so the clean `Mean` phenotype is unusable as
   an optimisation endpoint independently of being undocumented.

   On sg5, sg6 and sg7 the parent sits at the 61.78th, 21.85th and 27.41st percentile,
   495 / 1,012 / 940 mutants improve on it, and **all 65 members of the 5% tier beat the
   parent on all three**. Those are the readouts the Thean et al. library is actually
   for: making a PAM-relaxed Cas9 active on target sites the parent handles poorly.

The registered label is therefore **sg5 (PAM TGCCGT), log2 enrichment, read verbatim
from the pinned RAW table** - the readout with the largest parent-relative headroom
(0.8737 log2 units against sg6's 0.7409 and sg7's 0.5353), which is the clearest
optimisation problem of the three and is a property of the assay, declared before any
method ran. sg6 and sg7 are the pre-declared second and third reserve; `--readout`
prepares either without new code.

Reading the raw table rather than the clean one changes WHAT IS READ, not the pin:
`Data/DMS/Raw/SaCas9_raw.xlsx` is already a pinned entry of
`scripts/data/cross_protein_corpus.lock.json` at commit `c5dccd514a9a`. The lock needs
no new entry and none is added. Clean sg1-3 and raw sg1-3 were checked to agree at
**exactly 0.0** on all 1,296 rows, so nothing turns on the table either - only on the
provenance being the one table that carries all six columns.

FOUR-WAY LIBRARY VALIDATION (`validate`)
----------------------------------------
`docs/dataset-inventory.md` §1 cross-validated the six CombinGym landscapes' positions
and parent residues against three sources. The same check runs here against four, and
all four must agree before anything is written - a silent mis-mapping here is the
failure mode that made HER2 read the trastuzumab light chain:

  1. CombinGym's registry, `Data/Data_summary.xlsx`: `L887, N888, A889, N985, N986,
     L988, L989, R991`, i.e. positions (887, 888, 889, 985, 986, 988, 989, 991) and
     parent `LNANNLLR`;
  2. the positions that actually VARY across all 1,296 clean `genotype` strings;
  3. the wild-type FASTA `SaCas9(SaCas9-KKH).fasta` read at those 1-based positions -
     and, as a closed loop, the single `n_mut = 0` row is byte-identical to the FASTA;
  4. the RAW table's own column header `887-888-889-985-986-988-989-991` and its parent
     row `L-N-A-N-N-L-L-R`.

Plus four invariants: `n_mut` equals Hamming-distance-to-parent on all 1,296 rows; no
row varies at any of the other 1,045 sequence positions; the per-site alphabets multiply
to 1,296 which equals the distinct genotype count which equals the row count; and the
clean `Mean` reproduces `mean(sg1, sg2, sg3)` while differing from the mean over all six
on every row, which is `corpus-audit.md` §5.2 re-derived here.

PARSER VALIDATION, before any of it
-----------------------------------
`validate` first reproduces the two pre-registered reference figures from the committed
TrpB and GB1 tables - TrpB shallow 168/2,242 = 7.4933% and GB1 shallow 72/2,167 =
3.3226%, depth 0 excluded from both shells and retained in the universe. Genotypes are
canonicalised through `veritas.calibration.genotypes`, which fills the parent residue at
omitted positions and ASSERTS equal length; nothing here zips or truncates, and no sixth
canonicalisation is written.

CONVENTIONS, matched to the committed CR9114 and CR6261 files
-------------------------------------------------------------
  columns           mutant,DMS_score,DMS_score_bin,stage
  mutant            ProteinGym notation over the 8 library positions, wild-type
                    positions OMITTED, tokens in ascending position order
  wild-type row     dropped (the parent is not a library member) -> 1,295 rows
  DMS_score         the named readout verbatim, repr(float) so it round-trips
  DMS_score_bin     1 when DMS_score >= the parent's own value on that readout, else 0
  stage             'train' on every row
  row order         DMS_score descending, ties broken by `mutant` ascending

Usage
-----
    .venv\\Scripts\\python.exe -X utf8 -I scripts\\data\\prepare_sacas9.py validate
    .venv\\Scripts\\python.exe -X utf8 -I scripts\\data\\prepare_sacas9.py prepare
    .venv\\Scripts\\python.exe -X utf8 -I scripts\\data\\prepare_sacas9.py readouts
    .venv\\Scripts\\python.exe -X utf8 -I scripts\\data\\prepare_sacas9.py stats

`--corpus-dir` must point at the CORPUS DIRECTORY ITSELF - `<real repo root>\\downloads
\\combingym`, not the repository root, which still fails. It defaults to the real
repository root's copy via `veritas.paths`, because `downloads/` is gitignored and
therefore invisible from a git worktree.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import math
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from veritas import paths  # noqa: E402
from veritas.calibration.genotypes import (  # noqa: E402
    Combo, canonical_combo, depth, hamming,
)

#: The three agreeing sources, relative to the corpus directory.
CLEAN_REL = "Data/DMS/Clean/SaCas9_clean.xlsx"
RAW_REL = "Data/DMS/Raw/SaCas9_raw.xlsx"
FASTA_REL = "Data/WT Sequence/SaCas9(SaCas9-KKH).fasta"
SUMMARY_REL = "Data/Data_summary.xlsx"

#: Every sgRNA the raw table reports, in its own column order.
RAW_SGRNAS: Tuple[str, ...] = ("sg1", "sg2", "sg3", "sg5", "sg6", "sg7")
#: The three the clean workbook keeps (corpus-audit.md §5.2).
CLEAN_SGRNAS: Tuple[str, ...] = ("sg1", "sg2", "sg3")

#: The registered readout. See the module docstring: sg1/sg2/sg3 are unusable because
#: the parent already sits at the 99.6th-100th percentile on them.
REGISTERED_READOUT = "sg5"
#: Pre-declared second and third reserve, in order.
RESERVE_READOUTS: Tuple[str, ...] = ("sg6", "sg7")

#: `{readout: (registry key, prepared table name)}` for the three usable readouts.
PREPARED: Dict[str, Tuple[str, str]] = {
    "sg5": ("sacas9_sg5", "SACAS9_KKH_SG5_Thean_2022.csv"),
    "sg6": ("sacas9_sg6", "SACAS9_KKH_SG6_Thean_2022.csv"),
    "sg7": ("sacas9_sg7", "SACAS9_KKH_SG7_Thean_2022.csv"),
}

#: Pre-registered reference figures. If these do not reproduce, the parser is wrong and
#: nothing downstream is trustworthy (HANDOFF rule 2.2).
ANCHORS = {
    "trpb": ("TRPB_THEMA_Johnston_2024.csv", (183, 184, 227, 228), "VFVS",
             168, 2_242, 7.4933),
    "gb1": ("GB1_STRSG_Wu_2016.csv", (39, 40, 41, 54), "VDGV", 72, 2_167, 3.3226),
}

#: The protocol's own tier ladder and minimum tier size
#: (`scripts/freeze_cross_protein_corpus.py:TIER_LADDER`, `MIN_TIER_MEMBERS`).
TIER_LADDER: Tuple[float, ...] = (0.01, 0.02, 0.05, 0.10)
MIN_TIER_MEMBERS = 30


# ======================================================================================
# The pinned corpus is read through the freeze script's own stdlib .xlsx reader
# ======================================================================================
def _freeze_module():
    """`scripts/freeze_cross_protein_corpus.py`, imported rather than reimplemented.

    Its `read_xlsx` is the reader the freeze and the audit already use, so this script
    cannot disagree with them about what a cell holds, and no new dependency (pandas,
    openpyxl) enters a path that produces a frozen artefact.
    """
    script = paths.repo_path("scripts/freeze_cross_protein_corpus.py")
    spec = importlib.util.spec_from_file_location("_veritas_freeze_corpus", script)
    if spec is None or spec.loader is None:  # pragma: no cover - defensive
        raise ImportError(f"cannot load the freeze script at {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def default_corpus_dir() -> Path:
    """`<real repo root>/downloads/combingym`, naming every path on a failure.

    `downloads/` is gitignored, so a git worktree cannot see it and it must be resolved
    against the REAL repository root (HANDOFF rule 2.8). A bare "not found" here is what
    cost three agent runs.
    """
    candidates = [
        paths.repository_root() / "downloads" / "combingym",
        paths.checkout_root() / "downloads" / "combingym",
    ]
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    raise FileNotFoundError(
        "the pinned CombinGym corpus directory was not found. Checked, in order:\n"
        + "\n".join(f"  - {c}" for c in candidates)
        + "\nPass --corpus-dir pointing at the corpus directory ITSELF (the folder that "
          "holds Data/), not at the repository root."
    )


# ======================================================================================
# 0. Parser validation against the two pre-registered anchors
# ======================================================================================
def _shell_from_prepared(csv_name: str, positions: Sequence[int],
                         parent: str) -> Tuple[int, int]:
    """`(shallow hits, shallow variants)` for `1 <= depth <= 2` on a committed table.

    Depth 0 is excluded from the shell and retained in the universe, which is the
    convention the anchors were computed under.
    """
    table = paths.dataset_dir() / csv_name
    if not table.is_file():
        raise FileNotFoundError(f"{csv_name}: not at {table}")
    values: List[float] = []
    depths: List[int] = []
    with table.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            combo = canonical_combo(row["mutant"], tuple(positions), parent)
            if combo is None:
                continue
            d = depth(combo, parent)
            if d == 0:
                continue
            values.append(float(row["DMS_score"]))
            depths.append(d)
    n = len(values)
    tier = math.ceil(0.01 * n)
    threshold = sorted(values, reverse=True)[tier - 1]
    shallow = [i for i, d in enumerate(depths) if 1 <= d <= 2]
    hits = sum(1 for i in shallow if values[i] >= threshold)
    return hits, len(shallow)


def validate_parser() -> bool:
    print("parser validation (HANDOFF rule 2.2) - reproduce before trusting anything")
    ok = True
    for key, (name, positions, parent, e_hits, e_n, e_pct) in ANCHORS.items():
        hits, n = _shell_from_prepared(name, positions, parent)
        pct = 100.0 * hits / n if n else float("nan")
        good = (hits, n) == (e_hits, e_n) and abs(pct - e_pct) < 5e-4
        ok = ok and good
        print(f"  {'PASS' if good else 'FAIL'}  {key:4} shallow shell "
              f"{hits}/{n} = {pct:.4f}%   expected {e_hits}/{e_n} = {e_pct:.4f}%")
    return ok


# ======================================================================================
# 1. The four sources
# ======================================================================================
class Sources:
    """The four agreeing descriptions of the SaCas9 library, parsed and nothing more."""

    def __init__(self, corpus_dir: Path) -> None:
        fz = _freeze_module()
        self.corpus_dir = corpus_dir
        for rel in (CLEAN_REL, RAW_REL, FASTA_REL, SUMMARY_REL):
            if not (corpus_dir / rel).is_file():
                raise FileNotFoundError(f"missing pinned file: {corpus_dir / rel}")

        # ---- source 1: CombinGym's registry ------------------------------------------
        _s, _h, rows = fz.read_xlsx(corpus_dir / SUMMARY_REL, "Sheet1")
        row = next(r for r in rows if (r[0] or "").strip() == "SaCas9")
        tokens = [t.strip() for t in (row[6] or "").split(",")]
        self.registry_positions: Tuple[int, ...] = tuple(int(t[1:]) for t in tokens)
        self.registry_parent: str = "".join(t[0] for t in tokens)
        self.registry_sites = int(float(row[7]))
        self.registry_seq_len = int(float(row[5]))
        self.registry_variant_space = int(float(row[8]))
        self.registry_dms_rows = int(float(row[9]))

        # ---- source 3: the wild-type FASTA -------------------------------------------
        lines = (corpus_dir / FASTA_REL).read_text(encoding="utf-8").splitlines()
        self.fasta_header = lines[0]
        self.wild_type_sequence = "".join(
            l.strip() for l in lines[1:] if not l.startswith(">"))

        # ---- source 2: the clean workbook --------------------------------------------
        _s, header, rows = fz.read_xlsx(corpus_dir / CLEAN_REL, "Sheet1")
        self.clean_header = header
        at = {h.split(" ")[0]: i for i, h in enumerate(header)}
        self.clean_genotypes: List[str] = []
        self.clean_depths: List[int] = []
        self.clean_mean: List[float] = []
        self.clean_sg: Dict[str, List[float]] = {k: [] for k in CLEAN_SGRNAS}
        for r in rows:
            self.clean_genotypes.append((r[at["genotype"]] or "").strip())
            self.clean_depths.append(int(float(r[at["n_mut"]])))
            self.clean_mean.append(float(r[at["Mean"]]))
            for k in CLEAN_SGRNAS:
                self.clean_sg[k].append(float(r[at[k]]))

        # ---- source 4: the raw table -------------------------------------------------
        # Row 0 is a merged title and row 1 is the real header, so the header row is
        # LOCATED rather than assumed: a reader that took row 0 would silently find no
        # sgRNA column at all.
        _s, _h, rrows = fz.read_xlsx(corpus_dir / RAW_REL, "Sheet1")
        hdr = [(c or "").strip() for c in rrows[1]]
        if not any(h.startswith("sg") for h in hdr):
            raise ValueError(f"{RAW_REL}: no sgRNA column in the header row {hdr}")
        self.raw_header = hdr
        self.raw_positions: Tuple[int, ...] = tuple(int(x) for x in hdr[1].split("-"))
        self.pam: Dict[str, str] = {
            h.split(" ")[0]: h.split("PAM: ")[1].rstrip(")")
            for h in hdr if h.startswith("sg") and "PAM: " in h
        }
        raw_at = {h.split(" ")[0]: i for i, h in enumerate(hdr) if h.startswith("sg")}
        if tuple(sorted(raw_at)) != tuple(sorted(RAW_SGRNAS)):
            raise ValueError(
                f"{RAW_REL}: readouts {sorted(raw_at)} but this script declares "
                f"{sorted(RAW_SGRNAS)}")
        self.raw: Dict[Combo, Dict[str, float]] = {}
        for r in rrows[2:]:
            key = (r[1] or "").strip()
            if not key:
                continue
            self.raw[tuple(key.split("-"))] = {
                k: float(r[i]) for k, i in raw_at.items()}

    # ---- derived -------------------------------------------------------------------
    @property
    def positions(self) -> Tuple[int, ...]:
        return self.registry_positions

    @property
    def parent_residues(self) -> str:
        return self.registry_parent

    @property
    def parent_combo(self) -> Combo:
        return tuple(self.parent_residues)

    def combo_of(self, genotype: str) -> Combo:
        """One residue per library position. Always full length, asserted."""
        combo = tuple(genotype[p - 1] for p in self.positions)
        assert len(combo) == len(self.positions), "combo must be full library length"
        return combo

    @property
    def clean_combos(self) -> List[Combo]:
        return [self.combo_of(g) for g in self.clean_genotypes]

    @property
    def mutant_combos(self) -> List[Combo]:
        """The eligible universe: every genotype but the parent (benchmark rule)."""
        return [c for c in self.clean_combos if c != self.parent_combo]

    def mutant_string(self, combo: Combo) -> str:
        """ProteinGym notation, wild-type positions OMITTED, positions ascending."""
        return ":".join(
            f"{wt}{pos}{mut}"
            for pos, wt, mut in zip(self.positions, self.parent_residues, combo)
            if mut != wt
        )


# ======================================================================================
# 2. The four-way validation
# ======================================================================================
def validate_library(src: Sources) -> bool:
    ok = True

    def check(name: str, good: bool, detail: str) -> None:
        nonlocal ok
        ok = ok and good
        print(f"  {'PASS' if good else 'FAIL'}  {name:<46} {detail}")

    print(f"\nfour-way library validation ({src.corpus_dir})")
    parent_combo = src.parent_combo
    combos = src.clean_combos
    n = len(combos)
    wt = src.wild_type_sequence

    # --- 1 vs 3: registry against the wild-type FASTA ---------------------------------
    spelled = "".join(wt[p - 1] for p in src.positions)
    check("registry positions vs FASTA residues",
          spelled == src.parent_residues,
          f"FASTA spells {spelled!r}; registry declares {src.parent_residues!r}")
    check("FASTA length vs registry sequence length",
          len(wt) == src.registry_seq_len, f"{len(wt)} vs {src.registry_seq_len}")

    # --- 1 vs 2: registry against the positions that actually vary --------------------
    varying = tuple(sorted({
        i + 1 for g in src.clean_genotypes for i in range(len(g)) if g[i] != wt[i]
    }))
    check("varying sequence indices vs registry positions",
          varying == src.positions, f"{list(varying)}")
    outside = sum(
        1 for g in src.clean_genotypes
        if any(g[j] != wt[j] for j in range(len(wt)) if (j + 1) not in set(src.positions))
    )
    check("rows varying outside the library positions", outside == 0, f"{outside} of {n}")

    # --- 1 vs 4: registry against the raw table's own header and parent row -----------
    check("raw table position header vs registry",
          src.raw_positions == src.positions, f"{src.raw_header[1]!r}")
    check("raw table carries the parent row",
          parent_combo in src.raw, f"{'-'.join(parent_combo)}")
    check("raw genotype set == clean genotype set",
          set(src.raw) == set(combos),
          f"raw {len(src.raw)}, clean {len(set(combos))}")

    # --- closed loop: the depth-0 row is the FASTA, byte for byte ---------------------
    identical = [i for i, g in enumerate(src.clean_genotypes) if g == wt]
    check("exactly one row byte-identical to the FASTA",
          len(identical) == 1 and src.clean_depths[identical[0]] == 0,
          f"{len(identical)} row(s), n_mut {[src.clean_depths[i] for i in identical]}")

    # --- n_mut is Hamming-to-parent, on every row -------------------------------------
    bad = sum(1 for i, c in enumerate(combos)
              if hamming(c, parent_combo) != src.clean_depths[i])
    check("n_mut == Hamming-to-parent on every row", bad == 0, f"{bad} of {n} disagree")

    # --- the design space, exactly ----------------------------------------------------
    alphabets = tuple("".join(sorted({c[i] for c in combos}))
                      for i in range(len(src.positions)))
    product = 1
    for a in alphabets:
        product *= len(a)
    check("design-space product == distinct genotypes == rows",
          product == len(set(combos)) == n == src.registry_variant_space,
          f"{' * '.join(str(len(a)) for a in alphabets)} = {product}; "
          f"distinct {len(set(combos))}; rows {n}; registry {src.registry_variant_space}")
    for pos, wt_res, alpha in zip(src.positions, src.parent_residues, alphabets):
        print(f"            {wt_res}{pos}: {len(alpha)} residues -> {alpha}")

    # --- corpus-audit.md §5.2, re-derived ---------------------------------------------
    max_sg = max(abs(src.clean_sg[k][i] - src.raw[c][k])
                 for i, c in enumerate(combos) for k in CLEAN_SGRNAS)
    check("clean sg1-3 == raw sg1-3", max_sg == 0.0, f"max |delta| = {max_sg:.3e}")
    d3 = max(abs(src.clean_mean[i]
                 - sum(src.raw[c][k] for k in CLEAN_SGRNAS) / len(CLEAN_SGRNAS))
             for i, c in enumerate(combos))
    d6 = sum(1 for i, c in enumerate(combos)
             if src.clean_mean[i] != sum(src.raw[c].values()) / len(RAW_SGRNAS))
    check("clean Mean == mean(sg1,sg2,sg3)", d3 < 1e-15, f"max |delta| = {d3:.3e}")
    check("clean Mean != mean(all six) on every row", d6 == n,
          f"{d6} of {n} rows differ (corpus-audit.md §5.2)")

    # --- canonicalisation round-trips through the shared module -----------------------
    bad_round_trip = sum(
        1 for c in combos
        if canonical_combo(src.mutant_string(c), src.positions,
                           src.parent_residues) != c)
    check("mutant-string round trip via veritas.calibration.genotypes",
          bad_round_trip == 0, f"{bad_round_trip} of {n} fail")
    return ok


# ======================================================================================
# 3. readouts - the label decision, as auditable numbers
# ======================================================================================
def _spearman(x: Sequence[float], y: Sequence[float]) -> float:
    def rank(v: Sequence[float]) -> List[float]:
        order = sorted(range(len(v)), key=lambda i: v[i])
        out = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            average = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                out[order[k]] = average
            i = j + 1
        return out

    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def _greedy_diverse(hits: Sequence[Tuple[float, Combo]], d: int) -> int:
    kept: List[Combo] = []
    for _v, combo in sorted(hits, key=lambda p: (-p[0], p[1])):
        if all(hamming(combo, k) >= d for k in kept):
            kept.append(combo)
    return len(kept)


def ladder_tier(n: int) -> Tuple[float, int]:
    """The protocol's own tier: the shallowest ladder rung holding >= 30 members."""
    for fraction in TIER_LADDER:
        size = math.ceil(fraction * n)
        if size >= MIN_TIER_MEMBERS:
            return fraction, size
    fraction = TIER_LADDER[-1]
    return fraction, math.ceil(fraction * n)


def readouts(src: Sources) -> None:
    """Every number the label decision rests on. Assay structure, not performance."""
    mutants = src.mutant_combos
    n = len(mutants)
    depths = {c: depth(c, src.parent_residues) for c in mutants}
    fraction, tier = ladder_tier(n)
    print(f"\neligible universe (parent dropped) = {n}; protocol tier ladder picks "
          f"{fraction:.0%} = {tier} members (>= {MIN_TIER_MEMBERS})")

    print("\nparent-relative headroom - why sg1/sg2/sg3 and any mean over them are out")
    head = (f"  {'sgRNA':6} {'PAM':8} {'parent':>8} {'max':>8} {'headroom':>9} "
            f"{'n>parent':>9} {'parent pctile':>14} {'distinct':>9}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for k in RAW_SGRNAS:
        values = [src.raw[c][k] for c in mutants]
        parent_value = src.raw[src.parent_combo][k]
        above = sum(1 for v in values if v > parent_value)
        print(f"  {k:6} {src.pam.get(k, '?'):8} {parent_value:8.4f} {max(values):8.4f} "
              f"{max(values) - parent_value:9.4f} {above:9d} "
              f"{100.0 * (1 - above / n):13.2f}% {len(set(values)):9d}")

    print("\nrank agreement (Spearman rho over the eligible universe)")
    series = {k: [src.raw[c][k] for c in mutants] for k in RAW_SGRNAS}
    series["mean(sg1-3)"] = [
        sum(src.raw[c][k] for k in CLEAN_SGRNAS) / len(CLEAN_SGRNAS) for c in mutants]
    series["mean(all 6)"] = [
        sum(src.raw[c].values()) / len(RAW_SGRNAS) for c in mutants]
    names = list(series)
    print("  " + " " * 12 + " ".join(f"{k:>12}" for k in names))
    for a in names:
        print(f"  {a:12}" + " ".join(f"{_spearman(series[a], series[b]):12.3f}"
                                    for b in names))

    print(f"\nendpoint geometry at every ladder rung, per readout")
    head2 = (f"  {'sgRNA':6} {'tier':>5} {'raw':>4} {'threshold':>11} {'d>=2':>5} "
             f"{'removed@2':>10} {'d>=3':>5} {'removed@3':>10} {'mean depth':>10} "
             f"{'hits>parent':>12}")
    print(head2)
    print("  " + "-" * (len(head2) - 2))
    for k in RAW_SGRNAS:
        parent_value = src.raw[src.parent_combo][k]
        ordered = sorted((src.raw[c][k] for c in mutants), reverse=True)
        for fr in TIER_LADDER:
            size = math.ceil(fr * n)
            threshold = ordered[size - 1]
            hits = [(src.raw[c][k], c) for c in mutants if src.raw[c][k] >= threshold]
            c2 = _greedy_diverse(hits, 2)
            c3 = _greedy_diverse(hits, 3)
            print(f"  {k:6} {fr:5.0%} {len(hits):4d} {threshold:11.6f} {c2:5d} "
                  f"{1 - c2 / len(hits):10.2%} {c3:5d} {1 - c3 / len(hits):10.2%} "
                  f"{sum(depths[c] for _v, c in hits) / len(hits):10.3f} "
                  f"{sum(1 for v, _c in hits if v > parent_value):12d}")


# ======================================================================================
# 4. prepare
# ======================================================================================
def prepare(src: Sources, readout: str, out: Path) -> Path:
    parent_value = src.raw[src.parent_combo][readout]
    rows = [(src.raw[combo][readout], src.mutant_string(combo))
            for combo in src.mutant_combos]
    rows.sort(key=lambda p: (-p[0], p[1]))

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["mutant", "DMS_score", "DMS_score_bin", "stage"])
        for value, mutant in rows:
            writer.writerow([mutant, repr(value),
                             1 if value >= parent_value else 0, "train"])
    data = out.read_bytes()
    lf = data.replace(b"\r\n", b"\n")
    print(f"\nwrote {out}")
    print(f"  readout   {readout} (PAM {src.pam.get(readout)}), log2 enrichment, "
          f"read verbatim from {RAW_REL}")
    print(f"  rows      {len(rows)} (1,296 genotypes less the parent)")
    print(f"  parent    {parent_value!r}; DMS_score_bin = 1 at or above it "
          f"({sum(1 for v, _m in rows if v >= parent_value)} rows)")
    print(f"  bytes     {len(data):,}")
    print(f"  sha256    {hashlib.sha256(data).hexdigest()}")
    print(f"  sha256 LF {hashlib.sha256(lf).hexdigest()}")
    return out


# ======================================================================================
# 5. stats - the prepared table's own geometry, model-free
# ======================================================================================
def stats(table: Path, positions: Sequence[int], parent: str) -> None:
    combos: List[Combo] = []
    values: List[float] = []
    with table.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            combo = canonical_combo(row["mutant"], tuple(positions), parent, strict=True)
            assert combo is not None and len(combo) == len(positions)
            combos.append(combo)
            values.append(float(row["DMS_score"]))
    n = len(values)
    depths = [depth(c, parent) for c in combos]
    print(f"\n{table.name}: {n} eligible variants, depths {min(depths)}-{max(depths)}")
    print(f"  min {min(values):.6f}  max {max(values):.6f}  distinct {len(set(values))}  "
          f"ties at max {sum(1 for v in values if v == max(values))}  "
          f"ties at min {sum(1 for v in values if v == min(values))}")
    hist = Counter(depths)
    print("  depth histogram: " + ", ".join(f"{d}:{hist[d]}" for d in sorted(hist)))

    head = (f"  {'tier':>6} {'threshold':>11} {'raw':>5} {'d>=2':>5} {'d>=3':>5} "
            f"{'removed@2':>10} {'mean depth':>10}")
    print(head)
    print("  " + "-" * (len(head) - 2))
    for fraction in TIER_LADDER:
        tier = math.ceil(fraction * n)
        threshold = sorted(values, reverse=True)[tier - 1]
        hits = [(v, c) for v, c in zip(values, combos) if v >= threshold]
        hit_depths = [d for v, d in zip(values, depths) if v >= threshold]
        c2 = _greedy_diverse(hits, 2)
        print(f"  {fraction:6.0%} {threshold:11.6f} {len(hits):5d} {c2:5d} "
              f"{_greedy_diverse(hits, 3):5d} {1 - c2 / len(hits):10.2%} "
              f"{sum(hit_depths) / len(hit_depths):10.3f}")

    cut = len(positions) // 2
    shallow = [i for i, d in enumerate(depths) if 1 <= d <= cut]
    fraction, tier = ladder_tier(n)
    threshold = sorted(values, reverse=True)[tier - 1]
    shallow_hits = sum(1 for i in shallow if values[i] >= threshold)
    rate = shallow_hits / len(shallow)
    print(f"  primary shell 1 <= depth <= {cut}: {shallow_hits}/{len(shallow)} = "
          f"{rate:.4%} at the {fraction:.0%} tier, enrichment "
          f"{rate / (tier / n):.4f}x (resolvable at the >=500 floor: "
          f"{len(shallow) >= 500})")


# ======================================================================================
# 6. endpoint - the gate's verdict at the registered protocol and at the amendment
# ======================================================================================
#: The amendment, pre-registered in `docs/prereg-reserved-sacas9.md` §4. The budget is a
#: uniform 0.4x scaling of the registered 25 + 4 x 25 = 125, which preserves both the
#: 20% initial fraction and the five-round structure.
REGISTERED_BUDGET = 50
REGISTERED_INITIAL = 10
REGISTERED_ROUNDS = 4
REGISTERED_DISTANCE = 2


def endpoint(table: Path, positions: Sequence[int], parent: str) -> bool:
    """`check_endpoint_viability` at the registered protocol AND at the amendment.

    Both verdicts are printed because the first one is the finding: the registered
    125-measurement budget cannot discriminate on a 1,295-variant universe, which is
    what forces the amendment rather than a preference for it. The rule is imported from
    `veritas.preflight`, never reimplemented.
    """
    from veritas.preflight import check_endpoint_viability

    combos: List[Combo] = []
    values: List[float] = []
    with table.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            combo = canonical_combo(row["mutant"], tuple(positions), parent, strict=True)
            assert combo is not None and len(combo) == len(positions)
            combos.append(combo)
            values.append(float(row["DMS_score"]))
    n = len(values)

    def at(fraction: float) -> Tuple[float, int, int]:
        size = math.ceil(fraction * n)
        threshold = sorted(values, reverse=True)[size - 1]
        hits = [(v, c) for v, c in zip(values, combos) if v >= threshold]
        return threshold, len(hits), _greedy_diverse(hits, REGISTERED_DISTANCE)

    print(f"\n{table.name}: universe {n}")
    print("\n--- the REGISTERED protocol, unchanged: 1% tier, budget 125 ---")
    threshold, raw, ceiling = at(0.01)
    registered = check_endpoint_viability(
        library="sacas9_sg5", raw=raw, ceiling=ceiling,
        distance=REGISTERED_DISTANCE, budget=125, tier_size=raw)
    print(f"  tier threshold {threshold!r}")
    print(registered.render())

    print(f"\n--- the AMENDMENT: {ladder_tier(n)[0]:.0%} tier (the protocol's own ladder), "
          f"budget {REGISTERED_BUDGET} = {REGISTERED_INITIAL} + "
          f"{REGISTERED_ROUNDS} x {REGISTERED_INITIAL} ---")
    fraction, _size = ladder_tier(n)
    threshold, raw, ceiling = at(fraction)
    amended = check_endpoint_viability(
        library="sacas9_sg5", raw=raw, ceiling=ceiling,
        distance=REGISTERED_DISTANCE, budget=REGISTERED_BUDGET, tier_size=raw)
    print(f"  tier threshold {threshold!r}")
    print(amended.render())
    print(f"\n  budget / tier     = {REGISTERED_BUDGET}/{raw} = "
          f"{REGISTERED_BUDGET / raw:.1%}")
    print(f"  budget / universe = {REGISTERED_BUDGET}/{n} = {REGISTERED_BUDGET / n:.2%}")
    print(f"  uniform-random expectation at this budget = "
          f"{REGISTERED_BUDGET * raw / n:.3f} raw hits")
    print(f"\nregistered protocol: {registered.status}; amendment: {amended.status}")
    return amended.passed


# ======================================================================================
# CLI
# ======================================================================================
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command",
                        choices=("validate", "readouts", "prepare", "stats", "endpoint"))
    parser.add_argument("--readout", default=REGISTERED_READOUT, choices=sorted(PREPARED),
                        help="Which named sgRNA readout to prepare "
                             "(default: %(default)s, the registered one).")
    parser.add_argument("--corpus-dir", type=Path, default=None,
                        help="The CORPUS DIRECTORY ITSELF (the folder holding Data/), "
                             "not the repository root.")
    parser.add_argument("--out", type=Path, default=None,
                        help="Prepared table (default: scripts/data/<readout table>).")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    _key, name = PREPARED[args.readout]
    out = args.out if args.out is not None else paths.dataset_dir() / name

    corpus_dir = (paths.explicit_path(args.corpus_dir, what="--corpus-dir")
                  if args.corpus_dir is not None else default_corpus_dir())

    if args.command in ("stats", "endpoint"):
        if not out.is_file():
            print(f"no prepared table at {out}; run `prepare` first", file=sys.stderr)
            return 2
        src = Sources(corpus_dir)
        if args.command == "stats":
            stats(out, src.positions, src.parent_residues)
            return 0
        return 0 if endpoint(out, src.positions, src.parent_residues) else 1

    if not validate_parser():
        print("\nthe parser does not reproduce the pre-registered anchors; refusing to "
              "read or write anything (HANDOFF rule 2.2)", file=sys.stderr)
        return 1

    src = Sources(corpus_dir)
    if not validate_library(src):
        print("\nthe four sources do not agree on the library definition; refusing to "
              "write a prepared table", file=sys.stderr)
        return 1
    print("\nall four sources agree on the library definition.")

    if args.command == "readouts":
        readouts(src)
    elif args.command == "prepare":
        prepare(src, args.readout, out)
        stats(out, src.positions, src.parent_residues)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
