"""Drop `mutated_sequence` from the new dataset CSVs, preserving the WT sequences first.

The six new landscapes total 43.2 MB, almost all of it the `mutated_sequence` column -
CreiLOV alone is 29.4 MB. TrpB's committed CSV omits that column entirely, so dropping
it is convention-compliant.

But the column is also the only place the full wild-type sequences appear. Any sequence
can be reconstructed from `mutant` plus the wild type, so the wild type is extracted to
a small JSON first and the column then removed. Nothing is lost; it is just stored once
instead of once per row.

Recovering the wild type: take any row, read its own `mutated_sequence`, and reverse its
substitutions using each token's stated wild-type residue (`V183A` says position 183 is
`V` in the parent). Every row must agree, which is the correctness check - if two rows
imply different parents, the file is malformed and nothing is written for it.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SRC = Path(r"C:\Users\stoof\OneDrive\Documents\GitHub\VeritasBio\.claude\worktrees"
           r"\agent-af3c5e2352cf9f840\scripts\data")
OUT = Path(r"C:\Users\stoof\AppData\Local\Temp\claude"
           r"\c--Users-stoof-OneDrive-Documents-GitHub-VeritasBio"
           r"\0691478f-a671-445d-b8dc-a7025937cb24\scratchpad\slim")

FILES = [
    "CREILOV_CHLRE_Chen_2023.csv",
    "CR9114_H1_Phillips_2021.csv",
    "EQFP611_BLUE_Poelwijk_2019.csv",
    "EQFP611_RED_Poelwijk_2019.csv",
    "CR6261_H1_Phillips_2021.csv",
    "CR6261_H9_Phillips_2021.csv",
]

OUT.mkdir(parents=True, exist_ok=True)
wt_record: Dict[str, Dict[str, object]] = {}
csv.field_size_limit(10_000_000)


def tokens(mutant: str) -> Optional[List[Tuple[int, str, str]]]:
    s = mutant.strip()
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


for name in FILES:
    src = SRC / name
    if not src.exists():
        print(f"  ! missing {name}", file=sys.stderr)
        continue

    parent: Optional[str] = None
    disagreements = 0
    offset: Optional[int] = None
    rows_in = 0

    with open(src, newline="", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        cols = rdr.fieldnames or []
        if "mutated_sequence" not in cols:
            print(f"  - {name}: no mutated_sequence column, skipping")
            continue
        for r in rdr:
            rows_in += 1
            seq = (r.get("mutated_sequence") or "").strip()
            toks = tokens(str(r.get("mutant", "")))
            if not seq or toks is None:
                continue
            # Determine the index convention once, from a row with substitutions.
            if offset is None and toks:
                for cand in (1, 0):
                    if all(0 <= p - cand < len(seq) and seq[p - cand] == mut
                           for p, _w, mut in toks):
                        offset = cand
                        break
                if offset is None:
                    continue
            if offset is None:
                continue
            chars = list(seq)
            ok = True
            for p, w, mut in toks:
                i = p - offset
                if not (0 <= i < len(chars)) or chars[i] != mut:
                    ok = False
                    break
                chars[i] = w           # revert to the parent residue
            if not ok:
                continue
            implied = "".join(chars)
            if parent is None:
                parent = implied
            elif implied != parent:
                disagreements += 1

    if parent is None:
        print(f"  ! {name}: could not recover a wild-type sequence", file=sys.stderr)
        continue
    if disagreements:
        print(f"  ! {name}: {disagreements} rows imply a different parent - NOT written",
              file=sys.stderr)
        continue

    wt_record[name] = {"wild_type_sequence": parent,
                       "length": len(parent),
                       "position_offset": offset,
                       "rows": rows_in}

    # Rewrite without the column.
    dst = OUT / name
    with open(src, newline="", encoding="utf-8") as fin, \
            open(dst, "w", newline="", encoding="utf-8") as fout:
        rdr = csv.DictReader(fin)
        keep = [c for c in (rdr.fieldnames or []) if c != "mutated_sequence"]
        wtr = csv.DictWriter(fout, fieldnames=keep)
        wtr.writeheader()
        for r in rdr:
            wtr.writerow({k: r.get(k, "") for k in keep})

    before, after = src.stat().st_size, dst.stat().st_size
    print(f"  {name:<36} {before/1024/1024:>7.1f} MB -> {after/1024/1024:>6.1f} MB"
          f"   WT len {len(parent)}, offset {offset}, {rows_in:,} rows")

(OUT / "library_wt_sequences.json").write_text(
    json.dumps(wt_record, indent=2, sort_keys=True), encoding="utf-8")
total_before = sum((SRC / n).stat().st_size for n in FILES if (SRC / n).exists())
total_after = sum((OUT / n).stat().st_size for n in FILES if (OUT / n).exists())
print(f"\n  total {total_before/1024/1024:.1f} MB -> {total_after/1024/1024:.1f} MB")
print(f"  wrote {OUT / 'library_wt_sequences.json'}")
