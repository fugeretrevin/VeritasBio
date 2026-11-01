"""Where do the top-1% hits live, by mutation depth, versus the generator's depth-4 ceiling?

Mutation strings omit wild-type positions (HANDOFF rule 2.1), so depth is the count of
listed substitutions. Asserted below rather than assumed: no listed substitution may be
a synonymous wt->wt, and max depth may not exceed the library's position count.

Absent labels are excluded per HANDOFF rule 2.3, and *which* landscapes have them is read
from `veritas.calibration.landscapes.REGISTRY` rather than hardcoded here. HER2 is the
only landscape that declares one (`unlabelled_placeholder=0.0`, 25,114 rows): those rows
are absent labels, not measurements of zero, and including them reports a 36,391-row
universe with a 364-variant top-1% tier instead of the canonical 11,277 / 113.

That mattered: this script previously carried no exclusion at all, so it printed 364 while
`docs/day1-reachability.md` §4 - the document it is cited to reproduce - printed 113. The
prose had been corrected and the script had not. Deriving the rule from the registry is
what stops the two drifting again; GB1's 29,477 zero rows are NOT excluded, because no
rule mandates it and the registry records no placeholder for GB1.
"""
import collections
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from veritas.calibration.landscapes import REGISTRY  # noqa: E402

DATA = Path("scripts/data")
CEILING = 4

#: (label, registry key, parent functionality for the measured phenotype).
LANDSCAPES = [
    ("TrpB",         "trpb",         "functional"),
    ("GB1",          "gb1",          "functional"),
    ("CreiLOV",      "creilov",      "functional"),
    ("eqFP611-blue", "eqfp611_blue", "functional"),
    ("eqFP611-red",  "eqfp611_red",  "NON-functional"),
    ("CR9114-H1",    "cr9114_h1",    "NON-functional"),
    ("CR6261-H1",    "cr6261_h1",    "NON-functional"),
    ("CR6261-H9",    "cr6261_h9",    "NON-functional"),
    ("HER2",         "her2",         "functional"),
]


def load(label, key):
    """[(depth, score)], distinct positions, and the count of absent labels dropped.

    One loader, used by both tables below. They were separate loops once, which is how a
    correction reached one and not the other.
    """
    spec = REGISTRY[key]
    path = Path(spec.path)
    if not path.exists():
        path = DATA / Path(spec.path).name
    if not path.exists():
        return None, None, None, path

    placeholder = getattr(spec, "unlabelled_placeholder", None)
    rows, positions, dropped = [], set(), 0
    with path.open(encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh):
            mut = r["mutant"].strip()
            try:
                score = float(r["DMS_score"])
            except (TypeError, ValueError, KeyError):
                continue
            # Rule 2.3: an exact placeholder is an absent label, not a measurement.
            if placeholder is not None and score == placeholder:
                dropped += 1
                continue
            subs = [s for s in mut.split(":") if s]
            # Guard the depth definition: reject a synonymous entry, which would mean
            # wild-type positions are NOT omitted and depth != len(subs).
            for s in subs:
                if len(s) >= 3 and s[0] == s[-1]:
                    raise SystemExit(
                        f"{label}: synonymous substitution {s!r} in {mut!r} - wild-type "
                        "positions are not omitted; depth is not the count")
                if len(s) >= 3:
                    positions.add(s[1:-1])
            rows.append((len(subs), score))
    return rows, positions, dropped, path


def top1_threshold(rows):
    n = len(rows)
    k = max(1, round(n * 0.01))
    return sorted((s for _, s in rows), reverse=True)[k - 1]


loaded = {}
print(f"{'landscape':<14} {'parent':<15} {'N':>7} {'pos':>4} {'top1%':>6} "
      f"{'reachable(<=4)':>14} {'hits<=4':>8} {'hits>4':>7} {'% reachable':>12}")
print("-" * 104)

for label, key, parent in LANDSCAPES:
    rows, positions, dropped, path = load(label, key)
    if rows is None:
        print(f"{label:<14} MISSING {path}")
        continue
    if not rows:
        print(f"{label:<14} no usable rows in {path}")
        continue
    loaded[label] = (rows, parent, dropped, key)

    thresh = top1_threshold(rows)
    hits = [d for d, s in rows if s >= thresh]
    npos = len(positions)
    maxd = max(d for d, _ in rows)
    if maxd > npos:
        raise SystemExit(f"{label}: max depth {maxd} exceeds {npos} positions")

    reach = sum(1 for d, _ in rows if d <= CEILING)
    h_in = sum(1 for d in hits if d <= CEILING)
    h_out = len(hits) - h_in
    pct = 100.0 * h_in / len(hits) if hits else 0.0
    print(f"{label:<14} {parent:<15} {len(rows):>7} {npos:>4} {len(hits):>6} "
          f"{reach:>14} {h_in:>8} {h_out:>7} {pct:>11.1f}%")

for label, (_rows, _parent, dropped, key) in loaded.items():
    if dropped:
        print(f"\nnote: {label}: {dropped:,} row(s) dropped as absent labels "
              f"(rule 2.3, placeholder {REGISTRY[key].unlabelled_placeholder})")

print()
print("Hit-depth histogram for the deep landscapes (top-1% only):")
for label, key, parent in LANDSCAPES:
    if label not in loaded:
        continue
    rows, parent, _dropped, _key = loaded[label]
    if max(d for d, _ in rows) <= CEILING:
        continue
    thresh = top1_threshold(rows)
    hist = collections.Counter(d for d, s in rows if s >= thresh)
    tot = sum(hist.values())
    bars = " ".join(f"d{d}:{hist[d]}" for d in sorted(hist))
    print(f"  {label:<14} ({parent:<14}) n_hits={tot:<5} {bars}")
