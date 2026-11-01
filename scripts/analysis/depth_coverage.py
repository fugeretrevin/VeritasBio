"""Is the genotype loss in each prepared table DEPTH-CORRELATED?

`corpus.derived_table` flags seven landscapes that read a CombinGym "clean" workbook
while a raw table exists at the pin. The motivating defect (docs/cr9114-reprep.md) is
that CR9114's clean workbook dropped 16,384 genotypes with attrition that reached 100%
at depth 16, removing the matured antibody the landscape is about.

That failure mode is fatal to a claim about depth reachability specifically: if coverage
falls with depth, an engine that proposes deep variants is penalised by the TABLE, not by
the landscape. So this measures observed-vs-expected coverage per depth for every
landscape in the sweep.

Expected counts come from the design's generating function: a site with k alleles
contributes (1 + (k-1)x), since there is one way to carry the wild-type residue there and
k-1 ways to carry a substitution. The product's x^d coefficient is the number of depth-d
genotypes in the complete design.
"""
import csv
import collections
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import benchmark_panoptes_isolated as bench  # noqa: E402

SWEEP = ["creilov", "eqfp611_red", "eqfp611_blue", "gb1", "cr6261_h1", "cr6261_h9",
         "cr9114_h1", "trpb", "her2"]


def expected_by_depth(spec) -> list:
    """x^d coefficients of prod_sites (1 + (alleles-1) x)."""
    poly = [1]
    for slot in range(spec.n_sites):
        alt = len(spec.alphabet_for(slot)) - 1
        nxt = [0] * (len(poly) + 1)
        for d, c in enumerate(poly):
            nxt[d] += c              # wild type at this site
            nxt[d + 1] += c * alt    # a substitution at this site
        poly = nxt
    return poly


print(f"{'landscape':<14}{'depth':>6}{'observed':>10}{'expected':>10}{'coverage':>10}")
print("-" * 52)
verdicts = []
for key in SWEEP:
    spec = bench.LIBRARIES[key]
    observed = collections.Counter()
    for row in csv.DictReader(open(spec.dataset)):
        mutant = row["mutant"].strip()
        try:
            float(row["DMS_score"])
        except (TypeError, ValueError):
            continue
        observed[len([s for s in mutant.split(":") if s])] += 1

    expect = expected_by_depth(spec)
    coverages = []
    for d in range(len(expect)):
        if expect[d] == 0:
            continue
        cov = observed[d] / expect[d]
        coverages.append((d, observed[d], expect[d], cov))

    for d, obs, exp, cov in coverages:
        flag = "  <-- ZERO" if obs == 0 else ""
        print(f"{key if d == coverages[0][0] else '':<14}{d:>6}{obs:>10}{exp:>10}{cov:>9.1%}{flag}")

    # Depth-correlated? Compare coverage in the shallow half against the deep half,
    # weighted by expected count, and report the worst single depth.
    deep = [c for c in coverages if c[0] > spec.n_sites / 2]
    shallow = [c for c in coverages if 0 < c[0] <= spec.n_sites / 2]

    def wmean(rows):
        num = sum(o for _, o, _, _ in rows)
        den = sum(e for _, _, e, _ in rows)
        return num / den if den else float("nan")

    # Depth 0 is the wild type, whose mutation string would be empty; every one of these
    # tables omits that row by construction, so including it reports a spurious 0%.
    scored = [c for c in coverages if c[0] > 0]
    worst = min(scored, key=lambda c: c[3]) if scored else (0, 0, 0, float("nan"))
    verdicts.append((key, wmean(shallow), wmean(deep), worst))
    print()

print("=" * 72)
print(f"{'landscape':<14}{'shallow cov':>13}{'deep cov':>11}{'worst depth':>13}{'  verdict'}")
print("-" * 72)
for key, s, d, worst in verdicts:
    spec = bench.LIBRARIES[key]
    total_expected = sum(expected_by_depth(spec))
    drop = (s - d) if (s == s and d == d) else 0.0
    if total_expected > 50 * spec.published_measured_count if spec.published_measured_count else False:
        verdict = "SPARSE SAMPLE, not a complete design - coverage is not meaningful"
    elif worst[1] == 0:
        verdict = f"ZERO at depth {worst[0]}  **the CR9114 failure mode**"
    elif abs(drop) < 0.05:
        verdict = "flat - not depth-correlated"
    elif drop > 0:
        verdict = f"deep half {drop:.1%} lower - DEPTH-CORRELATED"
    else:
        verdict = f"deep half {-drop:.1%} higher - inverse"
    print(f"{key:<14}{s:>12.1%}{d:>11.1%}{worst[3]:>12.1%}  {verdict}")
