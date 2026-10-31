"""Measurement calibration for the cross-protein study.

Three defects that would silently invalidate the study, each found on a real run and
each fixed here in label space only - nothing in this package reads or changes a
selector, a surrogate or a campaign:

* `qualification` - the frontier policy's `--qualify-threshold` is an absolute fitness
  on landscapes whose scales are incomparable, so it gates nothing outside TrpB.
  Expresses the bar as a quantile of the eligible universe instead.
* `endpoint` - the pre-registered primary endpoint (diversity-adjusted yield at
  d >= 2) is degenerate where the diversity adjustment removes nothing, which is the
  case on HER2. Precondition check plus a per-target verdict.
* `genotype_sets` - whether two readouts of one protein really share a genotype set,
  which is what a within-genotype-set paired contrast requires (protocol O-11).

`landscapes` holds the label-only loader the three share. `report` prints every table.
"""
from __future__ import annotations

from veritas.calibration.endpoint import (
    DEGENERATE,
    NEAR_DEGENERATE,
    USABLE,
    diversity_ceiling,
    endpoint_verdict,
    greedy_diverse,
)
from veritas.calibration.genotype_sets import (
    assess_shared_set,
    compare_genotype_sets,
    depth_profile,
)
from veritas.calibration.genotypes import (
    GenotypeParseError,
    canonical_combo,
    depth,
    hamming,
    min_pairwise_hamming,
    pair_distance_histogram,
    pairs_closer_than,
    parse_mutations,
)
from veritas.calibration.landscapes import REGISTRY, Landscape, LandscapeSpec, load_landscape
from veritas.calibration.qualification import (
    absolute_for_quantile,
    calibrate,
    equivalent_quantile,
    p_qualify_in_draw,
    qualifying_count,
    qualifying_fraction,
)

__all__ = [
    "DEGENERATE",
    "NEAR_DEGENERATE",
    "REGISTRY",
    "USABLE",
    "GenotypeParseError",
    "Landscape",
    "LandscapeSpec",
    "absolute_for_quantile",
    "assess_shared_set",
    "calibrate",
    "canonical_combo",
    "compare_genotype_sets",
    "depth",
    "depth_profile",
    "diversity_ceiling",
    "endpoint_verdict",
    "equivalent_quantile",
    "greedy_diverse",
    "hamming",
    "load_landscape",
    "min_pairwise_hamming",
    "p_qualify_in_draw",
    "pair_distance_histogram",
    "pairs_closer_than",
    "parse_mutations",
    "qualifying_count",
    "qualifying_fraction",
]
