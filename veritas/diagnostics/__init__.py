"""
veritas.diagnostics - run diagnostics that attribute a failure to a specific stage.

`coverage` answers "did candidate GENERATION ever put the answer on the ballot?",
which is the question a terminal-variance or data-invariance guard cannot answer.
`trace` records the per-candidate scoring detail that makes a selection round
auditable, bounded from the start. `retrospective` is the only place a true
experimental outcome is allowed to meet a coverage metric, and it can only be
reached after the live accumulator has been sealed.
"""

from veritas.diagnostics.coverage import (  # noqa: F401
    CollapseThresholds,
    CoverageAccumulator,
    CoverageDiagnosis,
    CoverageReport,
    GenerationCollapseWarning,
    RoundCoverage,
    assert_generation_coverage,
    classify_coverage,
)
from veritas.diagnostics.depth_returns import (  # noqa: F401
    DepthReturnLedger,
    DepthReturnRow,
    DepthReturnTable,
    MeasuredSelection,
    OutcomeDisclosureError,
    write_depth_returns_csv,
)
from veritas.diagnostics.genotype import (  # noqa: F401
    GenotypeParseError,
    LibrarySpec,
    infer_library_from_variant_ids,
)
from veritas.diagnostics.retrospective import (  # noqa: F401
    DepthShellProfile,
    RetrospectiveAccessError,
    RetrospectiveCoverageReport,
    RetrospectiveFitnessAccess,
    analyse_reachable_set,
    depth_shell_profile,
)
from veritas.diagnostics.trace import (  # noqa: F401
    LEGACY_SELECTION_TRACE_COLUMNS,
    BoundedCandidateTrace,
    CandidateTraceRecord,
    TraceBudget,
    write_candidate_trace_csv,
)

__all__ = [
    "LEGACY_SELECTION_TRACE_COLUMNS",
    "BoundedCandidateTrace",
    "CandidateTraceRecord",
    "CollapseThresholds",
    "CoverageAccumulator",
    "CoverageDiagnosis",
    "CoverageReport",
    "DepthReturnLedger",
    "DepthReturnRow",
    "DepthReturnTable",
    "DepthShellProfile",
    "GenerationCollapseWarning",
    "GenotypeParseError",
    "LibrarySpec",
    "MeasuredSelection",
    "OutcomeDisclosureError",
    "RetrospectiveAccessError",
    "RetrospectiveCoverageReport",
    "RetrospectiveFitnessAccess",
    "RoundCoverage",
    "TraceBudget",
    "analyse_reachable_set",
    "assert_generation_coverage",
    "classify_coverage",
    "depth_shell_profile",
    "infer_library_from_variant_ids",
    "write_candidate_trace_csv",
    "write_depth_returns_csv",
]
