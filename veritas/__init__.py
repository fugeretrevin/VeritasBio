"""
veritas - library-aware domain packages shared by the benchmark and the production loop.

Nothing in here imports a benchmark script, an orchestrator or a worker at module scope:
these modules are meant to be importable from an analysis notebook, a test or a
production run without dragging a campaign's dependencies in with them.

`veritas.datasets` holds the dataset concern: the canonical schema for a combinatorial
library and the ingestion path that turns a prepared table into validated records. It
deliberately imports nothing from `workers` or `controller`, so a new library can be
described, validated and ingested without touching optimisation code.

`veritas.diagnostics` holds the measurement concern: generation-coverage accounting, the
bounded candidate trace and the per-depth return ledger. Anything in there that reads
true fitness is retrospective only and is gated behind an explicit access token, so a
number from the landscape cannot flow back into a round that has yet to happen.

`veritas.paths` holds the environment concern: the real repository root (which is not the
worktree you may be standing in), the gitignored results roots, the dataset directory,
and the rule that a failure to resolve must name every path it checked. Three agent runs
were lost in one night to "it is not on disk" conclusions drawn from inside a worktree;
see `docs/agents/paths.md`.
"""

__all__ = ["datasets", "diagnostics", "paths"]
