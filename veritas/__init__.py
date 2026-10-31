"""
veritas - library-level packages extracted from the benchmark harnesses.

Nothing in here imports a benchmark script, an orchestrator or a worker at module
scope: these modules are meant to be importable from an analysis notebook, a test or
a production run without dragging a campaign's dependencies in with them.
"""
