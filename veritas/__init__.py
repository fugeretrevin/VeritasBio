"""
veritas - library-aware domain packages shared by the benchmark and the production loop.

`veritas.datasets` holds the dataset concern: the canonical schema for a combinatorial
library and the ingestion path that turns a prepared table into validated records. It
deliberately imports nothing from `workers` or `controller`, so a new library can be
described, validated and ingested without touching optimisation code.
"""

__all__ = ["datasets"]
