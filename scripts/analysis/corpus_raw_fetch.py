#!/usr/bin/env python3
"""Fetch the RAW CombinGym tables for every frozen landscape, at the pinned commit.

Why this exists. `scripts/freeze_cross_protein_corpus.py` pins the ten CLEAN workbooks
and exactly one raw table (`bnAbs_CR9114_raw.csv`, added after the clean CR9114 workbook
was found to drop a whole quadrant of its design space). The other nine clean tables have
never been diffed against their raw counterparts, so the class of defect found in CR9114
has never been ruled out anywhere else.

Every clean table DOES have a raw counterpart at the same pinned commit. The blob SHA-1s
below were read from GitHub's tree API for the pinned tree

    https://api.github.com/repos/sitonglab/CombinGym/git/trees/
        911f19c7ba82e112056f154d3e8b2073f70dc7d1?recursive=1

which is the same authority `freeze_cross_protein_corpus.py` cites for its own pins, and
`git hash-object` reproduces them. A download whose blob SHA-1 does not match is deleted,
not kept: an unverified table must never reach an audit.

This script does NOT touch the lock file or the freeze script. Extending the pin is the
coordinator's call; see docs/corpus-audit.md for the recommendation.

Usage
-----
    python scripts/analysis/corpus_raw_fetch.py fetch
    python scripts/analysis/corpus_raw_fetch.py verify
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_CORPUS_DIR = REPO_ROOT / "downloads" / "combingym"

COMMIT = "c5dccd514a9acfe1089a5f1555faea98afd78019"
TREE_SHA = "911f19c7ba82e112056f154d3e8b2073f70dc7d1"
RAW_BASE = f"https://raw.githubusercontent.com/sitonglab/CombinGym/{COMMIT}/"

# (repo path, git blob SHA-1 at the pinned tree, size in bytes).
RAW_FILES: Tuple[Tuple[str, str, int], ...] = (
    ("Data/DMS/Raw/CreiLOV_raw.xlsx",
     "300672c721f1fc1a5398f689c9c961dd6eaa9df0", 22_619_510),
    ("Data/DMS/Raw/GB1_raw.xlsx",
     "08b80e4367dcbe8ea8bfa70dff5d0908380c1581", 6_633_543),
    ("Data/DMS/Raw/HIV-1_raw.tsv",
     "4611f98454408f6edfdacb570be57e3b34147ec6", 68_024),
    ("Data/DMS/Raw/RhlA_raw.xlsx",
     "a8704f81527fbd089bac63847fb734c0fa2c87a8", 287_545),
    ("Data/DMS/Raw/SaCas9_raw.xlsx",
     "f7fafde24a570a11143452de5e808d933ea1cd35", 242_262),
    ("Data/DMS/Raw/SpCas9_raw.xlsx",
     "072d1c0064647fda885652c788dd82d899859776", 112_221),
    ("Data/DMS/Raw/bnAbs_CR6261_raw.csv",
     "1158a5dcfcc54a66467056ce096ecfbe13a06c36", 358_244),
    # Already pinned by the freeze script; listed so `verify` covers the whole raw set.
    ("Data/DMS/Raw/bnAbs_CR9114_raw.csv",
     "8ac7352f539a9b4d33330dc3764b0ff5e86a094d", 12_258_680),
    ("Data/DMS/Raw/eqFP611_raw.xlsx",
     "0028edb5b41202512bfd26ca37bb14be25dd4875", 737_089),
)


def git_blob_sha1(path: Path) -> str:
    """Reproduces `git hash-object <path>`: SHA-1 over b'blob <len>\\0' + content."""
    size = path.stat().st_size
    digest = hashlib.sha1()
    digest.update(b"blob %d\0" % size)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def cmd_fetch(args: argparse.Namespace) -> int:
    corpus_dir: Path = args.corpus_dir
    failures: List[str] = []
    for rel, blob, size in RAW_FILES:
        target = corpus_dir / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and git_blob_sha1(target) == blob:
            print(f"have  {rel}")
            continue
        url = RAW_BASE + urllib.request.quote(rel)
        print(f"get   {rel}  ({size:,} bytes)")
        try:
            with urllib.request.urlopen(url, timeout=args.timeout) as response:
                payload = response.read()
            target.write_bytes(payload)
        except Exception as exc:  # noqa: BLE001 - reported, never raised
            print(f"FAIL  {rel}: {type(exc).__name__}: {exc}")
            failures.append(rel)
            continue
        actual = git_blob_sha1(target)
        if actual != blob:
            # An unverified table must not reach the audit.
            target.unlink(missing_ok=True)
            print(f"FAIL  {rel}: blob sha1 {actual} != pinned {blob} (file deleted)")
            failures.append(rel)
    print(f"\n{len(RAW_FILES) - len(failures)}/{len(RAW_FILES)} raw files present "
          f"and matching the pin.")
    return 1 if failures else 0


def cmd_verify(args: argparse.Namespace) -> int:
    corpus_dir: Path = args.corpus_dir
    failures: List[str] = []
    for rel, blob, size in RAW_FILES:
        path = corpus_dir / rel
        if not path.exists():
            failures.append(f"{rel}: absent (checked {path})")
            continue
        actual_blob = git_blob_sha1(path)
        actual_size = path.stat().st_size
        if actual_blob != blob:
            failures.append(f"{rel}: blob sha1 {actual_blob} != {blob}")
            continue
        if actual_size != size:
            failures.append(f"{rel}: size {actual_size} != {size}")
            continue
        print(f"ok    {rel}  sha256={sha256_file(path)}")
    print()
    if failures:
        print(f"FAILED {len(failures)} of {len(RAW_FILES)} raw file(s):", file=sys.stderr)
        for item in failures:
            print(f"  - {item}", file=sys.stderr)
        return 1
    print(f"RAW CORPUS VERIFIED: {len(RAW_FILES)}/{len(RAW_FILES)} files match "
          f"CombinGym @ {COMMIT[:12]} (tree {TREE_SHA[:12]})")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--corpus-dir", type=Path, default=DEFAULT_CORPUS_DIR)
    subparsers = parser.add_subparsers(dest="command", required=True)
    fetch = subparsers.add_parser("fetch", help="Download the pinned raw tables.")
    fetch.add_argument("--timeout", type=int, default=600)
    fetch.set_defaults(func=cmd_fetch)
    verify = subparsers.add_parser("verify", help="Re-check the raw tables against the pin.")
    verify.set_defaults(func=cmd_verify)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
