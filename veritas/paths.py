"""
veritas.paths - where things actually are, and what was checked when they are not.

Three agent runs were wasted in one night by the same class of mistake, so this module
exists to make that mistake impossible to make quietly:

  1. A coverage agent reported "there is no GB1 selection trace anywhere on disk - I
     scanned all 153 run folders" and fell back to the wrong seed block. It was standing
     in a git WORKTREE. `benchmark_comparison_results/` is gitignored (.gitignore:246),
     so it exists only in the main checkout and a worktree contains no run folders at
     all. Nothing was missing; the agent was looking in a tree that can never hold it.

  2. A dataset agent could not find `scripts/analysis/README.md` and
     `scripts/analysis/diversity.py` - the reference implementation it had been told to
     reuse - because its worktree predated the commit that added them. The files were in
     the repository, just not in that checkout.

  3. A benchmark agent reported ALDE "not on disk anywhere, no interpreter, no torch" and
     SILENTLY DROPPED the strongest baseline from a 48-seed run. ALDE is not in this
     repository. It is a sibling checkout with its OWN interpreter, and the main `.venv`
     (python 3.14, no torch) proves nothing about it either way. The root cause
     generalises: a worktree sits at `<repo>/.claude/worktrees/<name>`, so ITS parent is
     the worktrees folder, not the directory that holds sibling checkouts. That search
     was two levels too deep, which is why `sibling_checkouts_dir()` below is defined
     against the real repository root and nothing else.

The rules this module enforces:

  * A worktree resolves to the REAL repository root for anything gitignored or external.
    `git rev-parse --git-common-dir` is how that root is found, with the `.git` FILE
    parsed directly as a fallback so the answer does not need a git binary.
  * Every failure NAMES THE PATHS IT CHECKED, in order. "not found; looked in A, B, C"
    is the entire point of this module. A bare "unavailable" is what let a requested
    baseline vanish from a published comparison.
  * Nothing is ever silently relocated. A resolution either returns a path and says
    which root it came from, or raises and lists the candidates. A write target may be
    absent (`must_exist=False`) but is still never moved somewhere else.
  * An explicitly supplied path is a CONSTRAINT, not a hint: `explicit_path()` resolves
    it or raises, and never falls through to a discovered location. Quietly using a
    different directory than the operator named, and then reporting the result as
    theirs, is the same failure as the silent UNAVAILABLE row.

Scope. This module owns the repository root, the gitignored results roots and the
dataset directory. It deliberately does NOT resolve external tools: ALDE's checkout and
interpreter are resolved by `benchmark_panoptes_isolated.resolve_alde_paths`, which
records every candidate it tried and treats an explicit-but-missing `--alde-path` as
fatal. One implementation, not two. What this module contributes to that problem is
`sibling_checkouts_dir()`, the root such a search must start from.

See `docs/agents/paths.md` for the same rules written for a human or an agent.
"""

from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

__all__ = [
    "PathResolutionError",
    "Resolution",
    "checkout_root",
    "dataset_dir",
    "dataset_file",
    "describe_environment",
    "explain_repo_path",
    "explicit_path",
    "in_worktree",
    "repo_path",
    "repository_root",
    "results_dir",
    "sibling_checkouts_dir",
]


# ==============================================================================================
# Errors
# ==============================================================================================
class PathResolutionError(RuntimeError):
    """
    A path could not be resolved, and the message lists every candidate that was checked.

    Never raised without those candidates. An agent that is told only "not found" goes
    on to scan, guess, or substitute something else; an agent that is handed the three
    paths that were tried either fixes the configuration or reports a real absence.
    """


# ==============================================================================================
# Roots: the checkout you are standing in, and the repository it belongs to
# ==============================================================================================
# Used only when there is no `.git` anywhere above the starting point - a tarball, a
# container copy, or a bare directory in a test. Presence of `.git` always wins.
REPO_MARKERS: Tuple[str, ...] = (
    "AGENTS.md",
    "benchmark_panoptes_isolated.py",
    "ARCHITECTURE.md",
)

# The gitignored output roots. Named here so a caller cannot invent a fourth spelling.
RESULTS_DIRS: Tuple[str, ...] = (
    "benchmark_comparison_results",
    "benchmark_e2e_results",
)

DATASET_RELATIVE = "scripts/data"

_GIT_TIMEOUT_SECONDS = 60


def _as_path(value: Any) -> Path:
    """
    A caller-supplied or environment-supplied path, expanded but NOT resolved.

    Not resolved deliberately: the string an operator typed into `ALDE_PATH` is what
    belongs in an error message, so they can see their own value echoed back rather than
    a normalised one they have to re-recognise.
    """
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else Path(os.path.abspath(path))


def _default_start() -> Path:
    """This module's own location, which is inside whatever checkout imported it."""
    return Path(__file__).resolve().parent.parent


def _candidate_chain(start: Optional[Path]) -> Tuple[Path, ...]:
    begin = (_as_path(start) if start is not None else _default_start())
    try:
        begin = begin.resolve()
    except OSError:  # pragma: no cover - unreadable mount
        pass
    if begin.is_file():
        begin = begin.parent
    return (begin, *begin.parents)


def _find_checkout_root(start: Optional[Path] = None) -> Tuple[Path, bool]:
    """
    The root of the tree you are standing in, and whether git controls it.

    Returns `(root, git_controlled)`. A worktree and a main checkout both come back as
    themselves here; turning a worktree into its real repository root is
    `repository_root()`'s job, kept separate so both answers are available.
    """
    chain = _candidate_chain(start)
    for directory in chain:
        if (directory / ".git").exists():
            return directory, True
    for directory in chain:
        if any((directory / marker).exists() for marker in REPO_MARKERS):
            return directory, False
    raise PathResolutionError(
        "no repository root above "
        f"{chain[0]}: checked each of these directories for a '.git' entry and then for "
        f"any of the markers {list(REPO_MARKERS)}:\n"
        + "\n".join(f"  {index}. {d}" for index, d in enumerate(chain, start=1))
        + "\nThis is not a VeritasBio checkout. Pass start= explicitly, or run from "
          "inside the repository."
    )


def checkout_root(start: Optional[Path] = None) -> Path:
    """The tree you are standing in: a worktree if you are in one, else the repository."""
    return _find_checkout_root(start)[0]


def in_worktree(start: Optional[Path] = None) -> bool:
    """True when the starting point is inside a linked git worktree."""
    root, git_controlled = _find_checkout_root(start)
    return bool(git_controlled and (root / ".git").is_file())


def _real_root_from_git(checkout: Path) -> Optional[Path]:
    """Asks git for the shared `.git` directory. None when git cannot answer."""
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "--git-common-dir"],
            cwd=str(checkout), capture_output=True, text=True,
            timeout=_GIT_TIMEOUT_SECONDS,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    raw = completed.stdout.strip()
    if not raw:
        return None
    common = Path(raw)
    if not common.is_absolute():
        common = (checkout / common)
    try:
        common = common.resolve()
    except OSError:  # pragma: no cover
        return None
    if common.name == ".git" and common.is_dir():
        return common.parent
    return None


def _real_root_from_git_file(checkout: Path) -> Optional[Path]:
    """
    Parses the worktree's `.git` FILE, so the answer survives git being unavailable.

    A worktree's `.git` is a one-line `gitdir: <main>/.git/worktrees/<name>` pointer, and
    that directory holds a `commondir` file pointing back at the shared `.git`.
    """
    git_entry = checkout / ".git"
    if not git_entry.is_file():
        return None
    try:
        text = git_entry.read_text(encoding="utf-8").strip()
    except OSError:  # pragma: no cover
        return None
    if not text.startswith("gitdir:"):
        return None
    gitdir = Path(text.split(":", 1)[1].strip()).expanduser()
    if not gitdir.is_absolute():
        gitdir = checkout / gitdir
    commondir_file = gitdir / "commondir"
    common: Optional[Path] = None
    if commondir_file.is_file():
        try:
            relative = commondir_file.read_text(encoding="utf-8").strip()
        except OSError:  # pragma: no cover
            relative = ""
        if relative:
            candidate = Path(relative)
            common = candidate if candidate.is_absolute() else gitdir / candidate
    if common is None:
        # `<main>/.git/worktrees/<name>` -> `<main>/.git`
        parents = list(gitdir.parents)
        if len(parents) >= 2 and parents[0].name == "worktrees":
            common = parents[1]
    if common is None:
        return None
    try:
        common = common.resolve()
    except OSError:  # pragma: no cover
        return None
    return common.parent if common.name == ".git" else None


def repository_root(start: Optional[Path] = None, allow_git: bool = True) -> Path:
    """
    The REAL repository root, which is not the worktree you may be standing in.

    Everything gitignored (the benchmark results roots) and everything external (ALDE)
    is resolved against this, because a worktree has no gitignored content and its
    sibling directory is `.claude/worktrees/`, where no external tool has ever lived.

    `allow_git=False` forces the `.git`-file path, which is also the fallback when git is
    not installed; both routes must agree.
    """
    checkout, git_controlled = _find_checkout_root(start)
    if not git_controlled:
        return checkout
    if (checkout / ".git").is_dir():
        return checkout  # a main checkout is its own repository root
    if allow_git:
        root = _real_root_from_git(checkout)
        if root is not None:
            return root
    root = _real_root_from_git_file(checkout)
    if root is not None:
        return root
    raise PathResolutionError(
        f"{checkout} has a '.git' file, so it is a linked worktree, but the real "
        "repository root could not be determined. Checked:\n"
        f"  1. git rev-parse --git-common-dir (cwd {checkout})"
        f"{' - SKIPPED, allow_git=False' if not allow_git else ''}\n"
        f"  2. the 'gitdir:' pointer in {checkout / '.git'} and the 'commondir' file "
        "beside it\n"
        "Do not fall back to the worktree: gitignored run folders and external tool "
        "checkouts are not there."
    )


# ==============================================================================================
# Repository-relative resolution
# ==============================================================================================
@dataclass(frozen=True)
class Resolution:
    """
    One resolved path, plus which root it came from and what else was tried.

    `source` is "checkout" when the path was found in the tree the caller is standing in
    and "repository" when it was found in the real repository root instead. That
    distinction is recorded rather than smoothed over: a file found only in the
    repository root means the caller's worktree is out of date, which is worth knowing
    before its absence is reported as a finding.
    """

    relative: str
    path: Path
    source: str
    candidates: Tuple[Path, ...]
    exists: bool

    def __fspath__(self) -> str:
        return str(self.path)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "relative": self.relative,
            "path": str(self.path),
            "source": self.source,
            "candidates": [str(c) for c in self.candidates],
            "exists": self.exists,
        }


def _dedupe(paths: Iterable[Path]) -> Tuple[Path, ...]:
    seen: Dict[str, Path] = {}
    for path in paths:
        key = str(path).lower() if os.name == "nt" else str(path)
        seen.setdefault(key, path)
    return tuple(seen.values())


def _checked_list(candidates: Sequence[Path]) -> str:
    return "\n".join(f"  {index}. {path}" for index, path in enumerate(candidates, start=1))


def explain_repo_path(
    relative: str,
    *,
    start: Optional[Path] = None,
    search: str = "auto",
    must_exist: bool = True,
    allow_git: bool = True,
    description: Optional[str] = None,
) -> Resolution:
    """
    Resolves a repository-relative path and reports exactly how it was resolved.

    `search` picks the candidate roots, in order:
      * "auto"       - the current checkout first, then the real repository root. The
                       second candidate is what rescues a worktree that predates the
                       commit adding the file (incident 2), and `source` says which one
                       answered, so the fallback is recorded and not invisible.
      * "repository" - the real repository root ONLY. Correct for anything gitignored:
                       a worktree cannot hold it, so looking there would only produce a
                       misleading absence (incident 1).
      * "checkout"   - the current tree ONLY. Correct when the question really is "what
                       is in MY tree", e.g. verifying the diff under review.

    With `must_exist=False` the first candidate is returned as a write target. It is
    never swapped for another location.
    """
    if search not in ("auto", "repository", "checkout"):
        raise ValueError(f"search must be auto, repository or checkout, not {search!r}")

    relative_clean = str(relative).replace("\\", "/").strip("/")
    local = checkout_root(start)
    real = repository_root(start, allow_git=allow_git)

    roots: Tuple[Tuple[str, Path], ...]
    if search == "checkout":
        roots = (("checkout", local),)
    elif search == "repository":
        roots = (("repository", real),)
    else:
        roots = (("checkout", local), ("repository", real))

    candidates = _dedupe(root / relative_clean for _, root in roots)
    for source, root in roots:
        candidate = root / relative_clean
        if candidate.exists():
            return Resolution(
                relative=relative_clean, path=candidate, source=source,
                candidates=candidates, exists=True,
            )

    if not must_exist:
        source, root = roots[0]
        return Resolution(
            relative=relative_clean, path=root / relative_clean, source=source,
            candidates=candidates, exists=False,
        )

    label = description or f"{relative_clean!r}"
    extra = ""
    if relative_clean.split("/")[0] in RESULTS_DIRS:
        extra = (
            f"\n{relative_clean.split('/')[0]}/ is GITIGNORED, so it exists only where a "
            "run actually wrote it - never in a git worktree, and never in a fresh "
            "clone. An empty result here means no run has written it in THIS repository, "
            "not that the data is missing."
        )
    elif search == "auto" and local != real:
        extra = (
            f"\nThis is a git worktree ({local}); the real repository root is {real}. "
            "Both were checked. If the file was added recently, this worktree may "
            "predate the commit - rebase it rather than concluding the file is absent."
        )
    raise PathResolutionError(
        f"{label} not found. Checked, in order:\n{_checked_list(candidates)}{extra}"
    )


def repo_path(
    relative: str,
    *,
    start: Optional[Path] = None,
    search: str = "auto",
    must_exist: bool = True,
    allow_git: bool = True,
    description: Optional[str] = None,
) -> Path:
    """`explain_repo_path`, when only the path is wanted. Raises with candidates listed."""
    return explain_repo_path(
        relative, start=start, search=search, must_exist=must_exist,
        allow_git=allow_git, description=description,
    ).path


def results_dir(
    name: str,
    *,
    start: Optional[Path] = None,
    must_exist: bool = True,
    allow_git: bool = True,
) -> Path:
    """
    A gitignored benchmark output root, resolved against the REAL repository root.

    This is incident 1 in one function: run folders are never in a worktree, so the
    worktree is not even a candidate.
    """
    if name not in RESULTS_DIRS:
        raise ValueError(f"unknown results directory {name!r}; expected one of {list(RESULTS_DIRS)}")
    return explain_repo_path(
        name, start=start, search="repository", must_exist=must_exist,
        allow_git=allow_git, description=f"the gitignored results root {name!r}",
    ).path


def dataset_dir(*, start: Optional[Path] = None, allow_git: bool = True) -> Path:
    """The prepared-dataset directory (`scripts/data`)."""
    return explain_repo_path(
        DATASET_RELATIVE, start=start, allow_git=allow_git,
        description=f"the dataset directory {DATASET_RELATIVE!r}",
    ).path


def dataset_file(name: str, *, start: Optional[Path] = None, allow_git: bool = True) -> Path:
    """
    One prepared dataset by file name, or an error naming the datasets that ARE present.

    Listing the directory's contents in the error turns "which library is this called?"
    into reading a name instead of a filesystem scan.
    """
    directory = dataset_dir(start=start, allow_git=allow_git)
    candidate = directory / name
    if candidate.is_file():
        return candidate
    present = sorted(p.name for p in directory.glob("*.csv"))
    raise PathResolutionError(
        f"dataset {name!r} not found. Checked:\n  1. {candidate}\n"
        f"{len(present)} prepared dataset(s) in {directory}:\n"
        + "\n".join(f"  - {p}" for p in present)
    )


# ==============================================================================================
# Beside the repository, and paths an operator named
# ==============================================================================================
def sibling_checkouts_dir(
    *,
    start: Optional[Path] = None,
    allow_git: bool = True,
) -> Path:
    """
    The directory that holds checkouts BESIDE this repository (where ALDE lives).

    Defined against the real repository root, never against the tree you are standing
    in, and this is the whole trap: a worktree sits at `<repo>/.claude/worktrees/<name>`,
    so `checkout_root().parent` is the worktrees folder and `.parent.parent` is `.claude`.
    A search built that way is two levels too deep and finds nothing, which is exactly
    how an overnight run concluded ALDE was "not on disk anywhere" while it sat at
    `<repo>/../ALDE` the whole time.

    Resolving a specific external tool is NOT this module's job - see
    `benchmark_panoptes_isolated.resolve_alde_paths`, which searches from here and
    records every candidate it tried.
    """
    return repository_root(start, allow_git=allow_git).parent


def explicit_path(
    value: Any,
    *,
    what: str,
    must_exist: bool = True,
    must_be_dir: Optional[bool] = None,
) -> Path:
    """
    A path an operator supplied: a CONSTRAINT, not a hint.

    It resolves to exactly that path or it raises. There is no discovery, no nearest
    match and no fallback, because falling through from a missing `--alde-path` to a
    discovered sibling clone would run a different checkout than the one that was named
    and then report the result as theirs. The same applies to a `--dataset` that does
    not exist, a `--imputed-ids` file that was renamed, and an `--output-dir` whose
    parent is gone: each of those has to stop the run where the mistake was made.

    `what` names the flag or variable, so the message points at the thing to fix.
    """
    if value is None or str(value).strip() == "":
        raise PathResolutionError(f"{what} was not supplied, and this call requires it.")
    path = _as_path(value)
    if not must_exist:
        return path
    exists = path.is_dir() if must_be_dir else (path.is_file() if must_be_dir is False else path.exists())
    if exists:
        return path
    expected = (
        "a directory" if must_be_dir
        else ("a file" if must_be_dir is False else "a file or directory")
    )
    raise PathResolutionError(
        f"{what} was set to {path}, which is not {expected}. Checked:\n  1. {path}\n"
        "Nothing else was tried: an explicitly supplied path is never replaced by a "
        "discovered one, because using a different location than the one named - and "
        "reporting the result as though it came from the named one - is worse than "
        "stopping here. Fix the value, or clear it so the caller's own default applies."
    )


# ==============================================================================================
# A record for manifests, and a CLI for whoever is about to conclude something is missing
# ==============================================================================================
def describe_environment(
    *,
    start: Optional[Path] = None,
    allow_git: bool = True,
) -> Dict[str, Any]:
    """
    Where this process resolved everything, as a JSON-safe record for a manifest.

    Worth stamping into a run: a later reader can see whether the run was executed from a
    worktree and which roots its inputs came from, instead of inferring it. `exists` on a
    results root is a fact about this machine, not a verdict - a gitignored directory
    that is absent means no run has written it here yet.
    """
    local = checkout_root(start)
    real = repository_root(start, allow_git=allow_git)
    record: Dict[str, Any] = {
        "checkout_root": str(local),
        "repository_root": str(real),
        "in_worktree": local != real,
        "sibling_checkouts_dir": str(real.parent),
        "results_dirs": {},
        "dataset_dir": None,
    }
    for name in RESULTS_DIRS:
        target = real / name
        record["results_dirs"][name] = {
            "path": str(target),
            "exists": target.is_dir(),
            "gitignored": True,
        }
    try:
        record["dataset_dir"] = str(dataset_dir(start=start, allow_git=allow_git))
    except PathResolutionError as exc:
        record["dataset_dir_error"] = str(exc)
    return record


def main() -> int:
    """
    `python -m veritas.paths` - prints what resolves and what does not.

    Run this BEFORE concluding that something is missing from disk. It is the cheapest
    way to find out that you are standing in a worktree.
    """
    record = describe_environment()
    print(json.dumps(record, indent=2))
    return 0 if record.get("dataset_dir") else 1


if __name__ == "__main__":
    raise SystemExit(main())
