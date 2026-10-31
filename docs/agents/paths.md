# Resolving Paths: the standing rule

**Before you report that something is not on disk, run this:**

```powershell
.venv\Scripts\python.exe -X utf8 -m veritas.paths
```

It prints the checkout you are standing in, the **real** repository root, whether those
two differ, the directory that holds sibling checkouts, and each gitignored results root
with whether it exists. One command, and the three failures below do not happen.

---

## The rules

1. **Resolve gitignored and external paths against the REAL repository root, never
   against the worktree.** A linked worktree has no gitignored content at all, and the
   directory beside it is `.claude/worktrees/`, where nothing external has ever lived.
   `veritas.paths.repository_root()` gives the real root from anywhere (it asks
   `git rev-parse --git-common-dir`, and falls back to parsing the worktree's `.git`
   file, so it works without a git binary). `checkout_root()` gives the tree you are
   standing in, when that is genuinely the question.

2. **Name every path you checked when resolution fails.** "Not found" sends the next
   reader on a filesystem scan; "not found, checked A, B, C" gets the configuration
   fixed. Every error from `veritas.paths` lists its candidates in order, and
   `resolve_alde_paths` in `benchmark_panoptes_isolated.py` prints
   `Looked in: --alde-path: not set; ALDE_PATH: <path> (no such directory); sibling of
   the repository (<root>): ...`. Write errors that way.

3. **Never silently substitute a different location.** Resolve, or report precisely what
   was missing. If a lookup has an ordered list of candidates, say which one answered -
   `explain_repo_path()` returns `source` for exactly this reason.

4. **An explicitly supplied path is a constraint, not a hint.** A `--dataset`,
   `--alde-path` or `ALDE_PATH` that does not exist is FATAL. It must never fall through
   to a discovered location, because that runs something other than what the operator
   named and then reports the result as theirs. Use
   `veritas.paths.explicit_path(value, what="--dataset")`.

5. **An absent gitignored directory is not a finding.** It means no run has written it
   in this repository yet. It is never evidence about what a run produced elsewhere.

6. **Another interpreter's contents tell you nothing about this one.** This repository's
   `.venv` is python 3.14 with no torch. Probing it says nothing about any external
   tool, and importing a package here is not a test of whether that tool can run.

---

## Worked example 1: "no GB1 selection trace anywhere on disk"

A coverage-diagnostics agent reported:

> there is no GB1 selection trace anywhere on disk - I scanned all 153 run folders

and then fell back to a different seed block, producing a diagnosis of the wrong data.

**What was actually true.** `benchmark_comparison_results/` is gitignored
(`.gitignore:246`), so it exists only where a run wrote it: the main checkout. The agent
was working inside a git worktree, where that directory cannot exist. The 153 folders it
scanned were whatever else it found; the GB1 traces were in the main checkout the whole
time.

**The rule.** Results roots resolve against the real repository root:

```python
from veritas import paths

runs = paths.results_dir("benchmark_comparison_results")      # real root, always
runs = paths.results_dir("benchmark_e2e_results", must_exist=False)  # a write target
```

If it is absent, the error says so *and* says it is gitignored, so "absent here" is never
mistaken for "this data does not exist". Falling back to another seed block is the part
that turned an environment problem into a published wrong answer: **a fallback that
changes which data you analysed is not a fallback, it is a different experiment.**

---

## Worked example 2: the reference implementation that "did not exist"

A dataset agent was told to reuse `scripts/analysis/README.md` and
`scripts/analysis/diversity.py`. It reported both missing and wrote its own.

**What was actually true.** Its worktree had been created *before* the commit that added
those files. A worktree is a checkout of one commit; it cannot see what a later commit
added. The files were in the repository, just not in that tree.

**The rule.** For repository content, check your tree and then the real root, and record
which answered:

```python
resolution = paths.explain_repo_path("scripts/analysis/diversity.py")
resolution.path     # the file
resolution.source   # "checkout" or "repository"
```

`source == "repository"` means **your worktree is behind** - rebase it rather than
treating the gap as a finding. When neither has it, the error lists both candidates and
says the worktree may predate the commit.

---

## Worked example 3: the baseline that was "not on disk anywhere"

A benchmark agent reported ALDE "not on disk anywhere, no interpreter, no torch" and
**silently dropped it from a 48-seed run**. ALDE is our strongest baseline. Two complete
comparisons were published from a report with one `UNAVAILABLE` row in it before anyone
read the row.

**What was actually true.** ALDE is not in this repository. It is a separate checkout:

| | |
|---|---|
| checkout | `C:\Users\stoof\OneDrive\Documents\GitHub\ALDE` @ `d0b3593dd179` |
| interpreter | `C:\Users\stoof\OneDrive\Documents\GitHub\.venv-alde\Scripts\python.exe` |
| its packages | torch 2.1.1+cpu, gpytorch 1.11, botorch 0.9.4, numpy 1.26.4 |
| this repo's `.venv` | python 3.14, **no torch** - which proves nothing about ALDE |

The search was rooted at the worktree. A worktree sits at
`<repo>/.claude/worktrees/<name>`, so its parent is the worktrees folder and its
grandparent is `.claude`: **two levels too deep.** The directory that holds sibling
checkouts is the parent of the *real* repository root, which is what
`paths.sibling_checkouts_dir()` returns from anywhere.

**The rule.** Resolving ALDE is already done - do not write it again.
`benchmark_panoptes_isolated.resolve_alde_paths()` searches `--alde-path`, then
`ALDE_PATH`, then `ALDE/` under the sibling directories of both the worktree and the real
checkout; it records every candidate and its outcome, marks a discovered path as
discovered rather than letting it pass for one the operator pinned, and a requested
baseline that cannot run is now fatal rather than a row in a report.

Three distinct findings, three different fixes, and they must never collapse into one
word:

| finding | fix |
|---|---|
| checkout missing | clone it, or set `ALDE_PATH` |
| interpreter missing | create its environment (`alde.yml` pins python 3.11) |
| dependency missing | install into *that* environment, not this one |

And the defect found in the first version of that fix, which is rule 4 above: an
explicit `--alde-path` that does not exist must stop the run. Falling through to a
discovered sibling would run a different checkout than the operator named and report it
as theirs - the same class of failure as the silent `UNAVAILABLE` row.

---

## The API, briefly

```python
from veritas import paths

paths.repository_root()          # the REAL root, from a worktree or the main checkout
paths.checkout_root()            # the tree you are standing in
paths.in_worktree()              # True inside a linked worktree
paths.sibling_checkouts_dir()    # where checkouts beside this repo live (ALDE's parent)

paths.results_dir("benchmark_comparison_results")   # gitignored -> real root only
paths.dataset_dir()                                 # scripts/data
paths.dataset_file("GB1_STRSG_Wu_2016.csv")         # error lists the datasets present

paths.repo_path("scripts/analysis/diversity.py")    # checkout first, then real root
paths.explain_repo_path(...)                        # same, plus source and candidates
paths.explicit_path(args.dataset, what="--dataset") # a constraint, never a hint

paths.describe_environment()     # JSON-safe record; stamp it into a run manifest
```

Everything raises `paths.PathResolutionError`, and every message lists what was checked.

**Out of scope on purpose.** `veritas.paths` does not resolve external tools. ALDE has
one implementation (`resolve_alde_paths`) and a second would drift from it. What this
module contributes to that problem is the root such a search must start from.

---

## If you are about to write "X does not exist"

- [ ] Did you run `python -m veritas.paths`?
- [ ] Are you in a worktree? (`in_worktree()`)
- [ ] Is X gitignored? Then it is only where a run wrote it - look in the real root.
- [ ] Was X added recently? Your worktree may predate the commit.
- [ ] Is X outside the repository (another checkout, another interpreter)?
- [ ] Does your report list every path you checked?
- [ ] Are you about to continue with *different* data instead? Stop. Say what is
      missing, and let the run fail on the real cause.
