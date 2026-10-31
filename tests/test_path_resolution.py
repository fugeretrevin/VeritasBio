"""
tests/test_path_resolution.py - the environment resolution contract in `veritas.paths`.

These tests fix the behaviour that three separate agent failures in one night came from,
each of which wasted a whole run:

  * a coverage agent concluded "there is no GB1 selection trace anywhere on disk - I
    scanned all 153 run folders" and fell back to the wrong seed block. It was reading a
    git WORKTREE, and `benchmark_comparison_results/` is gitignored, so the directory
    does not exist inside a worktree at all;
  * a dataset agent could not find `scripts/analysis/README.md`, the reference
    implementation it had been told to reuse, because its worktree predated the commit
    that added it;
  * a benchmark agent reported ALDE "not on disk anywhere, no interpreter, no torch" and
    SILENTLY DROPPED the strongest baseline from a 48-seed run. ALDE is not in this repo:
    it is a sibling checkout with its own interpreter, and the main `.venv` (python 3.14,
    no torch) says nothing about it. The root cause was a search rooted at the worktree,
    whose parent is `<repo>/.claude/worktrees` - two levels short of the directory that
    holds sibling checkouts.

So: resolution is worktree-aware and resolves to the REAL repository root; every failure
NAMES THE PATHS IT CHECKED; a path an operator supplied is a constraint rather than a
hint; and nothing is ever silently relocated.

Resolving ALDE itself is NOT tested here, because it is not this module's job:
`benchmark_panoptes_isolated.resolve_alde_paths` owns it and has its own tests. What is
tested here is the root that such a search has to start from.
"""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from veritas.paths import PathResolutionError

THIS_REPO = Path(__file__).resolve().parent.parent


def _git(*args: str, cwd: Path) -> str:
    completed = subprocess.run(
        ["git", *args], cwd=str(cwd), capture_output=True, text=True, timeout=120,
    )
    if completed.returncode != 0:
        raise AssertionError(f"git {' '.join(args)} failed: {completed.stderr.strip()}")
    return completed.stdout.strip()


class _SyntheticRepo:
    """A real git repository with a real linked worktree, built in a temp directory.

    A synthetic repository rather than this one: the test needs to assert what happens
    when the worktree's checkout is OLDER than the main one, and it must not depend on
    the agent worktrees that happen to exist under `.claude/worktrees` today.
    """

    def __init__(self, root: Path) -> None:
        self.main = root / "main"
        self.main.mkdir(parents=True)
        _git("init", "-b", "trunk", cwd=self.main)
        _git("config", "user.email", "test@example.invalid", cwd=self.main)
        _git("config", "user.name", "Test", cwd=self.main)
        (self.main / "AGENTS.md").write_text("marker\n", encoding="utf-8")
        (self.main / ".gitignore").write_text("benchmark_comparison_results/\n", encoding="utf-8")
        _git("add", ".", cwd=self.main)
        _git("commit", "-m", "first", cwd=self.main)

        # The worktree is created from the FIRST commit, then the main checkout moves on.
        # That is incident 2: the file the agent was told to reuse exists in the repo but
        # not in the tree it happened to be standing in.
        self.worktree = root / "wt"
        _git("worktree", "add", "--detach", str(self.worktree), "HEAD", cwd=self.main)

        analysis = self.main / "scripts" / "analysis"
        analysis.mkdir(parents=True)
        (analysis / "README.md").write_text("reference implementation\n", encoding="utf-8")
        _git("add", ".", cwd=self.main)
        _git("commit", "-m", "second", cwd=self.main)

        # Gitignored, so it exists only in the main checkout - never in the worktree.
        (self.main / "benchmark_comparison_results" / "run_0001").mkdir(parents=True)


class TestMainCheckoutResolution(unittest.TestCase):
    """Resolution from the main checkout, which is the easy case and must stay easy."""

    def test_repository_root_is_the_checkout_that_holds_the_markers(self) -> None:
        from veritas.paths import repository_root

        root = repository_root()
        self.assertEqual(root, THIS_REPO)
        self.assertTrue((root / "AGENTS.md").is_file())
        self.assertTrue((root / "benchmark_panoptes_isolated.py").is_file())

    def test_a_nested_start_resolves_to_the_same_root(self) -> None:
        from veritas.paths import repository_root

        self.assertEqual(repository_root(start=THIS_REPO / "tests"), THIS_REPO)
        self.assertEqual(repository_root(start=THIS_REPO / "veritas" / "datasets"), THIS_REPO)

    def test_the_main_checkout_is_not_reported_as_a_worktree(self) -> None:
        from veritas.paths import checkout_root, in_worktree

        self.assertFalse(in_worktree(start=THIS_REPO))
        self.assertEqual(checkout_root(start=THIS_REPO), THIS_REPO)

    def test_a_tracked_file_resolves_and_reports_its_source(self) -> None:
        from veritas.paths import explain_repo_path, repo_path

        self.assertTrue(repo_path("AGENTS.md", start=THIS_REPO).is_file())
        resolution = explain_repo_path("scripts/data", start=THIS_REPO)
        self.assertEqual(resolution.source, "checkout")
        self.assertEqual(resolution.path, THIS_REPO / "scripts" / "data")

    def test_the_dataset_directory_and_a_known_dataset_resolve(self) -> None:
        from veritas.paths import dataset_dir, dataset_file

        self.assertTrue(dataset_dir(start=THIS_REPO).is_dir())
        self.assertTrue(dataset_file("GB1_STRSG_Wu_2016.csv", start=THIS_REPO).is_file())


class TestWorktreeResolution(unittest.TestCase):
    """The worktree case, which is where all three incidents actually happened."""

    def setUp(self) -> None:
        self._tmp = TemporaryDirectory(prefix="veritas_paths_")
        self.addCleanup(self._cleanup)
        self.repo = _SyntheticRepo(Path(self._tmp.name).resolve())

    def _cleanup(self) -> None:
        # Git keeps files read-only on Windows often enough to break rmtree.
        for path in Path(self._tmp.name).rglob("*"):
            try:
                os.chmod(path, 0o700)
            except OSError:
                pass
        try:
            self._tmp.cleanup()
        except OSError:
            pass

    def test_a_worktree_resolves_to_the_real_repository_root(self) -> None:
        from veritas.paths import checkout_root, in_worktree, repository_root

        self.assertTrue(in_worktree(start=self.repo.worktree))
        self.assertEqual(checkout_root(start=self.repo.worktree), self.repo.worktree)
        self.assertEqual(repository_root(start=self.repo.worktree), self.repo.main)

    def test_the_same_answer_without_shelling_out_to_git(self) -> None:
        # The `.git` FILE in a worktree points at the main repository, so the answer must
        # not depend on a git binary being on PATH.
        from veritas.paths import repository_root

        self.assertEqual(
            repository_root(start=self.repo.worktree, allow_git=False), self.repo.main
        )

    def test_a_gitignored_results_directory_resolves_into_the_real_root(self) -> None:
        # Incident 1: the run folders are not missing, they are in the main checkout.
        from veritas.paths import explain_repo_path, results_dir

        resolved = results_dir(
            "benchmark_comparison_results", start=self.repo.worktree
        )
        self.assertEqual(resolved, self.repo.main / "benchmark_comparison_results")
        self.assertTrue((resolved / "run_0001").is_dir())

        resolution = explain_repo_path(
            "benchmark_comparison_results", start=self.repo.worktree, search="repository"
        )
        self.assertEqual(resolution.source, "repository")
        self.assertNotIn(str(self.repo.worktree), str(resolution.path))

    def test_a_file_missing_from_the_worktree_is_found_in_the_real_root(self) -> None:
        # Incident 2: the worktree predates the commit that added the reference file.
        from veritas.paths import explain_repo_path, repo_path

        relative = "scripts/analysis/README.md"
        self.assertFalse((self.repo.worktree / relative).exists())

        resolution = explain_repo_path(relative, start=self.repo.worktree)
        self.assertEqual(resolution.source, "repository")
        self.assertEqual(resolution.path, self.repo.main / "scripts" / "analysis" / "README.md")
        # Both candidates are recorded, so the fallback is never invisible.
        self.assertEqual(len(resolution.candidates), 2)
        self.assertIn(self.repo.worktree / relative, resolution.candidates)
        self.assertEqual(
            repo_path(relative, start=self.repo.worktree).read_text(encoding="utf-8").strip(),
            "reference implementation",
        )

    def test_search_local_refuses_to_look_outside_the_worktree(self) -> None:
        from veritas.paths import repo_path

        with self.assertRaises(PathResolutionError) as caught:
            repo_path("scripts/analysis/README.md", start=self.repo.worktree, search="checkout")
        message = str(caught.exception)
        self.assertIn(str(self.repo.worktree), message)
        self.assertNotIn(str(self.repo.main / "scripts"), message)


class TestMissingTargetsNameWhatWasChecked(unittest.TestCase):
    """A bare 'unavailable' is what let a baseline get dropped. Every failure lists paths."""

    def test_a_missing_repo_path_lists_every_candidate(self) -> None:
        from veritas.paths import repo_path

        with self.assertRaises(PathResolutionError) as caught:
            repo_path("scripts/data/DOES_NOT_EXIST.csv", start=THIS_REPO)
        message = str(caught.exception)
        self.assertIn("DOES_NOT_EXIST.csv", message)
        self.assertIn(str(THIS_REPO / "scripts" / "data" / "DOES_NOT_EXIST.csv"), message)
        self.assertIn("checked", message.lower())

    def test_a_missing_dataset_names_the_dataset_directory_and_what_is_there(self) -> None:
        from veritas.paths import dataset_file

        with self.assertRaises(PathResolutionError) as caught:
            dataset_file("NO_SUCH_LIBRARY.csv", start=THIS_REPO)
        message = str(caught.exception)
        self.assertIn("NO_SUCH_LIBRARY.csv", message)
        self.assertIn("scripts", message)
        # It names what IS present, so the next step is reading a name rather than a scan.
        self.assertIn("GB1_STRSG_Wu_2016.csv", message)

    def test_a_missing_results_directory_says_it_is_gitignored(self) -> None:
        from veritas.paths import results_dir

        with TemporaryDirectory(prefix="veritas_paths_bare_") as raw:
            bare = Path(raw).resolve()
            (bare / "AGENTS.md").write_text("marker\n", encoding="utf-8")
            with self.assertRaises(PathResolutionError) as caught:
                results_dir("benchmark_comparison_results", start=bare, allow_git=False)
        message = str(caught.exception)
        self.assertIn("benchmark_comparison_results", message)
        self.assertIn("gitignored", message.lower())
        self.assertIn(str(bare / "benchmark_comparison_results"), message)

    def test_a_write_target_may_be_absent_but_is_never_relocated(self) -> None:
        from veritas.paths import results_dir

        with TemporaryDirectory(prefix="veritas_paths_bare_") as raw:
            bare = Path(raw).resolve()
            (bare / "AGENTS.md").write_text("marker\n", encoding="utf-8")
            target = results_dir(
                "benchmark_e2e_results", start=bare, must_exist=False, allow_git=False
            )
        self.assertEqual(target, bare / "benchmark_e2e_results")


class TestSiblingCheckoutsRoot(unittest.TestCase):
    """
    The directory that holds checkouts BESIDE the repository, which is where an external
    tool search has to start.

    This is incident 3's root cause reduced to one assertion: a worktree's own parent is
    `<repo>/.claude/worktrees`, two levels short of the answer, and a search rooted there
    reports a tool that is present as absent.
    """

    def test_from_the_main_checkout_it_is_the_parent_of_the_repository(self) -> None:
        from veritas.paths import sibling_checkouts_dir

        self.assertEqual(sibling_checkouts_dir(start=THIS_REPO), THIS_REPO.parent)

    def test_from_a_worktree_it_is_still_the_parent_of_the_real_repository(self) -> None:
        from veritas.paths import sibling_checkouts_dir

        worktrees = THIS_REPO / ".claude" / "worktrees"
        existing = (
            [p for p in worktrees.glob("*") if (p / ".git").exists()]
            if worktrees.is_dir() else []
        )
        if not existing:
            self.skipTest("no agent worktree present to resolve from")
        worktree = existing[0]
        resolved = sibling_checkouts_dir(start=worktree)
        self.assertEqual(resolved, THIS_REPO.parent)
        # The two wrong answers the overnight search actually produced.
        self.assertNotEqual(resolved, worktree.parent)
        self.assertNotEqual(resolved, worktree.parent.parent)
        self.assertNotIn("worktrees", str(resolved))


class TestAnExplicitPathIsAConstraint(unittest.TestCase):
    """
    A path an operator supplied resolves or raises. It is never replaced by a discovered
    one: running a different location than the one named, and reporting the result as
    though it came from the named one, is the same class of failure as the silent
    UNAVAILABLE row that dropped a baseline from a published comparison.
    """

    def test_an_existing_explicit_path_is_returned_unchanged(self) -> None:
        from veritas.paths import explicit_path

        resolved = explicit_path(THIS_REPO / "AGENTS.md", what="--manifest", must_be_dir=False)
        self.assertEqual(resolved, THIS_REPO / "AGENTS.md")

    def test_a_missing_explicit_path_raises_and_says_nothing_else_was_tried(self) -> None:
        from veritas.paths import explicit_path

        target = THIS_REPO / "scripts" / "data" / "NOT_A_DATASET.csv"
        with self.assertRaises(PathResolutionError) as caught:
            explicit_path(target, what="--dataset", must_be_dir=False)
        message = str(caught.exception)
        self.assertIn("--dataset", message)
        self.assertIn(str(target), message)
        self.assertIn("Nothing else was tried", message)

    def test_a_directory_supplied_where_a_file_was_required_is_refused(self) -> None:
        from veritas.paths import explicit_path

        with self.assertRaises(PathResolutionError) as caught:
            explicit_path(THIS_REPO / "scripts", what="--imputed-ids", must_be_dir=False)
        self.assertIn("not a file", str(caught.exception))

    def test_an_empty_value_is_refused_by_name(self) -> None:
        from veritas.paths import explicit_path

        with self.assertRaises(PathResolutionError) as caught:
            explicit_path("", what="ALDE_PATH")
        self.assertIn("ALDE_PATH", str(caught.exception))

    def test_a_write_target_need_not_exist_but_is_not_relocated(self) -> None:
        from veritas.paths import explicit_path

        target = THIS_REPO / "benchmark_e2e_results" / "not_yet"
        self.assertEqual(explicit_path(target, what="--output-dir", must_exist=False), target)


class TestEnvironmentRecord(unittest.TestCase):
    """The manifest record, which is what lets a later reader see where a run stood."""

    def test_it_is_json_serialisable_and_names_both_roots(self) -> None:
        import json

        from veritas.paths import describe_environment

        record = describe_environment(start=THIS_REPO)
        encoded = json.loads(json.dumps(record))
        self.assertEqual(encoded["repository_root"], str(THIS_REPO))
        self.assertFalse(encoded["in_worktree"])
        self.assertEqual(encoded["sibling_checkouts_dir"], str(THIS_REPO.parent))
        for name in ("benchmark_comparison_results", "benchmark_e2e_results"):
            self.assertTrue(encoded["results_dirs"][name]["gitignored"])
            self.assertTrue(encoded["results_dirs"][name]["path"].endswith(name))

    def test_the_record_points_at_the_real_root_from_a_worktree(self) -> None:
        from veritas.paths import describe_environment

        worktrees = THIS_REPO / ".claude" / "worktrees"
        existing = (
            [p for p in worktrees.glob("*") if (p / ".git").exists()]
            if worktrees.is_dir() else []
        )
        if not existing:
            self.skipTest("no agent worktree present to resolve from")
        record = describe_environment(start=existing[0])
        self.assertTrue(record["in_worktree"])
        self.assertEqual(record["repository_root"], str(THIS_REPO))
        # The gitignored run folders live in the real root, never in the worktree.
        self.assertNotIn(
            "worktrees", record["results_dirs"]["benchmark_comparison_results"]["path"]
        )

