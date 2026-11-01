"""`veritas check` must refuse, and must never pass by accident.

Why this exists. The command computes three ceilings - SCOPE, REACH, VIABILITY - and
composes them into one verdict with an exit code. Every defect this project has recorded
in that shape was a FALSE PASS: a metric that read normal while the thing it measured
had not happened. So these tests concentrate on the ways a verdict could be too
generous, and on the arithmetic the verdict rests on.

Three traps are pinned deliberately, each because the equivalent mistake has been made
in this repository already:

* A SKIPPED check must not count as a passed one. `--skip-viability` exists for speed
  and the veto removes between 37% and 100% of a hit set, so a plan cleared without it
  is not cleared. The composed verdict must be UNKNOWN and the exit code non-zero.
* The STRUCTURE audit must only ever downgrade. A bad input to the viability check can
  make the veto pass everything - `filter_spatially_viable_variants` returns its input
  unchanged when no position maps - so a structure problem that RAISED a verdict would
  be the exact failure the gate exists to catch.
* The depth arithmetic is checked against HANDOFF rule 2.2's registered anchors rather
  than against itself. Rule 2.2 exists because there are already five genotype
  canonicalisations in this repository and one of them silently inverted a published
  result. `check_reach` reads the landscape's own `depths` and `hit_mask`, so the
  anchors must reproduce exactly; if they ever do not, the reach number is wrong and
  every verdict built on it is wrong with it.

What these tests do NOT assert is a veto pass rate or a hit count for a landscape whose
structure or threshold could legitimately be re-derived. Two tests here have previously
pinned a defect in place by asserting a literal - `one_hot_width == n_sites * 20`, which
required every library to be saturated, and `disabled == ["dpp", "prior"]`, which
required the manifest to ignore module-level defaults. The anchors below are different:
they are registered constants that the project has committed to reproducing.
"""
from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from veritas import check as vcheck
from veritas.calibration.landscapes import REGISTRY, load_landscape

#: HANDOFF rule 2.2, shallow shell 1 <= depth <= 2: {landscape: hits}.
SHALLOW_SHELL_HITS = {"trpb": 168, "gb1": 72}


class TestStructureGrading(unittest.TestCase):
    """The audit of the viability check's own input.

    The three bands exist because three real libraries sit in three different states and
    grading them identically would be its own dishonesty: TrpB's 8EH0 maps 4 of 4
    positions and disagrees at one (184: P against the library's F) - the right protein,
    the right numbering, one substitution. HER2's 1N8Z resolves to chain A, the
    trastuzumab LIGHT chain, while that library's positions are CDR-H3 on the HEAVY
    chain - 9 of 10 disagree. CreiLOV's 8QI9 is a non-aligning homolog - 11 of 11.
    """

    def test_full_agreement_is_clean(self) -> None:
        self.assertEqual(vcheck.grade_structure(10, 0), vcheck.STRUCTURE_OK)

    def test_a_minority_disagreement_is_a_caveat_not_a_refusal(self) -> None:
        """TrpB's shape: one of four. Still the right molecule."""
        self.assertEqual(vcheck.grade_structure(4, 1), vcheck.STRUCTURE_CAVEAT)

    def test_a_majority_disagreement_is_untrusted(self) -> None:
        """HER2's shape (9 of 10) and CreiLOV's (11 of 11)."""
        self.assertEqual(vcheck.grade_structure(10, 9), vcheck.STRUCTURE_UNTRUSTED)
        self.assertEqual(vcheck.grade_structure(11, 11), vcheck.STRUCTURE_UNTRUSTED)

    def test_nothing_mapped_is_untrusted_and_not_clean(self) -> None:
        """The decisive case, and the counter-intuitive one.

        When no assayed position maps, `filter_spatially_viable_variants` returns its
        input unchanged and the veto passes 100% of everything. That is numerically
        indistinguishable from 'the veto found nothing to reject' and means the
        opposite. Grading it OK would turn a dead check into a passing one.
        """
        self.assertEqual(vcheck.grade_structure(0, 0), vcheck.STRUCTURE_UNTRUSTED)

    def test_the_bands_are_ordered(self) -> None:
        """Monotone in the mismatch rate, so a worse structure can never grade better."""
        grades = [vcheck.grade_structure(100, k) for k in (0, 20, 60, 100)]
        self.assertEqual(grades, [vcheck.STRUCTURE_OK, vcheck.STRUCTURE_CAVEAT,
                                  vcheck.STRUCTURE_UNTRUSTED, vcheck.STRUCTURE_UNTRUSTED])


class TestVerdictComposition(unittest.TestCase):
    def _report(self, *verdicts: str) -> vcheck.Report:
        plan = vcheck.Plan(positions=(1,), depth=1, source="test")
        checks = [vcheck.CheckOutcome(f"c{i}", v, "") for i, v in enumerate(verdicts)]
        return vcheck.Report(library="trpb", label="t", plan=plan, checks=checks)

    def test_all_pass_passes(self) -> None:
        report = self._report(vcheck.PASS, vcheck.PASS, vcheck.PASS)
        self.assertEqual(report.verdict, vcheck.PASS)
        self.assertEqual(report.exit_code, 0)

    def test_the_worst_check_decides(self) -> None:
        self.assertEqual(
            self._report(vcheck.PASS, vcheck.DEGRADED, vcheck.REFUSE).verdict,
            vcheck.REFUSE)

    def test_unknown_never_becomes_a_pass(self) -> None:
        """An unmeasured ceiling is not an open one.

        This is the test that matters most. A gate that exits 0 when it could not
        measure something is not a gate, and the project's own history is a sequence of
        metrics that read fine because the thing they measured never ran.
        """
        report = self._report(vcheck.PASS, vcheck.UNKNOWN)
        self.assertEqual(report.verdict, vcheck.UNKNOWN)
        self.assertNotEqual(report.exit_code, 0)

    def test_a_refusal_outranks_an_unknown(self) -> None:
        """A known-closed ceiling is a stronger statement than an unmeasured one."""
        self.assertEqual(self._report(vcheck.UNKNOWN, vcheck.REFUSE).verdict,
                         vcheck.REFUSE)

    def test_only_pass_exits_zero(self) -> None:
        for verdict in (vcheck.DEGRADED, vcheck.REFUSE, vcheck.UNKNOWN):
            with self.subTest(verdict=verdict):
                self.assertNotEqual(self._report(verdict).exit_code, 0)


class TestScopeCheck(unittest.TestCase):
    def test_a_position_the_library_does_not_vary_is_refused(self) -> None:
        plan = vcheck.Plan(positions=(183, 184, 999), depth=2, source="test")
        outcome = vcheck.check_scope(plan, REGISTRY["trpb"].positions)
        self.assertEqual(outcome.verdict, vcheck.DEGRADED)
        self.assertEqual(outcome.numbers["unmeasurable_positions"], [999])

    def test_nothing_measurable_is_a_refusal_not_a_degradation(self) -> None:
        plan = vcheck.Plan(positions=(900, 901), depth=2, source="test")
        self.assertEqual(vcheck.check_scope(plan, REGISTRY["trpb"].positions).verdict,
                         vcheck.REFUSE)

    def test_the_librarys_own_positions_pass(self) -> None:
        spec = REGISTRY["creilov"]
        plan = vcheck.Plan(positions=spec.positions, depth=3, source="test")
        self.assertEqual(vcheck.check_scope(plan, spec.positions).verdict, vcheck.PASS)


class TestReachCheck(unittest.TestCase):
    def test_the_shallow_shell_anchors_reproduce(self) -> None:
        """HANDOFF rule 2.2. If this fails, every reach number is wrong.

        `check_reach` must read the landscape's own `depths` and `hit_mask` rather than
        re-deriving substitution counts. These two anchors are the project's registered
        proof that a piece of genotype code agrees with the canonical one.
        """
        for key, want in SHALLOW_SHELL_HITS.items():
            with self.subTest(library=key):
                landscape = load_landscape(key)
                plan = vcheck.Plan(positions=landscape.spec.positions, depth=2,
                                   source="test")
                outcome = vcheck.check_reach(plan, landscape)
                self.assertEqual(outcome.numbers["n_reachable"], want)

    def test_reach_is_monotone_in_depth(self) -> None:
        """More reach cannot find fewer hits. Catches a sign or boundary error."""
        landscape = load_landscape("creilov")
        plan_positions = landscape.spec.positions
        counts = [
            vcheck.check_reach(
                vcheck.Plan(positions=plan_positions, depth=d, source="test"),
                landscape).numbers["n_reachable"]
            for d in range(1, 8)
        ]
        self.assertEqual(counts, sorted(counts))

    def test_full_depth_reaches_the_whole_hit_set(self) -> None:
        landscape = load_landscape("trpb")
        depth = int(np.asarray(landscape.depths).max())
        outcome = vcheck.check_reach(
            vcheck.Plan(positions=landscape.spec.positions, depth=depth, source="test"),
            landscape)
        self.assertEqual(outcome.numbers["n_reachable"], landscape.n_hits)
        self.assertEqual(outcome.verdict, vcheck.PASS)

    def test_a_shallow_plan_on_a_deep_landscape_is_refused(self) -> None:
        """CreiLOV's hits sit at depths 5-9; a depth-3 plan reaches 75 of 1,655.

        Asserted as a band rather than as the literal 75, because the hit set depends on
        the registered top-1% threshold and that is allowed to be re-derived. What must
        not change is the verdict.
        """
        landscape = load_landscape("creilov")
        outcome = vcheck.check_reach(
            vcheck.Plan(positions=landscape.spec.positions, depth=3, source="test"),
            landscape)
        self.assertEqual(outcome.verdict, vcheck.REFUSE)
        self.assertLess(outcome.numbers["pct_reachable"], 25.0)


class TestPlanParsing(unittest.TestCase):
    def test_a_raw_search_strategy_dump_is_accepted(self) -> None:
        """The planner's own serialisation, which is what a user has on disk.

        Requiring them to reshape it into this tool's format would invite exactly the
        transcription error the tool exists to prevent.
        """
        blob = {
            "hypothesis": "widen the pocket",
            "max_mutation_depth": 3,
            "target_clusters": [
                {"cluster_id": "pocket", "positions": [183, 184],
                 "rationale": "lines the substrate channel"},
                {"cluster_id": "distal", "positions": [227], "rationale": "allosteric"},
            ],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "plan.json"
            path.write_text(json.dumps(blob), encoding="utf-8")
            plan = vcheck.plan_from_json(path)
        self.assertEqual(plan.positions, (183, 184, 227))
        self.assertEqual(plan.depth, 3)
        self.assertIn("pocket", plan.rationales)

    def test_this_tools_own_json_output_round_trips(self) -> None:
        """So a report can be re-checked without re-asking the model."""
        payload = {"plan": {"positions": [39, 40], "max_mutation_depth": 2,
                            "source": "live:x", "hypothesis": "h"}}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "report.json"
            path.write_text(json.dumps(payload), encoding="utf-8")
            plan = vcheck.plan_from_json(path)
        self.assertEqual((plan.positions, plan.depth), ((39, 40), 2))

    def test_a_plan_without_a_depth_is_rejected_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "plan.json"
            path.write_text(json.dumps({"positions": [1, 2]}), encoding="utf-8")
            with self.assertRaises(ValueError):
                vcheck.plan_from_json(path)

    def test_a_plan_without_positions_is_rejected_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "plan.json"
            path.write_text(json.dumps({"max_mutation_depth": 2}), encoding="utf-8")
            with self.assertRaises(ValueError):
                vcheck.plan_from_json(path)


class TestSkippingIsNotPassing(unittest.TestCase):
    def test_skipping_the_veto_cannot_clear_a_plan(self) -> None:
        """`--skip-viability` is a speed switch, never a clearance.

        Run on TrpB at full depth, where SCOPE and REACH both pass, so the ONLY thing
        standing between this plan and a zero exit code is the skipped check. If that
        ever returns 0 the flag has become a way to launder a plan past the gate.
        """
        landscape = load_landscape("trpb")
        depth = int(np.asarray(landscape.depths).max())
        plan = vcheck.Plan(positions=landscape.spec.positions, depth=depth,
                           source="test")
        report = vcheck.run_check("trpb", plan, skip_viability=True)

        scope, reach = report.checks[0], report.checks[1]
        self.assertEqual((scope.verdict, reach.verdict), (vcheck.PASS, vcheck.PASS),
                         "this test only proves anything while scope and reach pass")
        self.assertEqual(report.verdict, vcheck.UNKNOWN)
        self.assertNotEqual(report.exit_code, 0)

    def test_an_unknown_library_is_an_error_not_a_verdict(self) -> None:
        plan = vcheck.Plan(positions=(1,), depth=1, source="test")
        with self.assertRaises(KeyError):
            vcheck.run_check("not_a_library", plan, skip_viability=True)


class TestCliSurface(unittest.TestCase):
    def test_depth_is_required_without_a_plan_source(self) -> None:
        self.assertEqual(vcheck.main(["--library", "trpb"]), 2)

    def test_every_registered_landscape_is_offered(self) -> None:
        """A library the benchmark measures but the gate cannot check is a gap."""
        action = {a.dest: a for a in vcheck.build_parser()._actions}["library"]
        self.assertEqual(set(action.choices), set(REGISTRY))

    def test_a_hand_specified_plan_renders_without_the_profiler(self) -> None:
        """The fast path must not require a PDB, a network call or a campaign."""
        self.assertEqual(
            vcheck.main(["--library", "gb1", "--depth", "1", "--skip-viability"]), 1)


if __name__ == "__main__":
    unittest.main()
