# VeritasBio — Agent Instructions & Skills

This repository uses automated AI agent engineering skills based on [Matt Pocock's Skills](https://github.com/mattpocock/skills) alongside domain-specific VeritasBio active learning and biophysics skills.

## Agent Skills Configuration

### Issue Tracker
Issues are tracked in GitHub Issues using the `gh` CLI. See [`docs/agents/issue-tracker.md`](docs/agents/issue-tracker.md).

### Triage Labels
Standard five-role triage vocabulary (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See [`docs/agents/triage-labels.md`](docs/agents/triage-labels.md).

### Domain Docs
Single-context domain architecture documentation layout. See [`docs/agents/domain.md`](docs/agents/domain.md) and [`ARCHITECTURE.md`](ARCHITECTURE.md).

---

## Available Skills Reference

### Engineering Workflows (Matt Pocock Skills)
- **`ask-matt`**: Ask which skill or flow fits your current engineering situation.
- **`wayfinder`**: Map big architectural decisions and plan the path through a problem.
- **`grill-with-docs`**: Deep interview to sharpen requirements and document findings in ADRs.
- **`to-spec`**: Write a comprehensive, agreed technical specification from an ideation thread.
- **`to-tickets`**: Slice a spec into tracer-bullet dependency-linked tickets.
- **`implement`**: Build a ticket or feature using test-driven development.
- **`implement-spec`**: Parallelized execution of a complete spec across subagents on an integration branch.
- **`tdd`**: Test-driven development (red-green-refactor loop) for clean, verifiable code.
- **`diagnosing-bugs`**: Methodical root-cause investigation with tight reproducing feedback loops.
- **`code-review`**: Standards and specification diff review.
- **`pr`**: Formats structured pull request descriptions with before/after verification evidence.
- **`retro`**: Post-build retrospective suggesting improvements to agent environment and checks.
- **`triage`**: Triages incoming issues and bug reports into actionable agent work.

### Productivity & Collaboration
- **`grill-me`**: Interactive interview to pressure-test design decisions.
- **`handoff`**: Prepares clean context bridges across sessions or subagents.
- **`teach`**: Educational explanations and walkthroughs.
- **`wait-what`**: Sanity checks assumptions when a plan feels confusing or counterintuitive.

### VeritasBio Domain Skills
- **`dev-orchestrator`**: Coordinates 6-step active learning and multi-track tasks.
- **`dev-quant-engine`**: Gaussian Process, DKL, GFlowNet, and Bayesian optimization.
- **`dev-cloud-worker`**: Modal A10G ESM-2 inference and mock protein evolution workers.
- **`dev-compliance-cro`**: Opentrons OT-2, Twist cDNA synthesis, Benchling LIMS, and FDA 21 CFR Part 11 Merkle proofs.
- **`dev-ui-viewer`**: 3Dmol WebGL molecular viewer and Plotly fitness landscape visualizations.
- **`impeccable`**: Frontend UX design, polishing, styling, and accessibility.

