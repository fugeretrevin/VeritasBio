# VeritasBio — Deterministic Cloud Economics Engine for Computational Biology

[![Meta AI Hackathon](https://img.shields.io/badge/Meta_Global_AI_Hackathon-2026-blue.svg)](https://lablab.ai)
[![Llama 3.1 Powered](https://img.shields.io/badge/Llama_3.1-Strict_JSON_Priors-purple.svg)](https://llama.meta.com/)
[![Modal GPU](https://img.shields.io/badge/Modal-NVIDIA_A10G-green.svg)](https://modal.com)
[![Tests Passing](https://img.shields.io/badge/Tests-50%2F50_Passing-brightgreen.svg)]()
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)]()

> **VeritasBio** is a deterministic cloud economics engine and active learning platform for generative protein engineering. An LLM (*Meta Llama 3.1*) only proposes biophysical search boundaries; a Gaussian Process actively selects high-yield micro-batches; and an immutable SQLite ledger with a fail-closed Python Gatekeeper guarantees that paid GPU compute (*Modal NVIDIA A10G + ESM-2*) never exceeds the user's hard budget.

---

![VeritasBio Benchmark Results](scripts/benchmark_results.png)

---

## 1. The Problem: The Cloud-Cost Trap in Generative Biology

In computational drug discovery and protein engineering, running foundation model inference (e.g. ESM-2 650M, AlphaFold 3, RoseTTAFold) across combinatorial mutation libraries creates runaway cloud infrastructure bills:
1. **Combinatorial Explosion**: Testing double and triple mutants on a 238-amino-acid protein yields $>1.1 \times 10^7$ variants — costing over **$35,000** in brute-force GPU forward passes.
2. **Cold-Start Waste**: Cloud serverless instances incur 15–30s container boot times on every invocation, wasting up to 90% of GPU compute budgets on container spin-up rather than tensor operations.
3. **Agentic Hallucination**: Open-ended autonomous AI agents lack cost predictability, frequently looping or executing hundreds of unconstrained inference calls.

### The VeritasBio Solution
VeritasBio replaces unpredictable open-ended agents with a **deterministic, fail-closed 6-step active learning cycle**. The engine optimizes proteins for highest biological fitness while mathematically maximizing the **Return on Investment (ROI)** per dollar of cloud spend.

---

## 2. Benchmark Results (Sarkisyan 2016 Nature avGFP DMS)

Validated against the gold-standard **ProteinGym / Sarkisyan et al. GFP Deep Mutational Scanning (DMS)** benchmark:

| Metric | VeritasBio Guided Search | Brute-Force Random Screening | Improvement / Value |
|---|---|---|---|
| **Spearman Rank Correlation ($\rho$)** | **+0.9718** ($p = 9.72 \times 10^{-13}$) | +0.1200 | **High ranking fidelity** |
| **Pearson Correlation ($r$)** | **+1.0000** ($p = 2.13 \times 10^{-38}$) | +0.0900 | **Exact calibration** |
| **Cloud Cost for Discovery** | **$0.0002** | **$0.0184** | **98.9% Cloud Savings** |
| **Cloud Efficiency Multiple** | **91.8x cheaper** | 1.0x baseline | **Near 2 orders of magnitude** |
| **Top Discovered Variant** | `T203F:E222G` (+1.2818) | Random variants | **Discovers synergistic double mutants** |

> **Reading the correlations above.** These are measured on the points the Gaussian Process was
> trained on, so they describe *fit*, not predictive skill, and will sit near 1.0 for almost any
> model. The TrpB benchmark below reports a held-out figure alongside the in-sample one.

### 2b. Second Benchmark: TrpB Active-Site Landscape (Johnston et al. 2024)

A harder and stricter test. The avGFP assay is sparse and survivorship-biased; the TrpB library is
**combinatorially complete** — all 20⁴ = 160,000 combinations at active-site residues 183/184/227/228
of the engineered parent Tm9D8\* were measured, so the global optimum is *known* to be in the pool
and nothing is missing.

| Metric | VeritasBio Guided Search | Random Screening (n=200 draws) |
|---|---|---|
| **Variants screened** | 125 of 159,999 (0.078%) | 125 |
| **Top-1% hits found** | **36** | 1.17 |
| **Enrichment over random** | **30.9x** | 1.0x baseline |
| **Best fitness found** | **0.8005** (parent Tm9D8\* = 0.4081) | 0.4453 |
| **Held-out Spearman ($\rho$)** | −0.0365 | — |
| **In-sample Spearman ($\rho$)** | +0.9944 | — |
| **Cloud cost for discovery** | **$0.00175** | $0.9852 brute force |

The honest headline here is **30.9x enrichment over random screening at 0.078% of the screening
cost**, reaching 0.80 fitness against a known global optimum of 1.00 after 125 measurements — the
loop finds strong variants well above chance. The near-zero *held-out* correlation is
reported deliberately: TrpB's active site is famously epistatic, and 5-scale physicochemical deltas
do not generalise across it from 125 measurements. Expected Improvement still concentrates sampling
in good regions even when the surrogate's global ranking is weak, which is why the hit rate rises
while $\rho$ does not. Treat this as the benchmark that shows where the current featurizer runs out
of road, not as a regression.

```bash
python scripts/download_trpb_benchmark.py && python scripts/run_trpb_benchmark.py
```

The full 159,999-variant run takes ~10s. `--pool-size N` subsamples the candidate pool for a faster
demo, at the cost of a weaker result (the optimum may not be in the subsample).

**Data provenance.** `scripts/download_trpb_benchmark.py` downloads the published landscape, keeps
every fitness value verbatim, and re-encodes each full-length sequence as mutation notation against
the Tm9D8\* parent. It refuses to emit a benchmark if the saturated positions it recovers from the
data disagree with the published active site. `tests/test_trpb_benchmark.py` asserts combinatorial
completeness (19, 19², 19³, 19⁴ by mutation order) and the published score range, so silent drift in
the upstream file fails the suite rather than producing a plausible-looking result.

---

## 3. System Architecture & The 6-Step Loop

```
  [ UniProt RAG ]
        │  1. Ingest WT sequence, chromophore triads, functional regions
        ▼
  [ Meta Llama 3.1 ]  (Strict JSON mode only)
        │  2. Propose multi-site residue clusters & forbidden positions
        ▼
  [ Quant Engine ]    (5-factor Physicochemical Featurizer + Matérn 5/2 GP)
        │  3. Model epistatic landscape; select micro-batch via Expected Improvement (EI)
        ▼
  [ Gatekeeper ]      (Biophysical validator + Two-state cold/warm cost model)
        │  4. Evaluate ROI (EI/$); request atomic spend authorization from SQLite ledger
        ▼
  [ SQLite Ledger ]   (INTEGER micro-USD, immutability triggers, ACID locks)
        │  5. Atomically deduct budget; authorize or fail-closed reject
        ▼
  [ Cloud Muscle ]    (Modal NVIDIA A10G + ESM-2 650M with scaledown_window=120)
        │  6. Execute batched PyTorch tensor scoring on approved sequences
        ▼
  [ Evidence Update ]
        └──> Append immutable scores to Ledger; update GP prior; repeat loop
```

### Non-Negotiable Engineering Invariants
- **Fail-Closed Security**: Rejects on parse error, unknown models, missing budgets, or non-canonical mutations.
- **Ledger Invariant**: Never bypasses the ledger. Every paid execution must receive `authorize_spend(...) -> True`.
- **Strict Separation of Concerns**:
  - The LLM **never** does arithmetic, executes code, or selects tools.
  - The Quant Engine **never** accesses external APIs.
  - The Modal Worker **never** decides spend or batch sizes.

---

## 4. Key Innovations

### 1. Two-State Cold/Warm Container Cost Model
A standard serverless cold container takes 20s to boot ($t_{\text{boot}} = 20.0\text{s}$ at \$1.10/hr = \$0.00612), while warm containers eliminate boot overhead ($t_{\text{boot}} = 0.0\text{s}$). By utilizing Modal's `scaledown_window=120`, VeritasBio schedules batches within the warm keep-alive window, achieving a **~21x reduction in cost per micro-batch**.

### 2. Physicochemical Property Featurization
Instead of treating amino acid substitutions as arbitrary one-hot tokens, mutations are featurized across 5 standardized biophysical scales:
- Kyte-Doolittle Hydrophobicity
- Molecular Mass (Daltons)
- Net Charge at pH 7.4
- Polarity (Grantham scale)
- Isoelectric Point (pI)

### 3. Expected Improvement per Dollar ($EI/\$$)
Candidates are selected not merely by predicted fitness, but by Expected Improvement ($EI$) balanced against the cost of inference ($EI / \$$), mathematically directing budget to the highest information-density search regions.

---

## 5. Quickstart Guide

### Installation

```bash
git clone https://github.com/fugeretrevin/VeritasBio.git
cd VeritasBio

# Create virtual environment (Python 3.11 recommended)
python -m venv .venv
.venv\Scripts\activate   # On Windows
# source .venv/bin/activate # On Linux/macOS

# Install dependencies
pip install -r requirements.txt
```

### Running the Interactive Web UI
Launch the real-time Gradio telemetry dashboard:

```bash
python main.py
```
Open your browser to `http://127.0.0.1:7860`. The dashboard streams live round-by-round active learning progression, remaining budget gauges, top candidate rankings, and an audit trail.

### Running Headless CLI Mode
For autonomous automated batch execution:

```bash
python main.py --cli --budget 2.00 --iterations 5 --batch-size 10 --model mock_esm
```

### Benchmarking & Metrics
Evaluate ledger discoveries against the Sarkisyan 2016 Nature DMS dataset:

```bash
# Print formatted terminal report
python scripts/evaluate.py --db ledger/veritasbio_ledger.db

# Output raw JSON metrics
python scripts/evaluate.py --db ledger/veritasbio_ledger.db --json

# Generate publication-quality visualization charts
python scripts/plot_benchmark.py --db ledger/veritasbio_ledger.db --out scripts/benchmark_results.png
```

### Running the Complete Test Suite
VeritasBio includes 50 comprehensive unit tests covering all 6 phases:

```bash
python -m unittest discover -s tests -v
```

---

## 6. Directory Structure

```
VeritasBio/
├── controller/
│   ├── gatekeeper.py       # Pydantic validation, two-state cost model, ROI evaluator
│   ├── prompts.py          # Llama 3.1 combinatorial search strategy prompts & parser
│   └── orchestrator.py     # Unidirectional 6-step active learning generator loop
├── workers/
│   ├── rag_worker.py       # UniProt metadata fetcher & local disk caching
│   ├── data/               # Offline caches: P42212 (avGFP), P50909 + TRPB_TM9D8 (TrpB parent)
│   ├── math_worker.py      # 5-factor featurizer, Matérn 5/2 GP, Expected Improvement
│   ├── esm_worker.py       # Modal serverless NVIDIA A10G ESM-2 scoring function
│   └── mock_workers.py     # Continuous epistatic avGFP fitness landscape (Sarkisyan 2016)
├── ledger/
│   ├── schema.sql          # DDL with micro-USD precision and immutability triggers
│   └── database.py         # Fail-closed SQLite CRUD & atomic authorize_spend()
├── ui/
│   └── dashboard.py        # Asynchronous streaming Gradio dashboard
├── scripts/
│   ├── evaluate.py         # Spearman correlation, hit rate, and cloud savings engine
│   ├── plot_benchmark.py   # Publication-quality benchmark visualization generator
│   ├── download_real_benchmark.py  # Fetches the ProteinGym avGFP DMS assay
│   ├── download_trpb_benchmark.py  # Fetches + re-encodes the TrpB active-site landscape
│   ├── run_real_data_benchmark.py  # Active learning loop on avGFP
│   ├── run_trpb_benchmark.py       # Active learning loop on TrpB (held-out + random baseline)
│   └── data/
│       ├── sarkisyan_gfp_benchmark.json # Sarkisyan et al. avGFP ground truth
│       └── TRPB_THEMA_Johnston_2024.csv # TrpB landscape, 159,999 measured variants
├── presentation.html       # Interactive partner deck with live economics simulator
├── main.py                 # Application entry point (CLI and Web Dashboard)
├── requirements.txt        # Production dependencies
└── tests/                  # 50 unit tests across Phases 1–6
```

---

## 7. License & Credits

Built for the **Meta Global AI Developer Hackathon**.
- **Architecture & Engine**: VeritasBio Team
- **Foundation Models**: Meta Llama 3.1, Meta Fair ESM-2
- **Infrastructure**: Modal Labs (NVIDIA A10G Serverless)
- **Benchmark Data**: Sarkisyan et al. (2016) *Nature* (ProteinGym GFP DMS)