---
name: pbcheck-map
description: Map of the pbcheck repository — what it is, where every module/script/doc lives, what is frozen and why, current verified numbers, and where to start. Load before walking the tree blind on any pbcheck task — continuing the Phase 0 measurement, touching gate_config/AMENDMENTS.md/PHASE0_SPEC.md, adding an amendment, running the synthetic gate, editing census_select/io_counts/naive_engine/moderated, or fixing a bug anywhere in src/pbcheck, scripts/, synthetic/, tests/ or docs/.
model: sonnet
effort: low
tools: Read, Glob, Grep
color: cyan
---

Verified against `c3c2556` on `main` (clean tree): `pytest -q` → **367 passed**, 0 failed, in ~20s
(18 files under `tests/`). `git log --oneline -20` shows the last five commits are Amendment 5 Part A,
two engine-affordability fixes, the stratum-list machine-checking pass, and Amendment 4 Part A —
i.e. today's work is the tip of the amendment chain, not a side branch of it.

## What this project is

pbcheck audits **pseudoreplication** in single-cell RNA-seq differential expression: naive per-cell
tests (the scanpy default) treat individual cells as independent replicates, but cells from one donor
are correlated, so the true unit of replication is the **donor**, and the naive approach massively
inflates false discoveries. **Phase 0** is a pre-registered *measurement* study — not a released tool —
that quantifies the size of that inflation on public CELLxGENE Census data by comparing naive per-cell
DE against donor-pseudobulk DE and a donor-permutation null. Nothing here is a scientific finding yet:
every number currently in the repo comes from synthetic oracles with known ground truth, and no real
Census stratum has been analysed.

## Where things live

### `src/pbcheck/` — the measurement engine (imported as `pbcheck`)

| module | owns |
|---|---|
| `design.py` | metadata-only auditor: donors/group, donor nesting, per-donor Cramér's V confound, imbalance — needs no counts |
| `gene_universe.py` | the single label-agnostic gene universe (`frozen_universe`) shared by both arms and every permutation, so the tested set can't shift under relabeling |
| `methods/de.py` | `DEResult`, the common return shape both arms produce |
| `methods/naive.py` | the naive per-cell arm — faithful scanpy `rank_genes_groups` Wilcoxon, the thing being audited |
| `methods/naive_engine.py` | the naive arm's donor-permutation null computed from per-donor sufficient statistics instead of re-running scanpy per permutation — see Traps |
| `methods/moderated.py` | the pseudobulk arm's **current** test: moderated eBayes (limma-style, Smyth 2004), selected by Amendment 2 |
| `methods/pseudobulk.py` | aggregation to donor pseudobulk (thin-donor filter, `MIN_CELLS=10`/`MIN_COUNTS=1000`) plus the **retired** DESeq2-Wald path (`deseq_from_pdata`), kept only so Amendment-1-era numbers stay reproducible |
| `mtc.py` | shared BH correction — `bh_both_arms` (paired, spec §5) and `bh_over_universe` (solo) |
| `permutation.py` | the donor-permutation null itself: `build_perms`, `labels_for`, `run_null` |
| `metrics.py` | genomic-inflation λ, the permutation floor, `NDegSeries` (carries its BH convention, read-only) |
| `gate_config.py` | every threshold in one place, tagged by provenance — see Frozen material below |
| `census_select.py` | `.obs`-only real-data stratum selection (spec §1): pinned Census version, `(dataset_id × cell_type)` strata, inclusion gate, confound pre-screen, candidate manifest. **Admits nothing** |
| `io_counts.py` | the counts gate (spec §9 item 2): integer-count assertion on values not dtype, frozen-universe sizing, fills the columns `census_select` left `pending`. **Admits nothing** |

`design.py`, `gene_universe.py`, `methods/`, `metrics.py`, `mtc.py`, `permutation.py` are the built
measurement engine. `census_select.py` + `io_counts.py` are the real-data harness's first two modules
(of the plan in spec §9); `controls.py`, `decision.py`, `report.py` are specified but **not built** —
see `pilot/README.md`'s "What remains".

### `scripts/` — drivers, one job each

| script | job |
|---|---|
| `synthetic_gate.py` | the calibration gate: runs the whole engine on synthetic oracles, prints/writes a GO-style readout. `python scripts/synthetic_gate.py` with no args is the protocol run |
| `pb_calibration_probe.py` | **frozen as evidence** — produced the Amendment 1/2 grid; deliberately not refactored, not imported by the shipped arm (see Traps) |
| `run_test_selection_grid.py` / `analyze_test_selection.py` | the 146-cell test-selection grid and its verdict-printer; the grid is committed at `pilot/testsel/`, frozen evidence for Amendment 2 |
| `census_candidates.py` | the two-pass streaming driver that runs `census_select` over the **whole** pinned Census (dispatched by `.github/workflows/census-candidates.yml`, manual only) |
| `freeze_stratum_list.py` | emits the frozen §1 stratum list from the committed candidate manifest; recomputes every figure it declares, `--check` verifies byte-identity |
| `fetch_preregistration_evidence.py` | re-fetches (does NOT need to be re-run — its output is pinned) the two external indexes the freeze reasons about |
| `check_upper_bound_claim.py` | reproducer for Amendment 4 Part A's Correction 1 probe table (named seeds, committed artifact) |
| `check_version_consistency.py` | guards `CITATION.cff` against drifting from `src/pbcheck/__init__.py`'s version |
| `proof_of_life.py` | **historical**, predates the frozen protocol; importing it warns — not comparable to anything current |

### `synthetic/oracles.py`

The correctness spec, not shipped runtime code (imported by tests via `sys.path` in
`tests/conftest.py`, deliberately outside `src/`). Builds synthetic `AnnData` with a donor random
effect; the NULL oracle (no true DE) and POSITIVE oracle (log2FC=1.0 injected into K=200 genes) are
what every gate number is computed on. `donor_sigma=0` is the falsification control.

### `tests/` — 18 files, 367 tests, all offline

No test needs the network or `cellxgene-census`; `pytest -m "not slow"` skips DESeq2-touching /
multi-permutation end-to-end tests. `tests/conftest.py` holds the shared oracle fixtures.
`tests/test_stratum_list_freeze.py`, `test_census_select.py`, `test_io_counts.py` and
`test_census_candidates.py` are the largest files — they exercise the real-data harness's admission
logic against synthetic `obs` frames, never a real Census read.

### `pilot/` — committed Phase 0 artifacts, no code

| dir | holds |
|---|---|
| `gate/` | `synthetic_gate_2026-08-15.json` — the committed calibration-gate run cited by README/CHANGELOG/Amendment 3 |
| `testsel/` | `summary.{csv,json}` — the frozen 146-cell test-selection grid (commit `72dec7b`) that Amendment 2's test choice rests on |
| `preregistration/` | the frozen §1 stratum list + the whole-Census candidate manifest it derives from + two pinned external indexes + the redacted proposal doc — see Frozen material |
| `upper_bound_check/` | Amendment 4 Part A Correction 1's reproducer artifact (named seeds); `replicates.jsonl` is a resume ledger, gitignored |
| `results/` | **empty by design** — gitignored except `.gitkeep`; see Traps |

### `docs/`

| doc | for |
|---|---|
| `PHASE0_SPEC.md` | the frozen pre-registration: thesis, decision rule, corrections A1–D5, the full methodology (§1–§10). Changes only via `AMENDMENTS.md` |
| `AMENDMENTS.md` | append-only, dated, numbered amendment log — currently **5 amendments** (Amendment 4 and 5 each have only a Part A so far; Part B is unrun for both) |
| `PREREGISTRATION_STRATUM_LIST.md` | the machine-checked §1 pre-registration act itself: fixes 12 datasets / 251 strata, not an amendment |
| `ENV_NOTES.md` | verified stack (Windows, Python 3.12.10) and upstream API drift (decoupler 2.x, PyDESeq2 0.5.4, scanpy 1.12.2) |
| `PILOT_FINDINGS.md` | **historical**, dated 2026-07-19, predates the frozen spec's final form. One of its claims is retracted by Amendment 1 (stated there, not edited here); read it as a dated snapshot, not current status |

`README.md` (repo root) carries the live status table and is kept in sync with the amendment log on
every amendment; `pilot/README.md` carries the fuller "what's built / what's gapped" narrative and the
decision rule's short form; `CHANGELOG.md` narrates the project's epochs; `CONTRIBUTING.md` states the
amendment-vs-engineering-change rule and dev setup.

## The frozen material — touching it without an amendment is a protocol violation

- **`docs/PHASE0_SPEC.md`** — the pre-registration. Frozen; a session may read it freely but must
  never edit it directly. Any change to a threshold, oracle, decision rule or dataset-selection
  criterion goes through `docs/AMENDMENTS.md` first.
- **`gate_config.PRE_REGISTERED`** — every frozen numeric threshold, in one place, each tagged with
  its spec section. `gate_config.INSTRUMENT_SANITY` is a separate block, explicitly **not**
  pre-registered (the gate script's own sanity checks) — don't confuse the two when reading it.
  `scripts/pb_calibration_probe.py` deliberately carries its **own copy** of these constants rather
  than importing `gate_config` — it is the independent reference instrument the Amendment 2 grid was
  measured with, and a reference that imports the thing it checks stops being a reference.
  `tests/test_gate_config.py` pins the two copies equal instead.
- **`pilot/gate/`, `pilot/testsel/`, `pilot/preregistration/`** — cited evidence.
  `.pre-commit-config.yaml` excludes all three from the fixing hooks (`end-of-file-fixer`,
  `trailing-whitespace`) on purpose: running the hooks once already silently rewrote the gate JSON
  before this exclusion existed (see CHANGELOG, 2026-08-16). `.gitattributes` additionally marks
  `pilot/preregistration/** -text` because that directory is pinned **by sha256** in
  `docs/PREREGISTRATION_STRATUM_LIST.md` and `scripts/freeze_stratum_list.py` refuses to run if the
  hash moves — line-ending normalization would silently break the pin. `pilot/gate/` and
  `pilot/testsel/` are deliberately **not** marked `-text`: they were committed with LF and are
  checked out with CRLF, and adding `-text` would make the CRLF working-tree bytes canonical and
  rewrite all three blobs on the next `git add` — protecting evidence by rewriting it is not
  protecting it (see `.gitattributes`'s own comment for the measured byte deltas).
- A session **may** edit: everything under `src/pbcheck/`, `scripts/` (other than treating
  `pb_calibration_probe.py` as live code — see Traps), `tests/`, and `synthetic/oracles.py`, and may
  add prose to `docs/AMENDMENTS.md` (append only) or `docs/PREREGISTRATION_STRATUM_LIST.md`'s own
  amendment-gated extensions. A session **must not** hand-edit `docs/PHASE0_SPEC.md`, silently rewrite
  anything under `pilot/`, or change a `gate_config.PRE_REGISTERED` value without a preceding
  `AMENDMENTS.md` entry.

## The amendment discipline

A change to a frozen threshold or a claim requires a **dated, numbered entry in `docs/AMENDMENTS.md`,
written and committed before the code that applies it** — including a "data visible at the time"
disclosure. This is not procedural decoration: **Amendment 1 exists because `synthetic_gate.py` once
quietly substituted an easier oracle (log2FC 1.5/K 150 for the pre-registered 1.0/200) and dropped two
binding criteria**, and every later amendment cites that episode as the reason the order matters.
`CONTRIBUTING.md`: "If your change makes `scripts/synthetic_gate.py` print different numbers with no
amendment explaining why, that is a protocol violation, not a refactor."

The log's own register includes **relaxations, named as such**. Amendment 3 *narrowed* the power
criterion's claimed domain (an operating envelope, not a lowered bar). **Amendment 5 Part A explicitly
relaxed a prohibition and said so in those words** — quoting it exactly: *"This entry narrows that
prohibition. Narrowing a prohibition is a relaxation, and it is called one here rather than described
as a clarification."* It permits `lambda_naive` and the naive arm's permutation FP floor to be reported
outside Amendment 3's operating envelope, at `min(n_A, n_B) ≥ 8` under measured per-stratum
calibration — while everything that touches the pseudobulk arm's power claim stays envelope-gated,
untouched. A future entry that widens what may be reported or eases a threshold must do the same: name
itself as a relaxation, not a clarification.

Amendment 4 and Amendment 5 both currently have **only a Part A** — criteria and mechanism declared,
committed *before* the run that will be judged by them. Neither has a Part B yet; `src/pbcheck/sigma_donor.py`
does not exist. Do not write Part B's numbers into a Part A commit, and do not let `gate_config` gain a
`SIGMA_GATE` block (Amendment 4's designated slot) before Part A's own criteria (V1–V11) are satisfied
by a real run.

## Current state (verified this session)

- `pytest -q` on `c3c2556`: **367 passed**, 0 failed (~20s). Coverage is report-only (no `--cov-fail-under`
  gate — deliberate, see `tests.yml`).
- The committed gate run (`pilot/gate/synthetic_gate_2026-08-15.json`) verdict: **`INSTRUMENT VALID
  WITHIN THE STATED OPERATING ENVELOPE`**. λ_pseudobulk 1.01 (band [0.9, 1.1]), pseudobulk perm-null FP
  rate 0.035 (7/200), λ_naive 54.57, naive floor 1162/1500 genes (77.4%), power 0.86 at the envelope
  boundary (`sigma_donor`=0.35, 8v8) against the required ≥0.60. **Power 0.60 at `sigma_donor`=0.5 is
  still 0.35 — unmet, not claimed.**
- Operating envelope (Amendment 3): minimum donors/group **4 / 8 / 13 / 23** at `sigma_donor` **0.2 /
  0.35 / 0.5 / 0.7**.
- The §1 stratum list is frozen and published (`docs/PREREGISTRATION_STRATUM_LIST.md`,
  `pilot/preregistration/stratum_list_2026-08-16.json` — reloaded and counted directly this session):
  **12 independent datasets, 251 stratum-contrasts** in the analysis set, plus **106** within-collection
  control strata over **5** sibling datasets (357 data rows total). **Every row's `admitted_to_sweep` is
  `False`** — verified directly from the JSON, all 251 + 106 rows.
- `sigma_donor` is **unanchored**: the estimator is specified (Amendment 4 Part A) but not validated,
  and Amendment 3's original claim that its mechanism (`sqrt(s0²)·ln 2`) is an *upper bound* on
  `donor_sigma` was itself retracted by Amendment 4 Part A Correction 1 — it is now an audit quantity
  of **unknown error sign**. **Nothing is admitted to the real sweep.**

## The traps a newcomer will hit

1. **The naive arm's permutation null must stay bit-exact against scanpy.** `methods/naive_engine.py`
   recovers the whole donor-permutation null from one ranking pass via per-donor sufficient statistics
   instead of re-running scanpy per permutation (a ~1000× speedup, load-bearing for the real sweep's
   cost). It is held to bit-for-bit agreement (`tests/test_naive_engine.py`, max ulp 0) against the slow
   `naive_engine="scanpy"` path. Touching the ranking/tie-correction/rank-sum arithmetic here is a
   correctness-of-a-pre-registered-statistic change, not an optimisation, and is held to that standard.
2. **The two BH conventions must never be pooled.** `BH_PAIRED` (`bh_both_arms`, both arms over one
   common tested set — the only cross-arm-comparable one) and `BH_SOLO` (`bh_over_universe`, the naive
   arm's own floor) answer different questions. `metrics.NDegSeries` is frozen and carries its
   convention; the old mixed `naive_ndeg` key was **removed, not deprecated**, specifically so a stale
   caller gets a `KeyError` instead of a silently-wrong number.
3. **`t-test_overestim_var` is not reproducible across machines, for any implementation** — scanpy's
   own `fast_array_utils.stats.mean_var` forms `E[x²]−E[x]²` in float32 storage precision, which loses
   all significant digits on log1p-normalised data (measured condition number ~1.9e14). This was found
   by CI, not locally. It doesn't touch any number the project currently reports (neither `run_null`
   nor the gate passes `method=`, so both always run Wilcoxon), but don't add an `rtol` to paper over a
   `t-test_overestim_var` mismatch — there is no correct tolerance to pick.
4. **The property suites (`tests/test_properties.py`, `tests/test_naive_engine.py`) are derandomised on
   purpose** (`derandomize=True`, Hypothesis). `.hypothesis/examples` is gitignored and Hypothesis
   reseeds from entropy per run by default, so a green local run means nothing about CI; derandomizing
   makes a CI failure reproduce locally by test name.
5. **`pilot/results/` is gitignored; `pilot/preregistration/` is not.** `pilot/results/` is where CI
   writes candidate manifests as *artifacts* — deliberately never committed, so a candidate list can't
   "pre-register itself by accident" just by surviving a CI run. `pilot/preregistration/` holds the
   evidence a **human** deliberately promoted into git with `docs/PREREGISTRATION_STRATUM_LIST.md`
   attached — a pre-registration whose evidence can expire with a CI artifact isn't auditable.
6. **The committed gate artifact's runtime figure is stale.** `pilot/gate/synthetic_gate_2026-08-15.json`
   was committed once (`44cb5ad`) and never rerun; its 511.8s runtime predates the naive-engine speedup
   landed the next day. Its *correctness* numbers still hold (all 187 leaf values reproduce exactly
   under the new engine, per `CHANGELOG.md`), but running `scripts/synthetic_gate.py` today takes
   seconds, not minutes — don't be alarmed if it finishes faster than the doc says, and don't assume a
   fast run is a broken one.
7. **`scripts/pb_calibration_probe.py` is frozen evidence, not a place to fix bugs going forward.** It
   generated the Amendment 1/2 grid and is deliberately not refactored or made to import `gate_config`;
   new pseudobulk-arm work belongs in `src/pbcheck/methods/moderated.py`.

## Where to start

**Continuing the science.** The single open thread every recent amendment closes on is `sigma_donor`
anchoring. The next concrete step is **Amendment 4 Part B**: build `src/pbcheck/sigma_donor.py` per the
estimator Part A specifies, run the confirmatory validation grid on the declared non-overlapping seed
range, and judge it against the pre-declared V1–V11 criteria — a criterion that fails is a result, not
something to loosen. Read Amendment 4 Part A in full first (`docs/AMENDMENTS.md`); it is long because
it pre-declares exactly what would make the estimator untrustworthy.

**Fixing a bug.** Start at `CONTRIBUTING.md` for the amendment-vs-engineering-change test, then
`pytest -q -m "not slow"` for the fast loop. If the fix doesn't move any number `scripts/synthetic_gate.py`
prints, it's a plain engineering change; if it does, it needs an `AMENDMENTS.md` entry written first
(see Amendment 1's and Amendment 2's own text for the pattern). Never touch a `gate_config.PRE_REGISTERED`
value, `docs/PHASE0_SPEC.md`, or anything under `pilot/` to make a test pass.
