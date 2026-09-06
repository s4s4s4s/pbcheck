# pbcheck

**An auditor of pseudoreplication in single-cell RNA-seq differential expression.**

> Status: **v0.1.0 - a released single-stratum audit tool, built on the engine of Phase 0, a
> pre-registered measurement study** whose decision rule was frozen in
> [`docs/PHASE0_SPEC.md`](https://github.com/s4s4s4s/pbcheck/blob/main/docs/PHASE0_SPEC.md) *before*
> the data was seen; every departure from it is dated and justified in
> [`docs/AMENDMENTS.md`](https://github.com/s4s4s4s/pbcheck/blob/main/docs/AMENDMENTS.md). The release runs outside
> that protocol and claims no Phase 0 result.

## Install

```bash
pip install pbcheck
```

Python 3.12 or newer is required.

## Quickstart

```bash
# 1. write a small synthetic dataset with donor structure, no download needed
pbcheck example example.h5ad

# 2. audit it (one line; works as written in bash, PowerShell and cmd)
pbcheck audit example.h5ad --donor donor --condition condition --test disease --ref ctrl --n-perm 50 --n-perm-pb 20

# 3. on your own data, check the design first without touching expression values
pbcheck audit yourfile.h5ad --donor <col> --condition <col> --test <a> --ref <b> --design-only
```

Each run writes `pbcheck_audit.json`, `pbcheck_report.md` and `pbcheck_report.html` into an
output directory. See [`docs/USAGE.md`](https://github.com/s4s4s4s/pbcheck/blob/main/docs/USAGE.md)
for the full CLI reference, exit codes, the JSON schema and what every number means.

The same audit is available from Python:

```python
import pbcheck

adata = pbcheck.example_adata()
settings = pbcheck.AuditSettings(donor_col="donor", condition_col="condition",
                                  test_level="disease", ref_level="ctrl")
payload = pbcheck.run_audit(adata, settings)
```

## What the numbers mean

pbcheck runs two tests on your data: a naive per-cell test (the common but incorrect default,
which treats every cell as an independent replicate) and a donor-pseudobulk test (which
aggregates to one profile per donor before testing, the correct unit of replication). Both are
also run many times with donor labels shuffled between conditions, so that by construction there
is no real signal to find. The **permutation floor** is how many genes a test still calls
significant under that shuffling (`*_floor_solo` / the pseudobulk floor); the **ratio of the
real-label result to that floor** (the `*_real_over_floor_*` fields) says how far above pure
noise the real result sits; **genomic inflation (lambda)** is a separate summary of the shape of
the whole real-label p-value distribution, not derived from the floor. The pseudobulk arm's floor
and false-positive rate serve as a **negative control**: they check that the pseudobulk arm itself
behaves as it should on this file. None of this certifies a published result.

Every audit ends in one of three statuses: `complete` (both arms ran), `naive_only` (the
pseudobulk arm could not run, for example because the count matrix was not raw counts, or too
few donor profiles survived aggregation), and `design_only` (pbcheck stopped before touching
expression data at all, for example because a donor was measured under both conditions, or
because there were too few donors). Full detail on each status is in
[`docs/USAGE.md`](https://github.com/s4s4s4s/pbcheck/blob/main/docs/USAGE.md).

The pseudobulk arm's calibration itself has a stated limit, quoted here exactly as pbcheck's own
report states it:

> The pseudobulk arm's power was established on synthetic oracles only inside the operating
> envelope declared in Amendment 3:
>
> sigma_donor 0.2: at least 4 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'not in the grid; Amendment 1 frontier only')
> sigma_donor 0.35: at least 8 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'ebayes power 0.793 at 8v8 (calibrated) -> n* <= 8')
> sigma_donor 0.5: at least 13 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'ebayes power 0.486 at 12v12, the largest n tested -> n* > 12')
> sigma_donor 0.7: at least 23 donors per group (power at least 0.6 at log2FC 1.0 in 200 genes; grid support 'ebayes power 0.003 at 8v8 -> n* far above 8')
>
> The arm's calibration was evaluated at one hard regime (sigma_donor 0.5, 8 against 8 donors) and
> nowhere else. pbcheck does not estimate sigma_donor for real data, so whether this stratum lies
> inside that envelope is not determined here.

## What pbcheck does not do

- No single combined score of how trustworthy a result is (see "What exists, and what does not"
  below for the exact wording).
- No Phase 0 real-data harness (`controls`, `decision`, the report module the spec reserves) is
  built; those stay "specified but not built".
- No GO/NO-GO verdict, no envelope-membership verdict, no `sigma_donor` estimate for real data.
  The report states the operating envelope and says whether this file lies inside it is not
  determined.
- No paired or repeated-measures designs: if a donor was measured under both conditions, pbcheck
  stops at the design audit rather than running either differential-expression arm on it. There
  is no covariate for a repeated measurement, and combining donor and condition into one column
  to work around this is the exact error the tool exists to catch.
- No claim about any real dataset or any publication of any kind. The demonstrations below are
  outside the pre-registered protocol and are labelled as such.
- No access to the 17 frozen datasets or 357 pre-registered strata; no Census access in the
  product path.
- No multi-stratum batch mode, no plots, no R, no DESeq2 test arm.
- No rounding, rescaling or other "repair" of a count matrix that fails the raw-integer check;
  such a matrix is dropped from the pseudobulk arm, never altered.

## The problem

A large fraction of published single-cell RNA-seq differential-expression (DE) results are contaminated by
**pseudoreplication**: individual cells are treated as independent replicates, when the true unit of replication
is the **donor / sample**. Because cells from one donor are correlated, per-cell tests (e.g. a Wilcoxon test
across cells, as in the default `scanpy` workflow) drastically inflate the false-discovery rate - reporting
hundreds to thousands of "differentially expressed" genes that do not replicate. The correct analysis aggregates
counts to one profile per donor per cell type (**pseudobulk**) and tests across donors.

This is well documented in the methods literature (Squair et al. 2021; Zimmerman et al. 2021; Murphy & Skene
2023), yet the naive per-cell approach remains widespread in applied papers. What appears to be missing - I am
not aware of one, in Python or R - is a tool that takes an *existing* analysis and estimates *how much* its DE
results are inflated by this error. Pointers to prior art are welcome and will be credited here.

## Phase 0 status

### Where this actually stands

| | |
|---|---|
| Measurement engine | built; tested in CI |
| Real CELLxGENE data | **not touched.** The stratum list is pre-registered (frozen 2026-08-16; verify the date in [`pilot/README.md`](https://github.com/s4s4s4s/pbcheck/blob/main/pilot/README.md)); no frozen stratum has been run. The v0.1.0 demonstrations use public datasets outside the frozen list |
| Every number in the Phase 0 sections | from synthetic oracles with known ground truth - instrument calibration, **not a result** |
| Pseudobulk arm's test | **moderated eBayes** (Amendment 2), replacing the DESeq2-Wald arm Amendment 1 found miscalibrated |
| Calibration (λ, FP rate) | **criteria met**, now at 200 permutations - the earlier marginal FP reading is resolved, see below |
| Power at the pre-registered effect size | **met, inside a stated envelope only** - 0.86 vs the required ≥ 0.60 at σ_donor = 0.35 with 8 donors/group. At σ_donor = 0.5 it is still **0.35: unmet, and not claimed** |
| Instrument validity | **established only within the operating envelope** declared by [Amendment 3](https://github.com/s4s4s4s/pbcheck/blob/main/docs/AMENDMENTS.md) - minimum donors/group 4 / 8 / 13 / 23 at σ_donor 0.2 / 0.35 / 0.5 / 0.7. Outside it this instrument makes no claim |
| Whether any *real* stratum falls inside that envelope | **unknown.** σ_donor is an unanchored knob of our own simulator - open since Amendment 1, still open |
| GO / NO-GO | not taken |

The Phase 0 measurement continues after this release and is not a gate for it; nothing in the
release changes a number `scripts/synthetic_gate.py` prints.

Reproduce the calibration run yourself - CPU only, no downloads, deterministic, about nine minutes:

```bash
python scripts/synthetic_gate.py
```

On a synthetic null with donor structure where the true number of DE genes is **exactly zero**, and a synthetic
positive at the pre-registered effect size (log2FC = 1.0, K = 200), default parameters, most recent run 511.8s
at 200 permutations ([full artifact](https://github.com/s4s4s4s/pbcheck/blob/main/pilot/gate/synthetic_gate_2026-08-15.json)):

| quantity | value | reads as |
|---|---|---|
| naive per-cell λ | **54.57** | grossly inflated |
| naive false-positive floor | **1162 / 1500 genes (77.4%)** | the error calls three quarters of the genome |
| donor-pseudobulk floor | **0** | correct FDR control under the null |
| pseudobulk λ | **1.01** | ✅ inside the pre-registered [0.9, 1.1] - was 1.25 under the retired DESeq2-Wald test |
| pseudobulk perm-null FP rate | **0.035** at target ≤ 0.05 | ✅ 7/200 permutations, Monte-Carlo SE 0.013, exact 95% CI [0.014, 0.071]. The previous 2/40 reading was compatible with a true rate of 0.12 (P = 0.13); at 200 permutations that is excluded (P = 2·10⁻⁵) |
| pseudobulk power (log2FC = 1.0, K = 200), **at σ_donor = 0.35** | **0.86** | ✅ above the required ≥ 0.60 - but read the envelope caveat below |
| the same power at σ_donor = 0.5 | **0.35** | ❌ below 0.60. Unchanged, unclaimed, and the reason the envelope exists |

The first two lines are why the tool is worth building. The calibration lines are why the test was replaced
(Amendment 2) and why replacing it worked - and the FP rate is no longer a marginal reading leaning on 40
permutations: Amendment 3 raised the count to 200 *before* the run, and the number resolved cleanly below α.

**The power line needs its caveat stated, not buried.** It is 0.86 because Amendment 3 re-scoped *where*
§8(c)'s threshold binds - to σ_donor = 0.35, the boundary of a declared operating envelope - not because
anything got easier. The pre-registered effect size (log2FC = 1.0, K = 200) and the ≥ 0.60 bar are
**unchanged**. At σ_donor = 0.5 the arm still delivers 0.35, no test in the committed selection grid clears
0.60 there at any donor count tested, and that failure stands on the record. What changed is that the
instrument now states the region it is valid in instead of being judged at an arbitrary simulator setting:

| σ_donor | 0.2 | 0.35 | 0.5 | 0.7 |
|---|---|---|---|---|
| minimum donors per group | 4 | 8 | 13 | 23 |

This **narrows** what pbcheck claims. A stratum near σ_donor ≈ 0.5 will need ≥ ~13 donors per group or must
be excluded - a stricter rule than existed before. And because σ_donor has never been anchored to real data,
**whether real strata land inside this envelope is still unknown**; Amendment 3 supplies the mechanism for
estimating it per stratum but explicitly not the anchor. The gate prints the envelope next to its verdict and
reports `INSTRUMENT VALID WITHIN THE STATED OPERATING ENVELOPE` - never an unqualified "valid".

### The measurement plan

pbcheck's value rests on an empirical claim: that across real public datasets, naive per-cell DE is *massively and
consistently* more inflated than a correct analysis. Before building the full tool we **measure** this over a
first pass of 8 - 12 datasets from the [CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/),
selected to span the outcome space rather than to cherry-pick (spec §1):

- **GO** - inflation is large and consistent → build the full auditor and publish the "map of false discoveries".
- **NO-GO** - inflation is small or erratic → reformat or pivot.

The pilot's measurement engine is checked first against **synthetic oracles with known ground truth** (a correct
auditor must report high inflation on a synthetic null and none on a synthetic positive) - only then is it trusted
on real data. Statistics are verified by hand against these oracles, not assumed.

## What exists, and what does not

Built and exercised by the test suite:

- **Design auditor** (metadata only): reads an `AnnData` `.obs` and flags the ingredients of pseudoreplication - donors per condition, nesting of donor within condition, batch⟂condition confounding measured **per donor**,
  group imbalance. Needs no counts, which is what lets pbcheck audit a published study from its metadata alone.
- **Both DE arms**: the naive per-cell test (the careless default, reproduced faithfully) and donor pseudobulk,
  corrected over one shared gene universe at one alpha so the two differ only by the test.
- **Donor-permutation null** and the inflation metrics (genomic-inflation λ, the false-positive floor).

- **Census stratum selection** (`census_select`, obs only, no network in the tests): opens the Census at the
  pinned version `2025-01-30` (never a mutable alias), forms `(dataset_id × cell_type)` strata with one
  `disease vs normal` contrast per disease term, applies the inclusion-gate items metadata can decide,
  pre-screens confounding per donor and emits a **candidate** manifest.

- **The counts gate** (`io_counts`, no network in the tests): the two inclusion-gate items only the count
  matrix can decide. Raw integrality is checked on the *values* rather than the dtype - Census raw is
  `float32` holding 0.0/1.0/7.0, and those are counts - with a normalised or log-transformed matrix
  **dropped, never rounded**, and only the stored entries of a sparse `X` scanned. The frozen universe is
  then sized with the arm's own aggregation and the C5 minimum applied, and the pending manifest columns
  (integer check, universe size, counts per cell, sequencing-depth bin) are filled. Where the load disagrees
  with the obs snapshot it was planned from, the row keeps its committed numbers and gains a discrepancy
  flag. It admits nothing either: `sigma_donor` and envelope membership stay pending on the anchor
  Amendment 3 leaves open, so every row still carries `admitted_to_sweep = False`.

- **The candidate run over the whole Census** (`scripts/census_candidates.py`, dispatched manually from the
  `census candidates` workflow): two passes, because one query does not fit the machine. A streaming pass
  folds ~50-65 M cells into per-`(dataset, cell type, disease, donor)` counts without ever materialising
  them, deciding only *which* datasets are worth reading; a second pass re-reads each of those and hands
  the per-cell frame to `census_select.screen_strata` unchanged. The coarse filter is conservative by
  construction - it may admit a dataset the gate then rejects, never the reverse - and the test suite
  cross-checks that against the real inclusion gate rather than against a hand-written expectation. Its
  output is a CI artifact (manifest JSON + CSV + run log), never a commit, and it admits nothing. **The run
  itself has not been made yet**; the driver and its job exist.

Built in v0.1.0, outside the protocol: the single-stratum audit (`pbcheck.audit`), its
Markdown/HTML report (`pbcheck.render`), the `pbcheck` CLI, the metadata-only mode and the
behaviour on non-integer counts (pseudobulk arm dropped, never rounded). Specified but not
built: the Phase 0 real-data harness (`controls`, `decision`, the spec section 9 `report`). No
`risk_score` exists and none is promised. See [`pilot/README.md`](https://github.com/s4s4s4s/pbcheck/blob/main/pilot/README.md), which also
lists the spec corrections that are named in module docstrings but not yet implemented, so that
nothing here reads as done when it is not.

## Demonstrations

Two committed demonstrations under `demo/` show pbcheck running against real public datasets
outside the pre-registered Phase 0 protocol, including one where the tool's own design audit
stops the run because the dataset uses a paired design. Both are labelled as demonstrations, not
as a Phase 0 measurement, and quote no result about the underlying publications. See
[`demo/README.md`](https://github.com/s4s4s4s/pbcheck/blob/main/demo/README.md) for the exact
commands, the datasets, their licences and the full disclaimer. `demo/` itself is produced by the
demonstration scripts of this release (`scripts/demo_kang2018.py`, `scripts/demo_two_arm.py`), not
committed by hand.

## Layout

```
src/pbcheck/          # the auditor engine (methods/) and the release: audit.py, cli.py, render/
pilot/                # Phase 0 artifacts: committed test-selection grid + gate runs, no code
scripts/              # the harness: synthetic gate, calibration probe, selection analyzer, Census candidate run
synthetic/            # synthetic-oracle generators with known ground truth
demo/                 # v0.1.0 demonstrations on real public data, outside the protocol
tests/                # unit + regression + property tests (oracles are the correctness spec)
docs/                 # PHASE0_SPEC.md (methodology), AMENDMENTS.md, ENV_NOTES.md, USAGE.md (CLI reference)
```

## Development install

Python 3.12 or newer is required: the verified stack (scanpy 1.12, anndata 0.13, pandas 3.0) does not resolve
below it.

```bash
python -m venv .venv

# POSIX
.venv/bin/python -m pip install -e ".[dev]"        # engine + oracles, no downloads
.venv/bin/python -m pip install -e ".[census]"     # add this to run the real-data audit

# Windows
.venv/Scripts/python -m pip install -e ".[dev]"
.venv/Scripts/python -m pip install -e ".[census]"
```

Reproduce the instrument's calibration run on synthetic data with known ground truth - CPU only, no downloads,
deterministic:

```bash
python scripts/synthetic_gate.py
```

## References

- Squair, J. W. et al. *Confronting false discoveries in single-cell differential expression.* Nat. Commun. (2021).
- Zimmerman, K. D. et al. *A practical solution to pseudoreplication bias in single-cell studies.* Nat. Commun. (2021).
- Murphy, A. E. & Skene, N. G. *A balanced measure shows superior performance of pseudobulk methods…* Nat. Commun. (2023).

## License

BSD-3-Clause. See [LICENSE](https://github.com/s4s4s4s/pbcheck/blob/main/LICENSE).

<!-- DOI badge added after the Zenodo release -->
