# pbcheck

**An auditor of pseudoreplication in single-cell RNA-seq differential expression.**

> Status: **Phase 0 — a pre-registered measurement study, not a released tool.** The decision rule was frozen
> in [`docs/PHASE0_SPEC.md`](docs/PHASE0_SPEC.md) *before* the data was seen; every departure from it is dated
> and justified in [`docs/AMENDMENTS.md`](docs/AMENDMENTS.md).

### Where this actually stands

| | |
|---|---|
| Measurement engine | built; 21 tests pass |
| Real CELLxGENE data | **not touched.** No dataset has been run, no stratum list pre-registered |
| Every number below | from synthetic oracles with known ground truth — instrument calibration, **not a finding** |
| Instrument validity | **not yet established.** The pseudobulk arm fails its own binding gate; see [Amendment 1](docs/AMENDMENTS.md) |
| GO / NO-GO | not taken |

Reproduce the calibration run yourself — CPU only, no downloads, deterministic, about five minutes:

```bash
python scripts/synthetic_gate.py
```

On a synthetic null with donor structure where the true number of DE genes is **exactly zero**:

| quantity | value | reads as |
|---|---|---|
| naive per-cell λ | **54.8** | grossly inflated |
| naive false-positive floor | **1166 / 1500 genes (77.8%)** | the error calls three quarters of the genome |
| donor-pseudobulk floor | **0** (mean 0.60) | correct FDR control under the null |
| pseudobulk λ | 1.25 | ⚠️ outside the pre-registered [0.9, 1.1] — the arm is not yet valid |

The first two lines are why the tool is worth building. The last line is why it is not finished: the pseudobulk
arm is the denominator of every inflation number, and until it is calibrated, none of them may be read as a
result. That is the honest state, and the gate script reports it as `INSTRUMENT NEEDS ATTENTION` rather than
quietly passing.

## The problem

A large fraction of published single-cell RNA-seq differential-expression (DE) results are contaminated by
**pseudoreplication**: individual cells are treated as independent replicates, when the true unit of replication
is the **donor / sample**. Because cells from one donor are correlated, per-cell tests (e.g. a Wilcoxon test
across cells, as in the default `scanpy` workflow) drastically inflate the false-discovery rate — reporting
hundreds to thousands of "differentially expressed" genes that do not replicate. The correct analysis aggregates
counts to one profile per donor per cell type (**pseudobulk**) and tests across donors.

This is well documented in the methods literature (Squair et al. 2021; Zimmerman et al. 2021; Murphy & Skene
2023), yet the naive per-cell approach remains widespread in applied papers. What appears to be missing — I am
not aware of one, in Python or R — is a tool that takes an *existing* analysis and estimates *how much* its DE
results are inflated by this error. Pointers to prior art are welcome and will be credited here.

## What exists, and what does not

Built and exercised by the test suite:

- **Design auditor** (metadata only): reads an `AnnData` `.obs` and flags the ingredients of pseudoreplication —
  donors per condition, nesting of donor within condition, batch⟂condition confounding measured **per donor**,
  group imbalance. Needs no counts, which is what lets pbcheck audit a published study from its metadata alone.
- **Both DE arms**: the naive per-cell test (the careless default, reproduced faithfully) and donor pseudobulk,
  corrected over one shared gene universe at one alpha so the two differ only by the test.
- **Donor-permutation null** and the inflation metrics (genomic-inflation λ, the false-positive floor).

Specified but **not** built — the real-data harness (`census_select`, `io_counts`, `controls`, `decision`,
`report`), the `risk_score`, the HTML report, and the no-raw-counts mode. See
[`pilot/README.md`](pilot/README.md), which also lists the spec corrections that are named in module
docstrings but not yet implemented, so that nothing here reads as done when it is not.

## The gate (Phase 0)

pbcheck's value rests on an empirical claim: that across real public datasets, naive per-cell DE is *massively and
consistently* more inflated than a correct analysis. Before building the full tool we **measure** this over a
first pass of 8–12 datasets from the [CELLxGENE Census](https://chanzuckerberg.github.io/cellxgene-census/),
selected to span the outcome space rather than to cherry-pick (spec §1):

- **GO** — inflation is large and consistent → build the full auditor and publish the "map of false discoveries".
- **NO-GO** — inflation is small or erratic → reformat or pivot.

The pilot's measurement engine is calibrated first on **synthetic oracles with known ground truth** (a correct
auditor must report high inflation on a synthetic null and none on a synthetic positive) — only then is it trusted
on real data. Statistics are verified by hand against these oracles, not assumed.

## Layout

```
src/pbcheck/          # the auditor engine (methods/)
pilot/                # Phase 0 harness: dataset selection, sweep, results
synthetic/            # synthetic-oracle generators with known ground truth
tests/                # unit + property tests (oracles are the correctness spec)
docs/                 # PHASE0_SPEC.md (methodology), ENV_NOTES.md
```

## Install (development)

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

Reproduce the instrument's calibration run on synthetic data with known ground truth — CPU only, no downloads,
deterministic:

```bash
python scripts/synthetic_gate.py
```

## References

- Squair, J. W. et al. *Confronting false discoveries in single-cell differential expression.* Nat. Commun. (2021).
- Zimmerman, K. D. et al. *A practical solution to pseudoreplication bias in single-cell studies.* Nat. Commun. (2021).
- Murphy, A. E. & Skene, N. G. *A balanced measure shows superior performance of pseudobulk methods…* Nat. Commun. (2023).

## License

BSD-3-Clause. See [LICENSE](LICENSE).
