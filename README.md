# pbcheck

**An auditor of pseudoreplication in single-cell RNA-seq differential expression.**

> Status: **Phase 0 — pilot audit.** This repository is currently a reproducible measurement study, not yet a
> released tool. The Phase 0 pilot decides whether the tool is worth building (see *The gate* below).

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

## What pbcheck will do

- **Design auditor** (metadata only): inspect an `AnnData` `.obs` and flag the ingredients of pseudoreplication —
  number of donors per condition, nesting of donor within condition, batch⟂condition confounding, group imbalance.
- **Inflation estimate**: compare a naive per-cell DE against a correct pseudobulk DE and against a
  **donor-level permutation null** (where the true number of DE genes is ≈ 0), and quantify the inflation.
- **Report**: a per-dataset `risk_score ∈ [0, 1]` and an HTML report.
- A mode that works **without raw counts**, auditing a published result table plus its cell metadata.

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
