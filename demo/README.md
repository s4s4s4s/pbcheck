# `demo/`: pbcheck v0.1.0 demonstrations

These are demonstrations outside the pre-registered Phase 0 protocol. Neither dataset is among the
17 frozen datasets or the 357 strata of `pilot/preregistration/stratum_list_2026-08-16.csv`; nothing
here is a Phase 0 measurement, no frozen data was opened, and no number here may be read as a
statement about the datasets' publications. Kang et al. 2018 is a paired design (each donor is
measured under both conditions); pbcheck stops at the design audit for such designs, and that stop
is the demonstration.

Both demonstrations run `pbcheck.audit.run_audit` and `pbcheck.render.write_outputs`, the same code
path the `pbcheck` command line uses, driven by `scripts/demo_kang2018.py` and
`scripts/demo_two_arm.py` (`scripts/demo_common.py` holds what the two share). Neither script runs
in CI; they are a developer's own session, and only their three output files per run are committed
here, never the input data (`data/` is git-ignored).

## Demonstration 1: the design gate on Kang et al. 2018 (`demo/kang2018_design_gate/`)

Command:

```
python scripts/demo_kang2018.py
```

which downloads the pertpy-prepared Kang et al. 2018 IFN-beta stimulation artifact (figshare API
file id 34464122; see the script's module docstring for the source note) and runs

```
pbcheck audit data/kang_2018.h5ad --donor replicate --condition label --test stim --ref ctrl \
    --celltype cell_type --celltype-value "CD4 T cells"
```

Chosen cell type: `CD4 T cells` (the largest of the file's 8 cell-type levels). Each of the file's 8
donors is measured under both `label` levels (`ctrl` and `stim`), so the design audit stops the run
before either differential-expression arm executes; that is what this demonstration shows.

Runtime: 0.42 s (the audit stage alone, from `pbcheck_audit.json`'s own `runtime_seconds`).
Machine: Windows-11-10.0.26200-SP0, Intel64 Family 6 Model 154 Stepping 3, 20 logical CPUs.

Citation: Kang, H., Subramaniam, M., Targ, S. et al. Multiplexed droplet single-cell
RNA-sequencing using natural genetic variation. Nat Biotechnol 36, 89-94 (2018).
doi:10.1038/nbt.4042

Licence: not confirmed by `scripts/demo_kang2018.py`; see the figshare item's own page for its
licence terms. Raw data is not redistributed by this repository.

## Demonstration 2: a two-arm audit on an unpaired dataset

Candidate order, fixed by the release plan: Stephenson et al. 2021 first, `example_reference`
(synthetic) only if that candidate cannot be used.

### Attempt 1: Stephenson et al. 2021, stopped before the design audit

Command:

```
python scripts/demo_two_arm.py --dataset stephenson2021 --out demo/stephenson2021_two_arm
```

The numbers in this subsection (the download size, the cell count, the value range) are the
developer's own observations from running this script once on a laptop on 2026-09-06, before the
attempt aborted; no console log or `pbcheck_audit.json` from that run is committed anywhere under
`demo/`, since the run never reached `run_audit`, so none of these numbers is independently
checkable from a committed artifact the way the two runtime numbers in this README are (those trace
to each committed `pbcheck_audit.json`'s own `runtime_seconds`). Treat them as unverifiable
developer-session notes, not as measurements this repository can reproduce a check against.

The script confirmed the download size (observed: 700,549,709 bytes, about 0.65 GiB, well under the
3 GB bound) from the figshare API endpoint recovered from pertpy's git history for the
`stephenson_2021_subsampled` artifact (file id 38171703; pertpy's current loader downloads the same
file from `https://exampledata.scverse.org/pertpy/stephenson_2021_subsampled.h5ad` instead, the same
figshare-to-scverse migration `scripts/demo_kang2018.py` documents for Kang 2018), downloaded and
verified it by sha256, and confirmed the `patient_id` (donor), `Status` (condition) and `cell_type`
columns exist.

Reason for stopping (exit code 3): `no_raw_counts_anywhere`. `adata.X` is `float32` but not
integer-valued (observed: log-normalised, values up to about 7.5 on 62,509 cells), the file has no
named `.layers` and no `.raw`, so `pbcheck.io_counts.check_integer_counts` failed on `X` and there
was no second matrix to check it on. Per the release plan, this stops the attempt without opening a
second real dataset on the auditor's own word; `example_reference` is the committed run instead.

Citation (recorded for completeness; no number from this attempt is used anywhere): Stephenson, E.,
Reynolds, G., Botting, R.A. et al. Single-cell multi-omics analysis of the immune response in
COVID-19. Nat Med 27, 904-916 (2021). doi:10.1038/s41591-021-01329-2. Licence: not confirmed by
`scripts/demo_two_arm.py`; no data from this dataset is redistributed by this repository, and none
of it is committed here.

### Committed run: `example_reference` (`demo/example_reference_two_arm/`)

`example_reference` is synthetic: `pbcheck.example.example_adata` at the `REFERENCE_SHAPE` the
`pbcheck example --shape reference` CLI command also builds (8,000 genes, 8 donors per group, 625
cells per donor), seed 0. It has one cell-type level, so the "largest cell type" selection is
trivial by construction; it carries no donor spanning both conditions, so the design audit was not
the reason this dataset was used.

Command:

```
python scripts/demo_two_arm.py --dataset example_reference --out demo/example_reference_two_arm
```

equivalent to

```
pbcheck example --shape reference --seed 0 --out data/example_reference.h5ad
pbcheck audit data/example_reference.h5ad --donor donor --condition condition --test disease \
    --ref ctrl --celltype cell_type --celltype-value example
```

Chosen cell type: `example` (the dataset's only level).

Runtime: 89.2 s (the audit stage alone, from `pbcheck_audit.json`'s own `runtime_seconds`); about
103 s wall clock for the whole script, including generating the synthetic data.
Machine: Windows-11-10.0.26200-SP0, Intel64 Family 6 Model 154 Stepping 3, 20 logical CPUs.

Licence: none; `example_reference` is generated by this repository's own code
(`src/pbcheck/example.py`), not a third-party dataset.

Defaults used by both demonstration 2 attempts: `--n-perm PRODUCT_N_PERM --n-perm-pb
PRODUCT_N_PERM_PB --seed 0` (`src/pbcheck/audit.py`'s own product constants; this file never types
their values, see the note against numeric tokens above).
