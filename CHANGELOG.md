# Changelog

All notable changes to this project are documented here. Format loosely follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). There are no tagged releases yet, so
everything lives under a single Unreleased section, grouped by the project's actual epochs of
work rather than by commit date — see `git log` for the literal timeline.

## [Unreleased]

### Initial engine + Phase 0 spec (2026-07-22)

- Added: project scaffolding (packaging, license, pinned environment), the frozen
  `docs/PHASE0_SPEC.md` pre-registration, and `docs/ENV_NOTES.md` recording the verified stack
  and upstream API drift.
- Added: the measurement engine — synthetic oracles with known ground truth, the frozen gene
  universe and shared BH correction, the two DE arms (naive per-cell, donor pseudobulk via
  DESeq2), the donor-permutation null and inflation metrics (λ, FP floor), and the obs-only
  design auditor.
- Added: the test suite, `scripts/synthetic_gate.py`, `scripts/pb_calibration_probe.py`, the
  README, the pilot status doc, and CI.
- Added: **Amendment 1** (`docs/AMENDMENTS.md`), recorded the same day the engine first ran end
  to end: the pseudobulk arm (then DESeq2-Wald) fails its binding validity gate (λ 1.21–1.25,
  permutation-null FP rate 0.35–0.50 against α = 0.05); the FP-rate criterion is restored as
  binding, the power oracle held at pre-registered parameters, and a statistics error in
  `PILOT_FINDINGS.md` retracted.
- Fixed: the design auditor's Cramér's V is measured per donor rather than weighted by cell
  count (`aecf8c5`), closing a confounding-detection blind spot Amendment 1's diagnosis exposed.

### Calibration probe and test-selection grid (2026-07-23)

- Added: the test-selection grid (`scripts/run_test_selection_grid.py`,
  `pilot/testsel/summary.{csv,json}`) sweeping test × `sigma_het` × (`sigma_donor`, `n_donors`)
  × generative arm, and `scripts/analyze_test_selection.py` to turn the grid into a verdict. The
  full 146-cell run is committed (`72dec7b`, `failed = 0`) — this is the evidence base Amendment
  2 later applies its selection rule to.
- Fixed: the naive arm now passes `tie_correct=True` / `pts=True` as the spec pins (§2); both
  arms corrected over one common tested set with the BH null reference pinned; the minimum
  gene-universe gate made enforced rather than decorative; Monte-Carlo error reported on the
  permutation floor; the real split checked against the permutation cell-count range.
- Changed: README rewritten to lead with actual state and numbers rather than aspirational
  status.

### 2026-08-15 remediation

One large pass closing gaps between what module docstrings claimed and what the code did,
culminating in Amendment 2.

- Fixed (R0, pre-Amendment-2 bug fixes): `build_perms` no longer hangs when available label sets
  are scarce but exceed `max_enumerate`; `min_donors` is threaded through to the
  `usable_for_pseudobulk` verdict rather than just its message text; `lambda_over_permutations`
  returns a consistent key set on the empty and normal paths; a relative `sys.path` hack removed
  from `test_permutation_mtc.py`; the blanket `Future`/`DeprecationWarning` filter replaced with
  one targeted, audited ignore; the dead `restrict_to` parameter removed; stale docstring claims
  (`min_cells`/`min_counts`, the calibration probe's status, the R1 rename) corrected to match
  the code.
- Added: **Amendment 2** (`docs/AMENDMENTS.md`), closing Amendment 1's open question and settling
  five further items claimed but not implemented. The pseudobulk arm's test is replaced with
  **moderated eBayes** (limma-style, Smyth 2004), selected by a pre-stated minimax-regret rule
  applied to the test-selection grid; the **paired BH** (`bh_both_arms`) is wired into `run_null`
  and the gate for the first time — an erratum, since every previously published number used
  per-arm BH over different tested sets; **B5** restored, its literal construction measured and
  shown to carry zero calibration information, and implemented instead as a machinery-only
  diagnostic; **C5** restored (frozen-universe `min_size` now enforced by every caller, the
  realised prior `d0`/shrinkage persisted as the moderated analog of `fitType`); **C3** retired
  with rationale (DESeq2-specific, no referent once DESeq2 leaves the arm; replaced by
  prior-strength disclosure and a null p-value uniformity check); **A2** deferred to Phase 1,
  honestly (a range check only; full cell-count stratification not implemented, with the
  rationale recorded); the pre-registered thin-donor filter (`min_cells`/`min_counts`)
  implemented manually, since decoupler 2.x dropped the parameters.
- Changed: docs truth sync — README and `pilot/README.md` updated to Amendment 2's actual gate
  numbers (λ_pseudobulk 1.02, permutation-null FP rate 0.05, power 0.35 against the required
  0.60); the Amendment 2 gate rerun recorded as a dated addendum rather than a rewrite of the
  amendment.
- Added (infrastructure): dependency floors pinned to the verified stack, with the `report`
  extra split out of runtime deps; the version single-sourced from `src/pbcheck/__init__.py`
  (`tool.hatch.version`), with a CI check that `CITATION.cff` doesn't drift from it; a curated
  public API (`__all__` / re-exports); the CI matrix extended to macOS, coverage made
  report-only with no threshold gate (deliberate — see CONTRIBUTING.md), and a build/twine
  sdist+wheel job added; `requirements.lock` regenerated with provenance and a
  lock-reproducibility CI job added; workflow permissions restricted and a real `ruff check`
  lint job added; `.pre-commit-config.yaml` added (ruff check plus hygiene hooks, with
  `docs/PHASE0_SPEC.md` and `docs/AMENDMENTS.md` excluded so their bytes are never mechanically
  touched); shared oracle fixtures added to `conftest.py`; a `slow` marker registered for
  DESeq2/multi-permutation end-to-end tests; Hypothesis property tests added for `build_perms`
  and BH monotonicity; self-tests added for the synthetic oracles; the grid driver's git side
  effects de-fanged (opt-in `--git-publish`, errors no longer swallowed).

### Amendment 3 — the validity gate passes, within a stated envelope (2026-08-15)

- Added: **Amendment 3** (`docs/AMENDMENTS.md`), written and committed before the code that applies
  it. **Change 1** re-scopes §8(c)'s power criterion from a point to a declared **operating
  envelope**: the pre-registered effect size (log2FC = 1.0, K = 200) and the ≥ 0.60 threshold are
  **unchanged**, and what changes is the region over which the instrument claims validity —
  calibration stays binding at the hard regime (`sigma_donor` = 0.5), power becomes binding at the
  envelope boundary (`sigma_donor` = 0.35, 8v8, where the committed grid measures `ebayes` at
  0.793), and the gate must report the envelope (minimum donors/group 4 / 8 / 13 / 23 at
  `sigma_donor` 0.2 / 0.35 / 0.5 / 0.7, from Amendment 1's frontier). This **narrows** the claimed
  validity domain rather than lowering a bar: power 0.60 at `sigma_donor` = 0.5 with 8 donors/group
  stays unmet at 0.35 and unclaimed, and strata near 0.5 will need ≥ ~13 donors/group or exclusion.
  Real-sweep stratum inclusion now also requires a per-stratum `sigma_donor` estimate; the amendment
  supplies the mechanism (the eBayes fit's own between-donor variance) and states plainly that the
  **anchor remains open**. **Change 2** raises the gate's permutation count 40 → 200, declared before
  the rerun that applies it.
- Changed: `gate_config` gains `CALIBRATION_EVAL_SIGMA`, `POWER_EVAL_SIGMA` and `OPERATING_ENVELOPE`
  (each row carrying its grid corroboration); `N_PERM` and `N_PERM_PB` both go 40 → 200, since
  `run_null` pairs `min(n_perm, n_perm_pb)` and raising one alone would have changed nothing.
  `scripts/synthetic_gate.py` runs the two arms at their two regimes, labels which is which
  everywhere, prints the envelope next to the verdict, and records each criterion with its
  evaluation regime and source document.
- Added: the gate run artifact `pilot/gate/synthetic_gate_2026-08-15.json` — **all six criteria pass**
  (511.8 s): λ_pseudobulk 1.01, permutation-null FP rate 0.035 (7/200, MC SE 0.013, exact 95% CI
  [0.014, 0.071]), pseudobulk floor 0 against a naive floor of 1162/1500 (77.4%) at λ_naive 54.57,
  and power 0.86 at the envelope boundary. Change 2 did what it was declared to do: the old 2/40
  reading was compatible with a true FP rate of 0.12 (P = 0.13), which 200 permutations exclude
  (P = 2·10⁻⁵).
- Changed: docs truth sync — README and `pilot/README.md` updated to the Amendment 3 gate numbers and
  to the envelope caveat, with the verdict stated as *valid within the stated operating envelope* and
  never as an unqualified "valid".
- Fixed: `naive.py` now recovers **both** expressed fractions the min-expression floor sensitivity
  check needs (§6/B5). scanpy 1.12.2's `rank_genes_groups_df` drops `pct_nz_reference` when an
  explicit `reference` is passed, so the reference fraction never reached the result table; but
  `uns['rank_genes_groups']['pts']` is a genes x groups frame that carries every level of the
  groupby. The arm reads that frame directly and reindexes it from `var_names` order onto the ranked
  gene order, emitting `pct_group` and `pct_reference` for both the Wilcoxon and t-test variants;
  the dataframe columns stay as a fallback. No amendment: the spec already required both fractions,
  so this is an implementation fix, not a protocol change. `tests/test_methods.py` now asserts both
  columns and checks them against per-condition nonzero rates computed from the matrix (which is
  what guards the reindex), replacing the test that recorded the gap.

### Real-data harness, module 1 — Census stratum selection (2026-08-15)

- Added: `pbcheck.census_select` (spec §9 item 1), the obs-only front of the real-data harness.
  `open_census()` refuses anything but the pinned `2025-01-30` and refuses the mutable aliases by
  name; `query_obs()` materialises the §1 columns, records the ones the Census schema does not
  expose instead of faking them, and picks up any library/pool identifier it does (D3); strata are
  `(dataset_id × cell_type)` with **one binary `disease vs normal` contrast per disease term**, never
  a pooled "any disease" arm; `apply_inclusion_gate()` drops thin donors (< 10 cells in that cell
  type — dropped, not merged) before applying the ≥ 3 donors/group rule, and rejects a donor that
  spans both conditions or a group carrying one donor; `confound_prescreen()` measures Cramér's V and
  perfect separation **per donor** — reusing `pbcheck.design`'s helpers rather than copying them —
  over assay, suspension type, tissue, a sequencing-depth bin and any pool id, **restricted to the
  levels that group two or more donors** (a level held by a single donor is a relabeling of
  `donor_id`, which is nested within condition by design, so scoring it would have excluded a
  40-donor stratum over a two-donor library coincidence), excluding on the three §1 names and
  tagging the rest; `emit_manifest()` writes the §1 manifest as JSON + CSV with the
  census version, package versions, the pre-registered config block and D4's excluded fraction in
  the header.
- Added: `tests/test_census_select.py` (50 tests, no network, synthetic obs frames) — including the
  pin test that no mutable census alias appears as an argument anywhere in the module, that every
  column the module cannot compute is written as `pending`, and two regression tests from the
  adversarial review: a singleton pool coincidence among 40 donors must tag and not exclude while a
  structural pool confound still excludes, and cells with a null `donor_id` must be dropped and
  counted — under pandas 3.0's `future.infer_string` a null survives `.astype(str)`, so such cells
  were passing the gate unseen although §1 item 3 requires `donor_id` present.
- Note: this module **admits nothing**. Inclusion-gate items 4 and 5 (integer counts, frozen-universe
  size) need `X` and are emitted as `pending` for `io_counts` / `gene_universe`; `sigma_donor` and
  envelope membership are `pending` on the anchor Amendment 3 explicitly leaves OPEN; every row
  carries `admitted_to_sweep = False`. It also publishes **no** covariate-affordability column:
  §1's C4 df rule lapsed with DESeq2 (Amendment 2 Change 1), the shipped moderated arm fits
  `~ 1 + x`, and a partial confound is tagged and neutralised by the permutation null, not
  modelled. Committing this manifest is not the §1 pre-registration of the
  stratum list.

### Real-data harness, module 2 — the counts gate (2026-08-15)

- Added: `pbcheck.io_counts` (spec §9 item 2), the `X`-reading half of the selection pipeline. It
  settles the two §1 inclusion-gate items obs cannot decide and fills the columns `census_select`
  had to leave `pending`. **Item 4** (`check_integer_counts` / `assert_integer_counts`): integrality
  is tested on the *values*, not on the dtype, because Census raw is routinely `float32` holding
  0.0/1.0/7.0 and rejecting a float dtype would throw away most of the Census; fractional values
  (a normalised or log1p matrix — §10 risk 6), non-finite values and negative values each fail with
  their own reason, and the response is a **DROP recorded in the manifest, never a rounding**. On a
  sparse matrix only the stored entries are scanned and nothing is ever densified, so the check
  costs the nonzeros rather than the 12 GB dense form of a real stratum. **Item 5 / C5**
  (`frozen_universe_check`): the universe is built with the arm's own code —
  `methods.pseudobulk.build_pseudobulk` (thin donors dropped, not merged, at the pre-registered
  `MIN_CELLS`/`MIN_COUNTS`) then `gene_universe.frozen_universe` — and `require_min_size` turns a
  small universe into a SKIP whose *measured size is still recorded*, since "too small" without the
  number cannot be audited. §3's "≥ 3 pseudosamples per group post-aggregation" is measured at the
  same time and reported as its own skip status: a donor can clear item 2's ≥ 10 cells and still
  fall below `min_counts` = 1000, which is only visible once counts are in hand.
- Added: the rest of the pending columns — `median_counts_per_cell_by_group` from `X` (only when it
  was pending; when `census_select` already read it from the Census's own `raw_sum`, the snapshot's
  value is kept and the X-derived one recorded beside it), and the §1 pre-screen's
  `sequencing_depth_bin`, computed by calling `census_select.confound_prescreen` with no other
  covariates rather than re-deriving the donor-level Cramér's V and the singleton-level rule a
  second time. `update_manifest_row(row, adata)` returns a **new** row and never mutates its input,
  `admitted_to_sweep` stays `False` on every row (Amendment 3's `sigma_donor` anchor is still OPEN),
  and `fitType` / `cell_type_ontology_depth` stay `pending` because they belong to the pseudobulk
  arm (C5 / Amendment 2 Change 4) and to `controls` (D5).
- Added: `reconcile_with_obs_snapshot()` — the obs snapshot and the load are **separate acts**, so a
  load with fewer cells or a missing donor produces a discrepancy **flag** and the row keeps the
  numbers it was committed with. Two differences are recognised rather than flagged: donors the obs
  gate itself recorded as dropped-thin (a load filters on dataset/cell type/disease, not on donor,
  so they come back), and a per-cell depth agreeing with the snapshot's `raw_sum` medians to 0.1%.
  The first allowance is **bounded and the bound is checked**: a donor the snapshot dropped at 4
  cells that arrives with 400 is not that decision reapplied but the two acts disagreeing about the
  fact §1 item 2 was applied to, so it is compared against `MIN_CELLS`, listed, flagged, and
  `agrees` is False. And when the integrality gate has failed the reconciliation runs with
  `compare_magnitudes=False`: donor and cell counts stay meaningful whatever the values are, but no
  loaded median is reported and the depth cross-check cannot fire — otherwise a ×0.5 rescale would
  be reported twice, once truthfully as a failed integer check and once as a spurious "depth
  discrepancy". The "nothing derived from a refused matrix" rule binds at every depth of the row,
  not only at its top level.
- Added: `load_stratum()` / `load_contrast()`, the only network-touching functions — a thin
  `cellxgene_census.get_anndata` wrapper with a lazy import, asking for the **raw** layer (§1) at
  `census_select.CENSUS_VERSION`, imported rather than restated so the pin has one home. It records
  the integrality verdict in `uns` and deliberately does **not** raise on a failure: §1's answer to
  a non-integer stratum is a drop in the manifest, and a loader that aborted the sweep could not
  report D4's excluded fraction. The obs-column argument is chosen by reading the installed
  signature (`obs_column_names` vs the older `column_names={"obs": …}`), the same API-drift
  discipline `docs/ENV_NOTES.md` records for decoupler and PyDESeq2.
- Added: `tests/test_io_counts.py` (52 tests, no network, no `cellxgene-census`) — integral floats
  pass and fractional values drop with their reason; a sparse matrix that raises on `toarray` proves
  the sparse path is never densified; a 50-gene stratum is a C5 SKIP that still reports 50; the
  manifest row is asserted unchanged byte for byte; the Census wrapper is driven against a stub
  module, including the older-client argument fallback and that a handle passed in is not closed.
  Four of them are regressions from the adversarial review of the first draft: no magnitude of a
  refused matrix survives anywhere inside `load_vs_obs_snapshot`, a ×0.5 rescale is not re-reported
  as a depth discrepancy, a donor that is no longer thin in the load is flagged instead of netted
  out, and one that is still thin stays unflagged (so the new check cannot over-fire).
- Note: filling items 4 and 5 removes two admission blockers and **admits nothing**. The extended
  manifest keeps `census_select`'s columns and adds four (`integer_check_detail`,
  `frozen_universe_detail`, `load_vs_obs_snapshot`, `counts_provenance`); `io_counts.emit_manifest`
  exists so those reach the CSV, since `census_select.emit_manifest` writes its own shorter field
  list, and its header carries a `counts_gate` block stating plainly that counts-gate exclusions are
  *not* part of the obs-only excluded fraction next to it (D4).

### Real-data harness — the Census candidate run (2026-08-16)

- Added: `scripts/census_candidates.py`, the driver that runs `census_select` over the **whole**
  pinned Census, and `.github/workflows/census-candidates.yml`, the manual job that runs it. The
  driver adds no method: every gate, threshold, column and exclusion reason is `census_select`'s,
  nothing in either module changed to accommodate it, and no command-line option can move a
  pre-registered number — the option list is asserted by a test, as is the absence of any
  `census_version` argument anywhere in the driver.
- Added: the two-pass streaming scheme, which exists because the naive shape does not fit the
  machine. After the §1 `value_filter` ~50-65 M cells remain, and `query_obs`'s
  `.concat().to_pandas()` peaks somewhere between 8 GB and 30+ GB depending on whether the Census
  returns `donor_id` dictionary-encoded or as Python objects — unknowable in advance on a 16 GB
  runner, and an OOM kill five hours in tells you nothing. **Pass 1** streams the four stratum
  columns straight off the SOMA reader with no `.concat()`, folding each Arrow batch into a
  `(dataset_id, cell_type, disease, donor_id) -> cells` counter that grows with the number of
  distinct keys rather than with the number of cells; addition is commutative, so the result cannot
  depend on where the Census cut its batches (pinned by a test that re-cuts and re-orders the
  stream — batch boundaries are not part of what the pinned version pins). **Pass 2** re-reads each
  surviving dataset with one scoped query and hands the **per-cell** frame to `screen_strata`
  unmodified — the inclusion gate counts *rows*, so feeding it the pass-1 aggregate would report
  every donor as holding a single cell and the whole Census as failing the gate. The coarse filter
  is deliberately conservative (§1 items 2 and 1 only; item 3's nesting rule is left to pass 2),
  and `tests/test_census_candidates.py` cross-checks that implication against the real
  `apply_inclusion_gate` on the same synthetic rows rather than against a hand-written
  expectation, because a filter that silently drops a usable stratum leaves no row, no reason and
  no count behind.
- Fixed (before it could bite): the scoped `value_filter` escapes literals as **Python** literals,
  not SQL ones. TileDB-SOMA parses a value filter with `ast.parse(expr, mode="eval")`, so the SQL
  escape `'O''Brien'` is implicit string concatenation — it parses, means `"OBrien"`, matches
  nothing, and the stratum would vanish from the manifest without an error anywhere. Cell-type
  labels are ontology free text and do carry apostrophes. Every scoped read is additionally checked
  for the keys it asked for and for emptiness where pass 1 counted cells, so a filter that stops
  applying fails loudly instead of producing a quietly short manifest.
- Added: a dataset above `--max-cells-per-dataset-query` (default 2 M) is **split, not skipped** —
  one query per cell type, which partitions the dataset's strata exactly, since a stratum is
  `(dataset_id × cell_type)`. A test asserts the split produces byte-identical manifest rows. The
  split stops there: a cut on `disease` would sever a contrast from its `normal` arm, and a cut on
  `donor_id` would sever a group from its donors. A single cell type over the cap is queried
  anyway — dropping a stratum for being large would bias the candidate list toward small ones.
- Fixed (adversarial review of the first draft): a run in which **no** dataset survived the coarse
  filter used to exit 0. Zero survivors across the whole Census is not an empty result — it is a
  moved threshold, a renamed stratum column or a `REFERENCE_LEVEL` that no longer matches the
  Census's `normal` — and the manifest it writes is well-formed, nearly empty and indistinguishable
  from a real one, so CI would have been green over a regression. `exit_code()` now fails that run
  loudly, and separately fails one in which every selected dataset failed, while a single flaky
  dataset out of sixty still exits 0 with its name in `notes.failed_datasets`: that run carries its
  result, and discarding five hours of work over one dropped connection would be its own bug.
- Fixed (second live dry run, run 31909947593): pass 1 was SIGTERMed 22 seconds in — exit 143, the
  hosted runner's memory-pressure kill rather than the kernel OOM killer's SIGKILL — before a
  single progress line. Two causes, both memory, both now fixed and both now visible. **The read
  buffers**: cellxgene-census's `DEFAULT_TILEDB_CONFIGURATION` sets `py.init_buffer_bytes` and
  `soma.init_buffer_bytes` to **1 GiB each**, allocated per column before the first batch arrives,
  and pass 1 asks for four string columns. The driver now opens the Census with a 128 MiB budget
  (the value cellxgene-census's own documentation uses in its override example), exposed as
  `--reader-buffer-mb` and as a workflow input so the next dispatch can go lower without a code
  change. `open_census` grew a `tiledb_config` passthrough for it — an engineering knob, never a
  data one: buffer sizes decide how many rows arrive per batch and nothing about which rows exist.
  **The conversion**: every batch was going through `.to_pandas()`, which turns four Arrow string
  columns into Python `object` arrays — one heap-allocated `str` per cell per column, the most
  expensive representation available, for data about to be reduced to group counts. Pass 1 now
  groups in Arrow (`Table.drop_null().group_by(keys).aggregate([(key, "count")])`) and materialises
  only the grouped keys. Dictionary-encoded columns are decoded to their **values** first, because
  Arrow assigns dictionary codes per batch — an accumulator keyed on codes would merge two
  different donors and inflate their cell counts silently, in the direction that manufactures
  candidates. A test pins the Arrow and pandas paths to identical output on data where a donor
  demonstrably carries different codes in different batches.
- Added (so the next failure explains itself): RSS is read from `/proc/self/status` and printed at
  start, after the **first** batch (with its row count and Arrow bytes, either side of the fold),
  every batch for the first five and every twentieth after that — calibrating cells/second should
  not require surviving twenty batches, which the killed run did not. A SIGTERM/SIGINT handler
  prints RSS and peak RSS before exiting, since a memory kill produces no traceback, and
  `resource.getrusage` peak RSS is printed on every exit path and recorded in the manifest notes.
- Fixed: the workflow's job summary never ran. GitHub executes `run:` as `bash -e -o pipefail`,
  and `set -uo pipefail` does not clear `-e`, so the driver failing in the pipe aborted the step
  before the `GITHUB_STEP_SUMMARY` block and the intended `exit "${status}"`. Reporting now lives
  in its own `if: always()` step, which is the only arrangement that survives a step that dies; the
  stack-provenance step tees into the same log, so the uploaded artifact carries the versions a run
  died with even when the driver never started.
- Added: `tests/test_census_candidates.py` (37 tests) — no network, and no `cellxgene-census`,
  including a test that the driver does not import it at module level, which is what keeps the
  suite runnable on the development machine. The Census access path is exercised against a SOMA
  stand-in that *applies* the equality clauses of the filter and can be told to fail a given
  dataset, so a mis-escaped literal comes back empty there exactly as it would on S3 and all three
  exit paths (partial failure, total failure, nothing survived) run end to end through `main()`.
- Fixed (first live dry run, run 31909378806): the §1 `value_filter` was being sent to `obs`
  verbatim, and the pinned Census refuses it — `SOMAError: 'Column organism does not exist in
  schema'`, raised on the reader constructor before a single cell was read. In schema 2.1.0 the
  organism is not an obs column but the **experiment**, `census["census_data"]["homo_sapiens"]`,
  which `query_obs` had always been opening correctly; the clause was therefore both unfilterable
  and redundant. `census_select` now separates the two things that were one string:
  `SPEC_VALUE_FILTER` is §1's text, kept verbatim for the manifest and asserted against the spec
  document, and `VALUE_FILTER` is what executes — the same text minus that one conjunct, pinned by
  a test that reconstructs it from the spec string rather than restating it. The manifest header
  carries **both**, plus `organism_realized_by`, so a reader comparing the applied filter against
  the protocol finds the missing conjunct accounted for instead of absent. The queried population
  is unchanged (every cell in the human experiment is *Homo sapiens*), no threshold or scope moved,
  and this is an executable translation of pre-registered text rather than a protocol change, so it
  carries no amendment. `io_counts.stratum_value_filter` inherits the fix through the same constant
  — its `get_anndata` wrapper was already passing `organism` as the API's own argument, which is
  the correct half of this and the reason nothing else had to move.
- Fixed (same root cause, the reason CI was green over it): the SOMA stand-ins in
  `tests/test_census_select.py` and `tests/test_census_candidates.py` never resolved a filter's
  names against their schema, so they accepted a filter the live parser rejects — a rubber stamp,
  not a stand-in. Both now parse the filter with `ast` and raise `Column X does not exist in
  schema` for any name outside the schema, exactly where tiledbsoma raises it, and their schema is
  the synthetic frame's columns **plus** the real Census columns those frames do not materialise
  (`is_primary_data`, and deliberately not `organism`). Three regression tests: the executable
  filter names no `organism` and is the spec text minus exactly one conjunct; the old text sent to
  a stand-in without that column fails as it failed in production; and every filter the driver
  builds — pass 1's and pass 2's scoped one — resolves entirely against an obs schema.
- Note: the run is manual (`workflow_dispatch`, `dry_run` defaulting to true), not cancellable by
  a second dispatch, capped at 330 minutes against the runner's hard 6-hour kill, and its output is
  an **artifact, never a commit** — `pilot/results/` is gitignored, and a candidate list that
  reached git by way of a CI job would have pre-registered itself by accident. Every row still
  carries `admitted_to_sweep = False`, and the manifest header records that D4's excluded fractions
  are computed over the datasets pass 2 actually read, not over the whole Census.

### The §1 stratum list is pre-registered (2026-08-16)

- Added: `docs/PREREGISTRATION_STRATUM_LIST.md` — the act §1 demands ("Pre-register the stratum
  list before computing any metric"), performed before any metric was computed on any of these
  strata. **Not an amendment**: it changes no threshold and supersedes no section. It fixes
  **12 independent datasets** and the **251 stratum-contrasts** they carry, and from this commit
  every change to either goes through `docs/AMENDMENTS.md`.
- Added: the evidence, committed in full rather than referenced —
  `pilot/preregistration/census_candidates_run31910799023_2026-08-15.{json,csv}`, the whole-Census
  candidate manifest of CI run **31910799023** (sha256 `33f8a800…`, 6 630 446 bytes, and
  `09eb110d…`, 4 513 660 bytes; `generated_utc` 2026-08-15T22:18:37Z, Census `2025-01-30`).
  `pilot/results/` stays gitignored — a manifest that reached git through a CI job "would have
  pre-registered itself by accident" — so the copy lives under `pilot/preregistration/`, put there
  by a deliberate human act with the document attached. A pre-registration whose evidence expires
  with a CI artifact is not auditable.
- Added: **the selection rule is deterministic and is code, not a table.** Every manifest row with
  `gate_status == "candidate"` belonging to one of the twelve datasets is in the analysis set;
  every other row is out. The twelve dataset ids are the whole of the judgement and the 251 strata
  are its arithmetic consequence — no stratum, cell type, disease term or donor-count tier is
  chosen by hand (§10 risk 13). The per-dataset "recommended strata" tables of the 2026-08-16
  proposal document were reading aids and are recorded as **not** the frozen set, so the smaller
  list cannot later be mistaken for the pre-registration.
- Added: `scripts/freeze_stratum_list.py`, which emits
  `pilot/preregistration/stratum_list_2026-08-16.{json,csv}`. It refuses to run unless the source
  artifact's size, sha256 and header stamps all match the pinned values (the hash guards the file,
  the stamps guard the constant), and it **recomputes every figure it declares** — per-dataset
  stratum counts and donor ceilings, the total, the Layer B subsets, the bin occupancy — aborting
  on any disagreement rather than silently re-deriving. The output carries no generation timestamp
  and no package versions, so regeneration is byte-identical on every platform; `--check` verifies
  the committed artifact against a fresh run.
- Added: **the cells-per-donor bins (D1), pre-registered numerically for the first time.** The spec
  references "the pre-registered bins" three times (decision rule item 2, §1 (iii), §7 item 3) and
  never states them, which strictly made §1 (iii) unsatisfiable. They are now
  `[10,30) [30,100) [100,300) [300,1000) [1000,3000) [3000,inf)` — half-decade log bins anchored at
  the inclusion gate's own 10-cells-per-donor floor, chosen independently of these data so no edge
  can have been fitted. All six are occupied over the frozen set's 502 group medians (range
  11.0 … 6671.5). This closes a gap in §1; it is not a change to it.
- Added: **Layer B, the pre-declared truncation, at all four tiers of the operating envelope.**
  Survivors of the twelve: 11 datasets / 227 strata at `sigma_donor` ≈ 0.2, **7 / 150** at ≈ 0.35,
  **5 / 94** at ≈ 0.5, **3 / 30** at ≈ 0.7. **Three of the four are below §1's own "8–12 datasets"
  floor** — every tier except the most optimistic — and the failure starts at ≈ 0.35, which is
  `gate_config.POWER_EVAL_SIGMA`, the envelope boundary Amendment 3 Change 1(b) makes binding and
  the instrument's own nominal operating point. The study in its pre-registered form is not
  executable there, recorded now, quoting Amendment 3, as "a live outcome of this study, not a
  failure mode to be designed around". Declaring the surviving subset before the anchor exists is
  what removes the freedom to pick a convenient one afterwards. The verdict column is
  `below_spec_dataset_floor`, computed from the manifest and parsed back out of the document by the
  tests, never typed.
- Added: **what the manifest could have supported instead, stated before the anchor exists.** Of
  the 68 candidate-bearing datasets, **62 / 33 / 21 / 12** hold a stratum at ≥ 4v4 / ≥ 8v8 /
  ≥ 13v13 / ≥ 23v23, and only **3** hold an exactly-3v3 stratum — none of which clears 13v13, so
  §1 (iii)'s mandatory 3v3 anchor costs a dataset slot at the hard tiers. Clustered the way D2
  requires, those 68 datasets sit in **50** collections and the tiers hold **46 / 25 / 15 / 10**;
  the frozen list retains **11 / 7 / 5 / 3**. **The truncation is therefore a consequence of
  selecting for §1 (iii)'s coverage axes rather than for donor counts, not a property of the public
  data at this Census pin** — worse for the study on feasibility, better on integrity, and
  impossible to reassemble the list on donor counts later and call it the original plan. An earlier
  draft argued this by constructing a rival twelve; the construction double-counted two collections
  and leaned on the one dataset that owns all 12 `excluded_confound` rows, so it was deleted in
  favour of the counts above, which are recomputed on every freeze.
- Added: **every same-collection sibling is frozen, and the set is computed rather than typed.** An
  earlier draft named two; there are **five** — SEA-AD MTG (18 candidate strata), Yoshida Airway
  (25), CAREBANK (11), KPMP scRNA v1.5 (43) and Emphysema immune (9). All **106** are emitted under
  `within_collection_control_rows` with `role = "within_collection_control"`, where the 251 carry
  `role = "analysis_set"`; the CSV twin holds both blocks, **357** data rows, told apart by that
  column. Collection membership comes from the pinned release table, so the set is recomputed on
  every freeze and an unnamed sibling aborts it — the manifest carries no collection column, which
  is exactly why three went unnoticed. A result from a control is reported as a within-collection
  control, never enters the D2 denominator, and promoting one to an independent dataset is an
  amendment. Naming two of five had left **79 runnable, unlisted strata** selectable after the
  fact, two of whose datasets clear every envelope tier.
- Added: `pilot/preregistration/stratum_list_proposal_2026-08-16.redacted.md`, the proposal the
  choice of twelve was made from — committed because §3.2, §4.2 and §8 cited it for the per-dataset
  rationale, the rejected datasets, **five named reserves**, the third 5′ candidate and the Mathys
  search method, and it existed only in a scratch directory: the one act of discretion in the freeze
  was justified by a file no reader could open. It is the reasoning, **not** part of the binding
  act; where it disagrees with the document, the document governs, and every established
  discrepancy is enumerated in §10. Both hashes are recorded in §2 — the copy as circulated and the
  redacted copy as committed — and the redaction is confined to six occurrences of one absolute
  filesystem path, replaced by `<REPO>` and enforced by a test. The reserves matter most: §9 item 7
  forbids substituting a replacement for a stratum that fails the counts gate, and pinning the list
  makes such a substitution detectable.
- Added: two external indices are **pinned and committed** rather than read live —
  `discover_index_2026-08-16.json` (2216 records) and `census_release_datasets_2025-01-30.json`
  (the release's own 1573-row dataset table). Collection membership and the assay / suspension /
  tissue / DOI table are recomputed from these, so claims that were previously uncheckable — every
  defect of the first three review rounds lived among them — now abort the freeze on disagreement.
  Their disagreement with each other is measured and reported rather than assumed away: 1567 of
  1573 release datasets resolve in Discover, `dataset_version_id` matches for **0**, and
  `collection_doi` differs for 61.
- Fixed (a circulated figure, corrected rather than absorbed): the 2026-08-16 proposal states
  Rexach's envelope ceiling as `min(A,B) = 11`. **It is 10** — the dataset's best design is A = 11
  versus B = 10 and no Rexach control group exceeds 10 donors, so the 11 was `max(n_donors_A)` read
  as the ceiling. Nothing downstream moves (Rexach was outside the σ = 0.5 tier either way), and
  the error, its cause and its magnitude are on the record.
- Note: the freeze **admits nothing**. All four blockers stand on all 2190 source rows and all 251
  frozen ones (`integer_check`, `frozen_universe_size`, `sigma_donor_estimate`,
  `envelope_membership`), `pooled` is `unresolved` on 1197/1197 candidates so the
  donor-pseudobulk-is-calibrated claim is unavailable on this entire list (a property of the Census
  pin, not of the selection), and the counts gate may still shrink the list — a shrinkage is a
  reported outcome, never a re-selection. §8(d)'s Mathys 2019 anchor is recorded in its own section:
  absent from CELLxGENE Discover entirely (full 2216-dataset index searched 2026-08-16, zero hits
  on six needles over every field of every record) **and absent from the pinned release itself** —
  the 1573 datasets published under `cell-census/2025-01-30/h5ads/` were enumerated directly,
  because "absent from Discover ⇒ absent from Census" does not follow: **6 of those 1573 are no
  longer listed in Discover**, and each was identified from its own h5ad (two unrelated
  collections) before the conclusion was allowed to stand. `cellxgene-census` could not be
  installed to read `census_info/datasets` through SOMA: `tiledbsoma` publishes no Windows wheel at
  any version. So the Census path is closed; the Synapse path needs a ROSMAP DUA — applied for on
  2026-08-16, **not granted at the time of writing** — and a second loader that §9's plan does not
  contain. It does not block the freeze and the freeze does not weaken it.
- Added: `tests/test_stratum_list_freeze.py` (77 tests, no network) — the source hashes and the
  proposal's; the twelve ids and the 251 strata re-derived from the raw manifest rather than
  through the module being tested, in both directions (nothing extra, nothing quietly omitted);
  Rexach's ceiling of 10; every bin occupied and the bins tiling `[10, inf)` from the gate's own
  floor; **all four** Layer B tiers equal to the envelope arithmetic with their
  `below_spec_dataset_floor` verdicts re-derived; §6's two tables parsed back out of the document
  and compared cell by cell, so a hand-typed verdict cannot disagree with the artifact; the
  manifest-wide 62 / 33 / 21 / 12 and the counterfactual's 12 / 12 / 11 / 11 with its explicit
  witness list; the 27 control strata, their role field and their disjointness from the 251; no row
  admitted; neither sibling in the analysis set; byte-identical regeneration and both committed
  halves equal to a fresh run. Ten tests drive the guards themselves — a wrong stratum count, a
  wrong Layer B subset, a **dropped** Layer B tier, a wrong manifest tier census, a wrong
  counterfactual maximum, a wrong sibling count, a sibling promoted into the D2 denominator, an
  unlabelled disease arm, and a missing or tampered source CSV must each abort the freeze — so the
  verification cannot become decorative the way `gate_config`'s history warns.
- Fixed: `load_source` **skipped the CSV hash check when the file was absent**
  (`csv_path is not None and csv_path.exists()`), so deleting the CSV turned a failing hash check
  into a passing run. A named-and-absent evidence file now aborts, as does one of the wrong size.
  The same guard covers the newly pinned proposal document. `--check` now compares **both**
  committed halves rather than the JSON alone.
- Added: `.gitattributes` marking `pilot/preregistration/` as `-text`. The committed CSV came off a
  Linux runner with CRLF terminators and `core.autocrlf` would rewrite the JSON on a Windows
  checkout; either conversion moves a recorded hash and would make the byte-identity test pass on
  one CI leg and fail on another. `.pre-commit-config.yaml` excludes the same directory from
  `check-added-large-files` rather than raising the 500 KB cap for the whole tree.
  `pilot/gate/` and `pilot/testsel/` are deliberately **not** marked, and `.gitattributes` now says
  why: they were committed with LF and are checked out with CRLF, so `-text` would make the working
  tree canonical and the next `git add` would rewrite all three blobs (measured: 8909 → 9174,
  21811 → 21958 and 63081 → 65856 bytes). Protecting evidence by rewriting it is not protecting it.
- Fixed (found by running the repo's own hooks while adding the above): `pre-commit run
  --all-files` **rewrote `pilot/gate/synthetic_gate_2026-08-15.json`.** `json.dump` writes no
  trailing newline, `end-of-file-fixer` appended one, and that artifact is cited evidence — its
  numbers appear in Amendment 3's addendum, the README status table and this changelog. The hooks
  are configured but not installed as a git hook, so nothing had ever run them over it. The
  fixing hooks now exclude `pilot/gate/`, `pilot/testsel/` and `pilot/preregistration/` alongside
  the two frozen protocol documents, on the same principle: a mechanical pass must not be able to
  edit evidence. The modification was reverted; no artifact byte changed.

### Amendment 4 Part A — the `sigma_donor` estimator is specified and its criteria pre-declared (2026-08-16)

- Added: **Amendment 4, Part A** (`docs/AMENDMENTS.md`). Amendment 3 declared a per-stratum
  `sigma_donor` estimate to be required work and did not do it; this entry specifies the estimator,
  the operating-envelope membership rule, and — in the order this log exists to enforce — commits the
  numeric criteria **before** the run that will be judged by them. Part B, after the confirmatory
  grid, records the outcome and fixes the aggregation functional by the mechanical rule declared
  here. The gate functional is the one deliberately open slot, and it is a required argument with no
  default in code rather than a TODO.
- Fixed (**Correction 1 to Amendment 3 Change 1**): that entry states `sqrt(s0^2) * ln 2` is an
  **upper bound** on `donor_sigma`. It is not. The `+1` of `log2(CPM + 1)` attenuates the donor
  random effect along with everything else inside the logarithm, so the quantity is
  `a_bar * sqrt(sigma^2 + v)` with two distortions of opposite sign, and it **understates**
  `donor_sigma` whenever `v < sigma^2 (1/a_bar^2 - 1)`. The failure is not at low expression, as
  Amendment 3's wording suggests, but on clean, deep data in a wide gene universe: at the gate's own
  `ORACLE_SIM` point the quantity falls below the truth in half of realisations, and at a
  15 000-gene universe it reads 0.332 against a true 0.350 in **64 of 64** seeds. Understating sigma
  admits strata the pseudobulk arm is not valid for, which is the dangerous direction. The quantity
  is demoted to an audit column of unknown error sign and gates nothing.
- Added: `scripts/check_upper_bound_claim.py` and `pilot/upper_bound_check/`, so Correction 1 rests
  on a committed artifact with named seeds rather than on prose — the defect Amendment 3 charged its
  own predecessor with.
- Fixed: `census_select.PENDING_FIELDS["sigma_donor_estimate"]` and the module docstring, which
  repeated the retracted claim. Live code is corrected in place; the published dated snapshot
  (`docs/PILOT_FINDINGS.md`) is corrected by this log, and the unpublished one
  (`docs/PREREGISTRATION_STRATUM_LIST.md`) carries the retraction at the point of the claim because
  it ships in the same push.
- Fixed (performance, shipped arm): `ebayes_from_pdata` rebuilt `set(pdata.var_names)` **inside** the
  universe comprehension, making the restriction O(G²). Measured at G = 15 000: **127.3 s → 0.008 s**,
  a 16 400× speedup with a bit-identical result. It is called once per permutation, so at spec §4's
  200 paired permutations over 251 strata it would have cost roughly 74 days of the real sweep. It
  was invisible because the gate has only ever run at the simulator's 1500 genes. The same
  construction survives in the superseded DESeq2 arm and in the frozen calibration probe, both
  deliberately untouched and both off this arm's path.

### The permutation null becomes affordable, and the floor stops mixing conventions (2026-08-16)

- Fixed (latent, would have corrupted the real sweep's headline): `run_null` returned one
  `naive_ndeg` array whose meaning changed at index `n_paired` — paired-BH below it, the naive
  arm's own solo BH above. The two agreed only because `N_PERM == N_PERM_PB == 200`; at spec §4's
  pre-registered counts for the real sweep (`n_perm` 1000, `n_perm_pb` >= 200) four fifths of the
  array would have come from the other convention and the headline floor would have been a median
  over both at once. The mixed array no longer exists: `naive_ndeg_paired` and `naive_ndeg_solo`
  are separate, the key `naive_ndeg` is **removed rather than deprecated** so a stale caller gets a
  `KeyError` instead of a number, and `metrics.NDegSeries` carries the BH convention with the
  counts (read-only, so the other convention cannot be appended). `perm_floor` propagates the label
  into the artifact and refuses a contradicting one.
- Added: `src/pbcheck/methods/naive_engine.py` — the naive arm's donor-permutation null computed
  from **per-donor sufficient statistics** instead of re-running scanpy per permutation. Under
  donor relabeling the pooled cell set is unchanged, so per-cell normalisation, the per-gene ranks
  and the tie correction are all invariant, and the Wilcoxon rank sum of the test group is the sum
  of its donors' rank sums; `t-test_overestim_var` (spec §2's robustness variant) collapses the
  same way onto `n`, `Σx`, `Σx²`. One ranking pass per stratum, then each permutation is an
  addition over at most `n_donors` rows.
- **This is an optimisation of a pre-registered statistic, not a new statistic, and it is held to
  that standard.** Agreement with scanpy 1.12.2 is bit-exact — max ulp 0 on p-values, log2fc and
  the expressed fractions, across donor counts, group sizes, tie structures, unequal cells per
  donor, single-donor-dominated strata, and sparse against dense. The scanpy kernels are vendored
  in numpy rather than imported privately, with a test pinning the two bit-for-bit so a change
  upstream turns the suite red instead of moving our numbers. The slow path stays callable
  (`naive_engine="scanpy"`) and is the reference the tests compare against.
- Measured: the synthetic gate runs in **7.1 s against 511.8 s**, and its artifact is unchanged —
  all 187 committed leaf values reproduce exactly, nothing removed, the 16 additions being the new
  engine's provenance, the BH-convention labels and the solo floor. Per-permutation speedups of
  **765×** (60 000 cells × 1500 genes) and **1117×** (150 000 × 2000) were measured directly;
  extrapolating the setup constant to the frozen list's largest stratum (547 665 × 15 000) gives
  1151×. The naive arm's share of the real sweep at `n_perm` = 1000 falls from thousands of
  core-hours to hours.
- Fixed: the same O(G²) universe restriction repaired in the moderated arm earlier today
  (`set(pdata.var_names)` rebuilt inside a comprehension) also stood in
  `methods/pseudobulk.py`, the superseded DESeq2 arm. Hoisted, provably output-identical. The copy
  in `scripts/pb_calibration_probe.py` is deliberately untouched: it is the frozen reference
  instrument the Amendment 2 selection grid was measured with.
- Known limitation, recorded rather than worked around: the engine's setup pass converts CSR to
  CSC and peaks at two copies of the sparse matrix, which for the largest stratum is on the order
  of 20 GB. Blockwise assembly would fix it; it has not been done, because doing it without a
  profile of a real stratum would be guesswork.
- Fixed (found by CI, not locally): the engine's cross-check against scanpy failed on every CI leg
  for `t-test_overestim_var` while passing here. The engine was not at fault — it reproduces scanpy
  bit-for-bit on a given machine, and 7200 local configurations differ by exactly 0.0. **scanpy's
  own value moves between machines.** `fast_array_utils.stats.mean_var` forms the variance as
  `E[x²] − E[x]²` with `power(x, 2)` called without a `dtype`, so the squares are taken in float32
  storage precision; the relative error of the variance goes as `(eps32/2 + n·eps64)·E[x²]/var`. On
  log1p-normalised data that condition number is enormous by construction — measured 1.9e14 on an
  adversarial gene, where the formula returns 4.79e-06 against a true variance of 4.55e-13, i.e.
  **zero significant digits** — and it does not improve with more cells (1.9e14 at 2 cells per
  group, 9.9e13 at 120). So `t-test_overestim_var`, spec §2's robustness variant, is not
  reproducible across machines **for any implementation**, ours included. This does not touch a
  single number the project reports: neither `run_null` nor the gate ever passes `method`, so both
  run the pre-registered Wilcoxon, whose agreement stays bitwise on every leg. It is escalated to
  the amendment log rather than absorbed here.
- Changed (the fix, which is not a wider tolerance): the `rtol=1e-12` on p-values is **removed and
  not replaced by a larger constant**. The assertions move to layers where the arithmetic is
  well conditioned — the per-group sufficient statistics against `math.fsum` at `n·eps64`, and a
  bitwise comparison against a line-by-line transcription of scanpy's `t_test` given identical
  statistics, on deliberately unbalanced groups so the `nobs2 = ns_group` substitution is
  observable. The end-to-end comparison keeps a per-gene rounding envelope **derived** from the
  dtype analysis rather than a constant. A six-mutation table in the report shows each layer
  catches a defect the others miss (shifted statistics, dropped ddof, removed overestim
  substitution, Student for Welch, double-counted donor row, removed tie correction).
- Changed: `derandomize=True` on the property suites, with example budgets raised 30 → 120 and
  50 → 100. Hypothesis reseeds from entropy per run and `.hypothesis/examples` is gitignored, so
  every machine explored a different set and a green local run meant nothing about CI. Now a CI
  failure reproduces here by test name. Stated plainly because it is the real lesson: this defect
  was invisible locally through 7200 configurations and CI found it on the first attempt.

### Amendment 5 Part A — the §6 reportability split (2026-08-16)

- Added: **Amendment 5, Part A** (`docs/AMENDMENTS.md`), written before the `sigma_donor` anchor
  exists. It **NARROWS Amendment 3 Change 1's blanket prohibition** on reporting anything from a
  stratum outside the operating envelope, and says in those words that narrowing a prohibition is a
  relaxation. `lambda_naive` and the naive arm's permutation false-positive floor become reportable
  at `min(n_A, n_B) >= 8` on **measured per-stratum calibration**; `real_label_ratio`, `concordance`
  and any claim that a published finding is false stay inside the envelope, unchanged. The ground is
  that the envelope gates on donor count to guarantee the pseudobulk arm is *powered*, while the two
  map quantities are naive-arm properties against the donor-permutation null: injecting NaNs into
  the pseudobulk arm leaves `lambda_naive` at 26.160, and over a 5.75× range of donor counts it
  moves −4.7 % against +620 % for cells-per-donor and +2683 % for σ. No threshold moves: α, the λ
  band, POWER_TARGET, the oracle's effect size, the envelope table and decision rule item 1 are
  untouched. The rule is **stricter** than the envelope at σ ≈ 0.2 (150 strata / 7 datasets against
  227 / 11) and looser at 0.5–0.7.
- Added, in the same entry: A2 cell-count stratification is **un-deferred at ≥ 8v8 as a blocking
  condition** — the map is made of floors and the floor scales with cell count, and the measured
  confound is large (within-group cells-per-donor max/min median 18.0, max 524.5; 40 % of candidate
  ≥ 8v8 strata above 3 on the between-group total-cell ratio). Amendment 2 Change 6's infeasibility
  argument was about 3v3 and inverts above 8v8, where the smallest tier holds 12 870 assignments.
- Added: eight pre-declared abandonment triggers with numbers, including the one that stops
  publication outright if the A2-stratified floor departs from the unstratified floor by more than
  2 Monte-Carlo SE on more than 20 % of strata, and the one that withdraws the rule entirely if the
  anchor lands at σ ≤ 0.2, where it is pure narrowing.
- Disclosed (and it is uncomfortable): **Change 3's own admission criterion voids calibrated strata
  at roughly an order of magnitude above its declared rate.** The declared 0.582 % assumed the
  permutations of one stratum are independent, which is measurably false. Two independent probes at
  σ = 0.2, 8v8, `n_perm` = 200 read 15/100 and 3/40 — 25.8× and 12.9× the nominal, consistent under
  sampling but not settled to better than a factor of two — and 5/100 at σ = 0.5, where the *fresh*
  null reads nominal, so nothing transfers between the two nulls. Trigger 8 is set at 10 %, which
  the present evidence straddles from both sides rather than clears.
- Disclosed: the arm's fresh-null rejection rate at σ = 0.2 sits above nominal across three
  independent seed blocks (56/1000, 27/300, 70/1000; pooled 153/2300 = 0.0665, exact p = 2.9e-04),
  homogeneous across blocks (χ² = 4.64, p = 0.098) and separating across σ (χ² = 9.90, p = 0.0071).
  A draft of the entry called this a failure to reproduce; that was a category error — a p of 0.21
  is not evidence of absence — and the sentence is withdrawn on the record. σ = 0.2 is the one
  envelope row with no grid support at all, so an arm-side excess sitting exactly there earns more
  disclosure, not less.

### The repository starts carrying its own map and reviewers (2026-08-16)

- Added: `.claude/agents/` — the five profile agents (`repo-scout`, `mechanic`, `terminal-runner`,
  `deep-diagnostician`, `reviewer`), each pinning its own model and tool set. They live **in the
  repository** rather than in a machine's user directory, because that is the only way a cloud
  session gets them; until now every role had to be improvised from a generic agent with a briefing
  pasted into the prompt.
- Added: `.claude/skills/pbcheck-map/SKILL.md` — the map a session reads instead of walking the tree
  blind: what the project is, what each module, script and document owns, what is frozen and why
  editing it is a protocol violation rather than a style question, the amendment discipline, the
  current state with numbers re-derived from the artifacts, and the traps. Every path in it was
  confirmed to exist and every number recomputed (367 tests at `c3c2556`, 251 analysis-set strata
  over 12 datasets plus 106 controls, `admitted_to_sweep` false on all 357). It references code by
  module, function and document section, never by line number, so it cannot rot into a document that
  lies without looking broken.
- Noted in that map rather than silently repaired: `pilot/gate/synthetic_gate_2026-08-15.json` still
  reports the 511.8 s runtime it was committed with, which predates the engine work of the following
  day. Its correctness numbers reproduce exactly; only the timing is historical, and regenerating
  cited evidence to refresh a stopwatch would be the wrong trade.

### Evidence for the cell-count stratification's tolerance (2026-08-16)

- Added: `scripts/a2_feasibility.py` and `pilot/preregistration/a2_feasibility_2026-08-16.{json,csv}`,
  the evidence a forthcoming amendment will cite when it fixes the matching tolerance for cell-count
  stratification of the permutation null. Both source artifacts are sha256-pinned before parsing, the
  output carries no timestamp or environment, and `--check` verifies byte-identical regeneration.
- **The committed data cannot answer the question exactly, and the artifact says so rather than
  papering over it.** A stratum row records only `n_donors`, `n_cells`, `min`, `median` and `max`
  cells per donor per group, never the per-donor vector, so the exact distribution of permuted
  cell totals is not computable from what is committed. The artifact therefore brackets it between
  two extreme per-donor vectors consistent with those five statistics, documents both constructions
  in its header, and records the median drift each construction incurs when clamped to the recorded
  total — for donor-rich strata that drift reaches thousands of cells and is reported per stratum
  rather than hidden. The binding count remains a load-time fact.
- Below the enumeration cap the assignment totals are computed by exact dynamic programming over all
  designs, not sampled; above it by seeded Monte Carlo whose seed is derived from the stratum's own
  identifiers so it does not depend on iteration order, with the standard error carried in every
  cell. Means, variances and the real split's z are computed by the exact finite-population formula
  in both branches. The dynamic program and the moment formulas were checked against brute-force
  enumeration on a real and a synthetic stratum.
- **A re-derivation disagreed with the planning document that motivated it, and the disagreement is
  recorded rather than reconciled.** At a tolerance of 5 % of total cells the pessimistic edge
  projects 11 of the frozen 150 strata unusable, not 5, and the loss is not confined to one dataset;
  the optimistic edge agrees at 3, all in the same dataset. The consequence is substantive: the
  planning document argued for 5 % from a knee — 2.5 % → 5 % saving five strata and 5 % → 7.5 %
  saving one — and on the re-derived ladder the pessimistic edge runs 21 → 11 → 5 and keeps
  improving past 5 %. The tolerance argument has to be rebuilt on these numbers, and no tolerance is
  fixed by this commit.
