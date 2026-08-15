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
