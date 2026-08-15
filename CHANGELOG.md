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
