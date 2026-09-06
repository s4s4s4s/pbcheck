# pbcheck usage

This is the full reference for the `pbcheck` command and the `pbcheck.audit` Python API: every
flag, every exit code, the three statuses an audit can end in, the fields of the JSON output, and
what each number in the report means. For a fast start see the Quickstart in the project
[README](../README.md).

## Install

```bash
pip install pbcheck
```

Python 3.12 or newer is required. `pip install pbcheck[report]` is a separate, unrelated extra
(it pulls in the dependencies for a not-yet-built HTML report pipeline under `pbcheck.report`);
it is not needed to run `pbcheck audit`, which already writes its own Markdown and HTML output
with no extra install.

## CLI reference

```
pbcheck --version

pbcheck example OUT.h5ad [--seed N] [--shape {small,reference}]
    # --seed default 0; --shape default small
    # small: 600 genes, 4 donors per condition, 60 cells per donor
    # reference: 8000 genes, 8 donors per condition, 625 cells per donor

pbcheck audit FILE.h5ad --donor COL --condition COL --test LEVEL --ref LEVEL
    [--celltype COL --celltype-value VALUE]
    [--batch COL ...]              # repeatable, one or more batch columns
    [--counts-layer NAME]
    [--n-perm N]                   # default 1000 (pbcheck.audit.PRODUCT_N_PERM, a product value)
    [--n-perm-pb N]                # default 200 (pbcheck.audit.PRODUCT_N_PERM_PB, a product value)
    [--seed N]                     # default 0
    [--alpha F]                    # default 0.05 (gate_config.ALPHA, a pre-registered constant)
    [--top N]                      # default 25
    [--design-only]
    [--out DIR]                    # default ./pbcheck_out/<file stem>_<celltype-value or all>
    [--format {json,md,html,all}]  # default all
    [--overwrite]
    [--quiet]
```

`pbcheck example` writes a small synthetic `.h5ad` file with donor structure, so the quickstart
needs no download and no real data.

`pbcheck audit` takes one `.h5ad` file, the names of the donor and condition columns in `.obs`,
and the two condition levels to compare (`--test` is the level of interest, `--ref` is the
level it is compared against). It writes three files into the output directory:

- `pbcheck_audit.json`, machine-readable, schema `pbcheck-audit/1` (see below).
- `pbcheck_report.md`, the same content as human-readable Markdown.
- `pbcheck_report.html`, the same content again as a standalone HTML page (no external assets,
  no JavaScript).

`--celltype COL --celltype-value VALUE` restricts the audit to one cell type; both flags must be
given together or neither. `--batch COL` may be repeated to name one or more columns pbcheck
checks for confounding with condition. `--counts-layer NAME` names the `.layers` entry pbcheck
should treat as the raw count matrix; without it, pbcheck looks in `.X`, then `.layers["counts"]`,
then a raw layer, then `.raw.X`, and uses the first one that passes the raw-integer-count check.
`--design-only` skips the differential-expression steps entirely and writes a report with only
the metadata audit (see "Metadata-only mode" below).

### Python API

The same audit is available without the CLI:

```python
import pbcheck

adata = pbcheck.example_adata()
settings = pbcheck.AuditSettings(
    donor_col="donor", condition_col="condition", test_level="disease", ref_level="ctrl",
)
payload = pbcheck.run_audit(adata, settings)
```

`pbcheck.audit_h5ad(path, settings)` does the same starting from a file path.
`pbcheck.render.render_markdown`, `render_html` and `write_outputs` turn a payload into the same
three files the CLI writes. The keyword arguments of `AuditSettings` mirror the CLI flags above.

## Exit codes

| Code | Meaning |
|---|---|
| `0` | Outputs were written. This includes the `naive_only` and `design_only` statuses below: pbcheck stopping early because of what it found in your data is a successful run, not an error. |
| `2` | A usage or input problem: a bad command-line argument, a missing file, a missing or misspelled column name, `--celltype` given without `--celltype-value` or the other way round. A one-line message on stderr says what was wrong. |
| `1` | Any other, unexpected failure. A full traceback is printed; this is a bug report, not a diagnosis of your data. |

## The three statuses an audit can end in

Every audit ends in exactly one `status`, and `status_reason` says which check produced it. In
plain words, in the order pbcheck checks them:

| Status | `status_reason` | When | What it means |
|---|---|---|---|
| `design_only` | `design_only_requested` | `--design-only` was given | You asked pbcheck to only look at the metadata, so it did: donor structure, condition balance, batch confounding. No expression data was read. |
| `design_only` | `donor_spans_conditions` | a donor appears under both condition levels | This is a paired design (see "Paired designs" below). pbcheck will not run either differential-expression arm on it. |
| `design_only` | `too_few_donors` | fewer than 3 donors in one of the two groups | There are too few donors to build a permutation null that means anything, so pbcheck stops before touching expression data. |
| `naive_only` | `non_integer_counts` | the count matrix fails the raw-integer-count check | The pseudobulk arm needs raw counts and none were found; it is skipped rather than run on transformed numbers it cannot interpret correctly. The naive (per-cell) arm still runs, on the matrix as found. |
| `design_only` | `universe_too_small` | the shared gene universe is too small after filtering | Not enough genes survive the frozen-universe filter (or, on the fallback path for non-integer counts, its own filter) to say anything useful. |
| `naive_only` | `too_few_profiles_after_thin_filter` | too few donor-pseudobulk profiles remain in a group after aggregation | The naive arm still ran, but there are not enough donor profiles per group for the pseudobulk arm to be meaningful. |
| `complete` | `null` | none of the above fired | Both arms ran: the naive per-cell test and the donor-pseudobulk test, each against its own donor-permutation null. |

`status_reason` is `null` only for `complete`; every other status carries exactly one of the six
reason strings above, which is the whole of `pbcheck.audit_schema.STATUS_REASON_VALUES`.

`naive_only` and `design_only` are not failures; they are pbcheck refusing to compute a number
that would not mean what it looks like it means for that file.

## Metadata-only mode (`--design-only`)

Runs only the stratum preparation and the design audit; the report is written with the same three
files, but its differential-expression sections are marked "not run". This is the fast, counts-free
way to check whether a design is even eligible for the rest of the audit before spending time on a
full run.

## Paired designs

pbcheck does not support paired or repeated-measures designs, where the same donor is measured
under both conditions being compared (for example, before and after treatment in the same person).
There is no `--subject` covariate and no mixed-model arm. Concretely:

> if your donors were measured under both conditions, pbcheck stops at the design audit; do not
> work around it by combining donor and condition into one column - that makes two halves of one
> person look like two people, which is the error this tool audits

If your data has this shape, pbcheck reports `status = design_only`, `status_reason =
donor_spans_conditions`, and explains in the report why it stopped there.

## What each number means

- **Genomic inflation factor (lambda, `lambda_naive_class` / `lambda_pseudobulk_class`).**
  A summary of how much more (or less) significant a test's real results look than a well-behaved
  test would produce, as its own separate quantity, not derived from the floor or from the
  real-over-floor ratio. `in_band` means lambda falls inside the tool's accepted band; `above_band`
  means the test is calling more than it should; `below_band` means fewer. `null` means that arm
  did not run. The band itself (`gate_config.LAMBDA_BAND`) is the donor-pseudobulk arm's band;
  where the naive arm's class is reported, the report shows the same band next to it as a
  description only, not as the naive arm's own criterion.
- **The permutation floor (`naive_floor_solo`, and the pseudobulk floor).** pbcheck shuffles
  condition labels between donors many times, so that by construction there is no real signal to
  find, and counts how many genes each test still calls significant. That count, at a chosen
  quantile across the shuffles, is the floor: how many "hits" a test produces on pure noise, under
  this data's own donor structure. A test's real-label result is only informative once it is read
  against its own floor, not in isolation.
- **The negative control (the pseudobulk arm's floor and false-positive rate).** The donor-pseudobulk
  test is expected to keep its own floor near zero and its permutation false-positive rate near the
  chosen significance level. That behaviour is the check that the pseudobulk arm's own error rate
  stays near that level on this file; it does not certify the naive arm, and it does not certify any published
  result.
- **The envelope.** The pseudobulk arm's calibration and power were established on synthetic data
  only, and only inside a stated range of donor counts and per-donor variability. Quoting the exact
  wording pbcheck's own report uses for this:

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

  In short: pbcheck has no way to tell you whether your specific file sits inside or outside the
  range it was validated on. It states the range and stops there.
- **No `risk_score` exists and none is promised.** pbcheck does not compute or output a single
  combined score of how trustworthy a result is. It reports the measurements above and nothing
  that summarises them into a verdict.

## JSON output schema (`pbcheck-audit/1`)

`pbcheck_audit.json` validates against `pbcheck.audit_schema.validate`. Top-level keys:

| Key | Holds |
|---|---|
| `schema_version` | Always `"pbcheck-audit/1"`. |
| `pbcheck_version`, `generated_utc` | The tool version and the UTC timestamp of the run. |
| `runtime_seconds`, `runtime_by_stage_seconds` | Total wall time and a per-stage breakdown (load, prepare, design, counts, pseudobulk build, universe, real-label arms, permutation null, render). |
| `status`, `status_reason` | One of the three statuses above and the reason string from the table above. |
| `input` | What was loaded and dropped: cell and gene counts, cells dropped for each reason (wrong condition, wrong cell type, missing donor/condition/cell type), the column and level names given, and where the count matrix was found. |
| `counts_check` | The result of the raw-integer-count check: pass/fail, dtype, how many values were checked, how many were non-integer, negative or non-finite. `null` when counts were never read (`design_only`). |
| `design` | The full metadata audit: donors per group, nesting of donor within condition, cell counts per donor, imbalance, batch-confounding measurements, and the flags derived from them. |
| `settings` | Every setting the run used, split into the tool's own defaults (`settings.tool`) and the small set of values that come from the pre-registered protocol (`settings.protocol_constants`, taken from `pbcheck.gate_config`). |
| `universe` | The size and origin of the shared gene universe both arms were tested over, and the rule used to build it. |
| `real_label` | The naive and pseudobulk arms' results on the real condition labels: how many genes each called significant, the top genes, and the accounting of how the two arms' gene lists were made comparable. `null` if neither arm ran. |
| `permutation_null` | The donor-permutation null results for each arm that ran: inflation, floor, and (for the pseudobulk arm) the false-positive rate, plus how many permutations were requested versus actually distinct. `null` in `design_only`. |
| `readout` | The plain-language read-out derived from the fields above: inflation classes, floors, ratios of real result to floor, and the sentences shown in the report. |
| `caveats` | The list of caveat texts shown in this report, each with a stable identifier (see the "What each number means" section and the report itself for their wording). |
| `provenance` | Which recorded synthetic measurement pbcheck's settings were fixed against, which settings are pre-registered protocol constants, and the platform and package versions the run used. |

## Settings (`settings.tool`)

Three of these are product constants that are not part of the pre-registered protocol constant
list (`pbcheck.gate_config`); each equals the value the measurement engine used as its own default
before this release existed, so shipping them changes no number, but they are the tool's own
values and can change in an ordinary engineering release without an amendment.

| Key | Product constant | Meaning |
|---|---|---|
| `min_donors_per_group` | `pbcheck.audit.PRODUCT_MIN_DONORS_PER_GROUP` | Minimum donors per group below which the run stops at `status = design_only`, `status_reason = too_few_donors`. |
| `universe_min_total_count` | `pbcheck.audit.PRODUCT_UNIVERSE_MIN_TOTAL_COUNT` | Minimum total count across pseudobulk profiles for a gene to enter the frozen universe. |
| `universe_min_prop` | `pbcheck.audit.PRODUCT_UNIVERSE_MIN_PROP` | Minimum proportion of profiles a gene must be detected in to enter the frozen universe. |

The rest of `settings.tool` (`n_perm_requested`, `n_perm_pb_requested`, `seed`, `top_n`,
`design_only`, `naive_method`, `naive_engine`, `pseudobulk_method`, `trend`,
`fallback_universe_min_prop`, `fallback_universe_min_size`,
`min_profiles_per_group_after_thin_filter`) mirrors the CLI flags and the fallback rule used when
the count matrix fails the raw-integer check (see "Metadata-only mode" and the status table above).

## Runtime

Machine: Windows-11-10.0.26200-SP0, Intel64 Family 6 Model 154 Stepping 3, GenuineIntel, 20 CPUs.
One machine's measurement (`scripts/measure_audit_runtime.py`), not a guarantee for any other
machine.

| shape | n_perm | n_perm_pb | wall seconds | stage timings (s) | peak working set | status |
|---|---|---|---|---|---|---|
| reference | 1000 | 200 | 228.6 | counts=1.7, design=0.0, permutation_null=99.5, prepare=0.7, pseudobulk_build=4.6, real_label=47.6, universe=0.0 | 2036 MiB | complete |
| reference | 200 | 200 | 192.7 | counts=1.7, design=0.0, permutation_null=77.4, prepare=0.7, pseudobulk_build=4.3, real_label=57.6, universe=0.0 | 2092 MiB | complete |
| example | 1000 | 200 | 14.4 | counts=0.0, design=0.0, permutation_null=2.1, prepare=0.0, pseudobulk_build=0.1, real_label=1.0, universe=0.0 | 301 MiB | complete |
| example | 200 | 200 | 12.5 | counts=0.0, design=0.0, permutation_null=1.9, prepare=0.0, pseudobulk_build=0.1, real_label=0.9, universe=0.0 | 301 MiB | complete |
| gate | 1000 | 200 | 40.0 | counts=0.1, design=0.0, permutation_null=18.7, prepare=0.0, pseudobulk_build=0.3, real_label=3.2, universe=0.0 | 513 MiB | complete |
| gate | 200 | 200 | 34.6 | counts=0.1, design=0.0, permutation_null=14.3, prepare=0.1, pseudobulk_build=0.4, real_label=4.1, universe=0.0 | 512 MiB | complete |

## Product defaults

`--n-perm` and `--n-perm-pb` default to `pbcheck.audit.PRODUCT_N_PERM` (1000) and
`pbcheck.audit.PRODUCT_N_PERM_PB` (200), the product's own permutation counts, not the
pre-registered protocol's counts (`gate_config.N_PERM`, `gate_config.N_PERM_PB`, which pbcheck's
product code never reads). `--alpha` defaults to `gate_config.ALPHA` (0.05), a pre-registered
protocol constant reused as-is because pbcheck does not define its own significance threshold.
`pbcheck example` defaults to `--seed 0 --shape small`. The runtime table above is what these
defaults cost in wall time on one measured machine.
