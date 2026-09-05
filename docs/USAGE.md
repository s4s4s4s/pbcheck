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
    # small: 600 genes, 4 donors per condition, 60 cells per donor
    # reference: 8000 genes, 8 donors per condition, 625 cells per donor

pbcheck audit FILE.h5ad --donor COL --condition COL --test LEVEL --ref LEVEL
    [--celltype COL --celltype-value VALUE]
    [--batch COL ...]              # repeatable, one or more batch columns
    [--counts-layer NAME]
    [--n-perm N]                   # default: the product defaults (see the runtime table)
    [--n-perm-pb N]                # default: the product defaults (see the runtime table)
    [--seed N]                     # default 0
    [--alpha F]                    # default: the tool's own significance threshold
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

| Status | When | What it means |
|---|---|---|
| `design_only` | `--design-only` was given | You asked pbcheck to only look at the metadata, so it did: donor structure, condition balance, batch confounding. No expression data was read. |
| `design_only` | a donor appears under both condition levels | This is a paired design (see "Paired designs" below). pbcheck will not run either differential-expression arm on it. |
| `design_only` | fewer than 3 donors in one of the two groups | There are too few donors to build a permutation null that means anything, so pbcheck stops before touching expression data. |
| `naive_only` | the count matrix fails the raw-integer-count check | The pseudobulk arm needs raw counts and none were found; it is skipped rather than run on transformed numbers it cannot interpret correctly. The naive (per-cell) arm still runs, on the matrix as found. |
| `design_only` | the shared gene universe is too small after filtering | Not enough genes survive the frozen-universe filter (or, on the fallback path for non-integer counts, its own filter) to say anything useful. |
| `naive_only` | too few donor-pseudobulk profiles remain in a group after aggregation | The naive arm still ran, but there are not enough donor profiles per group for the pseudobulk arm to be meaningful. |
| `complete` | none of the above fired | Both arms ran: the naive per-cell test and the donor-pseudobulk test, each against its own donor-permutation null. |

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
  test would produce. `calibrated` means lambda falls in the tool's accepted band; `inflated` means
  the test is calling more than it should; `under` means fewer. `null` means that arm did not run.
- **The permutation floor (`naive_floor_solo`, and the pseudobulk floor).** pbcheck shuffles
  condition labels between donors many times, so that by construction there is no real signal to
  find, and counts how many genes each test still calls significant. That count, at a chosen
  quantile across the shuffles, is the floor: how many "hits" a test produces on pure noise, under
  this data's own donor structure. A test's real-label result is only informative once it is read
  against its own floor, not in isolation.
- **The negative control (the pseudobulk arm's floor and false-positive rate).** The donor-pseudobulk
  test is expected to keep its own floor near zero and its permutation false-positive rate near the
  chosen significance level. That behaviour is the check that the pseudobulk arm itself is well
  calibrated on this file; it does not certify the naive arm, and it does not certify any published
  result.
- **The envelope.** The pseudobulk arm's calibration and power were established on synthetic data
  only, and only inside a stated range of donor counts and per-donor variability. Quoting the exact
  wording pbcheck's own report uses for this:

  > The pseudobulk arm's calibration and power were established on synthetic oracles only inside
  > the operating envelope declared in Amendment 3 (minimum donors per group 4 / 8 / 13 / 23 at
  > sigma_donor 0.2 / 0.35 / 0.5 / 0.7). pbcheck does not estimate sigma_donor for real data, so
  > whether this stratum lies inside that envelope is not determined here.

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
| `settings` | Every setting the run used, split into the tool's own defaults and the small set of values that come from the pre-registered protocol (`protocol_constants`, taken from `pbcheck.gate_config`). |
| `universe` | The size and origin of the shared gene universe both arms were tested over, and the rule used to build it. |
| `real_label` | The naive and pseudobulk arms' results on the real condition labels: how many genes each called significant, the top genes, and the accounting of how the two arms' gene lists were made comparable. `null` if neither arm ran. |
| `permutation_null` | The donor-permutation null results for each arm that ran: inflation, floor, and (for the pseudobulk arm) the false-positive rate, plus how many permutations were requested versus actually distinct. `null` in `design_only`. |
| `readout` | The plain-language read-out derived from the fields above: inflation classes, floors, ratios of real result to floor, and the sentences shown in the report. |
| `caveats` | The list of caveat texts shown in this report, each with a stable identifier (see the "What each number means" section and the report itself for their wording). |
| `provenance` | What pbcheck was calibrated against, which settings are pre-registered protocol constants, and the platform and package versions the run used. |

## Runtime

<!-- filled by the integrator from WP1's measurement -->

## Product defaults

`--n-perm` and `--n-perm-pb` default to the product's own permutation counts, not to the
pre-registered protocol's counts; the exact numbers are set from a measurement done in another
part of this project and are listed in the runtime table above, not repeated here so that this
document does not go stale when that measurement is updated.
