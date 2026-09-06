# Contributing

**pbcheck is a pre-registered measurement study that also ships a released tool; the study's rules come first.**
[`docs/PHASE0_SPEC.md`](docs/PHASE0_SPEC.md) is the frozen protocol: it was written before any
data was seen, and it is not edited after the fact. Any change to the measurement protocol — a
test, a threshold, an oracle, a decision rule, which datasets get selected — requires a dated,
numbered entry in [`docs/AMENDMENTS.md`](docs/AMENDMENTS.md), **written and committed before the
code that applies it**, including a "data visible at the time" disclosure: what results were
already known when the amendment was written. See Amendment 1 and Amendment 2 for the pattern to
follow, including how to disclose an amendment written *after* the deciding data was seen (Amendment
2 does this openly rather than pretending otherwise).

Engineering changes — bug fixes, CI, refactors that don't move a single number the gate reports —
don't need an amendment, but they must keep the gate's numbers reproducible. If your change makes
`scripts/synthetic_gate.py` print different numbers with no amendment explaining why, that is a
protocol violation, not a refactor.

If you are unsure which category a change falls into, assume it needs an amendment and ask before
writing code.

## Status

Phase 0 is mid-remediation. See [`pilot/README.md`](pilot/README.md) for what is built, what is
gapped, and what is still pending on real data; `docs/AMENDMENTS.md` is the append-only record of
every departure from the frozen spec.

## Dev setup

Requires Python >= 3.12 (the verified stack — scanpy 1.12, anndata 0.13, pandas 3.0 — does not
resolve below it).

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[dev]"   # Windows; POSIX: .venv/bin/python
pre-commit install
```

## Running things

```bash
pytest -q                          # full suite
pytest -q -m "not slow"            # fast loop, skips DESeq2 / multi-permutation end-to-end tests
ruff check .                       # lint (no --fix; see .pre-commit-config.yaml for why)
python scripts/synthetic_gate.py   # the synthetic calibration gate — CPU only, no downloads
```

`pre-commit run --all-files` runs the same lint plus basic hygiene hooks (trailing whitespace,
end-of-file newline, YAML validity, a large-file guard). `docs/PHASE0_SPEC.md` and
`docs/AMENDMENTS.md` are excluded from the fixing hooks — their bytes are frozen protocol text,
not free for a mechanical pass to touch.

## Commit style

Imperative, single-line subject stating the *intent* of the change, not a log of which files
moved — e.g. "Switch the pseudobulk arm to moderated eBayes", not "update methods.py". One
logical change per commit. If a commit changes the measurement protocol, its subject should make
that obvious, and the corresponding `docs/AMENDMENTS.md` entry should already exist in an earlier
commit.

## What this repo deliberately does not have

No CODEOWNERS, issue/PR templates, code of conduct, coverage threshold gate, generated docs
site, or logging framework. This is a solo research repo with a released tool; that scaffolding
is a deliberate omission, not an oversight.

## Releases

Tag `vX.Y.Z` on `main`. The version is kept in sync in `src/pbcheck/__init__.py` and
`CITATION.cff` (`scripts/check_version_consistency.py --tag`). `release.yml` publishes to PyPI by
trusted publishing; the GitHub Release it creates triggers the Zenodo archive.

### Release procedure

Run both of these locally, as the last two gates before pushing the tag; neither runs the tag
push itself.

1. `python scripts/check_version_consistency.py --tag vX.Y.Z --require-date` - confirms
   `src/pbcheck/__init__.py`, `CITATION.cff`, `.zenodo.json`, the `CHANGELOG.md` heading and the
   tag itself all name the same version, and that the `CHANGELOG.md` release date is a real ISO
   date matching `CITATION.cff`'s `date-released`. A non-zero exit means do not push the tag yet.
2. `python scripts/protocol_safety_check.py --with-gate` - on Windows with Python 3.12, the
   platform and interpreter recorded in the committed gate artifact
   (`pilot/gate/synthetic_gate_2026-08-15.json`). It fails if product code reads a Phase 0
   permutation-count constant, contains a forbidden protocol string, or has touched a frozen path
   since the last tag.

Decision of 2026-09-05 (owner): shipping v0.1.0 is treated as an engineering change under the
amendment test because it changes no number `scripts/synthetic_gate.py` prints and touches no
frozen file; the product's wording is a protocol surface (`src/pbcheck/render/text.py`, tested
for forbidden claims), and any future product output that could be read as a Phase 0 conclusion
needs an amendment.
