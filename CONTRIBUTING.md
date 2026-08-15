# Contributing

**pbcheck is a pre-registered measurement study, not a conventional software project.**
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

No CODEOWNERS, issue/PR templates, code of conduct, coverage threshold gate, PyPI release
workflow, generated docs site, or logging framework. This is a solo pre-release research repo;
that scaffolding is a deliberate omission, not an oversight.
