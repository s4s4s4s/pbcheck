"""Tests for :mod:`pbcheck.cli`, the ``pbcheck`` console script.

Subprocess-free throughout: every test calls ``cli.main([...])`` directly and reads its return
code plus ``capsys``'s captured stdout/stderr, the same interface the real console script and
``python -m pbcheck`` both resolve to. The one file this module writes to disk per test session
is ``x.h5ad`` (the shared ``audit_shape_oracle`` draw, see ``tests/conftest.py``); everything else
lives under a test's own ``tmp_path``.
"""

from __future__ import annotations

import importlib.metadata
import json
import re
import runpy
import sys

import numpy as np
import pytest

from pbcheck import AuditInputError, AuditSettings, audit_h5ad, audit_schema, cli, run_audit
from pbcheck.render.sections import build_sections, prose_blocks
from pbcheck.render.text import FORBIDDEN_PATTERNS

# The oracle's own column names and levels (synthetic/oracles.py), matching tests/test_audit.py.
DONOR_COL = "donor"
CONDITION_COL = "condition"
CELLTYPE_COL = "cell_type"
CELLTYPE_VALUE = "T_cell"
TEST_LEVEL = "disease"
REF_LEVEL = "ctrl"

N_PERM = 20
N_PERM_PB = 10

#: Silences the "unused import" complaint for names only referenced to prove the re-export works
#: (test_package_reexports below imports them directly instead of through this module).
_ = (AuditInputError, AuditSettings, audit_h5ad, run_audit)


@pytest.fixture(scope="module")
def oracle_h5ad(tmp_path_factory, audit_shape_oracle):
    """``audit_shape_oracle.adata`` written to disk once per module, for every CLI test to read."""
    path = tmp_path_factory.mktemp("cli_oracle") / "x.h5ad"
    audit_shape_oracle.adata.write_h5ad(path)
    return path


def _base_audit_args(file, out) -> list[str]:
    return [
        "audit", str(file),
        "--donor", DONOR_COL,
        "--condition", CONDITION_COL,
        "--test", TEST_LEVEL,
        "--ref", REF_LEVEL,
        "--celltype", CELLTYPE_COL,
        "--celltype-value", CELLTYPE_VALUE,
        "--n-perm", str(N_PERM),
        "--n-perm-pb", str(N_PERM_PB),
        "--out", str(out),
    ]


def test_audit_writes_three_files_and_validates(oracle_h5ad, tmp_path, capsys):
    out = tmp_path / "out"
    code = cli.main(_base_audit_args(oracle_h5ad, out))
    captured = capsys.readouterr()
    assert code == 0, captured.err

    json_path = out / "pbcheck_audit.json"
    md_path = out / "pbcheck_report.md"
    html_path = out / "pbcheck_report.html"
    for path in (json_path, md_path, html_path):
        assert path.exists(), path
        assert str(path) in captured.out

    payload = json.loads(json_path.read_text(encoding="utf-8"))
    audit_schema.validate(payload)
    assert payload["status"] == "complete"


def test_default_out_dir_names_stem_and_celltype(oracle_h5ad, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    args = [
        "audit", str(oracle_h5ad),
        "--donor", DONOR_COL,
        "--condition", CONDITION_COL,
        "--test", TEST_LEVEL,
        "--ref", REF_LEVEL,
        "--celltype", CELLTYPE_COL,
        "--celltype-value", CELLTYPE_VALUE,
        "--n-perm", str(N_PERM),
        "--n-perm-pb", str(N_PERM_PB),
    ]
    code = cli.main(args)
    captured = capsys.readouterr()
    assert code == 0, captured.err

    expected_dir = tmp_path / "pbcheck_out" / f"{oracle_h5ad.stem}_{CELLTYPE_VALUE}"
    assert (expected_dir / "pbcheck_audit.json").exists()


def test_design_only_needs_no_counts(audit_shape_oracle, tmp_path, capsys):
    """``--design-only`` never looks at the counts matrix: float zeros in ``X`` still succeed."""
    adata = audit_shape_oracle.adata.copy()
    adata.X = np.zeros(adata.shape, dtype=np.float64)
    path = tmp_path / "zeros.h5ad"
    adata.write_h5ad(path)

    out = tmp_path / "out"
    code = cli.main([
        "audit", str(path),
        "--donor", DONOR_COL,
        "--condition", CONDITION_COL,
        "--test", TEST_LEVEL,
        "--ref", REF_LEVEL,
        "--design-only",
        "--out", str(out),
    ])
    captured = capsys.readouterr()
    assert code == 0, captured.err
    payload = json.loads((out / "pbcheck_audit.json").read_text(encoding="utf-8"))
    assert payload["status"] == "design_only"
    assert payload["status_reason"] == "design_only_requested"


def test_missing_column_exit_2_lists_columns(oracle_h5ad, tmp_path, capsys):
    out = tmp_path / "out"
    code = cli.main([
        "audit", str(oracle_h5ad),
        "--donor", "no_such_column",
        "--condition", CONDITION_COL,
        "--test", TEST_LEVEL,
        "--ref", REF_LEVEL,
        "--out", str(out),
    ])
    captured = capsys.readouterr()
    assert code == 2
    assert "no_such_column" in captured.err
    assert "available columns" in captured.err
    assert not out.exists()


def test_celltype_flags_pairing_exit_2(oracle_h5ad, tmp_path, capsys):
    out = tmp_path / "out"
    code = cli.main([
        "audit", str(oracle_h5ad),
        "--donor", DONOR_COL,
        "--condition", CONDITION_COL,
        "--test", TEST_LEVEL,
        "--ref", REF_LEVEL,
        "--celltype", CELLTYPE_COL,
        "--out", str(out),
    ])
    captured = capsys.readouterr()
    assert code == 2
    assert "--celltype" in captured.err and "--celltype-value" in captured.err


def test_celltype_value_unknown_lists_levels(oracle_h5ad, tmp_path, capsys):
    out = tmp_path / "out"
    code = cli.main([
        "audit", str(oracle_h5ad),
        "--donor", DONOR_COL,
        "--condition", CONDITION_COL,
        "--test", TEST_LEVEL,
        "--ref", REF_LEVEL,
        "--celltype", CELLTYPE_COL,
        "--celltype-value", "no_such_level",
        "--out", str(out),
    ])
    captured = capsys.readouterr()
    assert code == 2
    assert "no_such_level" in captured.err
    assert "available values" in captured.err


def test_format_json_only(oracle_h5ad, tmp_path, capsys):
    out = tmp_path / "out"
    code = cli.main(_base_audit_args(oracle_h5ad, out) + ["--format", "json"])
    captured = capsys.readouterr()
    assert code == 0, captured.err
    assert (out / "pbcheck_audit.json").exists()
    assert not (out / "pbcheck_report.md").exists()
    assert not (out / "pbcheck_report.html").exists()


def test_overwrite_guard(oracle_h5ad, tmp_path, capsys):
    out = tmp_path / "out"
    first = cli.main(_base_audit_args(oracle_h5ad, out))
    assert first == 0, capsys.readouterr().err

    second = cli.main(_base_audit_args(oracle_h5ad, out))
    captured = capsys.readouterr()
    assert second == 2
    assert "pbcheck_audit.json" in captured.err or "pbcheck_report" in captured.err

    third = cli.main(_base_audit_args(oracle_h5ad, out) + ["--overwrite"])
    assert third == 0, capsys.readouterr().err


def test_non_ascii_path_round_trip(audit_shape_oracle, tmp_path, capsys):
    non_ascii_dir = tmp_path / "данные"
    non_ascii_dir.mkdir()
    h5ad_path = non_ascii_dir / "x.h5ad"
    audit_shape_oracle.adata.write_h5ad(h5ad_path)
    out = non_ascii_dir / "выход"

    code = cli.main(_base_audit_args(h5ad_path, out))
    captured = capsys.readouterr()
    assert code == 0, captured.err
    assert (out / "pbcheck_audit.json").exists()


def test_version_flag(capsys):
    code = cli.main(["--version"])
    captured = capsys.readouterr()
    assert code == 0
    assert captured.out.strip() == "pbcheck 0.1.0"


def test_python_m_entry(monkeypatch, capsys):
    """``python -m pbcheck`` resolves to ``pbcheck.__main__``, which this exercises directly."""
    monkeypatch.setattr(sys, "argv", ["pbcheck", "--version"])
    with pytest.raises(SystemExit) as excinfo:
        runpy.run_module("pbcheck.__main__", run_name="__main__")
    assert excinfo.value.code == 0
    assert "pbcheck 0.1.0" in capsys.readouterr().out


def test_example_then_audit_quickstart(tmp_path, capsys):
    example_path = tmp_path / "quickstart.h5ad"
    code = cli.main(["example", str(example_path)])
    assert code == 0, capsys.readouterr().err
    assert example_path.exists()

    out = tmp_path / "quickstart_out"
    code = cli.main([
        "audit", str(example_path),
        "--donor", "donor",
        "--condition", "condition",
        "--test", "disease",
        "--ref", "ctrl",
        "--n-perm", "50",
        "--n-perm-pb", "20",
        "--out", str(out),
    ])
    captured = capsys.readouterr()
    assert code == 0, captured.err

    payload = json.loads((out / "pbcheck_audit.json").read_text(encoding="utf-8"))
    assert payload["status"] == "complete"
    assert payload["readout"]["lambda_naive_class"] == "inflated"


def test_cli_report_prose_has_no_forbidden_patterns(tmp_path):
    example_path = tmp_path / "quickstart.h5ad"
    cli.main(["example", str(example_path)])
    out = tmp_path / "quickstart_out"
    code = cli.main([
        "audit", str(example_path),
        "--donor", "donor",
        "--condition", "condition",
        "--test", "disease",
        "--ref", "ctrl",
        "--n-perm", "50",
        "--n-perm-pb", "20",
        "--out", str(out),
        "--quiet",
    ])
    assert code == 0

    payload = json.loads((out / "pbcheck_audit.json").read_text(encoding="utf-8"))
    sections = build_sections(payload)
    for text in prose_blocks(sections):
        for pattern in FORBIDDEN_PATTERNS:
            assert re.search(pattern, text) is None, (
                f"prose block matches forbidden pattern {pattern!r}: {text!r}"
            )


def test_package_reexports():
    from pbcheck import AuditSettings, audit_h5ad, example_adata, run_audit  # noqa: F401
    assert callable(example_adata)


def test_console_script_entry_point():
    try:
        importlib.metadata.distribution("pbcheck")
    except importlib.metadata.PackageNotFoundError:
        pytest.skip("the pbcheck distribution is not installed (pip install -e . not run)")
    entries = importlib.metadata.entry_points(group="console_scripts")
    names = {entry.name for entry in entries}
    assert "pbcheck" in names
