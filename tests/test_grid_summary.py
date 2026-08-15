"""The grid driver reads the per-cell JSON the probe actually writes.

``scripts/run_test_selection_grid.py`` is a pure aggregator: it re-keys the probe's per-cell
JSON into ``pilot/testsel/summary.{csv,json}``. Because it reads with ``dict.get``, a key it
gets wrong does not fail -- it writes ``None`` for every row and the loss is invisible until
someone tries to use the column. That happened: the driver asked for ``power_sd`` and
``n_power_reps`` while the probe writes ``power_sd_across_datasets`` / ``power_se`` /
``power_n_reps``, so all 146 rows of the grid committed at 72dec7b carry no Monte-Carlo error
on their power figures (docs/AMENDMENTS.md, Amendment 3). These tests make the next such
drift a red test rather than a silent column of nulls.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
DRIVER = REPO / "scripts" / "run_test_selection_grid.py"
PROBE = REPO / "scripts" / "pb_calibration_probe.py"


def _load_driver():
    if not DRIVER.exists():  # pragma: no cover
        pytest.skip("grid driver not present")
    spec = importlib.util.spec_from_file_location("_run_test_selection_grid", DRIVER)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _driver_top_level_keys() -> set[str]:
    """Every ``o.get("...")`` in the driver -- i.e. the top-level keys it expects the probe to
    write. Read from the source rather than by running ``summarise``, so a key that is only
    reachable on some branch still counts."""
    tree = ast.parse(DRIVER.read_text(encoding="utf-8"))
    keys = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute) and node.func.attr == "get"
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "o"
                and node.args and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)):
            keys.add(node.args[0].value)
    assert keys, "no o.get(...) calls found -- the driver was restructured, update this test"
    return keys


def _probe_written_keys() -> set[str]:
    """Every string key appearing in a dict literal or an ``out[...] = `` assignment in the
    probe. A superset of the top-level output keys, which is the safe direction: this test
    only fails when the driver asks for a name the probe writes nowhere at all."""
    tree = ast.parse(PROBE.read_text(encoding="utf-8"))
    keys = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Dict):
            keys.update(k.value for k in node.keys
                        if isinstance(k, ast.Constant) and isinstance(k.value, str))
        elif isinstance(node, ast.Subscript) and isinstance(node.slice, ast.Constant):
            if isinstance(node.slice.value, str):
                keys.add(node.slice.value)
    return keys


def test_every_key_the_driver_reads_is_a_key_the_probe_writes():
    if not PROBE.exists():  # pragma: no cover
        pytest.skip("probe script not present")
    missing = sorted(_driver_top_level_keys() - _probe_written_keys())
    assert not missing, (
        f"the driver reads {missing} but the probe never writes those names; "
        "these columns would be None in every row of the summary"
    )


def test_the_power_monte_carlo_error_survives_the_summary():
    """The specific regression: mean, SD, SE and n must all reach the summary row."""
    mod = _load_driver()
    cell = {
        "lambda_pb": 1.01,
        "gen_arm": "directnb", "test": "ebayes",
        "sigma_het": 0.4, "sigma_donor": 0.5,
        "n_donors_per_group": 8, "n_perm": 3000,
        "power_at_pre_registered": 0.42,
        "power_sd_across_datasets": 0.05,
        "power_se": 0.01,
        "power_n_reps": 25,
        "seconds_elapsed": 72.3,
        "test_diagnostics": {"prior_df_d0": 8.0,
                             "shrinkage_factor_d0_over_d0_plus_d": 0.36},
    }
    row = _summarise_one(mod, cell)
    assert row["power"] == 0.42
    assert row["power_sd"] == 0.05
    assert row["power_se"] == 0.01
    assert row["n_power_reps"] == 25
    assert "power_se" in mod.SUMMARY_FIELDS
    assert set(row) == set(mod.SUMMARY_FIELDS)


def test_a_single_replicate_cell_reports_no_monte_carlo_error():
    """``--n-power-reps 1`` leaves SD and SE genuinely None; that must stay distinguishable
    from the key-mismatch bug, which produced the same None for a 25-replicate cell."""
    mod = _load_driver()
    cell = {
        "lambda_pb": 1.0, "power_at_pre_registered": 0.42,
        "power_sd_across_datasets": None, "power_se": None, "power_n_reps": 1,
    }
    row = _summarise_one(mod, cell)
    assert row["power"] == 0.42
    assert row["power_sd"] is None and row["power_se"] is None
    assert row["n_power_reps"] == 1


def _summarise_one(mod, cell: dict) -> dict:
    """Run the driver's real ``summarise`` over a single synthetic per-cell JSON."""
    with tempfile.TemporaryDirectory() as td:
        out_dir = Path(td) / "cells"
        out_dir.mkdir()
        (out_dir / "A__cell.json").write_text(json.dumps(cell), encoding="utf-8")
        summary_dir = Path(td) / "summary"
        assert mod.summarise(out_dir, summary_dir) == 1
        rows = json.loads((summary_dir / "summary.json").read_text(encoding="utf-8"))
    return rows[0]
