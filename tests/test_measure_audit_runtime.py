"""Test for ``scripts/measure_audit_runtime.py``'s ``--quick`` path.

Only the quick path runs inline here: the full six-row sweep on the reference shape can take
many minutes and belongs in the one-off hand-off measurement, not in the regular test suite.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "measure_audit_runtime.py"

# One header row (shape/n_perm/... labels), one separator row, one data row.
_TABLE_ROW_RE = re.compile(r"^\| quick \| 5 \| 5 \|.*\|$", re.MULTILINE)


def test_quick_path_produces_table_with_one_data_row(tmp_path):
    out_path = tmp_path / "table.md"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--quick", "--out", str(out_path)],
        cwd=str(REPO_ROOT), capture_output=True, text=True, encoding="utf-8", timeout=300,
    )
    assert result.returncode == 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    assert out_path.exists()
    table = out_path.read_text(encoding="utf-8")

    data_rows = _TABLE_ROW_RE.findall(table)
    assert len(data_rows) == 1, table
    assert "Machine:" in table
    assert "| shape | n_perm | n_perm_pb |" in table
