"""Measure ``pbcheck.audit.run_audit`` wall-clock runtime and peak memory on fixed shapes.

This script exists to execute the runtime rule stated in the release plan (WP1, part 3): time
``run_audit`` on the reference shape at two permutation-count settings and, if the cheaper setting
finishes within the cap, adopt the richer one as the product default. It also times the example and
gate shapes at both settings for the hand-off table.

Each measurement runs in its own subprocess (``--worker`` mode of this same script) so that a run
which exceeds ``--cap`` seconds can be killed by the parent without losing the rest of the sweep.
The worker measures its own peak working set right before it exits (``GetProcessMemoryInfo`` on
Windows, ``resource.getrusage`` elsewhere) and writes wall seconds, the payload's own
per-stage timings, platform and CPU count to a small JSON result file that the parent reads back.

Usage::

    python scripts/measure_audit_runtime.py --out runtime_table.md
    python scripts/measure_audit_runtime.py --quick   # one tiny shape, for the test suite
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

#: Oracle draw shapes measured by the runtime rule. Each maps to ``synthetic.oracles.null_oracle``
#: keyword arguments; ``oracle_seed`` is the seed passed to ``null_oracle`` itself, independent of
#: ``AuditSettings.seed`` (the donor-permutation seed).
SHAPES = {
    # The runtime rule's reference shape (release plan, WP1 "Runtime rule" paragraph): 10k cells,
    # a realistic gene universe.
    "reference": {
        "oracle_seed": 1,
        "shape_kwargs": {"n_genes": 8000, "n_donors_per_group": 8, "n_cells_per_donor": 625},
    },
    # tests/conftest.py AUDIT_ORACLE_SHAPE, the smallest draw that clears all three floors.
    "example": {
        "oracle_seed": 7,
        "shape_kwargs": {"n_genes": 600, "n_donors_per_group": 4, "n_cells_per_donor": 60},
    },
    # pbcheck.gate_config.ORACLE_SIM, the shape scripts/synthetic_gate.py runs the calibration
    # sweep on (seed 1, matching that script's own default --seed).
    "gate": {
        "oracle_seed": 1,
        "shape_kwargs": None,  # filled from gate_config.ORACLE_SIM at call time (frozen module)
    },
}

#: (n_perm, n_perm_pb) settings measured by the runtime rule.
SETTINGS = [(1000, 200), (200, 200)]

#: Tiny shape/setting used by ``--quick`` (and by the test suite): fast enough to run inline.
QUICK_SHAPE = {"n_genes": 300, "n_donors_per_group": 3, "n_cells_per_donor": 30}
QUICK_SETTING = (5, 5)

DEFAULT_CAP_SECONDS = 900.0


def _gate_shape_kwargs() -> dict:
    """The gate shape's kwargs, read from the frozen ``gate_config.ORACLE_SIM`` at call time."""
    sys.path.insert(0, str(REPO_ROOT / "src"))
    from pbcheck import gate_config as gc

    sim = dict(gc.ORACLE_SIM)
    return {
        "n_genes": sim["n_genes"],
        "n_donors_per_group": sim["n_donors_per_group"],
        "n_cells_per_donor": sim["n_cells_per_donor"],
        "dispersion": sim["dispersion"],
        "donor_sigma": sim["donor_sigma"],
    }


# ---------------------------------------------------------------------------
# Worker: runs in a subprocess, measures itself, writes a result JSON.
# ---------------------------------------------------------------------------


def _peak_working_set_bytes() -> int | None:
    """This process's own peak working set (Windows) or maximum resident set size (POSIX)."""
    if platform.system() == "Windows":
        import ctypes
        from ctypes import wintypes

        class PROCESS_MEMORY_COUNTERS(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD),
                ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        counters = PROCESS_MEMORY_COUNTERS()
        counters.cb = ctypes.sizeof(PROCESS_MEMORY_COUNTERS)
        handle = ctypes.windll.kernel32.GetCurrentProcess()
        psapi = ctypes.WinDLL("psapi.dll")
        psapi.GetProcessMemoryInfo.argtypes = [
            wintypes.HANDLE, ctypes.POINTER(PROCESS_MEMORY_COUNTERS), wintypes.DWORD,
        ]
        psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
        ok = psapi.GetProcessMemoryInfo(handle, ctypes.byref(counters), counters.cb)
        if not ok:
            return None
        return int(counters.PeakWorkingSetSize)
    else:
        import resource

        # ru_maxrss is kilobytes on Linux, bytes on macOS; normalise to bytes assuming Linux since
        # that is what CI/production POSIX hosts here run.
        return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024


def _run_worker(spec_path: str, result_path: str) -> None:
    sys.path.insert(0, str(REPO_ROOT / "src"))
    sys.path.insert(0, str(REPO_ROOT / "synthetic"))
    from oracles import null_oracle  # noqa: PLC0415 (path set up just above)

    from pbcheck.audit import AuditSettings, run_audit  # noqa: PLC0415 (same deferred path setup)

    spec = json.loads(Path(spec_path).read_text(encoding="utf-8"))
    oracle = null_oracle(seed=spec["oracle_seed"], **spec["shape_kwargs"])
    settings = AuditSettings(
        donor_col="donor",
        condition_col="condition",
        test_level="disease",
        ref_level="ctrl",
        celltype_col="cell_type",
        celltype_value="T_cell",
        n_perm=spec["n_perm"],
        n_perm_pb=spec["n_perm_pb"],
        seed=0,
    )

    started = time.perf_counter()
    payload = run_audit(oracle.adata, settings)
    elapsed = time.perf_counter() - started

    result = {
        "elapsed_seconds": elapsed,
        "runtime_seconds": payload.get("runtime_seconds"),
        "runtime_by_stage_seconds": payload.get("runtime_by_stage_seconds", {}),
        "status": payload.get("status"),
        "peak_working_set_bytes": _peak_working_set_bytes(),
        "platform": platform.platform(),
        "cpu_count": os.cpu_count(),
    }
    Path(result_path).write_text(json.dumps(result), encoding="utf-8", newline="\n")


# ---------------------------------------------------------------------------
# Parent: launches one worker subprocess per measurement, applies the cap.
# ---------------------------------------------------------------------------


def _measure_one(shape_kwargs: dict, oracle_seed: int, n_perm: int, n_perm_pb: int,
                  cap_seconds: float, workdir: Path) -> dict:
    spec_path = workdir / "spec.json"
    result_path = workdir / "result.json"
    if result_path.exists():
        result_path.unlink()
    spec_path.write_text(
        json.dumps({
            "shape_kwargs": shape_kwargs,
            "oracle_seed": oracle_seed,
            "n_perm": n_perm,
            "n_perm_pb": n_perm_pb,
        }),
        encoding="utf-8",
        newline="\n",
    )

    cmd = [sys.executable, str(Path(__file__).resolve()), "--worker",
           "--spec", str(spec_path), "--result", str(result_path)]
    wall_started = time.perf_counter()
    try:
        proc = subprocess.run(
            cmd, cwd=str(REPO_ROOT), timeout=cap_seconds,
            capture_output=True, text=True, encoding="utf-8",
        )
    except subprocess.TimeoutExpired:
        return {"capped": True, "wall_seconds": cap_seconds}
    wall_seconds = time.perf_counter() - wall_started

    if proc.returncode != 0 or not result_path.exists():
        raise RuntimeError(
            f"measurement worker failed (returncode={proc.returncode}):\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    result = json.loads(result_path.read_text(encoding="utf-8"))
    result["capped"] = False
    result["wall_seconds"] = wall_seconds
    return result


def _format_bytes(n: int | None) -> str:
    if n is None:
        return "n/a"
    return f"{n / (1024 ** 2):.0f} MiB"


def _render_table(rows: list[dict], machine_note: str) -> str:
    lines = [
        f"Machine: {machine_note}",
        "",
        "| shape | n_perm | n_perm_pb | wall seconds | stage timings (s) | peak working set | status |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        if row.get("capped"):
            wall = f"> {row['wall_seconds']:.0f} (cap)"
            stages = "n/a (capped)"
            peak = "n/a (capped)"
            status = "capped"
        else:
            wall = f"{row['wall_seconds']:.1f}"
            stage_items = sorted(
                (k, v) for k, v in row["runtime_by_stage_seconds"].items() if v is not None
            )
            stages = ", ".join(f"{k}={v:.1f}" for k, v in stage_items)
            peak = _format_bytes(row.get("peak_working_set_bytes"))
            status = row.get("status", "n/a")
        lines.append(
            f"| {row['shape']} | {row['n_perm']} | {row['n_perm_pb']} | {wall} | {stages} "
            f"| {peak} | {status} |"
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--spec", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--result", type=Path, help=argparse.SUPPRESS)
    ap.add_argument("--cap", type=float, default=DEFAULT_CAP_SECONDS,
                     help="seconds allowed per measurement before it is killed and recorded as "
                          "'> cap' (default %(default)s)")
    ap.add_argument("--out", type=Path, default=None,
                     help="write the Markdown table to this path (also printed to stdout)")
    ap.add_argument("--quick", action="store_true",
                     help="run one tiny shape/setting instead of the full six-row sweep")
    args = ap.parse_args(argv)

    if args.worker:
        _run_worker(str(args.spec), str(args.result))
        return 0

    machine_note = f"{platform.platform()}, {platform.processor()}, {os.cpu_count()} CPUs"

    with tempfile.TemporaryDirectory(prefix="pbcheck_runtime_") as tmp:
        workdir = Path(tmp)
        rows: list[dict] = []

        if args.quick:
            result = _measure_one(QUICK_SHAPE, oracle_seed=1, n_perm=QUICK_SETTING[0],
                                   n_perm_pb=QUICK_SETTING[1], cap_seconds=args.cap,
                                   workdir=workdir)
            result["shape"] = "quick"
            result["n_perm"], result["n_perm_pb"] = QUICK_SETTING
            rows.append(result)
        else:
            shapes = dict(SHAPES)
            shapes["gate"] = {**shapes["gate"], "shape_kwargs": _gate_shape_kwargs()}
            for shape_name, shape_spec in shapes.items():
                for n_perm, n_perm_pb in SETTINGS:
                    result = _measure_one(
                        shape_spec["shape_kwargs"], shape_spec["oracle_seed"],
                        n_perm, n_perm_pb, cap_seconds=args.cap, workdir=workdir,
                    )
                    result["shape"] = shape_name
                    result["n_perm"], result["n_perm_pb"] = n_perm, n_perm_pb
                    rows.append(result)
                    outcome = "> cap" if result.get("capped") else f"{result['wall_seconds']:.1f}s"
                    print(f"measured {shape_name} n_perm={n_perm} n_perm_pb={n_perm_pb}: {outcome}",
                          file=sys.stderr)

        table = _render_table(rows, machine_note)

    print(table)
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(table, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
