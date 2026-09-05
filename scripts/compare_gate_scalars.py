#!/usr/bin/env python
"""Compare the headline scalars of a fresh synthetic-gate run against the frozen pilot artifact.

    python scripts/compare_gate_scalars.py NEW.json PILOT.json

The release checklist asks one question: did the gate numbers move? ``diff`` cannot answer it,
because a rerun legitimately adds keys the pilot run did not have (``bh_mode``,
``naive_perm_floor_solo``, ``naive_engine``) and rewrites timestamps and runtimes. This script
therefore compares exactly the scalars the gate's verdict rests on and ignores every other key on
both sides:

    lambda_naive, lambda_pseudobulk, pb_perm_null_fp_rate, naive_perm_floor.median_count,
    pb_perm_floor.median_count, power_sensitivity, power_empirical_fdr, verdict

Numbers are compared at the precision the pilot file records them at: the pilot value is read as a
decimal literal, the new value is rounded to that many decimal places, and the two must then be
equal. Recording a value as ``0.86`` is a statement about two decimal places; recording it as
``54.571780406237494`` is a statement about all seventeen, and both are honoured as written.

Exit status:

    0   every compared scalar is equal
    1   at least one differs, or the new file does not carry one of them
    2   wrong usage, or a file that is missing or is not JSON

Any inequality is a finding for the owner of the gate, never something this script waives.
Standard library only.
"""

from __future__ import annotations

import decimal
import json
import sys
from decimal import Decimal
from pathlib import Path

USAGE = "usage: compare_gate_scalars.py NEW.json PILOT.json"

#: The scalars the gate's verdict rests on, as key paths into the artifact.
COMPARED_PATHS: tuple[tuple[str, ...], ...] = (
    ("lambda_naive",),
    ("lambda_pseudobulk",),
    ("pb_perm_null_fp_rate",),
    ("naive_perm_floor", "median_count"),
    ("pb_perm_floor", "median_count"),
    ("power_sensitivity",),
    ("power_empirical_fdr",),
    ("verdict",),
)

#: Working precision for the decimal comparisons: above the 17 significant digits a repr'd double
#: needs, so rounding the new value to the pilot's precision never loses a digit of its own.
_PRECISION = 60

_MISSING = object()


def load(path: Path) -> dict:
    """Read a gate artifact, keeping every float as the decimal literal the file records."""
    try:
        with path.open(encoding="utf-8") as handle:
            payload = json.load(handle, parse_float=Decimal)
    except OSError as exc:
        sys.stderr.write(f"cannot read {path}: {exc}\n")
        raise SystemExit(2) from exc
    except json.JSONDecodeError as exc:
        sys.stderr.write(f"{path} is not valid JSON: {exc}\n")
        raise SystemExit(2) from exc
    if not isinstance(payload, dict):
        sys.stderr.write(f"{path} is not a JSON object\n")
        raise SystemExit(2)
    return payload


def dig(payload: dict, path: tuple[str, ...]):
    """The value at ``path``, or :data:`_MISSING` when any step of it is absent."""
    node = payload
    for key in path:
        if not isinstance(node, dict) or key not in node:
            return _MISSING
        node = node[key]
    return node


def equal_at_pilot_precision(new, pilot) -> bool:
    """Is ``new`` equal to ``pilot`` at the precision the pilot value is recorded with?"""
    if isinstance(pilot, bool) or isinstance(new, bool):
        return new is pilot
    if isinstance(pilot, (Decimal, int)) and isinstance(new, (Decimal, int)):
        pilot_dec = Decimal(pilot)
        new_dec = Decimal(new)
        exponent = pilot_dec.as_tuple().exponent
        if not isinstance(exponent, int):  # NaN or Infinity in the pilot file.
            return str(new_dec) == str(pilot_dec)
        with decimal.localcontext() as context:
            context.prec = _PRECISION
            quantum = Decimal(1).scaleb(exponent)
            return new_dec.quantize(quantum, rounding=decimal.ROUND_HALF_EVEN) == pilot_dec
    return type(new) is type(pilot) and new == pilot


def show(value) -> str:
    """One cell of the table."""
    return "(missing)" if value is _MISSING else str(value)


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2:
        sys.stderr.write(USAGE + "\n")
        return 2
    new_payload = load(Path(args[0]))
    pilot_payload = load(Path(args[1]))

    rows = []
    differences = 0
    for path in COMPARED_PATHS:
        name = ".".join(path)
        pilot_value = dig(pilot_payload, path)
        new_value = dig(new_payload, path)
        if pilot_value is _MISSING or new_value is _MISSING:
            same = False
        else:
            same = equal_at_pilot_precision(new_value, pilot_value)
        differences += 0 if same else 1
        rows.append((name, show(pilot_value), show(new_value), "same" if same else "DIFFERS"))

    widths = [max(len(row[column]) for row in rows) for column in range(4)]
    header = ("scalar", "pilot", "new", "verdict")
    widths = [max(width, len(header[column])) for column, width in enumerate(widths)]
    line = "  ".join("-" * width for width in widths)
    print("  ".join(header[column].ljust(widths[column]) for column in range(4)))
    print(line)
    for row in rows:
        print("  ".join(row[column].ljust(widths[column]) for column in range(4)))
    print(line)

    if differences:
        print(f"{differences} of {len(rows)} compared scalars moved: report to the gate's owner")
        return 1
    print(f"all {len(rows)} compared scalars are unchanged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
