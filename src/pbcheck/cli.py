"""The ``pbcheck`` command-line interface.

Three subcommands: ``pbcheck --version``, ``pbcheck example OUT.h5ad`` (the offline quickstart
generator, :mod:`pbcheck.example`) and ``pbcheck audit FILE.h5ad ...`` (the single-stratum audit,
:mod:`pbcheck.audit` plus :mod:`pbcheck.render`). This module owns argument parsing, exit codes and
the run's progress lines; every word of the audit's own findings comes from
:func:`pbcheck.render.sections.summary_lines`, never retyped here.

Exit codes: ``0`` when the three output files were written, whatever ``status`` the payload landed
in (``complete``, ``naive_only`` or ``design_only`` are all a completed run, not a failure); ``2``
for a usage or input error (a bad flag, a missing file, a bad ``.obs`` column or level,
``--celltype``/``--celltype-value`` given one without the other, an existing output file without
``--overwrite``); ``1`` for any other exception, printed as a traceback so a bug is diagnosable
rather than swallowed.
"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path

from pbcheck import __version__
from pbcheck.audit import PRODUCT_N_PERM, PRODUCT_N_PERM_PB, AuditInputError, AuditSettings, audit_h5ad
from pbcheck.example import REFERENCE_SHAPE, SMALL_SHAPE, example_adata
from pbcheck.gate_config import ALPHA
from pbcheck.render import write_outputs
from pbcheck.render.sections import summary_lines

#: Process exit codes (see the module docstring).
EXIT_OK = 0
EXIT_USAGE = 2
EXIT_ERROR = 1

#: ``--format`` choices and what each expands to for :func:`pbcheck.render.write_outputs`.
_FORMAT_CHOICES = ("json", "md", "html", "all")
_FORMAT_EXPANSION: dict[str, tuple[str, ...]] = {
    "json": ("json",),
    "md": ("md",),
    "html": ("html",),
    "all": ("json", "md", "html"),
}

#: The value substituted for ``--celltype-value`` in the default output directory name when no
#: cell type was selected (section 1.1 of the plan: "the file's own stratum name").
_DEFAULT_OUT_ALL = "all"


def _reconfigure_streams() -> None:
    """Force UTF-8 on stdout/stderr so gene symbols and non-ASCII paths print on any console.

    ``TextIOWrapper.reconfigure`` is a no-op on streams that are not one (``capsys`` in tests
    replaces them with objects that may not offer it), hence the ``getattr`` guard.
    """
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name)
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(encoding="utf-8")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pbcheck",
        description="An auditor of pseudoreplication in single-cell RNA-seq differential expression.",
    )
    parser.add_argument("--version", action="store_true", help="print the version and exit")

    subparsers = parser.add_subparsers(dest="command")

    example = subparsers.add_parser(
        "example", help="write a small synthetic .h5ad file for the quickstart"
    )
    example.add_argument("out", type=Path, help="path of the .h5ad file to write")
    example.add_argument("--seed", type=int, default=0)
    example.add_argument(
        "--shape", choices=("small", "reference"), default="small",
        help="small: 600 genes, 4v4 donors, 60 cells/donor; "
             "reference: 8000 genes, 8v8 donors, 625 cells/donor",
    )

    audit = subparsers.add_parser("audit", help="audit one stratum of one .h5ad file")
    audit.add_argument("file", type=Path, help="the .h5ad file to audit")
    audit.add_argument("--donor", dest="donor_col", required=True, metavar="COL")
    audit.add_argument("--condition", dest="condition_col", required=True, metavar="COL")
    audit.add_argument("--test", dest="test_level", required=True, metavar="LEVEL")
    audit.add_argument("--ref", dest="ref_level", required=True, metavar="LEVEL")
    audit.add_argument("--celltype", dest="celltype_col", default=None, metavar="COL")
    audit.add_argument("--celltype-value", dest="celltype_value", default=None, metavar="VALUE")
    audit.add_argument(
        "--batch", dest="batch_cols", action="append", default=[], metavar="COL",
        help="repeatable",
    )
    audit.add_argument("--counts-layer", dest="counts_layer", default=None, metavar="NAME")
    audit.add_argument("--n-perm", dest="n_perm", type=int, default=PRODUCT_N_PERM)
    audit.add_argument("--n-perm-pb", dest="n_perm_pb", type=int, default=PRODUCT_N_PERM_PB)
    audit.add_argument("--seed", dest="seed", type=int, default=0)
    audit.add_argument("--alpha", dest="alpha", type=float, default=ALPHA)
    audit.add_argument("--top", dest="top_n", type=int, default=25)
    audit.add_argument("--design-only", dest="design_only", action="store_true")
    audit.add_argument("--out", dest="out", type=Path, default=None)
    audit.add_argument("--format", dest="format", choices=_FORMAT_CHOICES, default="all")
    audit.add_argument("--overwrite", dest="overwrite", action="store_true")
    audit.add_argument("--quiet", dest="quiet", action="store_true")

    return parser


def _default_out_dir(file: Path, celltype_value: str | None) -> Path:
    """``./pbcheck_out/<file stem>_<celltype value or "all">`` (the plan's default)."""
    suffix = celltype_value if celltype_value is not None else _DEFAULT_OUT_ALL
    return Path("pbcheck_out") / f"{file.stem}_{suffix}"


def _print_stage_lines(payload: dict, file: Path) -> None:
    """The run's progress lines, each carrying the real elapsed seconds of the stage it names.

    Built from ``payload["runtime_by_stage_seconds"]``, the same numbers the JSON output carries,
    so nothing printed here is a second measurement that could disagree with the file. The naive
    and pseudobulk real-label calls are timed together as one stage (``audit._run_both_arms``
    opens a single ``real_label`` block around both): "Naive arm" carries that measurement and
    "Pseudobulk arm" is printed as a bare announcement rather than attributed a fabricated second
    number.
    """
    stages = payload["runtime_by_stage_seconds"]

    def _emit(label: str, key: str) -> None:
        seconds = stages.get(key)
        if seconds is not None:
            print(f"{label} ... {seconds:.1f} s")

    load_seconds = stages.get("load")
    if load_seconds is not None:
        print(f"Loading {file} ... {load_seconds:.1f} s")
    _emit("Design audit", "design")
    _emit("Counts check", "counts")
    _emit("Pseudobulk profiles", "pseudobulk_build")

    real_label_seconds = stages.get("real_label")
    if real_label_seconds is not None:
        print(f"Naive arm ... {real_label_seconds:.1f} s")
        if payload["status"] == "complete":
            print("Pseudobulk arm")

    permutation_seconds = stages.get("permutation_null")
    if permutation_seconds is not None:
        n_naive = payload["readout"]["n_perm_naive_requested"]
        n_pb = payload["readout"]["n_perm_pb_requested"]
        print(
            f"Permutation null ({n_naive} naive / {n_pb} pseudobulk requested) "
            f"... {permutation_seconds:.1f} s"
        )


def _run_example(args: argparse.Namespace) -> int:
    shape = SMALL_SHAPE if args.shape == "small" else REFERENCE_SHAPE
    adata = example_adata(args.seed, **shape)
    out = args.out
    if out.parent != Path("."):
        out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out)
    print(f"Wrote {out} ({adata.n_obs} cells, {adata.n_vars} genes, shape={args.shape})")
    return EXIT_OK


def _run_audit(args: argparse.Namespace) -> int:
    if (args.celltype_col is None) != (args.celltype_value is None):
        print(
            "error: --celltype and --celltype-value must be given together", file=sys.stderr
        )
        return EXIT_USAGE

    settings = AuditSettings(
        donor_col=args.donor_col,
        condition_col=args.condition_col,
        test_level=args.test_level,
        ref_level=args.ref_level,
        celltype_col=args.celltype_col,
        celltype_value=args.celltype_value,
        batch_cols=tuple(args.batch_cols),
        counts_layer=args.counts_layer,
        n_perm=args.n_perm,
        n_perm_pb=args.n_perm_pb,
        seed=args.seed,
        alpha=args.alpha,
        design_only=args.design_only,
        top_n=args.top_n,
    )

    try:
        payload = audit_h5ad(args.file, settings)
    except AuditInputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception:
        traceback.print_exc()
        return EXIT_ERROR

    if not args.quiet:
        _print_stage_lines(payload, args.file)

    out_dir = args.out if args.out is not None else _default_out_dir(args.file, args.celltype_value)
    try:
        written = write_outputs(
            payload, out_dir, formats=_FORMAT_EXPANSION[args.format], overwrite=args.overwrite
        )
    except FileExistsError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception:
        traceback.print_exc()
        return EXIT_ERROR

    if not args.quiet:
        for line in summary_lines(payload):
            print(line)
        for path in written.values():
            print(str(path))
    return EXIT_OK


def main(argv: list[str] | None = None) -> int:
    """Parse ``argv`` (``sys.argv[1:]`` when ``None``) and run the requested subcommand."""
    _reconfigure_streams()
    parser = _build_parser()
    try:
        args = parser.parse_args(argv)
    except SystemExit as exc:
        return EXIT_USAGE if exc.code is None else int(exc.code)

    if args.version:
        print(f"pbcheck {__version__}")
        return EXIT_OK
    if args.command == "example":
        return _run_example(args)
    if args.command == "audit":
        return _run_audit(args)

    parser.print_usage(sys.stderr)
    return EXIT_USAGE


if __name__ == "__main__":
    sys.exit(main())
