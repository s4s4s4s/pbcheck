"""The gate's thresholds, in one place, with each one's provenance attached.

**Everything in :data:`PRE_REGISTERED` is frozen by ``docs/PHASE0_SPEC.md``. Changing any of these
values is a protocol change and requires a dated, numbered entry in ``docs/AMENDMENTS.md``, written
and committed BEFORE the code that applies it.** That is not a style preference: Amendment 1 exists
because ``scripts/synthetic_gate.py`` had quietly substituted an easier oracle (log2FC 1.5 / K 150
in place of the pre-registered 1.0 / 200) and dropped two binding criteria. Thresholds scattered as
literals across scripts is the condition under which that happens and goes unnoticed, so they live
here instead, each next to the spec section that pins it.

The second block, :data:`INSTRUMENT_SANITY`, is **not pre-registered** and is labelled as such. Those
are the synthetic-gate script's own sanity checks on the *measurement instrument* — "does the naive
arm read as inflated at all on data where we know it must" — not criteria from the decision rule.
They are separated so that no future reader can mistake a convenience threshold for a frozen one.

Note on ``scripts/pb_calibration_probe.py``: it carries its own copy of these constants and is
deliberately **not** made to import this module. The probe is the independent reference instrument
that the Amendment 2 selection grid was measured with and that ``tests/test_moderated.py`` validates
the shipped arm against; a reference that imports the thing it is checking is not a reference.
``tests/test_gate_config.py`` pins the two copies to be equal instead, which gives the single source
of truth without the circularity.
"""

from __future__ import annotations

from types import MappingProxyType

# ---------------------------------------------------------------------------
# PRE-REGISTERED — frozen by docs/PHASE0_SPEC.md. Amendment required to change.
# ---------------------------------------------------------------------------

#: Significance level, everywhere, including permutation-null metrics (spec §5 item 6).
ALPHA = 0.05

#: Genomic-inflation band a calibrated arm must sit in (decision rule item 1; §8(a)).
LAMBDA_BAND = (0.9, 1.1)

#: Pseudobulk sensitivity required on the synthetic-positive oracle (decision rule item 1; §8(c)).
POWER_TARGET = 0.60

#: The synthetic-positive oracle's effect size and gene count (§8(c)). Substituting an easier
#: oracle in gate code is prohibited by Amendment 1 Change 2.
ORACLE_LOG2FC = 1.0
ORACLE_K = 200

#: Minimum frozen-universe size before the pseudobulk arm may run (inclusion gate item 5; C5).
MIN_UNIVERSE_SIZE = 200

#: Thin-donor thresholds; donors below them are dropped, not merged (§1 item 2, §3).
#: Implemented manually as of Amendment 2 Change 7 — decoupler 2.x removed the parameters.
MIN_CELLS = 10
MIN_COUNTS = 1000

#: Naive inflation required across independent real datasets (decision rule item 2 / bands). These
#: bind the REAL sweep, not the synthetic instrument check below.
LAMBDA_NAIVE_GO = 2.0
LAMBDA_NAIVE_STRONG_GO = 3.0

#: Strong-GO band: naive floor at least this multiple of the BH complete-null expectation (~0).
NAIVE_FLOOR_OVER_NULL_STRONG_GO = 10

PRE_REGISTERED = MappingProxyType({
    "alpha": ALPHA,
    "lambda_band": LAMBDA_BAND,
    "power_target": POWER_TARGET,
    "oracle_log2fc": ORACLE_LOG2FC,
    "oracle_k": ORACLE_K,
    "min_universe_size": MIN_UNIVERSE_SIZE,
    "min_cells": MIN_CELLS,
    "min_counts": MIN_COUNTS,
    "lambda_naive_go": LAMBDA_NAIVE_GO,
    "lambda_naive_strong_go": LAMBDA_NAIVE_STRONG_GO,
    "naive_floor_over_null_strong_go": NAIVE_FLOOR_OVER_NULL_STRONG_GO,
})

# ---------------------------------------------------------------------------
# NOT PRE-REGISTERED — the synthetic gate's own instrument checks.
#
# These answer "is the instrument reading at all?" on data whose truth we set, and they are
# deliberately loose. They are NOT the decision rule and must never be cited as if they were. The
# spec's own naive threshold for the real sweep is LAMBDA_NAIVE_GO = 2.0 across independent
# datasets at matched cells-per-donor, which is a different and stricter claim.
# ---------------------------------------------------------------------------

#: On a synthetic null with donor structure the naive arm must at least register as inflated.
INSTRUMENT_LAMBDA_NAIVE_MIN = 1.5

#: ...and its false-positive floor must be a substantial fraction of the frozen universe.
INSTRUMENT_NAIVE_FLOOR_FRAC_MIN = 0.30

#: ...and it must exceed the pseudobulk floor by this factor (guarding against a floor of 1 vs 0
#: reading as a triumph).
INSTRUMENT_FLOOR_RATIO_MIN = 10

#: The synthetic oracle's operating point: 8v8 donors is the regime the spec says carries the
#: headline (permutations approximately orthogonal to truth — §4/A1). ``donor_sigma`` is a FREE KNOB
#: of the simulator and is **not** anchored to real data; both Amendment 1 and Amendment 2 close on
#: that as the outstanding threat to every power number computed here.
ORACLE_SIM = MappingProxyType({
    "n_genes": 1500,
    "n_donors_per_group": 8,
    "n_cells_per_donor": 250,
    "dispersion": 0.2,
    "donor_sigma": 0.5,
})

#: Permutation counts for the synthetic gate. Not pre-registered (spec §4 pins n_perm = 1000 and
#: n_perm_pb >= 200 for the REAL sweep); the synthetic gate runs fewer and reports its Monte-Carlo
#: error so the resolution limit is visible rather than assumed away.
N_PERM = 40
N_PERM_PB = 40

INSTRUMENT_SANITY = MappingProxyType({
    "instrument_lambda_naive_min": INSTRUMENT_LAMBDA_NAIVE_MIN,
    "instrument_naive_floor_frac_min": INSTRUMENT_NAIVE_FLOOR_FRAC_MIN,
    "instrument_floor_ratio_min": INSTRUMENT_FLOOR_RATIO_MIN,
    "n_perm": N_PERM,
    "n_perm_pb": N_PERM_PB,
    "oracle_sim": dict(ORACLE_SIM),
})


def manifest() -> dict:
    """Both blocks, tagged by provenance, for embedding in a run artifact.

    An artifact that records only the numbers cannot be audited against the spec later; one that
    records which numbers were frozen and which were not, can.
    """
    return {
        "pre_registered": dict(PRE_REGISTERED),
        "pre_registered_source": "docs/PHASE0_SPEC.md (frozen; changes go through docs/AMENDMENTS.md)",
        "instrument_sanity": dict(INSTRUMENT_SANITY),
        "instrument_sanity_source": "scripts/synthetic_gate.py — NOT pre-registered, not the decision rule",
    }
