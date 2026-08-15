"""The gate's thresholds, in one place, with each one's provenance attached.

**Everything in :data:`PRE_REGISTERED` is frozen by ``docs/PHASE0_SPEC.md`` and the amendments to
it. Changing any of these values is a protocol change and requires a dated, numbered entry in
``docs/AMENDMENTS.md``, written and committed BEFORE the code that applies it.** That is not a style
preference: Amendment 1 exists
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
#: oracle in gate code is prohibited by Amendment 1 Change 2. **Amendment 3 does not touch these**
#: — it re-scopes the *region* over which POWER_TARGET binds, never the effect it binds at.
ORACLE_LOG2FC = 1.0
ORACLE_K = 200

#: The two regimes at which the two halves of the validity gate are evaluated (Amendment 3
#: Change 1). Before Amendment 3 both halves ran at a single ``donor_sigma`` that the frozen spec
#: never pinned — so the power criterion was being judged at an arbitrary simulator setting where
#: the committed grid shows NO test reaches POWER_TARGET at any donor count it tested. Amendment 3
#: separates them:
#:
#: * **Calibration** (λ band, perm-null FP rate) is evaluated at the HARD regime, ``sigma_donor``
#:   0.5 — the worst donor variance in the grid's fully-populated tier and the conservative choice
#:   for a false-positive criterion, since the null gets harder as donor variance rises. Unchanged
#:   from every prior run; nothing about calibration is relaxed.
#: * **Power** is evaluated at the ENVELOPE BOUNDARY, ``sigma_donor`` 0.35 — the lowest-σ /
#:   smallest-n point at which the selected ``ebayes`` arm is grid-shown to clear POWER_TARGET
#:   (measured 0.793 at 8v8, calibrated: FP 0.043, λ 1.010), coinciding with Amendment 1's analytic
#:   n*(0.35) = 8.
#:
#: This NARROWS what the instrument claims, and lowers no bar: power 0.60 at ``sigma_donor`` 0.5
#: with 8 donors/group remains UNMET (measured 0.35) and is NOT claimed by anything here.
CALIBRATION_EVAL_SIGMA = 0.5
POWER_EVAL_SIGMA = 0.35

#: The declared operating envelope (Amendment 3 Change 1): the pseudobulk arm is valid only for
#: strata whose (``sigma_donor``, donors-per-group) lies inside it, i.e. where the selected test's
#: power at the UNCHANGED pre-registered oracle is ≥ POWER_TARGET. Outside it the arm is not
#: declared valid and no stratum result may be reported as a pbcheck measurement.
#:
#: ``min_donors_per_group`` is Amendment 1's power frontier (derived, then validated numerically to
#: |error| < 0.033); ``grid_support`` records what the committed 146-cell grid
#: (``pilot/testsel/summary.json``, ``72dec7b``) independently says about the same point, so the
#: envelope can be audited against frozen data rather than taken on the derivation's word.
#:
#: Every number here is SYNTHETIC. ``sigma_donor`` is a free knob of our simulator and is still not
#: anchored to any real dataset — Amendment 1's closing concern, restated by Amendments 2 and 3 and
#: still OPEN. Whether any real stratum falls inside this envelope is unknown.
OPERATING_ENVELOPE = (
    MappingProxyType({
        "sigma_donor": 0.2,
        "min_donors_per_group": 4,
        "grid_support": "not in the grid; Amendment 1 frontier only",
    }),
    MappingProxyType({
        "sigma_donor": 0.35,
        "min_donors_per_group": 8,
        "grid_support": "ebayes power 0.793 at 8v8 (calibrated) -> n* <= 8",
    }),
    MappingProxyType({
        "sigma_donor": 0.5,
        "min_donors_per_group": 13,
        "grid_support": "ebayes power 0.486 at 12v12, the largest n tested -> n* > 12",
    }),
    MappingProxyType({
        "sigma_donor": 0.7,
        "min_donors_per_group": 23,
        "grid_support": "ebayes power 0.003 at 8v8 -> n* far above 8",
    }),
)

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
    "calibration_eval_sigma": CALIBRATION_EVAL_SIGMA,
    "power_eval_sigma": POWER_EVAL_SIGMA,
    "operating_envelope": tuple(dict(row) for row in OPERATING_ENVELOPE),
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
#: of the simulator and is **not** anchored to real data; Amendments 1, 2 and 3 all close on that as
#: the outstanding threat to every power number computed here.
#:
#: This is the CALIBRATION regime (``donor_sigma`` = :data:`CALIBRATION_EVAL_SIGMA`, the hard point).
#: Since Amendment 3 Change 1 the power oracle is run at :data:`POWER_EVAL_SIGMA` instead, and the
#: gate labels which arm ran at which sigma rather than presenting one number for both.
ORACLE_SIM = MappingProxyType({
    "n_genes": 1500,
    "n_donors_per_group": 8,
    "n_cells_per_donor": 250,
    "dispersion": 0.2,
    "donor_sigma": CALIBRATION_EVAL_SIGMA,
})

#: Permutation counts for the synthetic gate. Not pre-registered (spec §4 pins n_perm = 1000 and
#: n_perm_pb >= 200 for the REAL sweep); the synthetic gate runs fewer and reports its Monte-Carlo
#: error so the resolution limit is visible rather than assumed away.
#:
#: **40 -> 200 by Amendment 3 Change 2**, declared before the rerun that applies it. At 40 paired
#: permutations the FP criterion's Monte-Carlo SE at p = 0.05 is sqrt(.05*.95/40) = 0.034, so the
#: 2/40 reading of the 2026-08-15 run could not be told apart from a true rate of ~0.12; at 200 it
#: is 0.015, which resolves it. The moderated arm makes this affordable (~5 ms/fit against the
#: retired DESeq2-Wald arm's ~2.8 s, where 200 permutations would have cost hours).
#:
#: **Both constants move together, and must.** ``permutation.run_null`` pairs
#: ``min(n_perm, n_perm_pb)`` permutations, so raising ``N_PERM_PB`` alone would leave the paired
#: count — and therefore the FP rate and its MC SE — pinned at ``N_PERM``, producing a run that
#: looked like 200 permutations and was not.
N_PERM = 200
N_PERM_PB = 200

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
