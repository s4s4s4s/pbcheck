"""Design auditor — the metadata-only half of pbcheck.

Before running any DE, most pseudoreplication risk is already visible in ``.obs``: how many donors
back each condition, whether donors nest cleanly within condition, whether a batch is perfectly
confounded with condition, and how imbalanced the groups are. This module reads only ``.obs`` and
emits a structured :class:`DesignReport` with human-readable flags. It never needs raw counts, which
is what lets pbcheck audit a published study from its cell metadata alone.

The thresholds encode the field's rule of thumb (Squair et al. 2021; Zimmerman et al. 2021): fewer
than ~3 donors per group makes donor-level inference essentially impossible, so a per-cell test there
is not merely inflated — it is uninterpretable.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class DesignReport:
    condition_col: str
    donor_col: str
    n_cells: int
    groups: dict[str, int]  # condition level -> n cells
    donors_per_group: dict[str, int]  # condition level -> n distinct donors
    cells_per_donor: dict[str, float]  # summary stats
    donor_nests_in_condition: bool  # each donor belongs to exactly one condition
    min_donors_per_group: int
    imbalance_ratio: float  # max/min group cell count
    batch_confounded: dict[str, float] = field(default_factory=dict)  # batch col -> donor-level Cramer's V
    batch_separates_condition: dict[str, bool] = field(default_factory=dict)  # batch col -> perfect separation
    flags: list[str] = field(default_factory=list)
    min_donors: int = 3  # threshold the verdict below is judged against; set by audit_design

    @property
    def usable_for_pseudobulk(self) -> bool:
        """Whether a donor-level pseudobulk DE is even meaningful for this design.

        The threshold defaults to 3 donors per group, matching the spec's inclusion gate (§1, item
        1) — not 2 — but tracks whatever ``min_donors`` the caller passed to :func:`audit_design`,
        so a stricter caller-chosen threshold actually changes the verdict and not just the flag
        text (it used to hardcode 3 here regardless of ``min_donors``). A design perfectly separated
        by a batch/assay/pool covariate is excluded outright: its inflation cannot be attributed to
        pseudoreplication rather than to the covariate.
        """
        return (
            self.donor_nests_in_condition
            and self.min_donors_per_group >= self.min_donors
            and not any(self.batch_separates_condition.values())
        )


def _cramers_v(a: pd.Series, b: pd.Series) -> float:
    """Association between two categoricals in [0, 1]; 1.0 == perfectly confounded.

    The caller is responsible for passing series at the correct unit of observation. For
    condition-vs-batch confounding that unit is the DONOR, never the cell — see
    :func:`_donor_level_table`.
    """
    tab = pd.crosstab(a, b)
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    chi2 = _chi2(tab.to_numpy())
    n = tab.to_numpy().sum()
    r, k = tab.shape
    denom = n * (min(r, k) - 1)
    return float(np.sqrt(chi2 / denom)) if denom > 0 else float("nan")


def _donor_level_table(obs: pd.DataFrame, donor_col: str, cols: list[str]) -> pd.DataFrame:
    """One row per (donor, *cols) combination — each donor counted once, not once per cell.

    Condition, batch, assay, suspension and pool identity are donor-level attributes. Measuring
    their association with condition over CELLS lets a donor with many cells dominate the estimate,
    which is the very cells-as-replicates fallacy this package exists to detect: a batch that is
    100% one condition can score near zero simply because the donors carrying it are small.
    Deduplicating first makes every donor contribute equally.
    """
    return obs[[donor_col, *cols]].astype(str).drop_duplicates()


def _perfectly_separates(a: pd.Series, b: pd.Series) -> bool:
    """True if knowing ``a`` determines ``b`` — every level of a maps to exactly one level of b.

    Checked explicitly rather than inferred from Cramer's V: V can sit below 1 on an unbalanced
    table that is nonetheless perfectly separating, and a perfectly separated design makes the
    inflation measurement uninterpretable rather than merely noisy.
    """
    tab = pd.crosstab(a, b)
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return False
    return bool(((tab > 0).sum(axis=1) == 1).all())


def _chi2(obs: np.ndarray) -> float:
    row = obs.sum(1, keepdims=True)
    col = obs.sum(0, keepdims=True)
    exp = row @ col / obs.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(exp > 0, (obs - exp) ** 2 / exp, 0.0)
    return float(terms.sum())


def audit_design(
    adata,
    *,
    condition_col: str = "condition",
    donor_col: str = "donor",
    batch_cols: list[str] | None = None,
    min_donors: int = 3,
) -> DesignReport:
    """Inspect ``adata.obs`` for the ingredients of pseudoreplication. Reads no counts."""
    obs = adata.obs
    for col in (condition_col, donor_col):
        if col not in obs.columns:
            raise ValueError(f"obs is missing required column '{col}'")

    groups = obs[condition_col].value_counts().to_dict()
    groups = {str(k): int(v) for k, v in groups.items()}

    donor_condition = obs[[donor_col, condition_col]].astype(str).drop_duplicates()
    donors_per_cond = donor_condition.groupby(condition_col)[donor_col].nunique().to_dict()
    donors_per_cond = {str(k): int(v) for k, v in donors_per_cond.items()}

    # A donor "nests" if it appears under exactly one condition.
    donor_ncond = donor_condition.groupby(donor_col)[condition_col].nunique()
    donor_nests = bool((donor_ncond == 1).all())

    cpd = obs.groupby(donor_col, observed=True).size()
    cells_per_donor = {
        "min": int(cpd.min()), "median": float(cpd.median()),
        "max": int(cpd.max()), "n_donors": int(cpd.size),
    }

    min_dpg = min(donors_per_cond.values()) if donors_per_cond else 0
    imbalance = (max(groups.values()) / max(min(groups.values()), 1)) if groups else float("nan")

    # Confounding is assessed at the DONOR level (see _donor_level_table): batch, assay, suspension
    # and pool are donor-constant attributes, and weighting them by cell count hides real confounds.
    batch_conf: dict[str, float] = {}
    batch_separates: dict[str, bool] = {}
    for bc in batch_cols or []:
        if bc in obs.columns:
            dtab = _donor_level_table(obs, donor_col, [bc, condition_col])
            batch_conf[bc] = round(_cramers_v(dtab[bc], dtab[condition_col]), 3)
            batch_separates[bc] = _perfectly_separates(dtab[bc], dtab[condition_col])

    flags: list[str] = []
    if not donor_nests:
        flags.append("donor spans multiple conditions — not a standard case/control design")
    if min_dpg < min_donors:
        flags.append(
            f"only {min_dpg} donor(s) in the smallest group (< {min_donors}); donor-level "
            "inference is unreliable and per-cell DE here is uninterpretable, not just inflated"
        )
    if len(groups) >= 2 and imbalance >= 5:
        flags.append(f"group cell-count imbalance {imbalance:.1f}x")
    for bc, v in batch_conf.items():
        if batch_separates.get(bc):
            flags.append(
                f"batch '{bc}' PERFECTLY separates condition at the donor level — the inflation "
                "measurement is uninterpretable for this design, not merely confounded"
            )
        elif not np.isnan(v) and v >= 0.8:
            flags.append(f"batch '{bc}' is nearly confounded with condition (donor-level Cramer's V={v})")

    return DesignReport(
        condition_col=condition_col,
        donor_col=donor_col,
        n_cells=int(obs.shape[0]),
        groups=groups,
        donors_per_group=donors_per_cond,
        cells_per_donor=cells_per_donor,
        donor_nests_in_condition=donor_nests,
        min_donors_per_group=min_dpg,
        imbalance_ratio=float(imbalance),
        batch_confounded=batch_conf,
        batch_separates_condition=batch_separates,
        flags=flags,
        min_donors=min_donors,
    )
