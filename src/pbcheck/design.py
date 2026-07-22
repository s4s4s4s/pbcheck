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
    batch_confounded: dict[str, float] = field(default_factory=dict)  # batch col -> Cramer's V
    flags: list[str] = field(default_factory=list)

    @property
    def usable_for_pseudobulk(self) -> bool:
        """Whether a donor-level pseudobulk DE is even meaningful for this design."""
        return self.donor_nests_in_condition and self.min_donors_per_group >= 2


def _cramers_v(a: pd.Series, b: pd.Series) -> float:
    """Association between two categoricals in [0, 1]; 1.0 == perfectly confounded."""
    tab = pd.crosstab(a, b)
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return float("nan")
    chi2 = _chi2(tab.to_numpy())
    n = tab.to_numpy().sum()
    r, k = tab.shape
    denom = n * (min(r, k) - 1)
    return float(np.sqrt(chi2 / denom)) if denom > 0 else float("nan")


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

    batch_conf = {}
    for bc in batch_cols or []:
        if bc in obs.columns:
            batch_conf[bc] = round(_cramers_v(obs[bc].astype(str), obs[condition_col].astype(str)), 3)

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
        if not np.isnan(v) and v >= 0.8:
            flags.append(f"batch '{bc}' is nearly confounded with condition (Cramer's V={v})")

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
        flags=flags,
    )
