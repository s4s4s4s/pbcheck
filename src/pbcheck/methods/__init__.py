"""DE methods used by the pbcheck auditor: the naive per-cell test and the correct pseudobulk test."""

from pbcheck.methods.de import DEResult
from pbcheck.methods.naive import naive_de
from pbcheck.methods.pseudobulk import pseudobulk_de

__all__ = ["DEResult", "naive_de", "pseudobulk_de"]
