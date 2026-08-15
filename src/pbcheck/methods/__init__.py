"""DE methods used by the pbcheck auditor: the naive per-cell test and the correct pseudobulk test.

The pseudobulk arm's test is moderated eBayes as of Amendment 2 Change 1; ``deseq_from_pdata`` is
superseded and retained only to reproduce Amendment-1-era numbers.
"""

from pbcheck.methods.de import DEResult
from pbcheck.methods.moderated import ebayes_from_pdata
from pbcheck.methods.naive import naive_de
from pbcheck.methods.pseudobulk import pseudobulk_de

__all__ = ["DEResult", "ebayes_from_pdata", "naive_de", "pseudobulk_de"]
