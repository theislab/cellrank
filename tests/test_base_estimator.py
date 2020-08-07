from anndata import AnnData
from cellrank.tools import GPCCA
from cellrank.tools._constants import Direction, _transition


class TestLoad:
    def test_read_key_same_adata(self, adata_cflare: AnnData):
        g = GPCCA(adata_cflare, obsp_key=f"{_transition(Direction.FORWARD)}")

        assert g.adata is adata_cflare
