from typing import Mapping, Sequence

from typing_extensions import Literal


class MarkerGenes:
    """Helper class for :class:`cellrank.external.kernels.WOT` to get proliferation/apoptosis genes."""

    # fmt: off
    _proliferation_markers: Mapping[str, Sequence[str]] = {
        # 97
        "human": ('ANLN', 'ANP32E', 'ATAD2', 'AURKA', 'AURKB', 'BIRC5', 'BLM',
                  'BRIP1', 'BUB1', 'CASP8AP2', 'CBX5', 'CCNB2', 'CCNE2', 'CDC20',
                  'CDC25C', 'CDC45', 'CDC6', 'CDCA2', 'CDCA3', 'CDCA7', 'CDCA8',
                  'CDK1', 'CENPA', 'CENPE', 'CENPF', 'CHAF1B', 'CKAP2', 'CKAP2L',
                  'CKAP5', 'CKS1B', 'CKS2', 'CLSPN', 'CTCF', 'DLGAP5', 'DSCC1',
                  'DTL', 'E2F8', 'ECT2', 'EXO1', 'FAM64A', 'FEN1', 'G2E3', 'GAS2L3',
                  'GINS2', 'GMNN', 'GTSE1', 'HELLS', 'HJURP', 'HMGB2', 'HMMR', 'HN1',
                  'KIF11', 'KIF20B', 'KIF23', 'KIF2C', 'LBR', 'MCM2', 'MCM4', 'MCM5',
                  'MCM6', 'MKI67', 'MLF1IP', 'MSH2', 'NASP', 'NCAPD2', 'NDC80',
                  'NEK2', 'NUF2', 'NUSAP1', 'PCNA', 'POLA1', 'POLD3', 'PRIM1',
                  'PSRC1', 'RAD51', 'RAD51AP1', 'RANGAP1', 'RFC2', 'RPA2', 'RRM1',
                  'RRM2', 'SLBP', 'SMC4', 'TACC3', 'TIPIN', 'TMPO', 'TOP2A', 'TPX2',
                  'TTK', 'TUBB4B', 'TYMS', 'UBE2C', 'UBR7', 'UHRF1', 'UNG', 'USP1',
                  'WDR76'),
        "mouse": (),
    }
    _apoptosis_markers: Mapping[str, Sequence[str]] = {
        # 161
        "human":
            ('ADD1', 'AIFM3', 'ANKH', 'ANXA1', 'APP', 'ATF3', 'AVPR1A', 'BAX',
             'BCAP31', 'BCL10', 'BCL2L1', 'BCL2L10', 'BCL2L11', 'BCL2L2', 'BGN',
             'BID', 'BIK', 'BIRC3', 'BMF', 'BMP2', 'BNIP3L', 'BRCA1', 'BTG2',
             'BTG3', 'CASP1', 'CASP2', 'CASP3', 'CASP4', 'CASP6', 'CASP7',
             'CASP8', 'CASP9', 'CAV1', 'CCNA1', 'CCND1', 'CCND2', 'CD14', 'CD2',
             'CD38', 'CD44', 'CD69', 'CDC25B', 'CDK2', 'CDKN1A', 'CDKN1B',
             'CFLAR', 'CLU', 'CREBBP', 'CTH', 'CTNNB1', 'CYLD', 'DAP', 'DAP3',
             'DCN', 'DDIT3', 'DFFA', 'DIABLO', 'DNAJA1', 'DNAJC3', 'DNM1L',
             'DPYD', 'EBP', 'EGR3', 'EMP1', 'ENO2', 'ERBB2', 'ERBB3', 'EREG',
             'ETF1', 'F2', 'F2R', 'FAS', 'FASLG', 'FDXR', 'FEZ1', 'GADD45A',
             'GADD45B', 'GCH1', 'GNA15', 'GPX1', 'GPX3', 'GPX4', 'GSN', 'GSR',
             'GSTM1', 'GUCY2D', 'H1-0', 'HGF', 'HMGB2', 'HMOX1', 'HSPB1',
             'IER3', 'IFITM3', 'IFNB1', 'IFNGR1', 'IGF2R', 'IGFBP6', 'IL18',
             'IL1A', 'IL1B', 'IL6', 'IRF1', 'ISG20', 'JUN', 'KRT18', 'LEF1',
             'LGALS3', 'LMNA', 'LUM', 'MADD', 'MCL1', 'MGMT', 'MMP2', 'NEDD9',
             'NEFH', 'PAK1', 'PDCD4', 'PDGFRB', 'PEA15', 'PLAT', 'PLCB2',
             'PLPPR4', 'PMAIP1', 'PPP2R5B', 'PPP3R1', 'PPT1', 'PRF1', 'PSEN1',
             'PSEN2', 'PTK2', 'RARA', 'RELA', 'RETSAT', 'RHOB', 'RHOT2',
             'RNASEL', 'ROCK1', 'SAT1', 'SATB1', 'SC5D', 'SLC20A1', 'SMAD7',
             'SOD1', 'SOD2', 'SPTAN1', 'SQSTM1', 'TAP1', 'TGFB2', 'TGFBR3',
             'TIMP1', 'TIMP2', 'TIMP3', 'TNF', 'TNFRSF12A', 'TNFSF10', 'TOP2A',
             'TSPO', 'TXNIP', 'VDAC2', 'WEE1', 'XIAP'),
        "mouse": (),
    }
    # fmt: on

    @classmethod
    def proliferation_markers(
        cls, organism: Literal["human", "mouse"]
    ) -> Sequence[str]:
        """Get proliferation markers for ``organism``."""
        try:
            return cls._proliferation_markers[organism]
        except KeyError:
            raise NotImplementedError(
                f"Proliferation markers for `{organism}` are not yet implemented."
            ) from None

    @classmethod
    def apoptosis_markers(cls, organism: Literal["human", "mouse"]) -> Sequence[str]:
        """Get apoptosis markers for ``organism``."""
        try:
            return cls._apoptosis_markers[organism]
        except KeyError:
            raise NotImplementedError(
                f"Apoptosis markers for `{organism}` are not yet implemented."
            ) from None
