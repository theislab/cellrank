from typing import Mapping, Sequence
from typing_extensions import Literal


class MarkerGenes:
    """Helper class for :class:`cellrank.external.kernels.WOT` to get proliferation/apoptosis genes."""

    # fmt: off
    _proliferation_markers: Mapping[str, Sequence[str]] = {
        # 97
        "human": ("ANLN", "ANP32E", "ATAD2", "AURKA", "AURKB", "BIRC5", "BLM",
                  "BRIP1", "BUB1", "CASP8AP2", "CBX5", "CCNB2", "CCNE2", "CDC20",
                  "CDC25C", "CDC45", "CDC6", "CDCA2", "CDCA3", "CDCA7", "CDCA8",
                  "CDK1", "CENPA", "CENPE", "CENPF", "CHAF1B", "CKAP2", "CKAP2L",
                  "CKAP5", "CKS1B", "CKS2", "CLSPN", "CTCF", "DLGAP5", "DSCC1",
                  "DTL", "E2F8", "ECT2", "EXO1", "FAM64A", "FEN1", "G2E3", "GAS2L3",
                  "GINS2", "GMNN", "GTSE1", "HELLS", "HJURP", "HMGB2", "HMMR", "HN1",
                  "KIF11", "KIF20B", "KIF23", "KIF2C", "LBR", "MCM2", "MCM4", "MCM5",
                  "MCM6", "MKI67", "MLF1IP", "MSH2", "NASP", "NCAPD2", "NDC80",
                  "NEK2", "NUF2", "NUSAP1", "PCNA", "POLA1", "POLD3", "PRIM1",
                  "PSRC1", "RAD51", "RAD51AP1", "RANGAP1", "RFC2", "RPA2", "RRM1",
                  "RRM2", "SLBP", "SMC4", "TACC3", "TIPIN", "TMPO", "TOP2A", "TPX2",
                  "TTK", "TUBB4B", "TYMS", "UBE2C", "UBR7", "UHRF1", "UNG", "USP1",
                  "WDR76"),
        # 98
        "mouse": ("Mcm4", "Smc4", "Gtse1", "Ttk", "Rangap1", "Ccnb2", "Cenpa",
                  "Cenpe", "Cdca8", "Ckap2", "Rad51", "Pcna", "Ube2c", "Lbr",
                  "Cenpf", "Birc5", "Dtl", "Dscc1", "Cbx5", "Usp1", "Hmmr", "Wdr76",
                  "Ung", "Hn1", "Cks2", "Kif20b", "Cdk1", "Slbp", "Aurkb", "Kif11",
                  "Cks1b", "Blm", "Msh2", "Gas2l3", "Tyms", "Hjurp", "Hells",
                  "Prim1", "Uhrf1", "Ndc80", "Mcm6", "Rrm1", "Mlf1ip", "Top2a",
                  "Hmgb2", "Ccne2", "G2e3", "Tmpo", "Nusap1", "Ncapd2", "Mcm2",
                  "Kif2c", "Cdca2", "Nasp", "Gmnn", "Cdc6", "Pold3", "Ckap2l",
                  "Fam64a", "Ubr7", "Fen1", "Bub1", "Brip1", "Atad2", "Psrc1",
                  "Rrm2", "Tipin", "Casp8ap2", "Tubb4b", "Kif23", "Exo1", "Rfc2",
                  "Pola1", "Mki67", "Tpx2", "Aurka", "Anln", "Chaf1b", "Hjurp",
                  "Tacc3", "Mcm5", "Anp32e", "Dlgap5", "Ect2", "Nuf2", "Cdc45",
                  "Ckap5", "Ctcf", "Clspn", "Cdca7", "Cdca3", "Rpa2", "Gins2",
                  "E2f8", "Cdc25c", "Nek2", "Cdc20", "Rad51ap1"),
    }
    _apoptosis_markers: Mapping[str, Sequence[str]] = {
        # 161
        "human":
            ("ADD1", "AIFM3", "ANKH", "ANXA1", "APP", "ATF3", "AVPR1A", "BAX",
             "BCAP31", "BCL10", "BCL2L1", "BCL2L10", "BCL2L11", "BCL2L2", "BGN",
             "BID", "BIK", "BIRC3", "BMF", "BMP2", "BNIP3L", "BRCA1", "BTG2",
             "BTG3", "CASP1", "CASP2", "CASP3", "CASP4", "CASP6", "CASP7",
             "CASP8", "CASP9", "CAV1", "CCNA1", "CCND1", "CCND2", "CD14", "CD2",
             "CD38", "CD44", "CD69", "CDC25B", "CDK2", "CDKN1A", "CDKN1B",
             "CFLAR", "CLU", "CREBBP", "CTH", "CTNNB1", "CYLD", "DAP", "DAP3",
             "DCN", "DDIT3", "DFFA", "DIABLO", "DNAJA1", "DNAJC3", "DNM1L",
             "DPYD", "EBP", "EGR3", "EMP1", "ENO2", "ERBB2", "ERBB3", "EREG",
             "ETF1", "F2", "F2R", "FAS", "FASLG", "FDXR", "FEZ1", "GADD45A",
             "GADD45B", "GCH1", "GNA15", "GPX1", "GPX3", "GPX4", "GSN", "GSR",
             "GSTM1", "GUCY2D", "H1-0", "HGF", "HMGB2", "HMOX1", "HSPB1",
             "IER3", "IFITM3", "IFNB1", "IFNGR1", "IGF2R", "IGFBP6", "IL18",
             "IL1A", "IL1B", "IL6", "IRF1", "ISG20", "JUN", "KRT18", "LEF1",
             "LGALS3", "LMNA", "LUM", "MADD", "MCL1", "MGMT", "MMP2", "NEDD9",
             "NEFH", "PAK1", "PDCD4", "PDGFRB", "PEA15", "PLAT", "PLCB2",
             "PLPPR4", "PMAIP1", "PPP2R5B", "PPP3R1", "PPT1", "PRF1", "PSEN1",
             "PSEN2", "PTK2", "RARA", "RELA", "RETSAT", "RHOB", "RHOT2",
             "RNASEL", "ROCK1", "SAT1", "SATB1", "SC5D", "SLC20A1", "SMAD7",
             "SOD1", "SOD2", "SPTAN1", "SQSTM1", "TAP1", "TGFB2", "TGFBR3",
             "TIMP1", "TIMP2", "TIMP3", "TNF", "TNFRSF12A", "TNFSF10", "TOP2A",
             "TSPO", "TXNIP", "VDAC2", "WEE1", "XIAP"),
        # 193
        "mouse": ("Ercc5", "Serpinb5", "Inhbb", "Steap3", "Btg2", "Phlda3", "Tnni1",
                  "Rgs16", "Ier5", "Slc19a2", "Adck3", "Ephx1", "Ptpn14", "Atf3",
                  "Notch1", "Rxra", "Ralgds", "Ak1", "Stom", "Ddb2", "Cd82", "Il1a",
                  "Pcna", "Bmp2", "Trib3", "Procr", "Blcap", "Ada", "Fgf13", "Irak1",
                  "Tspyl2", "Sat1", "Zmat3", "Hspa4l", "Slc7a11", "Tm4sf1", "Rap2b",
                  "Fbxw7", "S100a4", "S100a10", "Txnip", "Nhlh2", "Dnttip2", "Clca2",
                  "Wwp1", "Klf4", "Ikbkap", "Cdkn2a", "Cdkn2b", "Jun", "Slc35d1",
                  "Plk3", "Rnf19b", "Sfn", "Fuca1", "Epha2", "Wrap73", "Mxd4",
                  "Rchy1", "Iscu", "Triap1", "Prkab1", "Trafd1", "Pom121", "Pdgfa",
                  "Gadd45a", "Vamp8", "Retsat", "Tprkb", "Tgfa", "Mxd1", "Sec61a1",
                  "Xpc", "Ccnd2", "H2afj", "Ldhb", "Lrmp", "Tm7sf3", "Tgfb1",
                  "Sertad3", "Cebpa", "Klk8", "Bax", "Ppp1r15a", "Rpl18", "Aen",
                  "Rrp8", "Ccp110", "Nupr1", "Ptpre", "Hras", "Eps8l2", "Ctsd",
                  "Cd81", "Perp", "Rps12", "Tpd52l1", "Sesn1", "Foxo3", "Ddit4",
                  "Zfp365", "Prmt2", "Mknk2", "Dram1", "Apaf1", "Btg1", "Mdm2",
                  "Ddit3", "Gls2", "Dgka", "Cdkn2aip", "Hmox1", "Rrad", "Cdh13",
                  "Osgin1", "Cgrrf1", "Abhd4", "Kif13b", "Rb1", "Nudt15", "Tsc22d1",
                  "Casp1", "St14", "Ei24", "Vwa5a", "Zbtb16", "Rps27l", "Mapkapk3",
                  "Ip6k2", "Tcn2", "Lif", "Upp1", "Ccng1", "Cyfip2", "Gnb2l1",
                  "Hint1", "Gm2a", "Hist3h2a", "Alox8", "Trp53", "Tax1bp3", "Traf4",
                  "Cdk5r1", "Ppm1d", "Rad51c", "Tob1", "Krt17", "Hexim1", "Fdxr",
                  "Itgb4", "Sphk1", "Rhbdf2", "Baiap2", "Dcxr", "Hist1h1c", "Ninj1",
                  "Nol8", "F2r", "Ankra2", "Plk2", "Sdc1", "Gpx2", "Zfp36l1", "Fos",
                  "Ccnk", "Jag2", "Ndrg1", "Pmm1", "Plxnb2", "Vdr", "Csrnp2",
                  "Acvr1b", "Sp1", "Abat", "Socs1", "Abcc5", "Trp63", "Fam162a",
                  "App", "Rab40c", "Bak1", "Def6", "Cdkn1a", "Tap1", "Ier3", "Polh",
                  "Ccnd3", "Hbegf", "Hdac3", "Rad9a", "Ctsf", "Slc3a2", "Fas"),
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
