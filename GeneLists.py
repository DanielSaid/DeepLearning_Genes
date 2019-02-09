import numpy as np

EASY = ['ACACA', 'CPT2', 'CPT1A', 'NR1H4', 'NFKB2', 'RELA', 'RELB', 'PPARA', 'SCD', 'TGFB1']
DIFF = ['MAPK10', 'RXRG', 'PTEN', 'IL1A', 'ELOVL3', 'PMP2', 'GCKR', 'AKT3', 'MLXIPL', 'IL10']
CONN = ['AKT1', 'MTOR', 'RELA', 'NFKB2', 'TNF', 'PPARG', 'PTEN', 'ADIPOQ', 'PPARA', 'TGFB1']
COMBINED = list(np.unique(EASY + DIFF + CONN))
COMBINED_NO_TNF = list(np.unique(EASY + DIFF + ['AKT1', 'MTOR', 'RELA', 'NFKB2', 'PPARG', 'PTEN', 'ADIPOQ', 'PPARA', 'TGFB1']))

STEATOSIS = ['FABP4', 'ACACA', 'AKT1', 'AKT2', 'AKT3', 'PRKAA1', 'PRKAA2', 'ADIPOR1', 'ADIPOR2', 'ADIPOQ', 'BCL2A1', 'CPT2', 
				'CPT1A', 'CPT1C', 'CASP8', 'MLXIPL', 'FABP5', 'ELOVL3', 'FAS', 'FOXO1', 'NR1H4', 'RXRA', 'FASLG', 'FABP3', 'FABP7', 
				'PMP2', 'GCKR', 'IL1A', 'IL10', 'IRS1', 'IRS2', 'MAPK10', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'PPARA', 'PPARG', 'PTEN', 
				'RXRB', 'RXRG', 'SCD', 'SOCS3', 'SREBF1', 'TGFB1', 'TGFB2', 'TGFB3', 'TLR4', 'PNPLA3', 'MTOR']

# Genes are listed as they are called in human in vitro data. Rat version, according to excel sheet, is 
# the lower-case capitalized version of this, unless specified otherwise.

# ZNF767P is pseudogene, only present in human. Therefore removed from list.
GTX = [
	'CEACAM1', 'CLCN4', 'EML1', 'PWWP2B', 'UBE2E2', 'USP13', 'GMFG', 'PROSC', 'TTR', 'NR0B2', 'NAT8',
	'RBPMS', 'TBC1D9', 'SNX11', 'BCOR', 'ROBO2', 'DENND6B', 'APOM', 'NRIP3', 'PITHD1', 'AVEN', 'ZNRF3', 
	'BEAN1', 'SLC27A1', 'ANXA6', 'APOA4', 'BTD', 'EIF2D', 'AGFG1', 'NDUFA10', 'NFATC3', 'PLAA', 'FAN1', 'SLC40A1',
	'ANAPC5', 'MRPS5', 'GSTK1', 'HOGA1', 'FGA', 'SGK1', 'SLC6A4', 'SCRN2', 'CC2D1B'
]

# For rats, SLCO1B1 is Slco1b2, SLCO1B3 is Slco1b2, UGT2B4 is Ugt2b35 and CYP3A4 is Cyp3a18
CHOLESTASIS = [
	'HNF4A', 'SLC10A1', 'SLCO1B1', 'CYP7A1', 'CYP8B1', 'CYP27A1', 'CYP7B1', 'NR1H4', 'NR0B2', 
	'NR1I2', 'NR1I3', 'FGF19', 'ABCB11', 'SLC51A', 'ABCC3', 'UGT2B4', 'CYP3A4', 'SULT2A1'
]  # Dan: how to deal with SLCO1B3 when predicting to human (it maps to same gene as SLCO1B1)? I removed it!

NAFLD = [
	'PPARGC1A', 'IL6', 'SERPINE1', 'IL1B', 'STAT3', 'TCF7L2', 'CD14', 'IRS1', 'PNPLA3', 'PEMT', 
	'TM6SF2', 'SREBF1', 'HFE', 'SAMM50', 'FDFT1', 'NR1I2', 'PPARA', 'PPP1R3B', 'CHDH', 'LYPLAL1', 'SOD2', 'LEPR'
]

# For rats, FCN3 is Fibcd1
# Fibcd1 could not be found in rats vitro, so removed FCN3
# E2F1 could not be found in rats vitro, so removed
# MT2A could not be found in rats, so removed
CARCINOGENIC = [
	'GPC3', 'MDK', 'COL5A2', 'TP53BP2', 'XPO1', 'AFP', 'CCNA2', 'CCNE1', 'COL1A1', 'COL4A1', 'CTNNB1',
	'FBN1', 'FOXM1', 'STMN1', 'LGALS3BP', 'MARCKS', 'NME1', 'NRAS', 'PGK1', 'MAPK3', 'SMARCC1', 
	'COPS5', 'PEG10', 'HGFAC', 'IGFALS', 'LCAT', 'SLC22A1', 'ACADS', 'ACADVL', 'C9', 'DSG2', 'PLG', 'HAMP'
]


GTX_CAR = list(np.unique(GTX + CARCINOGENIC))
