import numpy as np
import pickle
import sys

sys.path.append('/Users/Daniel/Desktop/MaCSBio/mrp2-python/src/')
from typing import List, Dict, Callable, Tuple
import matplotlib.pyplot as plt
from SmallParser import read_data_xyz
from Types import Data, Header, ColInfo

EASY = ['ACACA', 'CPT2', 'CPT1A', 'NR1H4', 'NFKB2', 'RELA', 'RELB', 'PPARA', 'SCD', 'TGFB1']
DIFF = ['MAPK10', 'RXRG', 'PTEN', 'IL1A', 'ELOVL3', 'PMP2', 'GCKR', 'AKT3', 'MLXIPL', 'IL10']
CONN = ['AKT1', 'MTOR', 'RELA', 'NFKB2', 'TNF', 'PPARG', 'PTEN', 'ADIPOQ', 'PPARA', 'TGFB1']
COMBINED = list(np.unique(EASY + DIFF + CONN))
COMBINED_NO_TNF = list(
    np.unique(EASY + DIFF + ['AKT1', 'MTOR', 'RELA', 'NFKB2', 'PPARG', 'PTEN', 'ADIPOQ', 'PPARA', 'TGFB1']))
# Dan: the following can also be loaded from GeneLists.py
SMALL_SET_50 = ['FABP4', 'ACACA', 'AKT1', 'AKT2', 'AKT3', 'PRKAA1', 'PRKAA2', 'ADIPOR1', 'ADIPOR2', 'ADIPOQ', 'BCL2A1',
                'CPT2',
                'CPT1A', 'CPT1C', 'CASP8', 'MLXIPL', 'FABP5', 'ELOVL3', 'FAS', 'FOXO1', 'NR1H4', 'RXRA', 'FASLG',
                'FABP3', 'FABP7',
                'PMP2', 'GCKR', 'IL1A', 'IL10', 'IRS1', 'IRS2', 'MAPK10', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'PPARA',
                'PPARG', 'PTEN',
                'RXRB', 'RXRG', 'SCD', 'SOCS3', 'SREBF1', 'TGFB1', 'TGFB2', 'TGFB3', 'TLR4', 'PNPLA3', 'MTOR']

LOW_ERROR = ['FABP4', 'ACACA', 'AKT1', 'AKT2', 'AKT3', 'PRKAA1', 'PRKAA2', 'ADIPOR1', 'ADIPOR2', 'ADIPOQ', 'BCL2A1',
             'CPT2',
             'CPT1A', 'CPT1C', 'CASP8', 'MLXIPL', 'FABP5', 'ELOVL3', 'FAS', 'FOXO1', 'NR1H4', 'RXRA', 'FASLG',
             'FABP3', 'FABP7',
             'PMP2', 'GCKR', 'IL1A', 'IL10', 'IRS1', 'IRS2', 'MAPK10', 'NFKB1', 'NFKB2', 'RELA', 'RELB', 'PPARA',
             'PPARG', 'PTEN',
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
# GTX = [
# 	'CEACAM1', 'CLCN4', 'EML1', 'PWWP2B', 'UBE2E2', 'USP13', 'GMFG', 'PROSC', 'TTR', 'NR0B2', 'NAT8',
# 	'RBPMS', 'TBC1D9', 'SNX11', 'BCOR', 'ROBO2', 'DENND6B', 'NRIP3', 'PITHD1', 'AVEN', 'ZNRF3',
# 	'BEAN1', 'SLC27A1', 'BTD', 'EIF2D',  'NDUFA10', 'NFATC3', 'PLAA', 'FAN1', 'SLC40A1',
# 	 'MRPS5', 'GSTK1', 'HOGA1', 'FGA', 'SGK1', 'SLC6A4', 'SCRN2', 'CC2D1B'
# ]
# For rats, SLCO1B1 is Slco1b2, SLCO1B3 is Slco1b2, UGT2B4 is Ugt2b35 and CYP3A4 is Cyp3a18
CHOLESTASIS = [
    'HNF4A', 'SLC10A1', 'SLCO1B1', 'CYP7A1', 'CYP8B1', 'CYP27A1', 'CYP7B1', 'NR1H4', 'NR0B2',
    'NR1I2', 'NR1I3', 'FGF19', 'ABCB11', 'SLC51A', 'ABCC3', 'SLCO1B3', 'UGT2B4', 'CYP3A4', 'SULT2A1'
]

NAFLD = [
    'PPARGC1A', 'IL6', 'SERPINE1', 'IL1B', 'STAT3', 'TCF7L2', 'CD14', 'IRS1', 'PNPLA3', 'PEMT',
    'TM6SF2', 'SREBF1', 'HFE', 'SAMM50', 'FDFT1', 'NR1I2', 'PPARA', 'PPP1R3B', 'CHDH', 'LYPLAL1', 'SOD2', 'LEPR'
]

# For rats, FCN3 is Fibcd1
# Fibcd1 could not be found in rats vitro, so removed FCN3
# E2F1 could not be found in rats vitro, so removed
# MT2A could not be found in rats, so removed
# CARCINOGENIC = [
# 	'GPC3', 'MDK', 'COL5A2', 'TP53BP2', 'XPO1','CCNA2', 'CCNE1', 'COL1A1', 'COL4A1', 'CTNNB1',
# 	'FBN1', 'FOXM1', 'STMN1', 'LGALS3BP', 'MARCKS', 'NME1', 'NRAS', 'PGK1', 'MAPK3', 'SMARCC1',
# 	'COPS5', 'PEG10', 'HGFAC', 'IGFALS', 'LCAT', 'SLC22A1', 'C9', 'DSG2', 'PLG', 'HAMP'
# ]
CARCINOGENIC = [
    'GPC3', 'MDK', 'COL5A2', 'TP53BP2', 'XPO1', 'AFP', 'CCNA2', 'CCNE1', 'COL1A1', 'COL4A1', 'CTNNB1',
    'FBN1', 'FOXM1', 'STMN1', 'LGALS3BP', 'MARCKS', 'NME1', 'NRAS', 'PGK1', 'MAPK3', 'SMARCC1',
    'COPS5', 'PEG10', 'HGFAC', 'IGFALS', 'LCAT', 'SLC22A1', 'ACADS', 'ACADVL', 'C9', 'DSG2', 'PLG', 'HAMP'
]

GENERAL_47 = [
    'APAP', 'ADP', 'APL', 'AA', 'AM', 'ASA', 'AZP', 'BBr', 'BBZ', 'CBZ', 'CPZ', 'CIM',
    'CFB', 'CMA', 'CPA', 'DZP', 'DFNa', 'ET', 'FP', 'FT', 'GFZ', 'GBC', 'GF', 'HPL', 'HCB', 'IM', 'INAH',
    'KC', 'LBT', 'LS', 'MP', 'MTS', 'ANIT', 'NFT', 'OPZ', 'PH', 'PB', 'PhB', 'PHE', 'PTU', 'RIF', 'SS', 'TC',
    'TAA', 'TRZ', 'VPA', 'WY'
]

GTX_CAR = list(np.unique(GTX + CARCINOGENIC))


class ProcessData():
    def __init__(self):
        # self.source_x_r_vitr_original = None
        # self.source_y_r_vivo_original = None
        # self.target_x_h_vitr_original = None
        # self.source_x_r_norms = None
        # self.source_y_r_norms = None
        # self.target_x_h_norms = None
        # self.lables = None
        # self.all_usable_genes = None

        self.TargetX_train_ori = None
        self.SourceX_train_ori = None
        self.SourceY_train_ori = None
        self.TargetX_valid_ori = None
        self.SourceX_valid_ori = None
        self.SourceY_valid_ori = None
        self.TargetX_train_prepared = None
        self.SourceX_train_prepared = None
        self.SourceY_train_prepared = None
        self.TargetX_valid_prepared = None
        self.SourceX_valid_prepared = None
        self.SourceY_valid_prepared = None
        self.Labels_train = None
        self.Labels_valid = None
        self.X_norms = None
        self.Y_norms = None

    def correct_slopes(self, X: np.array, vivo=False) -> np.array:
        """ Change Learning data to shape of series instead of series """
        if not vivo:
            for i in range(0, X.shape[1], 3):
                X[:, i + 2] = X[:, i + 2] - X[:, i + 1]
                X[:, i + 1] = X[:, i + 1] - X[:, i + 0]
        else:
            for i in range(0, X.shape[1], 4):
                X[:, i + 3] = X[:, i + 3] - X[:, i + 2]
                X[:, i + 2] = X[:, i + 2] - X[:, i + 1]
                X[:, i + 1] = X[:, i + 1] - X[:, i + 0]
        return X

    def decorrect_slopes(self, X: np.array, vivo=False) -> np.array:
        """ Change slope of series back to original data """
        if not vivo:
            for i in range(0, X.shape[1], 3):
                X[:, i + 1] = X[:, i + 1] + X[:, i + 0]
                X[:, i + 2] = X[:, i + 2] + X[:, i + 1]
        else:
            for i in range(0, X.shape[1], 4):
                X[:, i + 1] = X[:, i + 1] + X[:, i + 0]
                X[:, i + 2] = X[:, i + 2] + X[:, i + 1]
                X[:, i + 3] = X[:, i + 3] + X[:, i + 2]
        return X

    def compute_x_norms(self, sourceX, targetX) -> Tuple:
        combined = np.append(sourceX, targetX, axis=0)
        maxX = np.max(combined, axis=0)
        minX = np.min(combined, axis=0)
        return (minX, maxX)

    def normalize_features(self, X: np.array) -> Tuple[np.array, Tuple]:
        """ Normalize using min/max per feature."""
        maxX = np.max(X, axis=0)
        minX = np.min(X, axis=0)
        X = (X - minX) / (maxX - minX)
        return X, (minX, maxX)

    def normalize_with_norms(self, X: np.array, norms: Tuple) -> np.array:
        """ Normalize using given norms. Tuple is (min, max) where min
        and max are either integers or np.arrays (when normalizing per feature)"""
        X = (X - norms[0]) / (norms[1] - norms[0])
        return X

    def denormalize_with_norms(self, X: np.array, norms: Tuple) -> np.array:
        """ Denormalize using given norms. Tuple is (min, max) where min
        and max are either integers or np.arrays (when normalizing per feature)"""
        X = (X * (norms[1] - norms[0])) + norms[0]
        return X

    def get_data_all_genes_v2(self, file_name, append_toxicity_labels, doubling=True, append_domain_label=True):
        # Load data
        # gene_list = SMALL_SET_50
        gene_list = GTX
        compound_list = GENERAL_47
        x_type = "rat_vitro"
        y_type = "human_vitro"
        z_type = "rat_vivo"

        X_RAT_VITRO, Z_HUMAN_VIVO, Y_RAT_VIVO, data_compounds, all_usable_genes, _, _ = read_data_xyz(compound_list,
                                                                                                      x_type=x_type,
                                                                                                      y_type=y_type,
                                                                                                      z_type=z_type,
                                                                                                      gene_list=gene_list,
                                                                                                      dataset="big")

        source_x_r_vitr = X_RAT_VITRO
        source_y_r_vivo = Y_RAT_VIVO
        target_x_h_vitr = Z_HUMAN_VIVO
        self.all_usable_genes = all_usable_genes
        dosages = np.array(['_C', '_L', '_M', '_H'])
        dosages = np.tile(np.repeat(dosages, 8), int(X_RAT_VITRO.shape[0] / 32))

        compound_lables = np.array([data_compounds[i] + dosages[i] for i in range(len(dosages))])
        self.labels = compound_lables

        # [print(l) for l in self.labels]
        # print(self.all_usable_genes)

        # remove doubling
        if not doubling:
            first_labels = self.labels[0::8]
            first_source_x_r_vitr = source_x_r_vitr[0::8]
            first_source_y_r_vivo = source_y_r_vivo[0::8]
            first_target_x_h_vitr = target_x_h_vitr[0::8]

            last_labels = self.labels[7::8]
            last_source_x_r_vitr = source_x_r_vitr[7::8]
            last_source_y_r_vivo = source_y_r_vivo[7::8]
            last_target_x_h_vitr = target_x_h_vitr[7::8]

            self.labels = self.labels[:(len(first_labels) + len(last_labels))]
            source_x_r_vitr = np.zeros(shape=(len(first_labels) + len(last_labels), source_x_r_vitr.shape[1]))
            source_y_r_vivo = np.zeros(shape=(len(first_labels) + len(last_labels), source_y_r_vivo.shape[1]))
            target_x_h_vitr = np.zeros(shape=(len(first_labels) + len(last_labels), target_x_h_vitr.shape[1]))

            self.labels[0::2] = first_labels
            self.labels[1::2] = last_labels
            source_x_r_vitr[0::2] = first_source_x_r_vitr
            source_y_r_vivo[0::2] = first_source_y_r_vivo
            target_x_h_vitr[0::2] = first_target_x_h_vitr
            source_x_r_vitr[1::2] = last_source_x_r_vitr
            source_y_r_vivo[1::2] = last_source_y_r_vivo
            target_x_h_vitr[1::2] = last_target_x_h_vitr

        print(source_x_r_vitr.shape)
        print(source_y_r_vivo.shape)
        print(target_x_h_vitr.shape)
        print(source_x_r_vitr.shape)
        print(source_y_r_vivo.shape)
        print(target_x_h_vitr.shape)

        shuffled_indices = np.random.permutation(np.arange(source_x_r_vitr.shape[0]))
        split = int(source_x_r_vitr.shape[0] * 0.8)

        self.TargetX_train_ori = target_x_h_vitr[shuffled_indices[:split]].astype(np.float32)
        self.SourceX_train_ori = source_x_r_vitr[shuffled_indices[:split]].astype(np.float32)
        self.SourceY_train_ori = source_y_r_vivo[shuffled_indices[:split]].astype(np.float32)
        self.TargetX_valid_ori = target_x_h_vitr[shuffled_indices[split:]].astype(np.float32)
        self.SourceX_valid_ori = source_x_r_vitr[shuffled_indices[split:]].astype(np.float32)
        self.SourceY_valid_ori = source_y_r_vivo[shuffled_indices[split:]].astype(np.float32)

        self.Labels_train = self.labels[shuffled_indices[:split]].copy()
        self.Labels_valid = self.labels[shuffled_indices[split:]].copy()
        short_labels_train = self.Labels_train
        short_labels_valid = self.Labels_valid

        if append_toxicity_labels:
            class_labels = pickle.load(open("../Data/labels.p", "rb"))
            train_toxicity_labels = np.zeros((self.SourceY_train_ori.shape[0], 2))
            valid_toxicity_labels = np.zeros((self.SourceY_valid_ori.shape[0], 2))

            for i in range(train_toxicity_labels.shape[0]):
                la = class_labels[short_labels_train[i].split('_')[0]]
                if 'GTX' in la:
                    train_toxicity_labels[i, 0] = la['GTX']
                else:
                    train_toxicity_labels[i, 0] = .5
                if 'C' in la:
                    train_toxicity_labels[i, 1] = la['C']
                else:
                    train_toxicity_labels[i, 1] = .5

            # for i in range(valid_toxicity_labels.shape[0]):
            #     la = class_labels[short_labels_valid[i].split('_')[0]]
            #     if 'GTX' in la:
            #         valid_toxicity_labels[i, 0] = 1
            #     if 'C' in la:
            #         valid_toxicity_labels[i, 1] = 1

            # print(train_toxicity_labels)
            # print(valid_toxicity_labels)

            self.SourceY_train_ori = np.concatenate((self.SourceY_train_ori, train_toxicity_labels), axis=1)
            self.SourceY_valid_ori = np.concatenate((self.SourceY_valid_ori, valid_toxicity_labels), axis=1)

        self.TargetX_train_prepared = self.TargetX_train_ori.copy()
        self.SourceX_train_prepared = self.SourceX_train_ori.copy()
        self.SourceY_train_prepared = self.SourceY_train_ori.copy()
        self.TargetX_valid_prepared = self.TargetX_valid_ori.copy()
        self.SourceX_valid_prepared = self.SourceX_valid_ori.copy()
        self.SourceY_valid_prepared = self.SourceY_valid_ori.copy()

        # Change Learning data to shape of series instead of series
        if append_toxicity_labels:
            self.TargetX_train_prepared = self.correct_slopes(self.TargetX_train_prepared)
            self.SourceX_train_prepared = self.correct_slopes(self.SourceX_train_prepared)
            self.SourceY_train_prepared[:, :-2] = self.correct_slopes(self.SourceY_train_prepared[:, :-2], vivo=True)
            self.TargetX_valid_prepared = self.correct_slopes(self.TargetX_valid_prepared)
            self.SourceX_valid_prepared = self.correct_slopes(self.SourceX_valid_prepared)
            self.SourceY_valid_prepared[:, :-2] = self.correct_slopes(self.SourceY_valid_prepared[:, :-2], vivo=True)
            # Normalising the data
            self.X_norms = self.compute_x_norms(self.SourceX_train_prepared, self.TargetX_train_prepared)

            self.TargetX_train_prepared = self.normalize_with_norms(self.TargetX_train_prepared, self.X_norms)
            self.SourceX_train_prepared = self.normalize_with_norms(self.SourceX_train_prepared, self.X_norms)
            self.TargetX_valid_prepared = self.normalize_with_norms(self.TargetX_valid_prepared, self.X_norms)
            self.SourceX_valid_prepared = self.normalize_with_norms(self.SourceX_valid_prepared, self.X_norms)

            self.SourceY_train_prepared[:, :-2], self.Y_norms = self.normalize_features(
                self.SourceY_train_prepared[:, :-2])
            self.SourceY_valid_prepared[:, :-2] = self.normalize_with_norms(self.SourceY_valid_prepared[:, :-2],
                                                                            self.Y_norms)
        else:
            self.TargetX_train_prepared = self.correct_slopes(self.TargetX_train_prepared)
            self.SourceX_train_prepared = self.correct_slopes(self.SourceX_train_prepared)
            self.SourceY_train_prepared = self.correct_slopes(self.SourceY_train_prepared, vivo=True)
            self.TargetX_valid_prepared = self.correct_slopes(self.TargetX_valid_prepared)
            self.SourceX_valid_prepared = self.correct_slopes(self.SourceX_valid_prepared)
            self.SourceY_valid_prepared = self.correct_slopes(self.SourceY_valid_prepared, vivo=True)

            # Normalising the data
            self.X_norms = self.compute_x_norms(self.SourceX_train_prepared, self.TargetX_train_prepared)

            self.TargetX_train_prepared = self.normalize_with_norms(self.TargetX_train_prepared, self.X_norms)
            self.SourceX_train_prepared = self.normalize_with_norms(self.SourceX_train_prepared, self.X_norms)
            self.TargetX_valid_prepared = self.normalize_with_norms(self.TargetX_valid_prepared, self.X_norms)
            self.SourceX_valid_prepared = self.normalize_with_norms(self.SourceX_valid_prepared, self.X_norms)

            self.SourceY_train_prepared, self.Y_norms = self.normalize_features(self.SourceY_train_prepared)
            self.SourceY_valid_prepared = self.normalize_with_norms(self.SourceY_valid_prepared, self.Y_norms)

        # append domain label: 0 for target and 1 for source
        if append_domain_label:
            self.TargetX_train_prepared = np.concatenate(
                (self.TargetX_train_prepared, np.zeros((self.TargetX_train_prepared.shape[0], 1))), axis=1)
            self.TargetX_valid_prepared = np.concatenate(
                (self.TargetX_valid_prepared, np.zeros((self.TargetX_valid_prepared.shape[0], 1))), axis=1)
            self.SourceX_train_prepared = np.concatenate(
                (self.SourceX_train_prepared, np.ones((self.SourceX_train_prepared.shape[0], 1))), axis=1)
            self.SourceX_valid_prepared = np.concatenate(
                (self.SourceX_valid_prepared, np.ones((self.SourceX_valid_prepared.shape[0], 1))), axis=1)

        print('final shapes')
        print('TargetX_train', self.TargetX_train_prepared.shape)
        print('TargetX_valid', self.TargetX_valid_prepared.shape)
        print('SourceX_train', self.SourceX_train_prepared.shape)
        print('SourceY_train', self.SourceY_train_prepared.shape)
        print('SourceX_valid', self.SourceX_valid_prepared.shape)
        print('SourceY_valid', self.SourceY_valid_prepared.shape)
        data = {
            'TargetX_train': self.TargetX_train_prepared,
            'TargetX_valid': self.TargetX_valid_prepared,
            'SourceX_train': self.SourceX_train_prepared,
            'SourceY_train': self.SourceY_train_prepared,
            'SourceX_valid': self.SourceX_valid_prepared,
            'SourceY_valid': self.SourceY_valid_prepared
        }

        text = {
            'short_labels_train': short_labels_train,
            'short_labels_valid': short_labels_valid,
            'gene_names': self.all_usable_genes,
            'Labels_valid': self.Labels_valid,
            'Labels_train': self.Labels_train
        }

        return data, text

    def load_data_for_loo(self, doubling=True):
        # Load data
        # gene_list = SMALL_SET_50
        gene_list = CARCINOGENIC
        compound_list = GENERAL_47

        x_type = "rat_vitro"
        y_type = "human_vitro"
        z_type = "rat_vivo"

        X_RAT_VITRO, Z_HUMAN_VIVO, Y_RAT_VIVO, data_compounds, all_usable_genes, _, _ = read_data_xyz(compound_list,
                                                                                                      x_type=x_type,
                                                                                                      y_type=y_type,
                                                                                                      z_type=z_type,
                                                                                                      gene_list=gene_list,
                                                                                                      dataset="big")

        source_x_r_vitr = X_RAT_VITRO
        source_y_r_vivo = Y_RAT_VIVO
        target_x_h_vitr = Z_HUMAN_VIVO
        self.all_usable_genes = all_usable_genes
        dosages = np.array(['_C', '_L', '_M', '_H'])
        dosages = np.tile(np.repeat(dosages, 8), int(X_RAT_VITRO.shape[0] / 32))

        compound_lables = np.array([data_compounds[i] + dosages[i] for i in range(len(dosages))])
        self.labels = compound_lables

        # [print(l) for l in self.labels]
        # print(self.all_usable_genes)

        # remove doubling
        if not doubling:
            first_labels = self.labels[0::8]
            first_source_x_r_vitr = source_x_r_vitr[0::8]
            first_source_y_r_vivo = source_y_r_vivo[0::8]
            first_target_x_h_vitr = target_x_h_vitr[0::8]

            last_labels = self.labels[7::8]
            last_source_x_r_vitr = source_x_r_vitr[7::8]
            last_source_y_r_vivo = source_y_r_vivo[7::8]
            last_target_x_h_vitr = target_x_h_vitr[7::8]

            self.labels = self.labels[:(len(first_labels) + len(last_labels))]
            source_x_r_vitr = np.zeros(shape=(len(first_labels) + len(last_labels), source_x_r_vitr.shape[1]))
            source_y_r_vivo = np.zeros(shape=(len(first_labels) + len(last_labels), source_y_r_vivo.shape[1]))
            target_x_h_vitr = np.zeros(shape=(len(first_labels) + len(last_labels), target_x_h_vitr.shape[1]))

            self.labels[0::2] = first_labels
            self.labels[1::2] = last_labels
            source_x_r_vitr[0::2] = first_source_x_r_vitr
            source_y_r_vivo[0::2] = first_source_y_r_vivo
            target_x_h_vitr[0::2] = first_target_x_h_vitr
            source_x_r_vitr[1::2] = last_source_x_r_vitr
            source_y_r_vivo[1::2] = last_source_y_r_vivo
            target_x_h_vitr[1::2] = last_target_x_h_vitr

        self.loaded_target_x_h_vitr = target_x_h_vitr
        self.loaded_source_x_r_vitr = source_x_r_vitr
        self.loaded_source_y_r_vivo = source_y_r_vivo

    def leave_one_out(self, i_one, append_toxicity_labels, doubling=True):
        # # Load data
        # gene_list = SMALL_SET_50
        # compound_list = GENERAL_47
        # x_type = "rat_vitro"
        # y_type = "human_vitro"
        # z_type = "rat_vivo"
        #
        # X_RAT_VITRO, Z_HUMAN_VIVO, Y_RAT_VIVO, data_compounds, all_usable_genes, _, _ = read_data_xyz(compound_list,
        #                                                                                               x_type=x_type,
        #                                                                                               y_type=y_type,
        #                                                                                               z_type=z_type,
        #                                                                                               gene_list=gene_list,
        #                                                                                               dataset="big")
        #
        # source_x_r_vitr = X_RAT_VITRO
        # source_y_r_vivo = Y_RAT_VIVO
        # target_x_h_vitr = Z_HUMAN_VIVO
        # self.all_usable_genes = all_usable_genes
        # dosages = np.array(['_C', '_L', '_M', '_H'])
        # dosages = np.tile(np.repeat(dosages, 8), int(X_RAT_VITRO.shape[0] / 32))
        #
        # compound_lables = np.array([data_compounds[i] + dosages[i] for i in range(len(dosages))])
        # self.labels = compound_lables
        #
        # # [print(l) for l in self.labels]
        # # print(self.all_usable_genes)
        #
        # # remove doubling
        # if not doubling:
        #     first_labels = self.labels[0::8]
        #     first_source_x_r_vitr = source_x_r_vitr[0::8]
        #     first_source_y_r_vivo = source_y_r_vivo[0::8]
        #     first_target_x_h_vitr = target_x_h_vitr[0::8]
        #
        #     last_labels = self.labels[7::8]
        #     last_source_x_r_vitr = source_x_r_vitr[7::8]
        #     last_source_y_r_vivo = source_y_r_vivo[7::8]
        #     last_target_x_h_vitr = target_x_h_vitr[7::8]
        #
        #     self.labels = self.labels[:(len(first_labels) + len(last_labels))]
        #     source_x_r_vitr = np.zeros(shape=(len(first_labels) + len(last_labels), source_x_r_vitr.shape[1]))
        #     source_y_r_vivo = np.zeros(shape=(len(first_labels) + len(last_labels), source_y_r_vivo.shape[1]))
        #     target_x_h_vitr = np.zeros(shape=(len(first_labels) + len(last_labels), target_x_h_vitr.shape[1]))
        #
        #     self.labels[0::2] = first_labels
        #     self.labels[1::2] = last_labels
        #     source_x_r_vitr[0::2] = first_source_x_r_vitr
        #     source_y_r_vivo[0::2] = first_source_y_r_vivo
        #     target_x_h_vitr[0::2] = first_target_x_h_vitr
        #     source_x_r_vitr[1::2] = last_source_x_r_vitr
        #     source_y_r_vivo[1::2] = last_source_y_r_vivo
        #     target_x_h_vitr[1::2] = last_target_x_h_vitr
        #
        # # print(source_x_r_vitr.shape)
        # # print(source_y_r_vivo.shape)
        # # print(target_x_h_vitr.shape)
        # # print(source_x_r_vitr.shape)
        # # print(source_y_r_vivo.shape)
        # # print(target_x_h_vitr.shape)
        #
        # # shuffled_indices = np.random.permutation(np.arange(source_x_r_vitr.shape[0]))
        # # split = int(source_x_r_vitr.shape[0] * 0.8)
        # # self.TargetX_train_ori = target_x_h_vitr[shuffled_indices[:split]].astype(np.float32)
        # # self.SourceX_train_ori = source_x_r_vitr[shuffled_indices[:split]].astype(np.float32)
        # # self.SourceY_train_ori = source_y_r_vivo[shuffled_indices[:split]].astype(np.float32)
        # # self.TargetX_valid_ori = target_x_h_vitr[shuffled_indices[split:]].astype(np.float32)
        # # self.SourceX_valid_ori = source_x_r_vitr[shuffled_indices[split:]].astype(np.float32)
        # # self.SourceY_valid_ori = source_y_r_vivo[shuffled_indices[split:]].astype(np.float32)
        # # self.Labels_train = self.labels[shuffled_indices[:split]].copy()
        # # self.Labels_valid = self.labels[shuffled_indices[split:]].copy()

        target_x_h_vitr = self.loaded_target_x_h_vitr.copy()
        source_x_r_vitr = self.loaded_source_x_r_vitr.copy()
        source_y_r_vivo = self.loaded_source_y_r_vivo.copy()

        all_indices = np.arange(source_x_r_vitr.shape[0])
        loo_start = i_one * 32
        loo_end = (i_one * 32) + 32
        if loo_end > len(all_indices):
            return None, None
        train_indices = np.delete(all_indices, range(loo_start, loo_end))
        valid_indices = all_indices[loo_start:loo_end]

        self.TargetX_train_ori = target_x_h_vitr[train_indices].astype(np.float32)
        self.SourceX_train_ori = source_x_r_vitr[train_indices].astype(np.float32)
        self.SourceY_train_ori = source_y_r_vivo[train_indices].astype(np.float32)
        self.TargetX_valid_ori = target_x_h_vitr[valid_indices].astype(np.float32)
        self.SourceX_valid_ori = source_x_r_vitr[valid_indices].astype(np.float32)
        self.SourceY_valid_ori = source_y_r_vivo[valid_indices].astype(np.float32)
        self.Labels_train = self.labels[train_indices].copy()
        self.Labels_valid = self.labels[valid_indices].copy()

        short_labels_train = self.Labels_train
        short_labels_valid = self.Labels_valid

        if append_toxicity_labels:
            class_labels = pickle.load(open("../Data/labels.p", "rb"))
            train_toxicity_labels = np.zeros((self.SourceY_train_ori.shape[0], 2))
            valid_toxicity_labels = np.zeros((self.SourceY_valid_ori.shape[0], 2))

            for i in range(train_toxicity_labels.shape[0]):
                la = class_labels[short_labels_train[i].split('_')[0]]
                if 'GTX' in la:
                    train_toxicity_labels[i, 0] = la['GTX']
                else:
                    train_toxicity_labels[i, 0] = .5
                if 'C' in la:
                    train_toxicity_labels[i, 1] = la['C']
                else:
                    train_toxicity_labels[i, 1] = .5

            for i in range(valid_toxicity_labels.shape[0]):
                la = class_labels[short_labels_valid[i].split('_')[0]]
                if 'GTX' in la:
                    valid_toxicity_labels[i, 0] = la['GTX']
                else:
                    valid_toxicity_labels[i, 0] = .5
                if 'C' in la:
                    valid_toxicity_labels[i, 1] = la['C']
                else:
                    valid_toxicity_labels[i, 1] = .5
            # for i in range(valid_toxicity_labels.shape[0]):
            #     la = class_labels[short_labels_valid[i].split('_')[0]]
            #     if 'GTX' in la:
            #         valid_toxicity_labels[i, 0] = 1
            #     if 'C' in la:
            #         valid_toxicity_labels[i, 1] = 1

            # print(train_toxicity_labels)
            # print(valid_toxicity_labels)

            self.SourceY_train_ori = np.concatenate((self.SourceY_train_ori, train_toxicity_labels), axis=1)
            self.SourceY_valid_ori = np.concatenate((self.SourceY_valid_ori, valid_toxicity_labels), axis=1)

        self.TargetX_train_prepared = self.TargetX_train_ori.copy()
        self.SourceX_train_prepared = self.SourceX_train_ori.copy()
        self.SourceY_train_prepared = self.SourceY_train_ori.copy()
        self.TargetX_valid_prepared = self.TargetX_valid_ori.copy()
        self.SourceX_valid_prepared = self.SourceX_valid_ori.copy()
        self.SourceY_valid_prepared = self.SourceY_valid_ori.copy()

        # Change Learning data to shape of series instead of series
        if append_toxicity_labels:
            self.TargetX_train_prepared = self.correct_slopes(self.TargetX_train_prepared)
            self.SourceX_train_prepared = self.correct_slopes(self.SourceX_train_prepared)
            self.SourceY_train_prepared[:, :-2] = self.correct_slopes(self.SourceY_train_prepared[:, :-2], vivo=True)
            self.TargetX_valid_prepared = self.correct_slopes(self.TargetX_valid_prepared)
            self.SourceX_valid_prepared = self.correct_slopes(self.SourceX_valid_prepared)
            self.SourceY_valid_prepared[:, :-2] = self.correct_slopes(self.SourceY_valid_prepared[:, :-2], vivo=True)
            # Normalising the data
            self.X_norms = self.compute_x_norms(self.SourceX_train_prepared, self.TargetX_train_prepared)

            self.TargetX_train_prepared = self.normalize_with_norms(self.TargetX_train_prepared, self.X_norms)
            self.SourceX_train_prepared = self.normalize_with_norms(self.SourceX_train_prepared, self.X_norms)
            self.TargetX_valid_prepared = self.normalize_with_norms(self.TargetX_valid_prepared, self.X_norms)
            self.SourceX_valid_prepared = self.normalize_with_norms(self.SourceX_valid_prepared, self.X_norms)

            self.SourceY_train_prepared[:, :-2], self.Y_norms = self.normalize_features(
                self.SourceY_train_prepared[:, :-2])
            self.SourceY_valid_prepared[:, :-2] = self.normalize_with_norms(self.SourceY_valid_prepared[:, :-2],
                                                                            self.Y_norms)
        else:
            self.TargetX_train_prepared = self.correct_slopes(self.TargetX_train_prepared)
            self.SourceX_train_prepared = self.correct_slopes(self.SourceX_train_prepared)
            self.SourceY_train_prepared = self.correct_slopes(self.SourceY_train_prepared, vivo=True)
            self.TargetX_valid_prepared = self.correct_slopes(self.TargetX_valid_prepared)
            self.SourceX_valid_prepared = self.correct_slopes(self.SourceX_valid_prepared)
            self.SourceY_valid_prepared = self.correct_slopes(self.SourceY_valid_prepared, vivo=True)

            # Normalising the data
            self.X_norms = self.compute_x_norms(self.SourceX_train_prepared, self.TargetX_train_prepared)

            self.TargetX_train_prepared = self.normalize_with_norms(self.TargetX_train_prepared, self.X_norms)
            self.SourceX_train_prepared = self.normalize_with_norms(self.SourceX_train_prepared, self.X_norms)
            self.TargetX_valid_prepared = self.normalize_with_norms(self.TargetX_valid_prepared, self.X_norms)
            self.SourceX_valid_prepared = self.normalize_with_norms(self.SourceX_valid_prepared, self.X_norms)

            self.SourceY_train_prepared, self.Y_norms = self.normalize_features(self.SourceY_train_prepared)
            self.SourceY_valid_prepared = self.normalize_with_norms(self.SourceY_valid_prepared, self.Y_norms)

        # append domain label: 0 for target and 1 for source
        self.TargetX_train_prepared = np.concatenate(
            (self.TargetX_train_prepared, np.zeros((self.TargetX_train_prepared.shape[0], 1))), axis=1)
        self.TargetX_valid_prepared = np.concatenate(
            (self.TargetX_valid_prepared, np.zeros((self.TargetX_valid_prepared.shape[0], 1))), axis=1)
        self.SourceX_train_prepared = np.concatenate(
            (self.SourceX_train_prepared, np.ones((self.SourceX_train_prepared.shape[0], 1))), axis=1)
        self.SourceX_valid_prepared = np.concatenate(
            (self.SourceX_valid_prepared, np.ones((self.SourceX_valid_prepared.shape[0], 1))), axis=1)

        # print('Final shapes')
        # print('\tTargetX_train', self.TargetX_train_prepared.shape)
        # print('\tTargetX_valid', self.TargetX_valid_prepared.shape)
        # print('\tSourceX_train', self.SourceX_train_prepared.shape)
        # print('\tSourceY_train', self.SourceY_train_prepared.shape)
        # print('\tSourceX_valid', self.SourceX_valid_prepared.shape)
        # print('\tSourceY_valid', self.SourceY_valid_prepared.shape)
        data = {
            'TargetX_train': self.TargetX_train_prepared,
            'TargetX_valid': self.TargetX_valid_prepared,
            'SourceX_train': self.SourceX_train_prepared,
            'SourceY_train': self.SourceY_train_prepared,
            'SourceX_valid': self.SourceX_valid_prepared,
            'SourceY_valid': self.SourceY_valid_prepared
        }

        text = {
            'short_labels_train': short_labels_train,
            'short_labels_valid': short_labels_valid,
            'gene_names': self.all_usable_genes,
            'Labels_valid': self.Labels_valid,
            'Labels_train': self.Labels_train
        }

        return data, text

    def correct_data(self, source_y_r, APPEND_TOXICITY_LABELS, target=False):
        # Denormalize all data
        if APPEND_TOXICITY_LABELS:
            # print(source_y_r.shape)
            # print(source_y_r[:,:-2].shape)
            if not target:
                source_y_r[:, :-2] = self.denormalize_with_norms(source_y_r[:, :-2], self.Y_norms)
            else:
                source_y_r[:, :-2] = self.denormalize_with_norms(source_y_r[:, :-2], self.Y_norms)
            source_y_r[:, :-2] = self.decorrect_slopes(source_y_r[:, :-2], vivo=True)
        else:
            if not target:
                source_y_r = self.denormalize_with_norms(source_y_r, self.Y_norms)
            else:
                source_y_r = self.denormalize_with_norms(source_y_r, self.Y_norms)

            # Decorrect for slopes, otherwise normalization won't be correct
            # X_train = self.decorrect_slopes(X_train)
            # X_valid = self.decorrect_slopes(X_valid)
            source_y_r = self.decorrect_slopes(source_y_r, vivo=True)

            # Normalize again, but using min/max over entire dataset instead of per feature.
            # This is to have a consistent error measurement across different experiments.

        return source_y_r

    def get_original(self):
        return {'target_x_train': self.TargetX_train_ori,
                'source_x_train': self.SourceX_train_ori,
                'source_y_train': self.SourceY_train_ori,
                'target_x_valid': self.TargetX_valid_ori,
                'source_x_valid': self.SourceX_valid_ori,
                'source_y_valid': self.SourceY_valid_ori}

    def get_prepared(self):
        return {'target_x_train': self.TargetX_train_prepared,
                'source_x_train': self.SourceX_train_prepared,
                'source_y_train': self.SourceY_train_prepared,
                'target_x_valid': self.TargetX_valid_prepared,
                'source_x_valid': self.SourceX_valid_prepared,
                'source_y_valid': self.SourceY_valid_prepared}

    def get_labels(self):
        return self.Labels_train, self.Labels_valid

    def get_norms(self):
        return self.X_norms, self.Y_norms

    def get_source_x_train_original(self):
        return self.SourceX_train_ori

    def get_source_y_train_original(self):
        return self.SourceY_train_ori

    def get_target_x_train_original(self):
        return self.TargetX_train_ori

    def get_source_x_valid_original(self):
        return self.SourceX_valid_ori

    def get_source_y_valid_original(self):
        return self.SourceY_valid_ori

    def get_target_x_valid_original(self):
        return self.TargetX_valid_ori

