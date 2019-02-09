#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import sys
import pickle
import hashlib
import itertools

import numpy as np

from Helper import add_once_to_list, create_gene_name_mapping
from Processing import get_hierarchical_data
from Types import Header, ColInfo, Data
from BigParser import BigDatasetParser
from SmallParser import SmallDatasetParser


def get_uda_data(rat_vitro: Data, rat_vivo: Data, human_vitro: Data):
    """
    Return X, Y, Z for unsupervised domain adaption.
    :param rat_vitro:
    :param rat_vivo:
    :param human_vitro:
    :return:
    """

    return None, None, None


def main_exportUdaData():
    """
    Reads all datasets + exports them as X, Y, Z for the UDA.
    """

    # TODO do non-hacky

    if len(sys.argv) != 4:
        print("Please specify 3 arguments for rat in vivo: [datafile] [rowfile (probes)] [columnfile (samples)]")
        sys.exit(1)

    parser = SmallDatasetParser()
    human_vitro, rat_vitro = parser.parse()

    gene_mapping = create_gene_name_mapping()
    genes_to_select = list(gene_mapping.keys())

    parser = BigDatasetParser()
    rat_vivo = parser.parse(*sys.argv[1:])

    if human_vitro.header.dosages != rat_vivo.header.dosages:
        print("ERROR: dosages NOK")
    if not set(human_vitro.header.compounds).issubset(set(rat_vivo.header.compounds)):
        print("ERROR: Compounds NOK")

    print("ignore following debug outputs")
    print(rat_vivo.header)
    print(len(rat_vivo.activations))

    for gene_name in genes_to_select:
        if gene_name not in rat_vivo.header.genes:
            print("ERROR: generated name '{}' seems to be wrong".format(gene_name))

    hierarchical_rat_vivo_wrong_names = get_hierarchical_data(rat_vivo, genes_to_select)
    hierarchical_rat_vivo = {gene_mapping[g]: hierarchical_rat_vivo_wrong_names[g] for g in
                             hierarchical_rat_vivo_wrong_names.keys()}

    # those are the genes we can use for all 3 datasets
    all_usable_genes = list(gene_mapping.values())

    hierarchical_rat_vitro = get_hierarchical_data(rat_vitro, all_usable_genes)
    hierarchical_human_vitro = get_hierarchical_data(human_vitro, all_usable_genes)

    X_RAT_VITRO = []
    Y_RAT_VIVO = []
    Z_HUMAN_VIVO = []
    labels = []

    # use data from this header (equal to human_vitro header)
    header = rat_vitro.header

    # for the three excluded compounds, some dosages are missing
    usable_compounds = [x for x in header.compounds if x not in ["ADP", "CPZ", "WY"]]

    # Data formatted as data[gene][compound][dosage][replicate][time]
    # name  index     index   index      index

    for compound_name in usable_compounds:
        compound_idx_vitro = header.compounds.index(compound_name)
        compound_idx_vivo = rat_vivo.header.compounds.index(compound_name)
        for d in range(len(header.dosages)):
            current_label = "Compound '{}' at dosage '{}'".format(compound_name, header.dosages[d])

            # add all time series for both replicates, then combine each pair
            activations_x_rat_vitro = []
            activations_y_rat_vivo = []
            activations_z_human_vitro = []
            for r in range(len(header.replicates)):
                # collect activations for each replicate, so cross product can be used later on
                act_x_rat_vitro = []
                act_y_rat_vivo = []
                act_z_human_vitro = []
                for gene in all_usable_genes:
                    act_x_rat_vitro.extend(hierarchical_rat_vitro[gene][compound_idx_vitro][d][r])
                    act_y_rat_vivo.extend(hierarchical_rat_vivo[gene][compound_idx_vivo][d][r])
                    act_z_human_vitro.extend(hierarchical_human_vitro[gene][compound_idx_vitro][d][r])
                activations_x_rat_vitro.append(act_x_rat_vitro)
                activations_y_rat_vivo.append(act_y_rat_vivo)
                activations_z_human_vitro.append(act_z_human_vitro)
            for x_rat_vitro, y_rat_vivo, z_human_vitro in itertools.product(activations_x_rat_vitro,
                                                                            activations_y_rat_vivo,
                                                                            activations_z_human_vitro):
                X_RAT_VITRO.append(x_rat_vitro)
                Y_RAT_VIVO.append(y_rat_vivo)
                Z_HUMAN_VIVO.append(z_human_vitro)
                labels.append(current_label)

    print("LEN X {}, LEN X[0] {}, LEN Y {}, LEN Y[0] {}, LEN Z {}, LEN Z[0] {}".format(
        len(X_RAT_VITRO),
        len(X_RAT_VITRO[0]),
        len(Y_RAT_VIVO),
        len(Y_RAT_VIVO[0]),
        len(Z_HUMAN_VIVO),
        len(Z_HUMAN_VIVO[0]))
    )

    with open("uda_xyz.p", "wb") as outfile:
        pickle.dump((X_RAT_VITRO, Y_RAT_VIVO, Z_HUMAN_VIVO, labels, all_usable_genes), outfile)

    print("Successfully exported X, Y, Z")


def main_exportCnnData():
    """
    Currently exports rat in vitro and human in vitro data for usage in convolutional neural networks, meaning
    each training sample is 2-dimensional, like:

    [[1 2 3], # timepoints of first gene
    [4,5,6]] # timepoints of second gene and so on
    :return:
    """

    parser = SmallDatasetParser()
    human_vitro, rat_vitro = parser.parse()

    if human_vitro.header != rat_vitro.header:
        raise ValueError("human and rat genes or header are not the same")

    genes = rat_vitro.header.genes

    rat_hierarchical = get_hierarchical_data(rat_vitro, genes)
    human_hierarchical = get_hierarchical_data(human_vitro, genes)
    header = human_vitro.header

    X_RAT = []
    Y_HUMAN = []
    labels = []

    # Data formatted as data[gene][compound][dosage][replicate][time]
    # name  index     index   index      index

    for c in range(len(header.compounds)):
        for d in range(len(header.dosages)):
            current_label = "Compound '{}' at dosage '{}'".format(header.compounds[c], header.dosages[d])
            # add all time series for both replicates, then combine each pair
            rat_activations = []
            human_activations = []
            for r in range(len(header.replicates)):
                # collect activations for each replicate, so cross product can be used later on
                rat_act_all_genes = []
                human_act_all_genes = []
                for gene in genes:
                    rat_act_all_genes.append(rat_hierarchical[gene][c][d][r])
                    human_act_all_genes.append(human_hierarchical[gene][c][d][r])
                rat_activations.append(rat_act_all_genes)
                human_activations.append(human_act_all_genes)
            for rat_act_repl, human_act_repl in itertools.product(rat_activations, human_activations):
                X_RAT.append(rat_act_repl)
                Y_HUMAN.append(human_act_repl)
                labels.append(current_label)

    with open("cnn_rat_human_vitro.p", "wb") as outfile:
        pickle.dump((X_RAT, Y_HUMAN, labels, header), outfile)

    print("Successfully exported X, Y")


if __name__ == '__main__':
    main_exportUdaData()
