#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import tensorflow as tf
from typing import List, Dict, Callable, Tuple

from Helper import create_gene_name_mapping

import numpy as np

from Types import Data


# Split X and Y set into training and test data, based on random split or exclude compound (loo-evaluation)
def split_train_test(X, Y, compounds, x_vivo, y_vivo, exclude_compound='', train_split=0.75):
    
    # If not leaving out one compound, randomly split in train and test according to ratio
    if exclude_compound == '':
        # print("Splitting randomly in", train_split, 1-train_split)
        shuffled_indices = np.random.permutation(np.arange(X.shape[0]))
        split = int(X.shape[0] * train_split)
        X_train = X[shuffled_indices[:split]].astype(np.float32)
        X_valid = X[shuffled_indices[split:]].astype(np.float32)
        Y_train = Y[shuffled_indices[:split]].astype(np.float32)
        Y_valid = Y[shuffled_indices[split:]].astype(np.float32)
        train_compounds = compounds[shuffled_indices[:split]]
        valid_compounds = compounds[shuffled_indices[split:]]
    else:
        # print("Excluding compound", exclude_compound)
        train_indices = np.where(compounds != exclude_compound)
        test_indices = np.where(compounds == exclude_compound)
        X_train = X[train_indices].astype(np.float32)
        X_valid = X[test_indices].astype(np.float32)
        Y_train = Y[train_indices].astype(np.float32)
        Y_valid = Y[test_indices].astype(np.float32)
        train_compounds = compounds[train_indices]
        valid_compounds = compounds[test_indices]

    # Change Learning data to shape of series instead of series
    X_train = correct_slopes(X_train, x_vivo)
    X_valid = correct_slopes(X_valid, x_vivo)
    Y_train = correct_slopes(Y_train, y_vivo)
    Y_valid = correct_slopes(Y_valid, y_vivo)

    X_train, norms_X = normalize_features(X_train)
    X_valid = normalize_with_norms(X_valid, norms_X)
    Y_train, norms_Y = normalize_features(Y_train)
    Y_valid = normalize_with_norms(Y_valid, norms_Y)

    return X_train, X_valid, Y_train, Y_valid, norms_X, norms_Y, train_compounds, valid_compounds


def get_all_gene_activations(data: Data, genes: List[str] = None) -> np.ndarray:
    """
    Get compounds*dosages*replicates arrays of length genes*timepoints, containing
    the gene activation data.
    :param data: data to generate the combinations from
    :param genes: genes to use, all if None
    :return:
    """
    header = data.header
    if not genes:
        genes = header.genes

    hierarchical = get_hierarchical_data(data, genes)

    nr_dosages = len(header.dosages)
    nr_compounds = len(header.compounds)
    nr_replicates = len(header.replicates)
    nr_genes = len(genes)
    nr_timepoints = len(header.times)

    x = np.ndarray((nr_dosages * nr_compounds * nr_replicates, nr_genes * nr_timepoints))

    idx_x = 0

    for idx_compound in range(len(header.compounds)):
        for idx_dosage in range(len(header.dosages)):
            for idx_replicate in range(len(header.replicates)):
                activations = []  # activations for exactly this combination
                for g in genes:
                    activations.extend(hierarchical[g][idx_compound][idx_dosage][idx_replicate])
                x[idx_x] = activations
                idx_x += 1

    x = x
    return x


def normalize_nparray(x: np.ndarray) -> np.ndarray:
    """Normalizes the array so all values lie between zero and one."""
    # Normalising the data
    max_x = max(map(max, x))
    min_x = min(map(min, x))
    return (np.array(x) - min_x) / (max_x - min_x), min_x, max_x


def multigenes_average(all_activations: List[np.array]) -> np.array:
    """Reduce multiple occurences of a single gene by calculating the arithmetic mean between all activations."""
    return np.mean(all_activations, 0)


def multigenes_median(all_activations: List[np.array]) -> np.array:
    """Reduce multiple occurences of a single gene by calculating the arithmetic median between all activations."""
    return np.median(all_activations, 0)


def get_hierarchical_data(data: Data, genes_to_use: List[str] = None,
                          multigene_reducer: Callable[[List[np.array]], np.array] = multigenes_median
                          ) -> Dict[str, List[List[List[List[float]]]]]:
    """
    Data formatted as data[gene][compound][dosage][replicate][time]
                           name  index     index   index      index
    Could be changed in the future, so don't rely on it.
    """

    activations = data.activations
    header = data.header

    if genes_to_use:
        considered_genes = list(set(header.genes).intersection(set(genes_to_use)))
        if len(considered_genes) != len(genes_to_use):
            print("WARN: {} genes specified, but only {} found in intersection".format(len(genes_to_use),
                                                                                       len(considered_genes)))
            print("Namely ", set(genes_to_use)-set(considered_genes), " is/are missing")
    else:
        considered_genes = list(set(header.genes))

    # TODO consider timepoints in case their order is not right in the big dataset, or verify for all
    # verified to not be a problem in: Rat In Vitro,

    # if anyone has a more readable version of this without loops, I'm in
    hierarchical = \
        {gene: [[[[]
                  for _ in header.replicates]
                 for _ in header.dosages]
                for _ in header.compounds]
         for gene in considered_genes}
    cols = header.columns

    for gene in [g for g in considered_genes if g is not None]:

        indices_for_gene = header.get_gene_indices(gene)
        gene_activs = [activations[idx] for idx in indices_for_gene]
        activ = multigene_reducer(gene_activs)

        for i in range(len(activ)):
            hierarchical[gene][cols[i].compound][cols[i].dosage][cols[i].replicate].append(activ[i])
    return hierarchical


def correct_slopes(X: np.array, vivo=False) -> np.array:
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

def decorrect_slopes(X: np.array, vivo=False) -> np.array:
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

def normalize_total(X: np.array) -> Tuple[np.array, Tuple]:
    """ Normalize using min/max of entire dataset"""
    maxX = np.max(X)
    minX = np.min(X)
    X = (X - minX) / (maxX - minX)
    return X, (minX, maxX)

def standardize_features(X: np.array) -> Tuple[np.array, Tuple]:
    """ Standardize per feature."""
    mean = np.mean(X, axis=0)
    corrected = X - mean
    std = np.std(corrected, axis=0)
    standardized = corrected / std
    return standardized, (mean, std)

def destandardize(X: np.array, stats: Tuple) -> np.array:
    """ Destandardize using (mean, std) per feature."""
    X *= stats[1]
    X += stats[0]
    return X

def standardize_with_stats(X: np.array, stats: Tuple) -> np.array:
    """ Standardize using (mean, std) per feature."""
    corrected = X - stats[0]
    standardized = corrected / stats[1]
    return standardized

def normalize_features(X: np.array) -> Tuple[np.array, Tuple]:
    """ Normalize using min/max per feature."""
    maxX = np.max(X, axis=0)
    minX = np.min(X, axis=0)
    X = (X - minX) / (maxX - minX)
    return X, (minX, maxX)

def normalize_with_norms(X: np.array, norms: Tuple) -> np.array:
    """ Normalize using given norms. Tuple is (min, max) where min 
    and max are either integers or np.arrays (when normalizing per feature)"""
    X = (X - norms[0]) / (norms[1] - norms[0])
    return X


def denormalize_with_norms(X: np.array, norms: Tuple) -> np.array:
    """ Denormalize using given norms. Tuple is (min, max) where min 
    and max are either integers or np.arrays (when normalizing per feature)"""
    X = (X * (norms[1] - norms[0])) + norms[0]
    return X

#get_all_gene_activations('/Users/Daniel/Desktop/MaCSBio/mrp2-python/src/Data/HumanData_51genes_1128samples.txt')
