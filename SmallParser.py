#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import itertools
from typing import List, Tuple

import numpy as np
import pickle
import os

from Helper import add_once_to_list, create_gene_name_mapping
from Processing import get_hierarchical_data
from Types import Data, Header, ColInfo
import CompoundLists  # Dan: added for testing read_data
import GeneLists


class SmallDatasetParser:
    """
    Allows to parse the small dataset for humans and rats alike.
    """

    def __init__(self):
        pass

    def parse(self, human_fname: str = "Data/HumanData_51genes_1128samples.txt",
              rat_fname: str = "Data/RatData_51genes_1128samples.txt") -> Tuple[Data, Data]:
        """
        Parses the small dataset with the human and rat gene expression data.
        """
        humandata = self.__parse_file(human_fname)
        ratdata = self.__parse_file(rat_fname)
        if humandata.header != ratdata.header:
            raise ValueError("human and rat data do not have the same header!")

        return humandata, ratdata

    def __parse_file(self, filename: str) -> Data:
        data = None
        with open(filename, "r") as file:
            genes = []
            activations = []
            header_str = file.readline()
            header = self.__parse_header(header_str)
            for line in file:
                split = line.strip().split("\t")
                gene = split[0]
                genes.append(gene)
                activations.append([float(ac) for ac in split[1:]])
                if len(activations[-1]) != len(header.columns):
                    raise ValueError(
                        "{} activations for gene {}, but header defined {}".format(len(activations[-1]), gene,
                                                                                   len(header.columns)))
            header.set_genes(genes)
            np_activations = np.array(activations)
            data = Data(header, np_activations)
        return data

    @staticmethod
    def __parse_header(header: str) -> Header:
        compounds = []
        dosages = []
        timepoints = []
        replicates = []

        columns = []  # describes entry in column index i for dataset

        split = header.strip().split("\t")

        # format is Compound_Dosage_Timepoint_Replicate
        for column in split:
            comp, dos, time, repl = column.split("_")
            comp_index = add_once_to_list(comp, compounds)

            # add dosages, timepoints, replicates if first compound
            if comp == compounds[0]:
                dos_index = add_once_to_list(dos, dosages)
                time_index = add_once_to_list(time, timepoints)
                repl_index = add_once_to_list(repl, replicates)
            elif dos not in dosages or time not in timepoints or repl not in replicates:
                print(Header(compounds, dosages, timepoints, replicates, []))
                raise ValueError("unknown dosage {}, time {} or replicate {} for {}".format(dos, time, repl, comp))
            else:
                dos_index = dosages.index(dos)
                time_index = timepoints.index(time)
                repl_index = replicates.index(repl)
            col_info = ColInfo(comp_index, dos_index, time_index, repl_index)
            columns.append(col_info)
        return Header(compounds, dosages, timepoints, replicates, columns)


def get_doubling_xy(first: Data, second: Data, genes_first=None, genes_second=None, compounds=None, max_replicates=2) -> Tuple[List[List[float]], List[List[float]]]:
    """
    Create doubled X and Y like Kurts parserDoubling, for which genes + headers have to be equal.
    :param first: e. g. human data
    :param second: e. g. rat data
    :param genes_first: genes to select from first dataset, default is all
    :param genes_second: genes to select from second dataset, default is same as genes_first
    :param compounds: default is all compounds in data. Realistically, should use CompoundLists.GENERAL_47
    :param max_replicates: how many replicates to use at most. For example, vivo has 5 replicates, but some are missing data.
                           the first three are safe to use, but using more than two creates complexity in later interpreting
                           the dataset (2x3x4 = 24 instances for each compound instead of 2x2x4 = 16)
    :return: X and Y
    """
    # Dan: this is the function that ACTUALLY reads the data; read_data only specifies which genes to read

    if first.header.dosages != second.header.dosages:
        raise ValueError("human and rat dosages are not the same")
    
    # use all genes if no special gene set indicated
    if genes_first is None:
        genes_first = first.header.genes

    if genes_second is None:
        genes_second = genes_first

    first_hierarchical = get_hierarchical_data(first, genes_first)  # Dan: contains selected genes
    second_hierarchical = get_hierarchical_data(second, genes_second)

    header = first.header
    times_first = len(first.header.times)
    times_second = len(second.header.times)

    if not compounds:
        compounds = header.compounds

    X = []
    Y = []

    # Data formatted as data[gene][compound][dosage][replicate][time]
    # name  index     index   index      index

    for compound_name in compounds:
        c1 = first.header.compounds.index(compound_name)
        c2 = second.header.compounds.index(compound_name)
        for d in range(len(header.dosages)):
            # add all time series for all replicates, then combine each pair
            first_activations = []
            second_activations = []

            # Use at most the first three replicates (for vivo); the others have missing data
            for r1 in range(min(len(first.header.replicates), max_replicates)):
                first_act = []
                for gene in genes_first:
                    first_act.extend(first_hierarchical[gene][c1][d][r1])
                    if len(first_hierarchical[gene][c1][d][r1]) < times_first:
                        print("Missing data! ", gene, compound_name, d, r1)
                        exit(0)

                first_activations.append(first_act)

            for r2 in range(min(len(second.header.replicates), max_replicates)):
                second_act = []
                for gene in genes_second:
                    second_act.extend(second_hierarchical[gene][c2][d][r2])
                    if len(second_hierarchical[gene][c2][d][r2]) < times_second:
                        print("Missing data! ", gene, compound_name, d, r2)
                        exit(0)

                second_activations.append(second_act)

            # combine replicates (i.e. doubling)
            for h, r in itertools.product(first_activations, second_activations):
                X.append(h)
                Y.append(r)

    return X, Y

# Dan: this function is for UDA
def get_doubling_xyz(first: Data, second: Data, third: Data, genes_first=None, genes_second=None, genes_third=None, 
    compounds=None, max_replicates=2) -> Tuple[List[List[float]], List[List[float]]]:
    
    # use all genes if no specfial gene set indicated
    if genes_first is None:
        genes_first = first.header.genes

    if genes_second is None:
        genes_second = genes_first

    if genes_third is None:
        genes_third = genes_second

    print(genes_first)
    print(genes_second)
    print(genes_third)

    first_hierarchical = get_hierarchical_data(first, genes_first)
    second_hierarchical = get_hierarchical_data(second, genes_second)
    third_hierarchical = get_hierarchical_data(third, genes_third)

    header = first.header
    times_first = len(first.header.times)
    times_second = len(second.header.times)
    times_third = len(third.header.times)

    if not compounds:
        compounds = header.compounds

    X = []
    Y = []
    Z = []

    # Data formatted as data[gene][compound][dosage][replicate][time]
    # name  index     index   index      index

    for compound_name in compounds:
        c1 = first.header.compounds.index(compound_name)
        c2 = second.header.compounds.index(compound_name)
        c3 = third.header.compounds.index(compound_name)
        for d in range(len(header.dosages)):
            # add all time series for all replicates, then combine each pair
            first_activations = []
            second_activations = []
            third_activations = []

            # Use at most the first three replicates (for vivo); the others have missing data
            for r1 in range(min(len(first.header.replicates), max_replicates)):
                first_act = []
                for gene in genes_first:
                    first_act.extend(first_hierarchical[gene][c1][d][r1])
                    if len(first_hierarchical[gene][c1][d][r1]) < times_first:
                        print("Missing data! ", gene, compound_name, d, r1)
                        exit(0)

                first_activations.append(first_act)

            for r2 in range(min(len(second.header.replicates), max_replicates)):
                second_act = []
                for gene in genes_second:
                    second_act.extend(second_hierarchical[gene][c2][d][r2])
                    if len(second_hierarchical[gene][c2][d][r2]) < times_second:
                        print("Missing data! ", gene, compound_name, d, r2)
                        exit(0)

                second_activations.append(second_act)

            for r3 in range(min(len(third.header.replicates), max_replicates)):
                third_act = []
                for gene in genes_third:
                    third_act.extend(third_hierarchical[gene][c3][d][r3])
                    if len(third_hierarchical[gene][c3][d][r3]) < times_third:
                        print("Missing data! ", gene, compound_name, d, r3)
                        exit(0)
                third_activations.append(third_act)

            # combine replicates (i.e. doubling)
            for h, r, r2 in itertools.product(first_activations, second_activations, third_activations):
                X.append(h)
                Y.append(r)
                Z.append(r2)

    return X, Y, Z


def get_x(data: Data, genes=None, compounds=None) -> np.ndarray:
    """
    Get gene activations per compound for given dataset and gene list
    :param data: e. g. human data
    :param genes: genes to select, default is all
    :return: X
    """

    # use all genes if no specfial gene set indicated
    if not genes:
        genes = data.header.genes


    hierarchical = get_hierarchical_data(data, genes)
    header = data.header

    # pickle.dump([hierarchical, header], open('small_human.p', 'wb'))
    # exit(0)

    X = []

    if not compounds:
        compounds = header.compounds

    # Data formatted as data[gene][compound][dosage][replicate][time]
    # name  index     index   index      index

    for compound_name in compounds:
        c = header.compounds.index(compound_name)
        for d in range(len(header.dosages)):
            for r in range(len(header.replicates)):
                # collect activations for each replicate, so cross product can be used later on
                act = []
                for gene in genes:
                    # print(gene)
                    for t in range(len(header.times)):
                        # print(compound_name, c, d, r, t, hierarchical[gene][c][d][r])
                        try:
                            act.append(hierarchical[gene][c][d][r][t])
                        except:
                            print("Fix failed miserably again!")
                            import sys
                            sys.exit(1337)
                X.append(np.array(act))
    return np.array(X)


# Read in data from files using big parser or small parser. Only supports X->Y for now (i.e. no 3-way and no autoencoders).
def read_data(compounds, x_type="rat_vitro", y_type="human_vitro", gene_list=None, dataset="big", numb_genes=50, genes_provided=None, orthologs=False, domain="both"):  # Dan: I added many parameters

    """ Dan: explanation of parameters:

    gene_list: string, e.g. STEATOSIS
    genes_provided: list
    domain: should be either "only_x", "only_y", or "both". However, I recommend always using "both" and just ignoring the undesired output if     necessary

    """
    if x_type == y_type:
        raise ValueError("AutoEncoder data not supported yet! x_type and y_type should be different.")

    if dataset == 'small':
        if x_type == "rat_vivo" or y_type == "rat_vivo":
            raise ValueError("Small dataset has no vivo data!")

        parser = SmallDatasetParser()
        humandata, ratdata = parser.parse()
        
        if x_type == 'human_vitro':
            print("Creating data for human_vitro -> rat_vitro using small dataset.")
            X, Y = get_doubling_xy(humandata, ratdata, gene_list)
        else:
            print("Creating data for rat_vitro -> human_vitro using small dataset.")
            X, Y = get_doubling_xy(ratdata, humandata, gene_list)

        genes_first = gene_list
        genes_second = gene_list

        # For autoencoder, 16 should be 8
        data_compounds = np.repeat(compounds, 16)

    else:
        print("Creating data for " + x_type + " -> " + y_type + ", using " + dataset + " dataset." )
        files = {
            'rat_vitro': os.path.join('Data', 'data_rat_vitro.p'),
            'rat_vivo': os.path.join('Data', 'data_rat_vivo.p'),
            'human_vitro': os.path.join('Data', 'data_human_vitro.p')
        }

        print("Loading pickle files")
        X_data = pickle.load(open(files[x_type], 'rb'))
        Y_data = pickle.load(open(files[y_type], 'rb'))

        # Randomly pick 30 genes from both domains (excluding those named None)
        if gene_list == 'random':
            if genes_provided == None:  # Dan: this creates random gene sets for both domains
                if orthologs == False:
                    print("Randomly selecting genes")
                    x_genes = np.unique([gene for gene in X_data.header.genes if not gene is None])
                    genes_first = np.random.choice(x_genes, numb_genes, replace=False).tolist()
                    y_genes = np.unique([gene for gene in Y_data.header.genes if not gene is None])
                    genes_second = np.random.choice(y_genes, numb_genes, replace=False).tolist()

                else:  # Dan: unnested orthologs (e.g. 20)
                    print("It's ortholog time!") 
                    data = np.genfromtxt("Data/ortholog_list.txt",dtype='str',delimiter=",")  # Dan: importing ortholog list
                    data = np.delete(data,0,0)
                    genes_first = np.random.choice(data[:,1], numb_genes, replace=False).tolist()  # Dan: choosing random genes
                    genes_second = []  # Dan: here comes a complicated way of choosing the correspoding human names
                    tempppppppppp = data[:,1].tolist()
                    for i in genes_first:
                        temp = tempppppppppp.index(i)
                        genes_second.append(data[temp,2])

            else:  # Dan: nested case
                if not orthologs:
                    print("Randomly selecting nested genes")
                    intersection = True  # Indicates whether any gene has been selected twice
                    x_genes = np.unique([gene for gene in X_data.header.genes if not gene is None])
                    y_genes = np.unique([gene for gene in Y_data.header.genes if not gene is None])
                    x_satisfied = False
                    y_satisfied = False
                    if domain == "only_x":
                        y_satisfied = True
                        genes_second = np.random.choice(y_genes, numb_genes, replace=False).tolist()
                    if domain == "only_y":
                        x_satisfied = True
                        genes_first = np.random.choice(x_genes, numb_genes, replace=False).tolist() # Dan: this is just to return something
                    while intersection == True:  # Dan: continues as long as there is some overlap
                        while x_satisfied == False:
                            genes_first = np.random.choice(x_genes, numb_genes, replace=False).tolist()
                            genes_intersec = list(set(genes_first).intersection(set(genes_provided)))
                            if len(genes_intersec) == 0:  # Dan: all newly selected genes must be unique
                                intersection = False
                                x_satisfied = True
                        while y_satisfied == False:
                            genes_second = np.random.choice(y_genes, numb_genes, replace=False).tolist()
                            genes_intersec = list(set(genes_second).intersection(set(genes_provided)))
                            if len(genes_intersec) == 0:
                                intersection = False
                                y_satisfied = True

                else: # Dan: nested ortholog case
                    print("Randomly selecting nested orthologs")
                    intersection = True
                    data = np.genfromtxt("Data/ortholog_list.txt",dtype='str',delimiter=",")  # Dan: importing ortholog list
                    data = np.delete(data,0,0)
                    #print("Amount of available orthologs: ", len(data))
                    if domain == "only_x":
                        genes_second = np.random.choice(y_genes, numb_genes, replace=False).tolist()  # Dan: just to return anything
                    if domain == "only_y":
                        genes_first = np.random.choice(x_genes, numb_genes, replace=False).tolist()
                    while intersection == True:  # Dan: continues as long as there is some overlap
                        genes_first = np.random.choice(data[:,1], numb_genes, replace=False).tolist()  # Dan: choosing random genes
                        genes_intersec = list(set(genes_first).intersection(set(genes_provided)))
                        if len(genes_intersec) == 0:  # Dan: newly selected genes must be unique
                            intersection = False
                            genes_second = []  # Dan: choosing the corresponding human names
                            temp1 = data[:,1].tolist()
                            for i in genes_first:
                                temp2 = temp1.index(i)
                                genes_second.append(data[temp2,2])

        else:  # Dan: list of desired genes provided
            gene_mapping = create_gene_name_mapping(gene_list,x_type)
            genes_first = list(gene_mapping.keys())
            gene_mapping = create_gene_name_mapping(gene_list,y_type)  # Dan: I added this
            genes_second = list(gene_mapping.keys())

        """print("Removed ADP and CPZ from compound list due to missing data")
        compounds.remove('ADP')
        compounds.remove('CPZ')"""
            # 2x2 doubling x4 dosages = 16 times the same compound in a row, unless using 3 replicates for vivo
        data_compounds = np.repeat(compounds, 16)
        #  else:
            # 2x2 doubling x4 dosages = 16 times the same compound in a row
           #  data_compounds = np.repeat(compounds, 16)
        print("Parsing data")
        X, Y = get_doubling_xy(X_data, Y_data, genes_first, genes_second, compounds)

    X = np.array(X)
    Y = np.array(Y)
    print("X shape: {}".format(X.shape))
    print("Y shape: {}".format(Y.shape))
    print("")

    return X, Y, data_compounds, genes_first, genes_second

# Takes orthologs as an input and selects additional ones to create a bigger 'nest'
def read_nested_orths(compounds, y_type="human_vitro", gene_list=None, numb_genes=50, genes_provided=None, orthologs=False, domain="both"):
    # TODO merge this with previous function (i.e. make the other one read orthologs)
    """ 
    gene_list: string, e.g. STEATOSIS
    genes_provided: list
    domain: can be either "only_x", "only_y", or "both" (only relevant for orthologs)
    """
    x_type = "rat_vitro"
    if x_type == y_type:
        raise ValueError("AutoEncoder data not supported yet! x_type and y_type should be different.")

    else:
        print("Creating data for " + x_type + " -> " + y_type + ", using big dataset." )
        files = {
            'rat_vitro': os.path.join('Data', 'data_rat_vitro.p'),
            'rat_vivo': os.path.join('Data', 'data_rat_vivo.p'),
            'human_vitro': os.path.join('Data', 'data_human_vitro.p')
        }

        print("Loading pickle files")
        X_data = pickle.load(open(files[x_type], 'rb'))
        Y_data = pickle.load(open(files[y_type], 'rb'))

        # Randomly pick genes from both domains (excluding those named None)
        if gene_list == 'random':
            if not orthologs:  # Dan: this ignores overlap in the selected genes for now
                print("Randomly selecting genes")
                x_genes = np.unique([gene for gene in X_data.header.genes if not gene is None])
                genes_first = np.random.choice(x_genes, numb_genes, replace=False).tolist()
                print("D1 genes:", genes_first)
                y_genes = np.unique([gene for gene in Y_data.header.genes if not gene is None])
                genes_second = np.random.choice(y_genes, numb_genes, replace=False).tolist()
                print("D2 genes:", genes_second)
            else:
                print("It's ortholog time!") 
                if genes_provided == None:
                    data = np.genfromtxt("Data/ortholog_list.txt",dtype='str',delimiter=",")  # Dan: importing ortholog list
                    data = np.delete(data,0,0)
                    genes_first = np.random.choice(data[:,1], numb_genes, replace=False).tolist()  # Dan: choosing random genes
                    genes_second = []  # Dan: here comes a complicated way of choosing the correspoding human names
                    tempppppppppp = data[:,1].tolist()
                    for i in genes_first:
                        temp = tempppppppppp.index(i)
                        genes_second.append(data[temp,2])

                else: # Dan: nested case
                    print("Randomly selecting nested orthologs")
                    intersection = True
                    data = np.genfromtxt("Data/ortholog_list.txt",dtype='str',delimiter=",")  # Dan: importing ortholog list
                    data = np.delete(data,0,0)
                    #print("Amount of available orthologs: ", len(data))
                    if domain == "only_x":
                        genes_second = np.random.choice(y_genes, numb_genes, replace=False).tolist()  # Dan: just to return anything
                    if domain == "only_y":
                        genes_first = np.random.choice(x_genes, numb_genes, replace=False).tolist()
                    while intersection == True:  # Dan: continues as long as there is some overlap
                        genes_first = np.random.choice(data[:,1], numb_genes, replace=False).tolist()  # Dan: choosing random genes
                        genes_intersec = list(set(genes_first).intersection(set(genes_provided)))
                        if len(genes_intersec) == 0:  # Dan: newly selected genes must be unique
                            intersection = False
                            genes_second = []  # Dan: choosing the corresponding human names
                            temp1 = data[:,1].tolist()
                            for i in genes_first:
                                temp2 = temp1.index(i)
                                genes_second.append(data[temp2,2])
        else:                    # Dan: gene list provided
            gene_mapping = create_gene_name_mapping(gene_list,x_type)
            genes_first = list(gene_mapping.keys())
            gene_mapping = create_gene_name_mapping(gene_list,y_type)  # Dan: I added this
            genes_second = list(gene_mapping.keys())

        """ Dan: make sure that correct amount of compounds are provided
        print("Removed ADP and CPZ from compound list due to missing data")
        compounds.remove('ADP')
        compounds.remove('CPZ')"""
        data_compounds = np.repeat(compounds, 16)
            # 2x2 doubling x4 dosages = 16 times the same compound in a row
         
        print("Parsing data")  # Dan: this is where the actual data is read (so far only the names of genes have been selected)
        X, Y = get_doubling_xy(X_data, Y_data, genes_first, genes_second, compounds)

    X = np.array(X)
    Y = np.array(Y)
    print("X shape: {}".format(X.shape))
    print("Y shape: {}".format(Y.shape))
    print("")

    return X, Y, data_compounds, genes_first, genes_second


# Read in data from files using big parser or small parser. Only supports X->Y for now (i.e. no 3-way and no autoencoders).
# Dan: this is for UDA
def read_data_xyz(compounds, x_type="rat_vitro", y_type="human_vitro", z_type="rat_vivo", gene_list=None, dataset="big"):
    files = {  # Dan: why doesn't it support 3-way? Isn't that what this is for?
        'rat_vitro': os.path.join('Data', 'data_rat_vitro.p'),
        'rat_vivo': os.path.join('Data', 'data_rat_vivo.p'),
        'human_vitro': os.path.join('Data', 'data_human_vitro.p')
    }

    print("Loading pickle files")
    X_data = pickle.load(open(files[x_type], 'rb'))
    Y_data = pickle.load(open(files[y_type], 'rb'))
    Z_data = pickle.load(open(files[z_type], 'rb'))

    # Randomly pick 30 genes from both domains (excluding those named None)
    if gene_list == 'random':
        print("Randomly selecting genes")
        x_genes = np.unique([gene for gene in X_data.header.genes if not gene is None])
        genes_first = np.random.choice(x_genes, 30, replace=False).tolist()
        print("D1 genes:", genes_first)
        y_genes = np.unique([gene for gene in Y_data.header.genes if not gene is None])
        genes_second = np.random.choice(y_genes, 30, replace=False).tolist()
        print("D2 genes:", genes_second)
        z_genes = np.unique([gene for gene in Z_data.header.genes if not gene is None])
        genes_third = np.random.choice(z_genes, 30, replace=False).tolist()
    else:  # Dan: what if gene_list stays None?
        gene_mapping = create_gene_name_mapping(gene_list)

        # Only human data uses uppercase gene names
        if x_type != 'human_vitro':
            genes_first = list(gene_mapping.keys())
        else:
            genes_first = list(gene_mapping.values())
        
        if y_type != 'human_vitro':
            genes_second = list(gene_mapping.keys())
        else:
            genes_second = list(gene_mapping.values())

        if z_type != 'human_vitro':
            genes_third = list(gene_mapping.keys())
        else:
            genes_third = list(gene_mapping.values())


    if x_type == 'rat_vivo' or y_type == 'rat_vivo' or z_type == 'rat_vivo':
        print("Removed ADP and CPZ from compound list due to missing data")
        if 'ADP' in compounds:
            compounds.remove('ADP')
        if 'CPZ' in compounds:
            compounds.remove('CPZ')

    data_compounds = np.repeat(compounds, 32)

    print("Parsing data")
    X, Y, Z = get_doubling_xyz(X_data, Y_data, Z_data, genes_first, genes_second, genes_third, compounds)

    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)

    print("X shape: {}".format(X.shape))
    print("Y shape: {}".format(Y.shape))
    print("Z shape: {}".format(Z.shape))
    print("")

    return X, Y, Z, data_compounds, genes_first, genes_second, genes_third


def main():  # Dan: this can be ignored as this script only 'helps' the RandomGenes.py file
    parser = SmallDatasetParser()
   # humandata, ratdata = parser.parse() Dan: I commented the following lines out

    #X, Y = get_doubling_xy(humandata, ratdata, ['FABP4', 'ACACA'])
    #print("CleanParser: len(X)={}, len(Y)={}".format(len(X), len(Y)))

    #print("Len X {}, Len Y {}".format(len(X), len(Y)))
    compounds = CompoundLists.GENERAL_47
    gene_list = GeneLists.GTX
    a,b,c,d,e,h,g = read_data_xyz(compounds, x_type="rat_vitro", y_type="human_vitro", z_type="rat_vivo", gene_list=gene_list, dataset="big")
    with open("STEATOSIS_xyz.p", 'wb') as f:
        pickle.dump([a,b,c,d,e,h,g], f)

if __name__ == '__main__':
    main()
