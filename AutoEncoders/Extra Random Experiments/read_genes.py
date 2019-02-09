import pickle
import numpy as np

a, b, genes, variances = pickle.load(open('../../Data/cnn_rat_human_vitro.p', 'rb')) # Dan: this replaces following line

#_ , _, genes, variances = pickle.load(open('data_Y1_human_vitro.p', 'rb'))


print(genes)
print(variances.mean())