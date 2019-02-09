#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from keras.layers import Input, Dense
from keras.models import Model
from matplotlib.mlab import PCA

from Processing import get_all_gene_activations, normalize_nparray
from SmallParser import SmallDatasetParser


# import tensorflow as tf


def main():
    smallParser = SmallDatasetParser()
    humandata, ratdata = smallParser.parse()

    geneSelection_easy = ['ACACA', 'CPT2', 'CPT1A', 'NR1H4', 'NFKB2', 'RELA', 'RELB', 'PPARA', 'SCD', 'TGFB1']
    x = get_all_gene_activations(humandata)
    x, min_x, max_x = normalize_nparray(x)

    pca(x)


def pca(x):
    x_pca = PCA(x)
    print("variance percentages from PCA: {}".format(x_pca.fracs))

    # http://blog.nextgenetics.net/?e=42
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    x = []
    y = []
    z = []
    for item in x_pca.Y:
        x.append(item[0])
        y.append(item[1])
        z.append(item[2])

    plt.close('all')  # close all latent plotting windows
    fig1 = plt.figure()  # Make a plotting figure
    ax = Axes3D(fig1)  # use the plotting figure to create a Axis3D object.
    pltData = [x, y, z]
    ax.scatter(pltData[0], pltData[1], pltData[2], 'bo')  # make a scatter plot of blue dots from the data

    # make simple, bare axis lines through space:
    xAxisLine = ((min(pltData[0]), max(pltData[0])), (0, 0),
                 (0, 0))  # 2 points make the x-axis line at the data extrema along x-axis
    ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')  # make a red line for the x-axis.
    yAxisLine = ((0, 0), (min(pltData[1]), max(pltData[1])),
                 (0, 0))  # 2 points make the y-axis line at the data extrema along y-axis
    ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')  # make a red line for the y-axis.
    zAxisLine = ((0, 0), (0, 0),
                 (min(pltData[2]), max(pltData[2])))  # 2 points make the z-axis line at the data extrema along z-axis
    ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')  # make a red line for the z-axis.

    # label the axes
    ax.set_xlabel("x-axis label")
    ax.set_ylabel("y-axis label")
    ax.set_zlabel("y-axis label")
    ax.set_title("The title of the plot")
    plt.show()  # show the plot


def autoenc(x):
    shape = x.shape
    encoding_dim = 20
    encoding_dim_2 = 10
    encoding_dim_3 = 7
    input = Input(shape=(x.shape[1],))
    # "encoded" is the encoded representation of the input
    encoded1 = Dense(encoding_dim, activation='relu')(input)
    encoded2 = Dense(encoding_dim_2, activation='relu')(encoded1)
    # encoded3 = Dense(encoding_dim_3, activation='relu')(encoded2)
    # decoded3 = Dense(encoding_dim_2, activation='relu')(encoded)
    decoded2 = Dense(encoding_dim, activation='relu')(encoded2)
    # forces sparse encoding
    # encoded = Dense(encoding_dim, activation='relu',
    #                activity_regularizer=regularizers.l1(10e-5))(input)
    # "decoded" is the lossy reconstruction of the input
    decoded1 = Dense(x.shape[1], activation='sigmoid')(decoded2)
    # this model maps an input to its reconstruction
    autoencoder = Model(input, decoded1)
    autoencoder.compile(optimizer='adam', loss='mean_squared_error')
    # split in training and validation data
    last_train_index = int(x.shape[0] * 0.75)
    np.random.shuffle(x)
    x_train = x[:last_train_index]
    x_test = x[last_train_index:]
    autoencoder.fit(x_train, x_train,
                    epochs=1000,
                    shuffle=True,
                    validation_data=(x_test, x_test),
                    # callbacks=[
                    #    EarlyStopping(monitor='val_loss', min_delta=0.0000000001, patience=100, verbose=0, mode='auto')]
                    )
    reconstructed_x_train = autoencoder.predict(x_train)
    reconstructed_x_test = autoencoder.predict(x_test)
    errors_train = np.absolute(x_train - reconstructed_x_train)
    errors_test = np.absolute(x_test - reconstructed_x_test)
    print("Average training error (not loss): {}".format(np.average(errors_train)))
    print("Average test error (not loss): {}".format(np.average(errors_test)))


if __name__ == '__main__':
    main()
