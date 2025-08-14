# -*- coding: utf-8 -*-

import numpy
import math
import random
import argparse
from data_generation import make_shift_scale_data
#import bibench
#import bibench.all as bb


def gettrend2(cols=100, rows=1000, clustrows=100, clustcols=10,
              clusts=3, noise=0, base_loc=0, base_scale=1,
              laprow=0, lapcol=0, shuffle=False):

    data, expected = make_shift_scale_data(
        ncols=cols, nrows=rows,
        nclustrows=clustrows, nclustcols=clustcols,
        nclusts=clusts, base_loc=base_loc, base_scale=base_scale,
        noverlap_rows=laprow, noverlap_cols=lapcol,
        shuffle=shuffle
    )

    return data, expected


def noise_data(cols=100, rows=1000, clustrows=30, clustcols=20,
               clusts=3, noise=0, shuffle=False, LR=0, LC=0):
    data, expected = gettrend2(cols=cols, rows=rows,
                               clustrows=clustrows, clustcols=clustcols,
                               clusts=clusts, noise=noise,
                               base_loc=0, base_scale=1,
                               laprow=LR, lapcol=LC, shuffle=shuffle)
    ncols = len(expected[0].cols)
    nrows = len(expected[0].rows)
    n_d = math.sqrt(noise)
    noise_level = n_d
    noise_rnum = int(noise_level * nrows)
    noise_cnum = int(noise_level * ncols)
    noise_rindex = []
    noise_cindex = []

    for ele in expected:
        rnoise = random.sample(ele.rows, noise_rnum)
        cnoise = random.sample(ele.cols, noise_cnum)
        noise_rindex.append(rnoise)
        noise_cindex.append(cnoise)

    for i in range(clusts):
        for ele2 in noise_rindex[i]:
            for ele3 in noise_cindex[i]:
                data[ele2][ele3] = data[ele2][ele3] + numpy.random.normal(0, 0.1)

    return data, expected


def write(expected, file_path):
    with open(file_path, 'w') as fo:
        for ele in expected:
            for ele2 in ele.rows:
                fo.write(str(ele2) + " ")
            fo.write("\n")
            for ele2 in ele.cols:
                fo.write(str(ele2) + " ")
            fo.write("\n\n")


def write_data(data, file_path):
    flag = 0
    with open(file_path, 'w') as fo:
        fo.write('o' + "\t")
        for row in data:
            if flag == 0:
                i = 0
                for ele in row:
                    fo.write("cond" + str(i) + "\t")
                    i += 1
                flag = 1
                fo.write("\n")
            x = flag - 1
            fo.write("gene" + str(x) + "\t")
            flag += 1
            for ele in row:
                fo.write(str(ele) + "\t")
            fo.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate bicluster data (Python 2 compatible)")
    parser.add_argument("--output_data", type=str, required=True,
                        help="Path to save the generated data matrix")
    parser.add_argument("--output_clusters", type=str, required=True,
                        help="Path to save the generated cluster indices")
    parser.add_argument("--cols", type=int, default=200,
                        help="Number of columns in the data matrix")
    parser.add_argument("--rows", type=int, default=200,
                        help="Number of rows in the data matrix")
    parser.add_argument("--clustrows", type=int, default=20,
                        help="Number of rows per bicluster")
    parser.add_argument("--clustcols", type=int, default=20,
                        help="Number of columns per bicluster")
    parser.add_argument("--clusts", type=int, default=3,
                        help="Number of biclusters to generate")
    parser.add_argument("--noise", type=float, default=0,
                        help="Noise level (0 means no noise)")
    parser.add_argument("--shuffle", action="store_true",
                        help="Shuffle rows and columns in the data matrix")

    args = parser.parse_args()

    data, clusters = noise_data(
        cols=args.cols,
        rows=args.rows,
        clustrows=args.clustrows,
        clustcols=args.clustcols,
        clusts=args.clusts,
        noise=args.noise,
        shuffle=args.shuffle
    )

    write_data(data, args.output_data)
    write(clusters, args.output_clusters)
    print("Data file saved to: {}".format(args.output_data))
    print("Cluster file saved to: {}".format(args.output_clusters))
