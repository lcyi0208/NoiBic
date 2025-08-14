#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

"""
Standalone synthetic bicluster data generator (legacy-compatible, no rpy2 required).

Implements the "shift-scale" bicluster model:
- background ~ dist(loc=background_loc, scale=background_scale)
- inside each bicluster: data = outer(rowscale, colbase) + rowshift (broadcasted)
- optional global noise and per-bicluster noise (Gaussian)
- optional row/column shuffling (preserving biclusters)

Only depends on numpy and Python stdlib.

Example:
    python synthetic_bicluster_legacy.py --nrows 300 --ncols 50 --nclusts 3 --noise 0.1 --seed 42
"""

import argparse
import random

import numpy as np


# ----------------------------- Utilities -----------------------------

def improper_normal(loc=0.0, scale=1.0, size=None):
    """
    Same as numpy.random.normal, but if scale is 0, returns just the mean 'loc'.
    """
    if scale == 0:
        if size is None:
            return loc
        else:
            return np.zeros(size) + loc
    return np.random.normal(loc, scale, size)


class Bicluster(object):
    def __init__(self, rows, cols, data=None):
        self.rows = list(rows)
        self.cols = list(cols)
        self.data = data

    def shape(self):
        return (len(self.rows), len(self.cols))

    def area(self):
        return len(self.rows) * len(self.cols)

    def __repr__(self):
        return "Bicluster(rows={}, cols={})".format(len(self.rows), len(self.cols))


# ----------------------- Core generation helpers ----------------------

def _add_bicluster_noise_(data, rowmatrix, colmatrix, stdevs):
    """
    Adds noise from Normal(0, stdevs[i]) to bicluster i in 'data'.
    rowmatrix: (nrows, nclusts) 0/1; colmatrix: (ncols, nclusts) 0/1
    """
    noisy = data.copy()
    for j in range(rowmatrix.shape[1]):
        rows_vec = rowmatrix[:, j]
        cols_vec = colmatrix[:, j]
        scale = float(stdevs[j]) if stdevs is not None else 0.0
        if scale > 0:
            rows_mask = rows_vec > 0
            cols_mask = cols_vec > 0
            shape = (int(rows_mask.sum()), int(cols_mask.sum()))
            block_noise = np.random.normal(loc=0.0, scale=scale, size=shape)
            noisy[np.ix_(rows_mask, cols_mask)] += block_noise
    return noisy


def add_noise(data, scale):
    """Adds Normal(0, scale) iid noise to the whole matrix."""
    return data + improper_normal(scale=scale, size=data.shape)


def _shuffle_(data, expected, new_rows=None, new_cols=None):
    """
    Shuffles rows/columns of the dataset while preserving biclusters.
    """
    nrows, ncols = data.shape
    if new_rows is None:
        new_rows = list(range(nrows))
        random.shuffle(new_rows)
    if new_cols is None:
        new_cols = list(range(ncols))
        random.shuffle(new_cols)

    shuffled_data = data[new_rows][:, new_cols]

    inv_row = dict((old, new) for new, old in enumerate(new_rows))
    inv_col = dict((old, new) for new, old in enumerate(new_cols))

    shuffled_biclusters = []
    for b in expected:
        new_b_rows = [inv_row[r] for r in b.rows]
        new_b_cols = [inv_col[c] for c in b.cols]
        shuffled_biclusters.append(Bicluster(new_b_rows, new_b_cols, shuffled_data))

    return shuffled_data, shuffled_biclusters


def _make_row_matrix_(nrows, nclusts, nclust_rows, noverlap_rows):
    """
    Row-membership matrix M (nrows x nclusts) where M[i, j] == 1 if row i in bicluster j.
    Biclusters are placed consecutively with optional row overlaps (fixed overlap size between neighbors).
    """
    if nclust_rows + (nclust_rows - noverlap_rows) * (nclusts - 1) > nrows:
        raise ValueError("Biclusters are too large to fit in the dataset with the requested overlaps.")
    M = np.zeros((nrows, nclusts), dtype=np.int8)
    for j in range(nclusts):
        start = (nclust_rows * j) - (noverlap_rows * j)
        stop = start + nclust_rows
        M[start:stop, j] = 1
    return M


def _make_col_matrix_(ncols, nclusts, nclust_cols, noverlap_cols):
    """Column-membership matrix analogous to _make_row_matrix_."""
    return _make_row_matrix_(ncols, nclusts, nclust_cols, noverlap_cols)


def _make_expected_biclusters_(row_matrix, col_matrix, data):
    """Return a list of Bicluster objects derived from the membership matrices."""
    biclusters = []
    nclusts = row_matrix.shape[1]
    for j in range(nclusts):
        row_line = row_matrix[:, j]
        col_line = col_matrix[:, j]
        rows = list(np.where(row_line > 0)[0])
        cols = list(np.where(col_line > 0)[0])
        biclusters.append(Bicluster(rows, cols, data))
    return biclusters


def _make_biclusters_(row_matrix, col_matrix, bicluster_colbase, bicluster_rowshift, bicluster_rowscale):
    """
    Build the full data matrix for all biclusters (no background yet) using:
        block = outer(rowscale, colbase) + rowshift (broadcasted down columns)
    Biclusters occupy rectangles indicated by row_matrix and col_matrix.
    """
    nrows, nclusts = row_matrix.shape
    ncols, nclusts_check = col_matrix.shape
    assert nclusts == nclusts_check

    data = np.zeros((nrows, ncols), dtype=float)
    for j in range(nclusts):
        row = row_matrix[:, j]
        col = col_matrix[:, j]
        base = bicluster_colbase[:, j]
        shift = bicluster_rowshift[:, j]
        scale = bicluster_rowscale[:, j]
        block = np.outer(scale, base) + shift.reshape(-1, 1)
        mask = np.outer(row, col) > 0
        data[mask] = block[mask[np.ix_(row > 0, col > 0)] if False else block == block]
        # 上面一行保持与原实现一致；实质就是复制对应矩形区域
        # 更直观也可写成：
        # rows_mask = row > 0; cols_mask = col > 0
        # data[np.ix_(rows_mask, cols_mask)] = block
    return data


def _do_both_overlap_(noverlap_rows, noverlap_cols):
    return (noverlap_rows > 0) and (noverlap_cols > 0)


def _general_model_(nrows, ncols,
                    nclusts, nclustrows, nclustcols,
                    noverlap_rows, noverlap_cols,
                    colbase, rowshift, rowscale,
                    background_generator):
    """
    Compose background + biclusters given parameters and membership matrices.
    """
    k = np.identity(nclusts)
    row_matrix = _make_row_matrix_(nrows, nclusts, nclustrows, noverlap_rows)
    col_matrix = _make_col_matrix_(ncols, nclusts, nclustcols, noverlap_cols)

    # Boolean layout of bicluster rectangles (avoid '@' for old Python)
    layout = np.dot(np.dot(row_matrix, k), col_matrix.T)
    bicluster_bools = layout > 0
    background_bools = ~bicluster_bools

    background = background_generator(nrows, ncols)

    if _do_both_overlap_(noverlap_rows, noverlap_cols):
        # both-overlap case
        biclusters = np.outer(rowscale, colbase) + rowshift.reshape(-1, 1)
    else:
        biclusters = _make_biclusters_(row_matrix, col_matrix, colbase, rowshift, rowscale)

    data = np.empty((nrows, ncols), dtype=float)
    data[bicluster_bools] = biclusters[bicluster_bools]
    data[background_bools] = background[background_bools]
    return data, row_matrix, col_matrix


def _set_defaults_(nrows, ncols, nclusts, nclustrows, nclustcols, bicluster_noise):
    if nclustrows is None:
        nclustrows = nrows // nclusts
        if nclustrows == nrows:
            nclustrows = nrows // 2
    if nclustcols is None:
        nclustcols = ncols // nclusts
        if nclustcols == ncols:
            nclustcols = ncols // 2
    if bicluster_noise is None:
        bicluster_noise = [0.0] * nclusts
    return nclustrows, nclustcols, bicluster_noise


def _make_data_(nrows, ncols, nclusts,
                noverlap_rows, noverlap_cols,
                colbase, rowshift, rowscale,
                noise, nclustrows=None, nclustcols=None, bicluster_noise=None,
                background_loc=0.0, background_scale=1.0, shuffle=False,
                dist=improper_normal):
    """Apply the model and add noise."""
    background_generator = lambda x, y: dist(loc=background_loc, scale=background_scale, size=(x, y))

    data, row_matrix, col_matrix = _general_model_(nrows, ncols, nclusts,
                                                   nclustrows, nclustcols,
                                                   noverlap_rows, noverlap_cols,
                                                   colbase, rowshift, rowscale,
                                                   background_generator)

    if noise and noise > 0:
        data = add_noise(data, scale=noise)

    # add extra noise per bicluster
    data = _add_bicluster_noise_(data, row_matrix, col_matrix, bicluster_noise)

    expected = _make_expected_biclusters_(row_matrix, col_matrix, data)

    if shuffle:
        data, expected = _shuffle_(data, expected)

    return data, expected


def make_shift_scale_data(nrows=300, ncols=50,
                          nclusts=3, nclustrows=None, nclustcols=None,
                          noverlap_rows=0, noverlap_cols=0,
                          background_loc=0.0, background_scale=1.0,
                          base_loc=0.0, base_scale=1.0,
                          shift_loc=0.0, shift_scale=1.0,
                          scale_loc=0.0, scale_scale=1.0,
                          noise=0.0, bicluster_noise=None,
                          shuffle=False,
                          dist=improper_normal):
    """
    Create a synthetic dataset with shift-scale biclusters.
    Returns:
        data: np.ndarray [nrows x ncols]
        expected: list of Bicluster
    """
    nclustrows, nclustcols, bicluster_noise = _set_defaults_(nrows, ncols, nclusts,
                                                             nclustrows, nclustcols,
                                                             bicluster_noise)

    # Draw parameters for biclusters
    if _do_both_overlap_(noverlap_rows, noverlap_cols):
        colbase = dist(loc=base_loc,  scale=base_scale,  size=ncols)
        rowshift = dist(loc=shift_loc, scale=shift_scale, size=nrows)
        rowscale = dist(loc=scale_loc, scale=scale_scale, size=nrows)
    else:
        colbase = dist(loc=base_loc,  scale=base_scale,  size=(nclustcols, nclusts))
        rowshift = dist(loc=shift_loc, scale=shift_scale, size=(nclustrows, nclusts))
        rowscale = dist(loc=scale_loc, scale=scale_scale, size=(nclustrows, nclusts))

    return _make_data_(nrows, ncols, nclusts,
                       noverlap_rows, noverlap_cols,
                       colbase, rowshift, rowscale,
                       noise, nclustrows, nclustcols, bicluster_noise,
                       background_loc, background_scale, shuffle, dist)


# ------------------------------ CLI ----------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Generate synthetic bicluster data (shift-scale model).")
    p.add_argument("--nrows", type=int, default=300)
    p.add_argument("--ncols", type=int, default=50)
    p.add_argument("--nclusts", type=int, default=3)
    p.add_argument("--nclustrows", type=int, default=None)
    p.add_argument("--nclustcols", type=int, default=None)
    p.add_argument("--noverlap_rows", type=int, default=0)
    p.add_argument("--noverlap_cols", type=int, default=0)

    p.add_argument("--background_loc", type=float, default=0.0)
    p.add_argument("--background_scale", type=float, default=1.0)

    p.add_argument("--base_loc", type=float, default=0.0)
    p.add_argument("--base_scale", type=float, default=1.0)
    p.add_argument("--shift_loc", type=float, default=0.0)
    p.add_argument("--shift_scale", type=float, default=1.0)
    p.add_argument("--scale_loc", type=float, default=0.0)
    p.add_argument("--scale_scale", type=float, default=1.0)

    p.add_argument("--noise", type=float, default=0.0, help="Global iid Gaussian noise stdev.")
    p.add_argument("--bicluster_noise", type=float, nargs='*', default=None,
                   help="Per-bicluster noise stdevs; defaults to zeros.")
    p.add_argument("--shuffle", action="store_true", help="Shuffle rows/cols after generation.")
    p.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    p.add_argument("--save_prefix", type=str, default="synthetic_legacy")

    return p.parse_args()


def main():
    args = parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)
        random.seed(args.seed)

    if args.bicluster_noise is None:
        bic_noise = None
    else:
        if len(args.bicluster_noise) == 1:
            bic_noise = [args.bicluster_noise[0]] * args.nclusts
        else:
            if len(args.bicluster_noise) != args.nclusts:
                raise ValueError("--bicluster_noise length must be 1 or equal to nclusts")
            bic_noise = args.bicluster_noise

    data, expected = make_shift_scale_data(
        nrows=args.nrows, ncols=args.ncols,
        nclusts=args.nclusts, nclustrows=args.nclustrows, nclustcols=args.nclustcols,
        noverlap_rows=args.noverlap_rows, noverlap_cols=args.noverlap_cols,
        background_loc=args.background_loc, background_scale=args.background_scale,
        base_loc=args.base_loc, base_scale=args.base_scale,
        shift_loc=args.shift_loc, shift_scale=args.shift_scale,
        scale_loc=args.scale_loc, scale_scale=args.scale_scale,
        noise=args.noise, bicluster_noise=bic_noise,
        shuffle=args.shuffle,
        dist=improper_normal
    )

    prefix = args.save_prefix
    np.save(prefix + "_data.npy", data)
    np.savetxt(prefix + "_data.tsv", data, fmt="%.6f", delimiter="\t")

    with open(prefix + "_biclusters.txt", "w") as f:
        for b in expected:
            f.write(" ".join(map(str, b.rows)) + "\n")
            f.write(" ".join(map(str, b.cols)) + "\n\n")

    nrows, ncols = data.shape
    Row = np.zeros((nrows, len(expected)), dtype=int)
    Col = np.zeros((ncols, len(expected)), dtype=int)
    for j, b in enumerate(expected):
        Row[b.rows, j] = 1
        Col[b.cols, j] = 1
    np.savetxt(prefix + "_row_membership.tsv", Row, fmt="%d", delimiter="\t")
    np.savetxt(prefix + "_col_membership.tsv", Col, fmt="%d", delimiter="\t")

    print("Data shape: {}".format(data.shape))
    for i, b in enumerate(expected):
        print("  Bicluster {}: rows={}, cols={}, area={}".format(i, len(b.rows), len(b.cols), b.area()))
    print("Saved: {}_data.npy/.tsv, {}_biclusters.txt, {}_row_membership.tsv, {}_col_membership.tsv"
          .format(prefix, prefix, prefix, prefix))


if __name__ == "__main__":
    main()
