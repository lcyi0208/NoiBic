# -*- coding: utf-8 -*-
"""
Generate BiBench synthetic datasets for biclustering experiments.

Output root:
    /home/longchaoyi/biclustering_tool/syndata_bibench

Each dataset writes:
    dataset_XX_data.tsv
    dataset_XX_bics.txt
    dataset_XX_bics.json
"""

from __future__ import print_function

import json
import os
import random

import numpy


ROOT = "/home/longchaoyi/biclustering_tool/syndata_bibench"
REPLICATES = 10
DEFAULT_BIC_NOISE = 0.10
NOISE_SD = 1
MIN_VALUE = -20.0
MAX_VALUE = 20.0


def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)


class Bicluster(object):
    def __init__(self, rows, cols, data=None):
        self.rows = rows
        self.cols = cols
        self.data = data


def clip(value):
    return max(MIN_VALUE, min(MAX_VALUE, value))


def clipped_gaussian(size):
    return numpy.clip(numpy.random.normal(0, 1, size), MIN_VALUE, MAX_VALUE)


def make_expected(rows, cols, bic_rows, bic_cols, bic_num,
                  overlap_rows=0, overlap_cols=0):
    if bic_rows + (bic_rows - overlap_rows) * (bic_num - 1) > rows:
        raise ValueError("biclusters are too large for dataset rows")
    if bic_cols + (bic_cols - overlap_cols) * (bic_num - 1) > cols:
        raise ValueError("biclusters are too large for dataset columns")

    expected = []
    for index in range(bic_num):
        row_start = index * (bic_rows - overlap_rows)
        col_start = index * (bic_cols - overlap_cols)
        bic_r = list(range(row_start, row_start + bic_rows))
        bic_c = list(range(col_start, col_start + bic_cols))
        expected.append(Bicluster(bic_r, bic_c))
    return expected


def increasing_values(size, low, high):
    values = numpy.random.uniform(low, high, size)
    values.sort()
    return values


def plant_constant(data, bic):
    value = random.uniform(5.0, 10.0)
    for row in bic.rows:
        for col in bic.cols:
            data[row][col] = value


def plant_additive(data, bic):
    """
    Column-trend additive pattern.

    Every row has the same ordered column relationship:
        value(row, col_j) = row_shift(row) + additive_factor(j)
    """
    factors = increasing_values(len(bic.cols), -8.0, 8.0)
    for row in bic.rows:
        row_shift = random.uniform(-4.0, 4.0)
        for offset, col in enumerate(bic.cols):
            data[row][col] = clip(row_shift + factors[offset])


def plant_multiplicative(data, bic):
    """
    Column-trend multiplicative pattern.

    Every row has the same ordered column relationship:
        value(row, col_j) = row_scale(row) * multiplicative_factor(j)
    """
    factors = increasing_values(len(bic.cols), 1.0, 8.0)
    for row in bic.rows:
        row_scale = random.uniform(0.5, 3.0)
        for offset, col in enumerate(bic.cols):
            data[row][col] = clip(row_scale * factors[offset])


def plant_shift_scale(data, bic):
    """
    Column-trend shift-scale pattern used as order-preserving data.

    Every row is a shifted and scaled version of the same increasing
    column profile, so column order is preserved within each row.
    """
    base = increasing_values(len(bic.cols), -8.0, 8.0)
    for row in bic.rows:
        row_shift = random.uniform(-4.0, 4.0)
        row_scale = random.uniform(0.5, 2.5)
        for offset, col in enumerate(bic.cols):
            data[row][col] = clip(row_shift + row_scale * base[offset])


def make_data(kind, rows, cols, bic_rows, bic_cols, bic_num,
              overlap_rows=0, overlap_cols=0, shuffle=False):
    data = clipped_gaussian((rows, cols))
    expected = make_expected(rows, cols, bic_rows, bic_cols, bic_num,
                             overlap_rows, overlap_cols)

    planters = {
        "constant": plant_constant,
        "additive": plant_additive,
        "multiplicative": plant_multiplicative,
        "shift_scale": plant_shift_scale,
    }
    if kind not in planters:
        raise ValueError("Unknown kind: {0}".format(kind))

    for bic in expected:
        planters[kind](data, bic)

    for bic in expected:
        bic.data = data

    if shuffle:
        data, expected = shuffle_data(data, expected)

    return data, expected


def shuffle_data(data, expected):
    row_order = list(range(data.shape[0]))
    col_order = list(range(data.shape[1]))
    random.shuffle(row_order)
    random.shuffle(col_order)

    shuffled = data[row_order].T[col_order].T
    shuffled_expected = []
    for bic in expected:
        new_rows = [row_order.index(r) for r in bic.rows]
        new_cols = [col_order.index(c) for c in bic.cols]
        shuffled_expected.append(Bicluster(new_rows, new_cols, shuffled))
    return shuffled, shuffled_expected


def add_bicluster_noise(data, expected, noise_ratio=DEFAULT_BIC_NOISE,
                        noise_sd=NOISE_SD):
    """
    Add N(0, noise_sd) noise to about noise_ratio of cells in each bicluster.

    The old script selected sqrt(noise_ratio) rows and sqrt(noise_ratio)
    columns, whose cross product is about noise_ratio of the bicluster cells.
    This version directly samples cells, so 10%, 20%, 30%, 40% are exact up to
    rounding and are easier to interpret.
    """
    noisy = numpy.array(data, copy=True)
    noisy_cells_by_bic = []

    for bic in expected:
        cells = [(r, c) for r in bic.rows for c in bic.cols]
        n_noisy = int(round(len(cells) * noise_ratio))
        selected = random.sample(cells, min(n_noisy, len(cells)))

        for row, col in selected:
            noisy[row][col] = clip(noisy[row][col] + numpy.random.normal(0, noise_sd))

        noisy_cells_by_bic.append(selected)

    for bic in expected:
        bic.data = noisy

    return noisy, expected, noisy_cells_by_bic


def write_data(data, file_path):
    with open(file_path, "w") as fo:
        fo.write("X")
        for col in range(data.shape[1]):
            fo.write("\tcond{0}".format(col))
        fo.write("\n")

        for row_index, row in enumerate(data):
            fo.write("gene{0}".format(row_index))
            for value in row:
                fo.write("\t{0}".format(value))
            fo.write("\n")


def write_bics_txt(expected, file_path):
    with open(file_path, "w") as fo:
        fo.write("Number of planted biclusters: {0}\n\n".format(len(expected)))
        for index, bic in enumerate(expected):
            fo.write("Bicluster #{0} ({1}, {2})\n".format(
                index, len(bic.rows), len(bic.cols)))
            fo.write("Rows=" + " ".join(str(r) for r in bic.rows) + "\n")
            fo.write("Columns=" + " ".join(str(c) for c in bic.cols) + "\n\n")


def write_bics_json(expected, file_path, meta, noisy_cells_by_bic):
    payload = dict(meta)
    payload["biclusters"] = {}

    for index, bic in enumerate(expected):
        payload["biclusters"][str(index)] = {
            "Rows": list(map(int, bic.rows)),
            "Columns": list(map(int, bic.cols)),
            "NoisyCells": [[int(r), int(c)] for r, c in noisy_cells_by_bic[index]],
        }

    with open(file_path, "w") as fo:
        json.dump(payload, fo, sort_keys=True)


def save_dataset(out_dir, name, kind, rows, cols, bic_rows, bic_cols, bic_num,
                 overlap_rows=0, overlap_cols=0,
                 bic_noise=DEFAULT_BIC_NOISE, shuffle=False):
    ensure_dir(out_dir)

    data, expected = make_data(
        kind=kind,
        rows=rows,
        cols=cols,
        bic_rows=bic_rows,
        bic_cols=bic_cols,
        bic_num=bic_num,
        overlap_rows=overlap_rows,
        overlap_cols=overlap_cols,
        shuffle=shuffle,
    )
    data, expected, noisy_cells = add_bicluster_noise(
        data, expected, noise_ratio=bic_noise, noise_sd=NOISE_SD)

    meta = {
        "kind": kind,
        "rows": rows,
        "cols": cols,
        "biclusterRows": bic_rows,
        "biclusterColumns": bic_cols,
        "biclusterNumber": bic_num,
        "overlapRows": overlap_rows,
        "overlapColumns": overlap_cols,
        "biclusterNoiseRatio": bic_noise,
        "biclusterNoiseDistribution": "N(0, 0.1)",
    }

    write_data(data, os.path.join(out_dir, name + "_data.tsv"))
    write_bics_txt(expected, os.path.join(out_dir, name + "_bics.txt"))
    write_bics_json(expected, os.path.join(out_dir, name + "_bics.json"),
                    meta, noisy_cells)


def repeat(out_dir, kind, rows, cols, bic_rows, bic_cols, bic_num,
           overlap_rows=0, overlap_cols=0, bic_noise=DEFAULT_BIC_NOISE):
    for replicate in range(1, REPLICATES + 1):
        name = "dataset_{0:02d}".format(replicate)
        save_dataset(
            out_dir=out_dir,
            name=name,
            kind=kind,
            rows=rows,
            cols=cols,
            bic_rows=bic_rows,
            bic_cols=bic_cols,
            bic_num=bic_num,
            overlap_rows=overlap_rows,
            overlap_cols=overlap_cols,
            bic_noise=bic_noise,
        )
        print("done", out_dir, name)


def main():
    random.seed()
    numpy.random.seed()
    ensure_dir(ROOT)

    # 1-3. Constant, additive, multiplicative. Default: 10000x100, 3 non-overlap bics.
    for kind in ["constant", "additive", "multiplicative"]:
        repeat(
            out_dir=os.path.join(ROOT, kind),
            kind=kind,
            rows=10000,
            cols=100,
            bic_rows=1000,
            bic_cols=20,
            bic_num=3,
        )

    # 4. Overlap 10%, 20%, 30%, 40% on rows and columns.
    # These are shift-scale column-trend order-preserving datasets.
    for ratio in [0.1, 0.2, 0.3, 0.4]:
        repeat(
            out_dir=os.path.join(ROOT, "overlap", "none", str(ratio)),
            kind="shift_scale",
            rows=10000,
            cols=100,
            bic_rows=1000,
            bic_cols=20,
            bic_num=3,
            overlap_rows=int(round(1000 * ratio)),
            overlap_cols=int(round(20 * ratio)),
        )

    # 5. Different matrix sizes.
    # These are shift-scale column-trend order-preserving datasets.
    size_settings = []
    for cols in [100, 500, 1000, 1500, 2000]:
        size_settings.append((1000, cols))
    for rows in [5000, 10000, 15000, 20000]:
        size_settings.append((rows, 100))

    for rows, cols in size_settings:
        bic_rows = max(100, int(round(rows * 0.10)))
        bic_cols = max(int(round(cols * 0.10)), 20)
        repeat(
            out_dir=os.path.join(ROOT, "bicluster",
                                 "rows{0}_cols{1}".format(rows, cols)),
            kind="shift_scale",
            rows=rows,
            cols=cols,
            bic_rows=bic_rows,
            bic_cols=bic_cols,
            bic_num=3,
        )

    # 6. Different bicluster counts.
    # These are shift-scale column-trend order-preserving datasets.
    for bic_num in [1, 3, 5]:
        repeat(
            out_dir=os.path.join(ROOT, "bicluster_num", str(bic_num)),
            kind="shift_scale",
            rows=10000,
            cols=100,
            bic_rows=1000,
            bic_cols=20,
            bic_num=bic_num,
        )

    # 7. Different bicluster noise ratios.
    # These are shift-scale column-trend order-preserving datasets.
    for noise in [0.1, 0.2, 0.3, 0.4,0.5]:
        repeat(
            out_dir=os.path.join(ROOT, "Noise", str(noise)),
            kind="shift_scale",
            rows=10000,
            cols=100,
            bic_rows=1000,
            bic_cols=20,
            bic_num=3,
            bic_noise=noise,
        )


if __name__ == "__main__":
    main()
