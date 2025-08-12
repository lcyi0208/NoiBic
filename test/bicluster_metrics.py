#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Compute list-wise bicluster similarity between two sets (G: expected, D: found).

Metrics:
- Jaccard-based Recovery = (1/|G|) * sum_{e in G} max_{f in D} J(e,f)
- Jaccard-based Relevance = (1/|D|) * sum_{f in D} max_{e in G} J(e,f)

Optional (asymmetric) recovery/relevance pair:
- recovery(e,f) = |e ∩ f| / |e|
- relevance(e,f) = |e ∩ f| / |f|

Input file format:
Each bicluster uses two lines:
  line 1: row indices (space-separated)
  line 2: col indices (space-separated)
Biclusters are separated by a blank line. Lines starting with '#' are ignored.
"""

import argparse
from collections import namedtuple
ListScore = namedtuple('ListScore', ['relevance', 'recovery'])


class Bicluster(object):
    def __init__(self, rows, cols):
        self.rows = frozenset(rows)
        self.cols = frozenset(cols)

    def area(self):
        return len(self.rows) * len(self.cols)

    def intersection(self, other):
        return Bicluster(self.rows & other.rows, self.cols & other.cols)

    def union(self, other):
        return Bicluster(self.rows | other.rows, self.cols | other.cols)

def parse_block(two_lines):
    r = list(map(int, two_lines[0].split()))
    c = list(map(int, two_lines[1].split()))
    return Bicluster(r, c)

def load_biclusters(path):
    blocks = []
    cur = []
    with open(path, "r") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                if not line and cur:
                    if len(cur) != 2:
                        raise ValueError("Malformed bicluster block: %s" % cur)
                    blocks.append(cur)
                    cur = []
                continue
            cur.append(line)
            if len(cur) == 2:
                blocks.append(cur)
                cur = []
        if cur:
            if len(cur) != 2:
                raise ValueError("Malformed bicluster block at EOF: %s" % cur)
            blocks.append(cur)
    return [parse_block(b) for b in blocks]

def area_union(b1, b2):
    return b1.area() + b2.area() - b1.intersection(b2).area()

def area_intersection(b1, b2):
    return b1.intersection(b2).area()

def jaccard(b1, b2):
    u = area_union(b1, b2)
    if u == 0:
        return 0.0
    return float(area_intersection(b1, b2)) / u

def recovery_pair(expected, found):
    a_e = expected.area()
    if a_e == 0:
        return 0.0
    return float(area_intersection(expected, found)) / a_e

def relevance_pair(expected, found):
    a_f = found.area()
    if a_f == 0:
        return 0.0
    return float(area_intersection(expected, found)) / a_f

def avg_of_maxes(A, B, score_fn):
    if not A or not B:
        raise ValueError("Empty bicluster list.")
    vals = []
    for a in A:
        best = max(score_fn(a, b) for b in B)
        vals.append(best)
    return sum(vals) / len(vals)

def jaccard_list_recovery_relevance(G, D):
    rec = avg_of_maxes(G, D, jaccard)
    rel = avg_of_maxes(D, G, jaccard)
    return rec, rel

def rr_list(G, D):
    rec = avg_of_maxes(G, D, recovery_pair)
    rel = avg_of_maxes(D, G, relevance_pair)
    return rec, rel

def build_pairwise_matrix(A, B, score_fn):
    return [[score_fn(a,b) for b in B] for a in A]

def jaccard_list(expected, found):
    rec, rel = jaccard_list_recovery_relevance(expected, found)
    return ListScore(relevance=rel, recovery=rec)
  
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--expected", "-g", required=True)
    ap.add_argument("--found", "-d", required=True)
    ap.add_argument("--metric", "-m", choices=["jaccard", "rr"], default="jaccard")
    ap.add_argument("--show-matrix", action="store_true")
    args = ap.parse_args()

    G = load_biclusters(args.expected)
    D = load_biclusters(args.found)

    if args.metric == "jaccard":
        rec, rel = jaccard_list_recovery_relevance(G, D)
        score_fn = jaccard
    else:
        rec, rel = rr_list(G, D)
        score_fn = recovery_pair

    print("# Metric:", args.metric,
      "Recovery (S(G,D))=%.8f" % rec,
      "Relevance (S(D,G))=%.8f" % rel)


    if args.show_matrix:
        print("\n# Pairwise score matrix (rows=expected, cols=found):")
        M = build_pairwise_matrix(G, D, score_fn)
        for row in M:
            print(" ".join("%.8f" % v for v in row))

if __name__ == "__main__":
    main()
