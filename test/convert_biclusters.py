#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Convert RecBic-like bicluster dump to 2-line-per-bicluster format.

Input blocks look like:
  PC_Genes [..]: gene40 gene42 ...
  NC_Genes [..]: gene41 gene47 ...
  Conds     [..]: cond40 cond41 ...

Output (for each bicluster):
  <row indices merged from PC_ + NC_>
  <column indices from Conds>
  <blank line>
"""

import re
import sys

gene_pat = re.compile(r'gene(\d+)')
cond_pat = re.compile(r'cond(\d+)')

def parse_nums(line, pat):
    """Extract all integers matched by pattern, dedup & sort."""
    nums = pat.findall(line)
    # py2: map returns list; ensure int conversion works
    nums = [int(x) for x in nums]
    nums = sorted(set(nums))
    return nums

def convert(in_fp, out_fp):
    cur_rows_pc = None
    cur_rows_nc = None
    cur_cols = None

    for raw in in_fp:
        line = raw.strip()
        if not line:
            continue

        if line.startswith('PC_Genes'):
            cur_rows_pc = parse_nums(line, gene_pat)
        elif line.startswith('NC_Genes'):
            cur_rows_nc = parse_nums(line, gene_pat)
        elif line.startswith('Conds'):
            cur_cols = parse_nums(line, cond_pat)

        # when we have Conds, we assume one bicluster block is complete
        if cur_cols is not None:
            rows = []
            if cur_rows_pc: rows.extend(cur_rows_pc)
            if cur_rows_nc: rows.extend(cur_rows_nc)
            # merge + dedup + sort
            rows = sorted(set(rows))

            # write two lines + blank line (Py2/3 safe)
            out_fp.write(' '.join([str(x) for x in rows]) + '\n')
            out_fp.write(' '.join([str(x) for x in cur_cols]) + '\n')
            out_fp.write('\n')

            # reset for next block
            cur_rows_pc = None
            cur_rows_nc = None
            cur_cols = None

def main():
    # Simple argv parsing to stay Py2-friendly:
    #   python convert_biclusters_legacy.py input.txt output.txt
    # or use "-" for stdin/stdout
    argv = sys.argv[1:]
    if len(argv) == 0 or len(argv) > 2:
        sys.stderr.write("Usage: python %s <input|-] [output|-]\n" % sys.argv[0])
        sys.exit(1)

    in_path = argv[0]
    out_path = argv[1] if len(argv) == 2 else '-'

    in_fp = sys.stdin if in_path == '-' else open(in_path, 'r')
    out_fp = sys.stdout if out_path == '-' else open(out_path, 'w')

    try:
        convert(in_fp, out_fp)
    finally:
        if in_fp is not sys.stdin:
            in_fp.close()
        if out_fp is not sys.stdout:
            out_fp.close()

if __name__ == '__main__':
    main()
