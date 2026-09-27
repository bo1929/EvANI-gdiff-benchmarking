#!/usr/bin/env python3
import csv
import sys

NAMES = [f"SE{i:03d}.fa" for i in range(1, 16)]
COL = "d"

dist_tsv, grid_tsv = sys.argv[1], sys.argv[2]

pairs = {}
with open(dist_tsv) as fh:
    rows = [ln for ln in fh if not ln.startswith("#")]
for row in csv.DictReader(rows, delimiter="\t"):
    a, b = sorted((row["genome_a"], row["genome_b"]))
    pairs[(a, b)] = float(row[COL])

with open(grid_tsv, "w") as fh:
    for a in NAMES:
        for b in NAMES:
            d = 0.0 if a == b else pairs.get(tuple(sorted((a, b))))
            ani = 100.0 * (1.0 - d) if d is not None else float("nan")
            fh.write(f"{a}\t{b}\t{ani:.4f}\n")
