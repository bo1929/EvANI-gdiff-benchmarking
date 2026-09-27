#!/usr/bin/env python3
import math
import sys

NAMES = [f"SE{i:03d}.fa" for i in range(1, 16)]

pairs_tsv, grid_tsv = sys.argv[1], sys.argv[2]

pairs = {}
with open(pairs_tsv) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 3 or f[0] == f[1]:
            continue
        try:
            v = float(f[2])
        except ValueError:
            continue
        if not math.isfinite(v):
            continue
        pairs[tuple(sorted((f[0], f[1])))] = v

with open(grid_tsv, "w") as fh:
    for a in NAMES:
        for b in NAMES:
            v = 100.0 if a == b else pairs.get(tuple(sorted((a, b))))
            ani = float(v) if v is not None else float("nan")
            fh.write(f"{a}\t{b}\t{ani:.4f}\n")
