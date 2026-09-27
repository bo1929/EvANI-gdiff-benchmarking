#!/usr/bin/env python3
import math
import sys

# usage: dashing2_pairs.py <mode> <k>
#   distance    --cmpout holds a Mash/Poisson distance -> ANI = 100*(1-d)
#   containment --cmpout holds a containment similarity C -> D = -ln(C)/k
MODE = sys.argv[1] if len(sys.argv) > 1 else "distance"
K = float(sys.argv[2]) if len(sys.argv) > 2 else 21.0


def to_ani(v):
    if v != v:
        return None
    if MODE == "containment":
        if v <= 0.0:
            return None
        return 100.0 * (1.0 + math.log(min(v, 1.0)) / K)
    if v < 0.0 or v > 1.0:
        return None
    return 100.0 * (1.0 - v)


lines = [ln.rstrip("\n") for ln in sys.stdin]

sources = None
for ln in lines:
    if ln.startswith("#Sources"):
        sources = ln.split("\t")[1:]
        break
if not sources:
    sys.exit(0)

names = [s.rsplit("/", 1)[-1] for s in sources]


def label(s):
    return s[:-7] + ".fa" if s.endswith("_dna.fa") else s


out = {}
for ln in lines:
    if not ln or ln.startswith("#"):
        continue
    f = ln.split("\t")
    row = f[0].rsplit("/", 1)[-1]
    if row not in names:
        continue
    i = names.index(row)
    for j in range(i + 1, len(names)):
        try:
            v = float(f[j + 1])
        except (ValueError, IndexError):
            continue
        ani = to_ani(v)
        if ani is None:
            continue
        a, b = sorted([label(names[i]), label(names[j])])
        out[(a, b)] = ani

for (a, b), ani in sorted(out.items()):
    print(f"{a}\t{b}\t{ani:.6f}")
