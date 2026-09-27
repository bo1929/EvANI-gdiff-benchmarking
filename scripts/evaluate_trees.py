#!/usr/bin/env python3
"""Build a 15x15 distance matrix per method/case, infer NJ*/BIONJ* (R
ape::njs / bionjs), and compare to the ground-truth tree (RF / FN / FP / nRF).

Reads ordered grids <EVANI>/<method>/<case>.tsv (`genome_a genome_b ani`,
nan = missing) and writes <OUT>/tree-eval.tsv + tree-eval-wide.tsv.
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from njstar_ape import nj_star, bionj_star
from treemetrics import split_metrics
from reference_tree import reference_tree, reference_leaf_set

METHODS = ["gdiff", "mash", "skani", "fastani", "dashing2"]


def iter_cases(data_root: Path):
    """Yield (study, rate, rep) for all benchmark cases."""
    for study in ["mutation", "duplication", "lgt"]:
        sdir = data_root / study
        if not sdir.is_dir():
            continue
        for rate_d in sorted(sdir.iterdir()):
            if not rate_d.is_dir():
                continue
            rate = rate_d.name
            for rep_d in sorted(rate_d.iterdir()):
                if not rep_d.is_dir():
                    continue
                rep = rep_d.name
                db = rep_d / "DB"
                n = len(list(db.glob("SE*_dna.fa"))) if db.is_dir() else 0
                if n >= 2:
                    yield study, rate, rep


def load_evani_grid(path: Path):
    """Load an EvANI/sample_tool-style ordered grid (a<TAB>b<TAB>ani).

    Returns {(a,b): d} with a<b and d = 1 - ANI/100; nan/inf -> missing.
    """
    dist = {}
    if not path.exists():
        return dist
    for ln in path.read_text().splitlines():
        if not ln.strip():
            continue
        fields = ln.split("\t")
        if len(fields) < 3:
            continue
        a, b, ani = fields[0], fields[1], fields[2]
        try:
            v = float(ani)
        except ValueError:
            continue
        if not math.isfinite(v):  # NaN/inf -> missing
            continue
        if a == b:
            continue
        a, b = sorted([a, b])
        dist[(a, b)] = 1.0 - v / 100.0
    return dist


def build_matrix(pairs, leaf_labels):
    """Return (matrix, n_missing) for the 15 leaves (SE001..SE015)."""
    n = len(leaf_labels)
    M = [[None] * n for _ in range(n)]
    nmissing = 0
    for i in range(n):
        M[i][i] = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            name = (f"{leaf_labels[i]}.fa", f"{leaf_labels[j]}.fa")
            a, b = sorted(name)
            d = pairs.get((a, b))
            if d is None:
                nmissing += 1
            else:
                M[i][j] = M[j][i] = d
    return M, nmissing


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=Path("data/simulated_dataset"))
    ap.add_argument("--evani", type=Path, default=Path("methods"),
                    help="grid root <EVANI>/<method>/<case>.tsv")
    ap.add_argument("--method-dir", action="append", default=[],
                    metavar="METHOD=DIR",
                    help="override the grid directory for one method, e.g. "
                         "fastani=methods/fastani-complete (repeatable)")
    ap.add_argument("--out", type=Path, default=Path("results"))
    args = ap.parse_args()

    grid_dirs = {m: args.evani / m for m in METHODS}
    for spec in args.method_dir:
        name, sep, d = spec.partition("=")
        if not sep or name not in METHODS:
            ap.error(f"--method-dir expects METHOD=DIR with METHOD in {METHODS}")
        grid_dirs[name] = Path(d)

    leaf_labels = reference_leaf_set()
    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for study, rate, rep in iter_cases(args.data):
        case = f"{study}_{rate}_{rep}"
        ref_nwk = reference_tree(study, rate)

        for method in METHODS:
            pairs = load_evani_grid(grid_dirs[method] / f"{case}.tsv")
            M, nmissing = build_matrix(pairs, leaf_labels)
            mp = {leaf_labels[i]: str(i) for i in range(len(leaf_labels))}

            import re
            def relabel(nwk):
                out = re.sub(r"\b(SE\d{3})\b", lambda mo: mp[mo.group(1)], nwk)
                return out
            ref_relabelled = relabel(ref_nwk)

            for builder, fn in [("NJ*", nj_star), ("BIONJ*", bionj_star)]:
                try:
                    est_nwk = fn(M)
                except Exception as e:  # noqa
                    msg = " ".join(str(e).split())  # keep the TSV single-line
                    rows.append(dict(study=study, rate=rate, rep=rep, case=case,
                                     method=method, builder=builder,
                                     n_missing=nmissing, RF="", FN="", FP="",
                                     nRF="", nFN="", error="no_tree: " + msg))
                    continue
                met = split_metrics(ref_relabelled, est_nwk, n_leaves=len(leaf_labels))
                rows.append(dict(study=study, rate=rate, rep=rep, case=case,
                                 method=method, builder=builder,
                                 n_missing=nmissing, **{k: met[k] for k in
                                                        ("RF", "FN", "FP", "nRF", "nFN")}))

    # write long table
    fields = ["study", "rate", "rep", "case", "method", "builder", "n_missing",
              "RF", "FN", "FP", "nRF", "nFN", "error"]
    with open(out_dir / "tree-eval.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # wide summary by (study, rate, method, builder)
    import collections
    agg = collections.defaultdict(list)
    for r in rows:
        key = (r["study"], r["rate"], r["method"], r["builder"])
        agg[key].append(r)

    with open(out_dir / "tree-eval-wide.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["study", "rate", "method", "builder", "n_cases", "n_failed",
                    "n_missing_mean", "RF_mean", "RF_median", "RF_sd",
                    "FN_mean", "FN_median", "nRF_mean", "nFN_mean", "nRF_median"])
        for (study, rate, method, builder), rs in sorted(agg.items()):
            def num(x):
                vals = [r[x] for r in rs if isinstance(r.get(x), (int, float))]
                return vals
            rfs = num("RF"); fns = num("FN"); nrfs = num("nRF"); nfns = num("nFN")
            nm = [r["n_missing"] for r in rs]
            n_failed = sum(1 for r in rs if r.get("error"))
            import statistics
            w.writerow([study, rate, method, builder, len(rs), n_failed,
                        f"{statistics.mean(nm):.1f}",
                        f"{statistics.mean(rfs):.2f}" if rfs else "",
                        f"{statistics.median(rfs):.1f}" if rfs else "",
                        f"{statistics.stdev(rfs):.2f}" if len(rfs) > 1 else "",
                        f"{statistics.mean(fns):.2f}" if fns else "",
                        f"{statistics.median(fns):.1f}" if fns else "",
                        f"{statistics.mean(nrfs):.4f}" if nrfs else "",
                        f"{statistics.mean(nfns):.4f}" if nfns else "",
                        f"{statistics.median(nrfs):.4f}" if nrfs else ""])
    print(f"wrote {out_dir/'tree-eval.tsv'} and {out_dir/'tree-eval-wide.tsv'}")
    print(f"rows: {len(rows)}")


if __name__ == "__main__":
    main()