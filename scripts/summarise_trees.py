#!/usr/bin/env python3
"""Summarise species-tree reconstruction accuracy and produce figures.

Reads results/tree-eval.tsv (long) + tree-eval-wide.tsv and writes results/summary-*.tsv
"""

from __future__ import annotations

import csv
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

HERE = Path(__file__).resolve().parent
RES = HERE.parent / "results"

BUILDER_ORDER = ["NJ*", "BIONJ*"]
STUDIES = {
    "mutation": ["5", "10", "25", "50", "100", "200", "300", "400", "500"],
    "duplication": ["0.0000", "0.0005", "0.0010", "0.0020"],
    "lgt": ["0.0001", "0.0005", "0.0010", "0.0020"],
}
STUDY_LABEL = {"mutation": "Mutation rate", "duplication": "Duplication rate", "lgt": "LGT rate"}
METHOD_ORDER = ["gdiff", "fastani", "mash", "skani", "dashing2"]
PAL = {"gdiff": "#E64B35", "fastani": "#4DBBD5", "mash": "#00A087", "skani": "#3C5488", "dashing2": "#F39B7F"}
METHOD_LABEL = {"gdiff": "gdiff", "fastani": "fastANI", "mash": "mash", "skani": "skani", "dashing2": "dashing2"}


def read_tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f


def mean(vs):
    vs = [x for x in vs if x is not None]
    return statistics.mean(vs) if vs else None


def _write(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    rows = read_tsv(RES / "tree-eval.tsv")
    ok = [r for r in rows if r.get("RF") not in ("", None)]

    # ---- overall ----
    overall = []
    for method in METHOD_ORDER:
        for builder in BUILDER_ORDER:
            sub = [r for r in ok if r["method"] == method and r["builder"] == builder]
            if not sub:
                continue
            overall.append(
                {
                    "method": method,
                    "builder": builder,
                    "n": len(sub),
                    "n_failed": sum(1 for r in rows if r["method"] == method and r["builder"] == builder and r.get("error")),
                    "FN_mean": mean([fnum(r["FN"]) for r in sub]),
                    "FN_median": statistics.median([fnum(r["FN"]) for r in sub]),
                    "RF_mean": mean([fnum(r["RF"]) for r in sub]),
                    "RF_median": statistics.median([fnum(r["RF"]) for r in sub]),
                    "nRF_mean": mean([fnum(r["nRF"]) for r in sub]),
                }
            )
    overall.sort(key=lambda r: (r["builder"], r["nRF_mean"]))
    _write(RES / "summary-overall.tsv", overall)
    print("=== overall (nRF_mean; 0 = perfect) ===")
    for r in overall:
        print(
            f"  {r['builder']:<6} {METHOD_LABEL[r['method']]:<9} "
            f"nRF={r['nRF_mean']:.4f}  FN={r['FN_mean']:.2f}  (n={r['n']}, failed={r['n_failed']})"
        )

    # ---- by study ----
    by_study = []
    for method in METHOD_ORDER:
        for builder in BUILDER_ORDER:
            for study in ["mutation", "duplication", "lgt"]:
                sub = [r for r in ok if r["method"] == method and r["builder"] == builder and r["study"] == study]
                if not sub:
                    continue
                by_study.append(
                    {
                        "study": study,
                        "method": method,
                        "builder": builder,
                        "n": len(sub),
                        "nRF_mean": mean([fnum(r["nRF"]) for r in sub]),
                        "FN_mean": mean([fnum(r["FN"]) for r in sub]),
                    }
                )
    _write(RES / "summary-by-study.tsv", by_study)
    print("\n=== by study (nRF_mean) ===")
    for r in by_study:
        print(f"  {r['study']:<10} {r['builder']:<5} {METHOD_LABEL[r['method']]:<9} " f"nRF={r['nRF_mean']:.4f}  FN={r['FN_mean']:.2f}")

    # ---- per-study, per-rate summary files ----
    global_rdata = {}
    for study, rates in STUDIES.items():
        by_rate = []
        for method in METHOD_ORDER:
            for builder in BUILDER_ORDER:
                for rate in rates:
                    sub = [r for r in ok if r["method"] == method and r["builder"] == builder and r["study"] == study and r["rate"] == rate]
                    nrfs = [fnum(r["nRF"]) for r in sub]
                    if not nrfs:
                        continue
                    by_rate.append(
                        {
                            "rate": rate,
                            "method": method,
                            "builder": builder,
                            "n": len(sub),
                            "nRF_mean": mean(nrfs),
                            "nRF_sd": statistics.stdev([x for x in nrfs if x is not None]) if len(nrfs) > 1 else None,
                            "FN_mean": mean([fnum(r["FN"]) for r in sub]),
                        }
                    )
        _write(RES / f"summary-{study}-rate.tsv", by_rate)


if __name__ == "__main__":
    main()
