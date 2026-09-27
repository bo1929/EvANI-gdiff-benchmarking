#!/usr/bin/env python3
"""NJ* and BIONJ* via R's ape::njs / ape::bionjs.

Distance convention: an (n x n) Python matrix where `M[i][j] is None` (or NaN)
means a missing pairwise distance.  ape codes missing as -1 internally; we pass
NA through `as.dist()` and let njs/bionjs convert it.
"""
from __future__ import annotations

import os
import re
import subprocess
import tempfile

__all__ = ["nj_star", "bionj_star", "infer_trees"]

FS = 15  # ape's njs/bionjs default and the value used throughout this repo

_RSCRIPT = "Rscript"

_cache: dict = {}


def _key(M):
    return tuple(
        tuple(None if (j <= i or M[i][j] is None or M[i][j] != M[i][j]) else float(M[i][j])
              for j in range(len(M)))
        for i in range(len(M))
    )


def _relabel(nwk: str, n: int) -> str:
    return re.sub(r"t(\d+)\b", lambda mo: str(int(mo.group(1)) - 1), nwk.strip())


def infer_trees(M, fs: int = FS):
    """Return (nj_star_newick, bionj_star_newick) for distance matrix M.

    Both trees are inferred in a single R process.  Leaves are relabelled to
    "0".."n-1" to match the rest of the pipeline.
    """
    n = len(M)
    if n < 3:
        raise ValueError("need >=3 taxa")
    key = (_key(M), int(fs))
    if key in _cache:
        return _cache[key]

    lines = [
        "suppressMessages(library(ape))",
        "options(warn=-1)",
        f"m <- matrix(NA_real_, {n}, {n})",
        "diag(m) <- 0",
    ]
    for i in range(n):
        for j in range(i + 1, n):
            v = M[i][j]
            if v is None or v != v:  # missing / NaN
                continue
            lines.append(f"m[{i + 1},{j + 1}] <- {float(v)!r}")
            lines.append(f"m[{j + 1},{i + 1}] <- {float(v)!r}")
    lines.append(f"rownames(m) <- paste0('t', 1:{n})")
    lines.append(f"colnames(m) <- paste0('t', 1:{n})")
    lines.append("d <- as.dist(m)")
    lines.append("for (meth in c('njs', 'bionjs')) {")
    lines.append(f"  tr <- tryCatch(get(meth)(d, fs={int(fs)}), error=function(e) NULL)")
    lines.append("  if (is.null(tr)) {")
    lines.append("    cat(meth, ':__ERR__', sep=''); cat('\\n')")
    lines.append("  } else {")
    lines.append("    cat(meth, ':', write.tree(tr), sep=''); cat('\\n')")
    lines.append("  }")
    lines.append("}")

    fh = tempfile.NamedTemporaryFile("w", suffix=".R", delete=False)
    try:
        fh.write("\n".join(lines))
        fh.close()
        out = subprocess.run([_RSCRIPT, fh.name], capture_output=True,
                             text=True, timeout=600)
    finally:
        try:
            os.unlink(fh.name)
        except OSError:
            pass

    if out.returncode != 0:
        raise RuntimeError("ape njs/bionjs failed: " + out.stderr.strip()[-400:])

    trees = {}
    for ln in out.stdout.splitlines():
        ln = ln.strip()
        for meth, outkey in (("njs", "njs"), ("bionjs", "bionjs")):
            pref = meth + ":"
            if ln.startswith(pref):
                body = ln[len(pref):].strip()
                if body and body != "__ERR__":
                    trees[outkey] = _relabel(body, n)
    if "njs" not in trees or "bionjs" not in trees:
        raise RuntimeError("ape njs/bionjs produced no tree: " + out.stdout[:200])

    result = (trees["njs"], trees["bionjs"])
    _cache[key] = result
    return result


def nj_star(D, s: int = FS):
    return infer_trees(D, s)[0]


def bionj_star(D, s: int = FS):
    return infer_trees(D, s)[1]
