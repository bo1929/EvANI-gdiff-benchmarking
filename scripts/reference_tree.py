#!/usr/bin/env python3
"""Reference (ground-truth) trees for the EvANI benchmark cases.

Layout:
  data/trees/simulated_tree_mutation_{5,10,25,50,100,200,400,500}.nwk
       -- one per mutation rate, all sharing the SAME topology (branch
          lengths differ).
  data/trees/simulated_tree_duplication_lgt.nwk
       -- the same topology again (used for both duplication and lgt studies).
  data/simulated_dataset/mutation/RealTree300.drw (ALF "RealTree" RSF format)
       -- the source tree for mutation rate 300 (no .nwk shipped).

We return, for every (study, rate, replicate) case, the reference newick
(string), i.e. the topology + branch lengths, with leaf labels SE001..SE015.
Mutation-300 is synthesized from the .drw file via a tiny RSF parser.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent

TREE_DIR = _ROOT / "data" / "trees"
REALTREE300 = _ROOT / "data" / "simulated_dataset" / "mutation" / "RealTree300.drw"

_MUT_RATES = ["5", "10", "25", "50", "100", "200", "300", "400", "500"]


# --------------------------------------------------------------------------- #
# RSF (Simple format used by ALF / "RealTree"):  Tree(..,Leaf('SE001',len,[])...)
# --------------------------------------------------------------------------- #
def _rsf_to_newick(text):
    s = text[text.index("Tree("):]
    i = 0
    n = len(s)

    def skip_ws(i):
        while i < n and s[i] in " \t\r\n":
            i += 1
        return i

    def parse(i):
        i = skip_ws(i)
        tok = "Tree" if s.startswith("Tree", i) else "Leaf"
        i = skip_ws(i + len(tok))
        assert s[i] == "("
        i += 1
        if tok == "Leaf":
            i = skip_ws(i)
            m = re.match(r"'([^']*)'", s[i:])
            name = m.group(1)
            i = skip_ws(i + m.end())
            assert s[i] == ","
            i += 1
            i = skip_ws(i)
            m = re.match(r"[0-9.]+", s[i:])
            length = float(m.group(0))
            i = skip_ws(i + m.end())
            assert s[i] == ","
            i += 1
            i = skip_ws(i)
            assert s[i:i + 2] == "[]"
            i += 2
            assert s[i] == ")"
            return {"type": "leaf", "name": name, "length": length, "children": ()}, i + 1
        # Tree( A, len, B, [] )
        c1, i = parse(i)
        i = skip_ws(i)
        assert s[i] == ","
        i += 1
        i = skip_ws(i)
        m = re.match(r"[0-9.]+", s[i:])
        length = float(m.group(0))
        i = skip_ws(i + m.end())
        assert s[i] == ","
        i += 1
        c2, i = parse(i)
        i = skip_ws(i)
        assert s[i] == ","
        i += 1
        i = skip_ws(i)
        assert s[i:i + 2] == "[]"
        i += 2
        assert s[i] == ")"
        return {"type": "node", "length": length, "children": (c1, c2)}, i + 1

    root, _ = parse(0)
    return _render_newick(root)


def _render_newick(node):
    if node["type"] == "leaf":
        return f"{node['name']}:{node['length']}"
    inner = ",".join(_render_newick(c) for c in node["children"])
    return f"({inner}):{node['length']}"


# --------------------------------------------------------------------------- #
# reference tree per case
# --------------------------------------------------------------------------- #
def reference_tree(study, rate):
    """Return ground-truth newick string for a case (study, rate)."""
    tree_dir = TREE_DIR
    if study == "mutation":
        if rate == "300":
            drw = REALTREE300
            return _rsf_to_newick(drw.read_text()) + ";"
        f = tree_dir / f"simulated_tree_mutation_{rate}.nwk"
        return f.read_text().strip()
    # duplication / lgt share one topology (same branch lengths were used in
    # the ALF simulations; the EvANI repo ships one tree for both)
    f = tree_dir / "simulated_tree_duplication_lgt.nwk"
    return f.read_text().strip()


def reference_leaf_set():
    return [f"SE{i:03d}" for i in range(1, 16)]


if __name__ == "__main__":
    for r in _MUT_RATES:
        t = reference_tree("mutation", r)
        import re as _re
        leaves = _re.findall(r"\bSE\d{3}\b", t)
        print(f"mutation {r}: {len(leaves)} leaves, len={len(t)}")
    print("dup:", len(re.findall(r"SE\d{3}", reference_tree("duplication", "0.0005"))))
    print("lgt:", len(re.findall(r"SE\d{3}", reference_tree("lgt", "0.0010"))))