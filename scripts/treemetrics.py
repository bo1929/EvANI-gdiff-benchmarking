#!/usr/bin/env python3
"""Robust newick parsing + tree split extraction + RF / FN / FP metrics.

Splits are canonicalized over the *unrooted* tree: each internal edge induces a
2-partition of the leaf label set; we store the smaller side as a frozenset
(the complement is redundant).  Trivial partitions (size 1 or n-1) are skipped.

Metrics between a reference and estimated tree over the same leaf set (both
binary):
  FN = reference splits missing from estimate
  FP = estimated splits missing from reference
  RF = FN + FP                     (Robinson-Foulds symmetric distance)
  nRF = RF / (2*(n-3))
  nFN = FN / (n-3)
"""
from __future__ import annotations

__all__ = ["tree_splits", "split_metrics"]


# --------------------------------------------------------------------------- #
# newick parsing
# --------------------------------------------------------------------------- #
class _Node:
    __slots__ = ("name", "length", "children", "parent", "_leafset")

    def __init__(self, name=None, length=None):
        self.name = name
        self.length = length
        self.children = []
        self.parent = None


def parse_newick(s):
    """Parse a newick string into a rooted _Node tree.  Returns root."""
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]
    s = s.strip()
    pos = [0]

    def read_name():
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in "(),:":
            pos[0] += 1
        return s[start:pos[0]].strip()

    def read_length():
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in "(),;":
                pos[0] += 1
            try:
                return float(s[start:pos[0]])
            except ValueError:
                return None
        return None

    def parse_node():
        node = _Node()
        if pos[0] < len(s) and s[pos[0]] == "(":
            pos[0] += 1  # consume (
            while True:
                child = parse_node()
                node.children.append(child)
                child.parent = node
                if pos[0] < len(s) and s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                break
            # expect )
            assert pos[0] < len(s) and s[pos[0]] == ")"
            pos[0] += 1
        else:
            node.name = read_name()
        node.length = read_length()
        return node

    root = parse_node()
    if pos[0] < len(s):
        print("warn: trailing", repr(s[pos[0]:]))
    return root


def iter_leaves(root):
    if root.name is not None:
        yield root.name
        return
    for c in root.children:
        yield from iter_leaves(c)


def tree_splits(newick):
    """Return set of nontrivial unrooted splits (smaller side): frozenset of
    labels."""
    root = parse_newick(newick)
    leaves = sorted(iter_leaves(root))
    n = len(leaves)

    # leaf-set under each node (as frozensets)
    def fill(node):
        if node.name is not None:
            node._leafset = frozenset({node.name})
        else:
            sets = [fill(c) for c in node.children]
            node._leafset = frozenset().union(*sets) if sets else frozenset()
        return node._leafset

    fill(root)
    splits = set()
    all_leaves = frozenset(leaves)
    # Canonicalize each split to the side that does NOT contain a fixed anchor
    # leaf.  This is unambiguous even for balanced splits (|side| == n/2),
    # where "keep the smaller side" would map the same split to two different
    # frozensets depending on the tree's orientation.
    anchor = leaves[0]

    def collect(node, parent_side=None):
        # every edge from node to its parent induces a split = node's leafset
        if node is root:
            for c in node.children:
                collect(c)
            return
        s = node._leafset
        size = len(s)
        if size != 1 and size != n - 1:
            canon = (all_leaves - s) if anchor in s else s
            splits.add(frozenset(canon))
        for c in node.children:
            collect(c)
        return

    collect(root)
    return splits


# --------------------------------------------------------------------------- #
# metrics
# --------------------------------------------------------------------------- #
def split_metrics(ref_newick, est_newick, n_leaves=None):
    """Compare reference vs estimated tree.  n_leaves: optional leaf count;
    if omitted, inferred from the label set of ref|est splits."""
    ref = tree_splits(ref_newick)
    est = tree_splits(est_newick)
    labels = set()
    for sp in ref:
        labels |= set(sp)
    for sp in est:
        labels |= set(sp)
    n = n_leaves if n_leaves is not None else len(labels)

    fn = len(ref - est)
    fp = len(est - ref)
    return {
        "n_splits_ref": len(ref),
        "n_splits_est": len(est),
        "FN": fn,
        "FP": fp,
        "RF": fn + fp,
        "nRF": (fn + fp) / (2 * (n - 3)) if n > 3 else 0.0,
        "nFN": fn / (n - 3) if n > 3 else 0.0,
    }


if __name__ == "__main__":
    import sys
    t1 = "((a,b),(c,d));"
    t2 = "((a,c),(b,d));"
    print("splits t1:", tree_splits(t1))
    print("splits t2:", tree_splits(t2))
    print("metrics:", split_metrics(t1, t2))