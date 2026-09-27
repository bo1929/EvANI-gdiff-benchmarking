#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1

DATA=data/simulated_dataset
METHODS=methods
FASTANI=methods/fastani-complete
OUT=results

python3 scripts/evaluate_trees.py \
    --data "$DATA" \
    --evani "$METHODS" \
    --method-dir "fastani=$FASTANI" \
    --out "$OUT"
