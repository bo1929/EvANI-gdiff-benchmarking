#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

PY=/Users/asapci/micromamba/envs/evani/bin/python3

"$PY" scripts/summarise_trees.py
