#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

GDIFF=./gdiff
DATA=data/simulated_dataset
OUT=methods/gdiff
RAW=$OUT/raw-dist

STUDIES="mutation duplication lgt"
SKETCH_ARGS="--frac 0.5 -k 23 -w 23 -l 333 --sample-size 1000"
DIST_ARGS="--hdist-th 3"
JOBS=12
THREADS=1

mkdir -p "$OUT" "$RAW"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

process() {
    local case=$1 db=$2 il="$TMP/$1.ilist" f b
    [ -s "$OUT/$case.tsv" ] && return 0

    : > "$il"
    for f in "$db"/SE*_dna.fa; do
        [ -f "$f" ] || continue
        b=$(basename "$f"); b=${b%_dna.fa}
        printf '%s.fa\t%s\n' "$b" "$f" >> "$il"
    done

    "$GDIFF" --num-threads "$THREADS" sketch --input-list "$il" $SKETCH_ARGS -o "$TMP/$case.gdiff"
    "$GDIFF" --num-threads "$THREADS" dist "$TMP/$case.gdiff" $DIST_ARGS -o "$RAW/$case.dist.tsv"
    python3 scripts/gdiff_dist_to_grid.py "$RAW/$case.dist.tsv" "$OUT/$case.tsv"
    rm -f "$TMP/$case.gdiff" "$il"
}

pids=(); names=(); rc=0
for study in $STUDIES; do
    for rate in "$DATA/$study"/*/; do
        for rep in "$rate"*/; do
            db="${rep}DB"
            [ -d "$db" ] || continue
            cid="$(basename "$study")_$(basename "$rate")_$(basename "$rep")"
            process "$cid" "$db" &
            pids+=($!); names+=("$cid")
            if [ ${#pids[@]} -ge "$JOBS" ]; then
                for i in "${!pids[@]}"; do
                    wait "${pids[$i]}" || { echo "failed: ${names[$i]}" >&2; rc=1; }
                done
                pids=(); names=()
            fi
        done
    done
done
for i in "${!pids[@]}"; do wait "${pids[$i]}" || { echo "failed: ${names[$i]}" >&2; rc=1; }; done

echo "grids: $(ls "$OUT"/*.tsv | wc -l | tr -d ' ') -> $OUT"
exit $rc
