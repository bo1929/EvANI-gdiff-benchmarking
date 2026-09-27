#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

DATA=data/simulated_dataset
OUT=methods
STUDIES="mutation duplication lgt"
METHODS="dashing2" # skani fastani mash 
JOBS=8
THREADS=4
FORCE=1

MASH=mash
SKANI=skani
FASTANI=fastANI
DASHING2=./dashing2-osx

MASH_K=21
MASH_S=10000
SKANI_ARGS="--slow --min-af 0"
FASTANI_FRAGLEN=1000
FASTANI_MINFRACTION=0.1
DASHING2_MODE=distance
DASHING2_K=21
DASHING2_ARGS="--mash-distance --full -L 14 -k $DASHING2_K"
# DASHING2_MODE=containment
# DASHING2_K=23
# DASHING2_ARGS="--symmetric-containment -S 2048 -k $DASHING2_K"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

emit_mash() {
    local db=$1 tag=$2
    local d="$TMP/mash-$tag" f b names=() i j x
    mkdir -p "$d"
    for f in "$db"/SE*_dna.fa; do
        [ -f "$f" ] || continue
        b=$(basename "$f"); b=${b%_dna.fa}
        "$MASH" sketch -k "$MASH_K" -s "$MASH_S" -o "$d/$b" "$f" >/dev/null 2>&1
        names+=("$b")
    done
    for ((i=0; i<${#names[@]}; i++)); do
        for ((j=i+1; j<${#names[@]}; j++)); do
            x=$("$MASH" dist "$d/${names[$i]}.msh" "$d/${names[$j]}.msh" 2>/dev/null | tail -1 | cut -f3)
            [ -n "$x" ] || continue
            awk -v a="${names[$i]}.fa" -v b="${names[$j]}.fa" -v x="$x" \
                'BEGIN{printf "%s\t%s\t%.6f\n", a, b, (1-x)*100}'
        done
    done
    rm -rf "$d"
}

emit_skani() {
    local db=$1 tag=$2
    local ql="$TMP/skani-$tag.ql" f
    for f in "$db"/SE*_dna.fa; do [ -f "$f" ] && echo "$f"; done > "$ql"
    "$SKANI" dist -t "$THREADS" $SKANI_ARGS --ql "$ql" --rl "$ql" 2>/dev/null | awk -F'\t' '
        function base(p,   a,n,s){ n=split(p,a,"/"); s=a[n]; sub(/_dna\.fa$/,"",s); return s }
        NR>1 { r=base($1); q=base($2); if(r==q) next; if(r<q) next; print r".fa\t"q".fa\t"$3 }'
}

emit_fastani() {
    local db=$1 tag=$2
    local ql="$TMP/fani-$tag.ql" rl="$TMP/fani-$tag.rl" out="$TMP/fani-$tag.out" f dbabs
    case $db in /*) dbabs=$db ;; *) dbabs="$PWD/$db" ;; esac
    : > "$ql"; : > "$rl"
    for f in "$dbabs"/SE*_dna.fa; do
        [ -f "$f" ] || continue
        echo "$f" >> "$ql"; echo "$f" >> "$rl"
    done
    "$FASTANI" --ql "$ql" --rl "$rl" -o "$out" --fragLen "$FASTANI_FRAGLEN" \
        --minFraction "$FASTANI_MINFRACTION" -t "$THREADS" >/dev/null 2>&1 || true
    awk -F'\t' '
        function base(p,   a,n,s){ n=split(p,a,"/"); s=a[n]; sub(/_dna\.fa$/,"",s); return s }
        { r=base($1); q=base($2); if(r==q) next; if(r<q) next; print r".fa\t"q".fa\t"$3 }' "$out" 2>/dev/null
}

emit_dashing2() {
    local db=$1 tag=$2
    local out="$TMP/dashing2-$tag.txt" f files=()
    for f in "$db"/SE*_dna.fa; do [ -f "$f" ] && files+=("$f"); done
    [ ${#files[@]} -ge 2 ] || return 0
    "$DASHING2" cmp $DASHING2_ARGS -p "$THREADS" --cmpout "$out" "${files[@]}" >/dev/null 2>&1
    python3 scripts/dashing2_pairs.py "$DASHING2_MODE" "$DASHING2_K" < "$out"
}

process() {
    local method=$1 cid=$2 db=$3 outdir
    case $method in
        fastani) outdir=$OUT/fastani-complete ;;
        *)       outdir=$OUT/$method ;;
    esac
    if [ "$FORCE" != 1 ] && [ -s "$outdir/$cid.tsv" ]; then return 0; fi
    mkdir -p "$outdir"
    local raw="$TMP/$method-$cid.pairs"
    emit_"$method" "$db" "$cid-$$-$RANDOM" > "$raw"
    python3 scripts/pairs_to_grid.py "$raw" "$outdir/$cid.tsv"
}

pids=(); names=(); rc=0
for method in $METHODS; do
    for study in $STUDIES; do
        for rate in "$DATA/$study"/*/; do
            for rep in "$rate"*/; do
                db="${rep}DB"
                [ -d "$db" ] || continue
                cid="$(basename "$study")_$(basename "$rate")_$(basename "$rep")"
                process "$method" "$cid" "$db" &
                pids+=($!); names+=("$method/$cid")
                if [ ${#pids[@]} -ge "$JOBS" ]; then
                    for i in "${!pids[@]}"; do
                        wait "${pids[$i]}" || { echo "failed: ${names[$i]}" >&2; rc=1; }
                    done
                    pids=(); names=()
                fi
            done
        done
    done
done
for i in "${!pids[@]}"; do wait "${pids[$i]}" || { echo "failed: ${names[$i]}" >&2; rc=1; }; done

for method in $METHODS; do
    d=$OUT/$method; [ "$method" = fastani ] && d=$OUT/fastani-complete
    echo "$method: $(ls "$d"/*.tsv 2>/dev/null | wc -l | tr -d ' ') grids -> $d"
done
exit $rc
