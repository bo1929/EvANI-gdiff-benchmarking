# reproduce-tree-inference

Reproduction of the EvANI species-tree inference results and figures using the
new `gdiff` binary. All commands are run from this directory:

```bash
./00_methods.sh                      # mash/skani/fastANI/dashing2 grids -> methods/
./01_gdiff.sh                        # gdiff grids        -> methods/gdiff/
./02_tree_eval.sh                    # NJ*/BIONJ* + RF/FN -> results/tree-eval.tsv
./03_summarise.sh                    # summaries + PDFs   -> results/
Rscript 04_evani-tree-inference.R    # bar chart          -> results/figs/
```

Each script keeps its configuration (paths, binaries, parameters, jobs) as
plain variables at the top. Steps 2-4 read what steps 0-1 wrote, so run them in
order. Step 0 is optional: the `methods/` grids are already present, and
`00_methods.sh` skips a grid that exists unless `FORCE=1`.

## Notes

- `fastani-complete` is the fastANI `--fragLen 1000 --minFraction 0.1` run. It
  is named "complete" in the parent repo because, unlike the documented
  `--fragLen 3000` grid (undefined ANI for 20/85 cases), it produces a complete
  matrix. The R figure uses this variant.
- The shipped `methods/` grids for mash, skani, fastANI and dashing2 were
  regenerated with `00_methods.sh` and match the copies in `../EvANI/` exactly.
- `data/simulated_dataset` is a symlink to the original dataset to avoid
  duplicating it; nothing in the parent directory is modified.

## Requirements

- new `gdiff` binary at `../gdiff` (`GDIFF` in `01_gdiff.sh`)
- `mash`, `skani`, `fastANI` and `dashing2` for `00_methods.sh` (paths at the
  top of that script); not needed if you keep the shipped grids
- `Rscript` with `ape`, plus `dplyr`, `ggplot2`, `vroom`, `latex2exp`,
  `cowplot`, `tidyr` for step 4
- python3 with `pandas`, `seaborn`, `matplotlib`, `numpy` for step 3
  (`PY` in `03_summarise.sh`); steps 0-2 need only the standard library
