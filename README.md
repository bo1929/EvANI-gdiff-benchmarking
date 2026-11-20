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

## References

`NJ*` (and `BIONJ*`) in step 2 are the *incomplete-distance* variants of
neighbor joining, not the standard NJ. They are needed because skani and
dashing2 leave pairwise distances undefined (NaN) for very divergent genomes,
whereas `ape::nj` / `ape::bionj` require a complete matrix. They are not
reimplemented here: `scripts/njstar_ape.py` calls the reference
implementations `ape::njs` / `ape::bionjs` with `fs = 15` (ape's default,
the paper's top-`s` pair parameter) and missing pairs passed as `NA`.

- **Criscuolo, A. & Gascuel, O. (2008)** Fast NJ-like algorithms to deal with
  incomplete distance matrices. *BMC Bioinformatics* **9**:166.
  [doi:10.1186/1471-2105-9-166](https://doi.org/10.1186/1471-2105-9-166)
  — introduces `NJ*`, `BIONJ*`, the `Q*` pair-selection criterion and the
  `s` parameter that `fs` controls.
- Saitou, N. & Nei, M. (1987) The neighbor-joining method: a new method for
  reconstructing phylogenetic trees. *Mol Biol Evol* **4**:406-425.
  [doi:10.1093/oxfordjournals.molbev.a040454](https://doi.org/10.1093/oxfordjournals.molbev.a040454)
  — the NJ algorithm that `NJ*` generalizes.
- Gascuel, O. (1997) BIONJ: an improved version of the NJ algorithm based on a
  simple model of sequence data. *Mol Biol Evol* **14**:685-695.
  [doi:10.1093/oxfordjournals.molbev.a025808](https://doi.org/10.1093/oxfordjournals.molbev.a025808)
  — the BIONJ algorithm that `BIONJ*` generalizes.
- Paradis, E. & Schliep, K. (2019) ape 5.0: an environment for modern
  phylogenetics and evolutionary analyses in R. *Bioinformatics* **35**:526-528.
  [doi:10.1093/bioinformatics/bty633](https://doi.org/10.1093/bioinformatics/bty633)
  — provides `njs()` / `bionjs()`.

