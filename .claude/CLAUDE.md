# CLAUDE.md

Snakemake port of the ENCODE4 cCRE pipeline (Moore lab, steps 0-8). Outputs MUST match the May 2022 V7
release byte for byte, except the deviations below. README.md holds the layout, idioms and style
conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run unit                                # differential unit tests (bedtools 2.30.0 oracle)
pixi run reference -m <step> -n 16           # compare data/results/v7 with the May 2022 release
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local <target>
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` is the rule env. Rule dependencies go there, never in the root.
- Its `oracle` env holds bedtools 2.30.0 alone. Its zlib pin conflicts with the UCSC tools.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.

## Fidelity

- Every difference from the reference MUST be reported with evidence before it is accepted.
- Keep deterministic bugs of the original. Record each under Quirks kept.
- Fix non-determinism and host dependence (set order, glob order, FMA). The fix becomes the new baseline.
- Logic lives in Python. External binaries stay only where output bytes depend on them.
- The reference tree belongs to another user and keeps changing. Check a file's mtime before trusting it.

### Deviations from the release

- D1. filter-tf-rpeaks.py joined a TF cluster's targets, biosamples and experiments in set() order. The
  port sorts them. Affects `TF-Clusters.txt`, `tmp.no-overlap` and `rTFBS-Summary.txt`.
- D2. Step 7 bgzips and tabix-indexes its BED, where the release wrote plain `.bed`. The bigBeds match.
- D3. Step 1 writes `Processed-DHSs/<signal>.<row>.bed`, not `output.<signal>.<row>`. Step 2 still
  concatenates them in the en_US order of the old names.
- D4. `REF`'s ATAC files come from a Sep 2023 run on glibc's FMA log, which reproduces them exactly. The
  port keeps the May 2022 non-FMA log. About 150 of 2.9M z-scores per file move by at most 9e-16. Signal
  and rank columns match. Steps 6-8 do not read ATAC.
- D5. `REF` hg38 `Individual-Files` mix two master lists. 1,675 files date from May 17-19, 4 from May 26,
  after a `DNase-List.txt` edit. The port follows the final list. `ENCFF346LEZ_ENCFF118OBT_ENCFF536VOI.bed`
  is stale. `ENCFF269VAY.bigBed` has no `.bed`. The port names the neural crest cell's file
  `ENCFF592CNZ_NA_ENCFF346DEM_ENCFF192FOK`, as the published scripts do. It equals `REF`'s
  `ENCFF592CNZ_ENCFF346DEM_ENCFF192FOK.bed`. GM12866's `NA.bed.gz` is empty. `REF` has no file for it.

### Release over the published scripts

- Step 8: the published calculate-center-distance.py keeps `|d| < 200`, an off-by-one against step 6's
  inclusive 200 bp window. It leaves 93 hg38 PLS without a gene. `PLS-Gene-List.hg38-V4.txt` (Dec 2023)
  keeps them. The port follows the release. The run-time copy of the script (2021) has no assignment mode.
- Step 3 (May 5) lacked ENCSR008LHT's peak file, which reached the mirror in June 2022. `tf_excluded`
  leaves it out of `tmp.summits`. Step 7 (May 19) used its peaks: SK-N-SH matches `REF` only with them.

### Quirks kept

1. Step 0 sweeps up to 1E-4942. Both genomes have FDR orders 4943 and 4944, which no threshold filters.
2. Step 0 prints an FDR without an exponent as `'%.2E'` before reading its order (0.009999 reads as
   order 2). A `+` exponent, as in 0, passes every threshold.
3. bedtools merge prints the min FDR with 10 significant digits. An FDR below the double range, such as
   1e-4944, prints as 0.
4. Step 1 sorts names with `sort -k4,4`, whose key includes the tab before field 4.
5. Step 2's `grep -v chrEBV` and `grep -v chrM` drop any line containing those strings.
6. make-region-accession.py keys earlier regions by chrom, start and end joined without separators.
   Distinct coordinates can collide.
7. Step 3 places a resized peak's summit at the new start plus the original offset.
8. Step 3's awk compares a non-numeric FRiP such as NA with 0.003 as text, so it passes.
9. Step 4's rank loop starts from `before=0`. A top z-score of exactly 0 would rank 0.
10. Two hg38 H3K4me3 rows list the file `NA`. Their z-score files are empty. Step 5 skips them. Step 7
    puts `NA` in their biosamples' file names.
11. In steps 6-7, `awk 'FNR==NR {x[$1];next} ...' first second` with an empty first file selects nothing,
    even under `!(key in x)`.
12. Step 7's No-DNase branch: an empty no1 or no2 falls back to the whole group or to no1. One cCRE can
    then get two rows.
13. Same branch: a mark whose above-1.64 lines miss every classified cCRE empties the file (quirk 11, twice).
14. Step 7's `blank.txt` line, `blank blank`, passes `$2 > 1.64` as text. It matches no cCRE.
15. calculate-center-distance.py's `elif` never fires. The first closest TSS bedtools prints wins.

### Host dependence removed

- Collation: `_bed` sorts with glibc's strxfrm under en_US.UTF-8. It checks at import that glibc still
  orders a probe set as the release did.
- log: `_zscore` ports glibc 2.39's non-FMA `log` to plain double operations. It matches non-FMA glibc on
  all 10,800,000 six-digit values in [1e-6, 1e6).
- Sums: numpy < 2.3 summed in 8192-element buffers. `_zscore` transliterates the buffered pairwise sum.

## Reference release (`REF`)

- `REF` is `/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry/V7-{hg38,mm10}`. Read only.
- `pixi run reference` takes `-m <step>` and `-k <genome>`. Steps: tf_list, dhs, rdhs, tf, anchors, signal,
  maxz, ccres, cts, pls. Modes absorb D1 (sorted lists) and D4 (z within 1e-14). `.gz` outputs stream decompressed.
- `data/manifests/2026-10-13/` holds sizes, mtimes and MD5s of `REF`. `md5sum -c --quiet md5.txt` there
  lists what changed since.
- `REF` edits after May 2022, which the test or the frozen inputs absorb:
  - hg38 `TF-List.Filtered.txt` (Feb 2023). Compare with columns 1-8 of `Filtered-Mod`.
  - mm10 `TF-List.Filtered.txt`: column 6 went from `,` to `:` (2023).
  - `TF-List.Filtered-Mod.txt` (Apr 2025) keeps the 2,509 hg38 experiments of May 2022, but flags 279
    biosamples TF=yes instead of 121.
  - mm10 `Master-Cell-List.txt` and `Individual-Files` were rebuilt in Aug 2025. The 2022 files are in `Old/`.
  - `Previous-Versions/cCREs.*.bed.gz` gained V4 rows (2025). mm10 mark lists changed in 2025, the mm10
    hotspot list in 2023. `REF/maxZ/` was rebuilt on GENCODE47 (2025).
- hg38 `REF` has no `tmp.bed`. The script's `rDHS-PreAltius.txt` is `rDHS-PreStam.txt`.

## Gotchas

- `data/` is a gitignored symlink to `/data/zusers/ramirezc/data/ccres-smk`.
- Inputs come from `freeze_input`, a checkpoint. A frozen file MUST match `config/inputs.b3`.
  A new input needs a new checksum line, written by hand after review.
- Shell blocks MUST call scripts as `{workflow.basedir:q}/scripts/x.py`.
- A shell command is one argv string, capped at 128 KiB. Rules over thousands of files MUST pass
  list files or directories, never `{input}`.
- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- Resources MUST be flat values. zlurm has no accounts. A runtime above the partition cap pends forever.
- The rule env has no PyYAML. The launcher env has it.
- `np.isin` on string arrays took 218 s per call. Use polars `is_in`.
- A pipe such as `pixi run lint | tail` hides the exit code. Check it before committing.
