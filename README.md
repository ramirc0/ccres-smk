# ccres-smk

A Snakemake port of the ENCODE4 cCRE pipeline
([weng-lab/ENCODE-cCREs, Version-4/cCRE-Pipeline](https://github.com/weng-lab/ENCODE-cCREs/tree/master/Version-4/cCRE-Pipeline),
steps 0-8) for hg38 and mm10. Outputs MUST match the May 2022 V7 release byte for byte, except for
the deviations listed below. Where the original depends on the host or on set order, the port
fixes it, and its output becomes the new baseline.

## Layout

```
.tests/integration/      toy inputs (one window per genome) for `pixi run test`
.tests/reference/        comparison of a run with the May 2022 release
.tests/unit/             differential tests against bedtools 2.30.0, GNU sort and glibc
config/
  config.yaml.template   run settings, per-genome constants, input sources
  inputs.b3              BLAKE3 sums of every frozen input
workflow/
  Snakefile              includes the rules + rule all
  rules/
    common.smk           config, list readers, ENCODE path resolution, wildcards
    inputs.smk           freeze_input checkpoint, download_encode
    dhs.smk              steps 0-2
    tf.smk               step 3
    signal.smk           anchors, steps 4-5
    classify.smk         steps 6-8
  scripts/               argparse scripts, each with build_parser() + main()
    _bed.py              exact BED text IO, en_US GNU sort, bedtools 2.30.0 interval ops
    _registry.py         representative-peak loop, region accessions
    _zscore.py           glibc's non-FMA log and numpy's buffered pairwise sum
    _cres.py             awk joins and the TSS/proximal/distal split of steps 6-7
    _encode.py           ENCODE file paths: local mirror, then downloads
  envs/pixi.toml         rule env (pixi workspace)
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      zlurm (partition per rule, no account)
pixi.toml                launcher env + test and lint tasks
pyproject.toml           ruff + numpydoc config
```

Data lives under `data/`, a gitignored symlink to the data volume:

```
data/inputs/<genome>/            frozen inputs
data/encode/<exp>/<file>.<ext>   ENCODE files missing from the local mirror
data/results/<run_id>/<genome>/  outputs, laid out like the release's V7-<genome>/
```

## Pipeline

| Step | Rule | Output under `data/results/<run_id>/<genome>/` | Original |
|---|---|---|---|
| 0 | `call_dhs` | `rDHS/DHSs/<peak>.{DHSs,Excluded}.bed` | 0_Call-DHSs.sh |
| 1 | `process_dhs` | `rDHS/Processed-DHSs/<signal>.<row>.bed` | 1_Process-DHSs.sh |
| 2 | `combine_dhs`, `create_rdhs` | `rDHS/<g>-DHS-{All,Filtered}.bed`, `rDHS/<g>-rDHS-*` | 2_Create-rDHSs.sh |
| 3 | `tf_list`, `tf_peaks`, `tf_rpeaks`, `tf_clusters` | `TF/` | 3_Curate-TF-Clusters.sh |
| - | `anchors` | `<g>-Anchors.bed` | none (rDHSs + rTFBSs, sorted) |
| 4 | `extend_anchors`, `signal_zscore` | `signal-output/<mark>/<exp>-<file>.txt` | 4_Calculate-Signal-Zscores.sh |
| 5 | `max_zscore` | `<g>-<mark>-maxZ.txt` | 5_Determine-Max-Zscores.sh |
| 6 | `classify_ccres` | `<g>-cCREs-Unfiltered.bed` | 6_Classify-cCREs.sh |
| 7 | `master_cell_list`, `cell_type_ccres` | `Cell-Type-Specific/` | 7_Call-Cell-Type-Specific.sh |
| 8 | `pls_genes` | `<g>-PLS-Gene-List.txt` | 8_Assign-PLS-Gene.sh |

Step 0 sweeps 4,941 FDR thresholds. The port merges only at thresholds where some position has
that FDR order, since elsewhere the merged set repeats. It also drops each merged interval once it
settles. An interval of 50-350 bp is a DHS, and one under 50 bp can only shrink, so neither changes
the output again. The largest hotspot file (72M positions) takes 35 s.

## Inputs

`freeze_input` copies each source in `config.yaml` to `data/inputs/<genome>/` and fails unless the
copy matches its line in `config/inputs.b3`. Sources can change upstream, so a changed source can
never reach a run unnoticed. Two inputs are derived:

- `cCREs-Previous.bed`: the V1-V3 rows of `Previous-Versions/cCREs.<genome>.bed.gz`. The file now
  also holds V4, the release this port reproduces.
- `TF-List.Cell-Type.txt`: the rows of `TF-List.Filtered-Mod.txt` whose biosample the release
  flagged TF=yes in `Master-Cell-List.txt` (121 hg38 biosamples). A 2025 edit of the list matches
  more biosamples, which would move them to step 7's TF branch.

ENCODE files resolve to `<encode_root>/<exp>/<file>.<ext>`, a local mirror. A file missing there
downloads from ENCODE's public S3 bucket to `data/encode/`. Its key and MD5 come from a snapshot of
the ENCODE file report (`encode_metadata` in the config).

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
pixi run snakemake -n -p --profile profiles/local  # dry run
pixi run snakemake --profile profiles/local        # local
pixi run snakemake --profile profiles/slurm        # zlurm
pixi run test                                      # all rules on the toy inputs
pixi run unit                                      # unit tests
```

## Fidelity

### Deviations from the release

- **D1. Set order.** filter-tf-rpeaks.py joined the targets, biosamples and experiments of a TF
  cluster in Python set() order, which changes from run to run. The port joins them sorted. This
  affects `TF-Clusters.txt`, `tmp.no-overlap` and `rTFBS-Summary.txt`.
- **D2. Step 7 product.** The port writes bgzipped, tabix-indexed BED. It writes no bigBed.
- **D3. Step 1 names.** `Processed-DHSs/<signal>.<row>.bed` replaces `output.<signal>.<row>`. Step 2
  still concatenates them in the en_US order of the old names.
- **D4. ATAC z-scores.** The release's ATAC files come from a September 2023 run whose glibc used
  its FMA log. Recomputing with the FMA log matches them exactly. The port uses the non-FMA log of
  the May 2022 run, which moves about 150 of 2.9M z-scores per file by at most 9e-16. Signal and
  rank columns stay identical. Steps 6-8 do not read ATAC.
- **D5. Step 7 file set.** The release's hg38 `Individual-Files` mix runs of two master cell lists.
  `DNase-List.txt` changed on May 26, 2022, and only some biosamples were rerun. The port follows
  the final list. Its two biosamples whose H3K4me3 file is `NA` (quirk 10) have no release file.

### Where the port follows the release over the published scripts

- Step 8: the published calculate-center-distance.py keeps a TSS whose start lies less than 200 bp
  from the PLS center. `PLS-Gene-List.hg38-V4.txt` also keeps the 93 hg38 pairs at exactly 200 bp.
  The port keeps them.
- Steps 3 and 7: ENCSR008LHT's peak file was not available to the May 2022 run, so its peaks are
  missing from `tmp.summits`. `tf_excluded` leaves it out.

### Quirks kept

The port reproduces these deterministic behaviors of the original. Fixing any of them changes
the output.

1. Step 0 sweeps FDR thresholds up to 1E-4942. Both genomes have positions of order 4943 and 4944,
   which no threshold filters. They stay in the last merged set.
2. Step 0 prints an FDR without an exponent as `'%.2E'` before reading its order, so rounding sets
   the order (0.009999 reads as order 2). An FDR whose text has a `+` exponent, such as 0, passes
   every threshold.
3. bedtools merge prints the minimum FDR with 10 significant digits. An FDR below the double range,
   such as 1e-4944, prints as 0.
4. Step 1 sorts names with `sort -k4,4`, whose key includes the tab before field 4.
5. Step 2's `grep -v chrEBV` and `grep -v chrM` drop any line containing those strings.
6. make-region-accession.py keys earlier regions by chrom, start and end concatenated without
   separators. Distinct coordinates can collide.
7. Step 3's `tmp.summits` places a summit at the new start plus the original offset. For a resized
   peak that point is not the summit.
8. Step 3's awk compares a non-numeric FRiP such as NA with 0.003 as text, so it passes.
9. Step 4's rank loop starts from `before=0`. A top z-score of exactly 0 would rank 0.
10. Two hg38 H3K4me3 rows list the file `NA`. Their z-score files are empty. Step 5 skips them, and
    step 7 puts `NA` in their biosamples' file names.
11. In steps 6-7, `awk 'FNR==NR {x[$1];next} ...' first second` with an empty first file selects
    nothing, even under `!(key in x)`.
12. In step 7's No-DNase branch, an empty no1 or no2 falls back to the whole group or to no1. One
    cCRE can then get two rows.
13. In the same branch, a mark whose above-1.64 lines miss every classified cCRE empties the
    biosample's file (quirk 11, twice).
14. Step 7's `blank.txt` placeholder line, `blank blank`, passes `$2 > 1.64` as text. It matches no
    cCRE.
15. calculate-center-distance.py meant to keep the nearest of several closest TSSs, but its `elif`
    never fires. The first TSS bedtools prints wins.

### Host dependence removed

- Collation. The original sorted and globbed under en_US.UTF-8. `_bed` sorts with glibc's strxfrm
  and checks at import that glibc still orders a probe set as the release did.
- log. On CPUs with FMA, glibc picks an FMA build of `log`, which rounds differently (D4).
  `_zscore` ports glibc 2.39's non-FMA `log` to plain double operations. It matches non-FMA glibc
  on all 10,800,000 six-digit values in [1e-6, 1e6).
- Sums. numpy < 2.3 summed in 8192-element buffers. numpy 2.3 dropped them. `_zscore` transliterates
  the buffered pairwise sum.

## Conventions

Follow the [Nextstrain Snakemake style guide][style-guide] and
[Snakemake's best-practices guide][best-practices]:

- No `run:` blocks; custom logic lives in `workflow/scripts/` called from `shell:`.
- Raw, triple-quoted `shell` blocks (`r"""`), one option per line.
- Every interpolation quoted with `:q`; config mapped into `shell` via `params:`.
- `input:`/`output:` are literal path strings; all paths relative.
- Logs tee'd to file and terminal (`exec &> >(tee {log:q})`).
- Every rule has a `log:`, `benchmark:` and `software:`.
- Scripts are called as `{workflow.basedir:q}/scripts/x.py`.
- Docstrings follow numpydoc. `pixi run lint` checks the scripts with ruff and `numpydoc lint`.

[style-guide]: https://docs.nextstrain.org/en/latest/reference/snakemake-style-guide.html
[best-practices]: https://snakemake.readthedocs.io/en/stable/snakefiles/best_practices.html

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake, the pixi software-deployment plugin and
the SLURM executor plugins. The rule env is a second pixi workspace, `workflow/envs/`. Its
`pixi.lock` pins every package, including UCSC tools 482, whose output bytes the release depends
on. Its `oracle` env holds bedtools 2.30.0, the version of the release, for the unit tests only.
An `onstart` hook installs the env before the first job, since parallel SLURM jobs would race to
install it.
