# ccres-smk

A Snakemake port of the ENCODE4 cCRE pipeline
([weng-lab/ENCODE-cCREs, Version-4/cCRE-Pipeline](https://github.com/weng-lab/ENCODE-cCREs/tree/master/Version-4/cCRE-Pipeline),
steps 0-8) for hg38 and mm10. Outputs MUST match the May 2022 V7 release byte for byte. Where the
original depends on the host or on set and glob order, the port fixes it, and its output becomes the
new baseline.

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
| 7 | `master_cell_list`, `frozen_master_cell_list`, `cell_type_ccres`, `cell_type_bigbed` | `Cell-Type-Specific/` | 7_Call-Cell-Type-Specific.sh |
| 8 | `pls_genes` | `<g>-PLS-Gene-List.txt` | 8_Assign-PLS-Gene.sh |

Step 0 sweeps 4,941 FDR thresholds. The port merges only at thresholds where some position has
that FDR order, since elsewhere the merged set repeats. It also drops each merged interval once it
settles. An interval of 50-350 bp is a DHS, and one under 50 bp can only shrink, so neither changes
the output again. The largest hotspot file (117M positions) takes under 90 s.

## Inputs

`freeze_input` copies each source in `config.yaml` to `data/inputs/<genome>/` and fails unless the
copy matches its line in `config/inputs.b3`. Sources can change upstream, so a changed source can
never reach a run unnoticed. Two inputs are derived:

- `cCREs-Previous.bed`: the V1-V3 rows of `Previous-Versions/cCREs.<genome>.bed.gz`. The file now
  also holds V4, the release this port reproduces.
- `TF-List.Cell-Type.txt`: the TF experiments of the biosamples the release flagged TF=yes in
  `Master-Cell-List.txt`. hg38 takes them from `TF-List.Filtered-Mod.txt` (121 biosamples). A 2025
  edit of that list matches more biosamples, which would move them to step 7's TF branch. mm10
  takes them from two May 2022 files that hold the list's columns (36 biosamples).

mm10's step 7 freezes the release's May 2022 `Master-Cell-List.txt`. The mm10 mark lists were
edited in 2025 and no longer rebuild it.

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
