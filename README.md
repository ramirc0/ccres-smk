# ccres-smk

A Snakemake port of the ENCODE4 cCRE pipeline
([weng-lab/ENCODE-cCREs, Version-4/cCRE-Pipeline](https://github.com/weng-lab/ENCODE-cCREs/tree/master/Version-4/cCRE-Pipeline),
steps 0-8) for hg38 and mm10. Outputs MUST match the May 2022 V7 release byte for byte. Where the
original depends on the host or on set and glob order, the port fixes it, and its output becomes the
new baseline.

## Layout

```
.tests/reference/        pytest comparison of a run with the May 2022 release
.tests/unit/             differential tests of the interval and collation code
config/
  config.yaml.template   run settings, per-genome constants, input sources
  inputs.b3              BLAKE3 sums of every frozen input
workflow/
  Snakefile              includes the rules + rule all
  rules/
    common.smk           config, list readers, ENCODE path resolution, wildcards
    inputs.smk           freeze_input checkpoint, download_encode
  scripts/               argparse scripts, each with build_parser() + main()
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

## Inputs

`freeze_input` copies each source in `config.yaml` to `data/inputs/<genome>/` and fails unless the
copy matches its line in `config/inputs.b3`. The sources sit in another user's tree and keep
changing, so a changed source can never reach a run unnoticed. Two inputs are derived:

- `cCREs-Previous.bed`: the V1-V3 rows of `Previous-Versions/cCREs.<genome>.bed.gz`. The file now
  also holds V4, the release this port reproduces.
- `TF-List.Cell-Type.txt`: rows of today's `TF-List.Filtered-Mod.txt` for the 121 biosamples the
  release flagged TF=yes in `Master-Cell-List.txt`.

ENCODE files resolve to the local mirror, `/data/projects/encode/data/<exp>/<file>.<ext>`. Missing
ones download from the public portal, with their MD5 checked against the portal record.

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
pixi run snakemake -n -p --profile profiles/local  # dry run
pixi run snakemake --profile profiles/local        # local
pixi run snakemake --profile profiles/slurm        # zlurm
pixi run unit                                      # unit tests
pixi run reference -m <step> -n 16                 # compare with the release
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
