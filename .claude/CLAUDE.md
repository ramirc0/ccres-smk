# CLAUDE.md

Snakemake port of the ENCODE4 cCRE pipeline (Moore lab, steps 0-8). Outputs MUST match the May 2022 V7
release byte for byte, except the deviations listed in README.md. README.md holds the layout, idioms
and style conventions. Follow it.

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
- Keep deterministic bugs of the original. Record each in the README quirks ledger.
- Fix non-determinism and host dependence (set order, glob order, FMA). The fix becomes the new baseline.
- Logic lives in Python. External binaries stay only where output bytes depend on them.
- The reference tree belongs to another user and keeps changing. Check a file's mtime before trusting it.

## Reference release (`REF`)

- `REF` is `/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/V7/Registry/V7-{hg38,mm10}`. Read only.
- `pixi run reference` takes `-m <step>` and `-k <genome>`. Steps: tf_list, dhs, rdhs, tf, anchors, signal,
  maxz, ccres, cts, pls. Comparison modes absorb D1 (sorted lists), D2 (gunzip) and D4 (z within 1e-14).
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
- hg38 `REF` has no `tmp.bed` or `tmp.no-overlap`. The script's `rDHS-PreAltius.txt` is `rDHS-PreStam.txt`.
- D5: 1,675 hg38 `Individual-Files` date from May 17-19, 4 from May 26. `ENCFF346LEZ_ENCFF118OBT_ENCFF536VOI.bed`
  is stale. `ENCFF269VAY.bigBed` has no `.bed`.

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
