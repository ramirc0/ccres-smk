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
