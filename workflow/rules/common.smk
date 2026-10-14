# Config, per-genome lists, ENCODE path resolution, and wildcards.

import os
import sys
from functools import cache

from snakemake.exceptions import WorkflowError

sys.path.insert(0, f"{workflow.basedir}/scripts")
from _encode import encode_path


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run, one subdirectory per genome
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
INPUTS = config["inputs_dir"]           # run-independent frozen inputs
ENCODE_ROOT = config["encode_root"]
DOWNLOAD_DIR = config["download_dir"]
GENOMES = list(config["genomes"])
MARKS = ["DNase", "H3K4me3", "H3K27ac", "CTCF", "ATAC"]

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


# The pixi plugin installs the env inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install it once before any job starts.
onstart:
    shell("pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked")


wildcard_constraints:
    genome="|".join(GENOMES),
    mark="|".join(MARKS),
    exp=r"ENCSR\d{3}[A-Z]{3}",
    acc=r"ENCFF\d{3}[A-Z]{3}",
    row=r"\d+",


def frozen(genome, name):
    """Return the path of a frozen input, which `freeze_input` writes.

    Parameters
    ----------
    genome : str
        Genome key in `config["genomes"]`.
    name : str
        Frozen input name, a key of `config["sources"][genome]`.

    Returns
    -------
    str
        Path `<inputs_dir>/<genome>/<name>`.
    """
    return f"{INPUTS}/{genome}/{name}"


def frozen_list(genome, name):
    """Return the path of a frozen list that input functions read.

    These lists set which jobs exist, so they must be frozen before the DAG
    is built: `snakemake freeze` writes them.

    Parameters
    ----------
    genome : str
        Genome key in `config["genomes"]`.
    name : str
        Frozen input name, a key of `config["sources"][genome]`.

    Returns
    -------
    str
        Path `<inputs_dir>/<genome>/<name>`.

    Raises
    ------
    WorkflowError
        If the list is missing.
    """
    path = frozen(genome, name)
    if not os.path.exists(path):
        raise WorkflowError(f"{path} is missing. Run `snakemake freeze` first.")
    return path


@cache
def read_rows(path):
    """Read a tab-separated list file, once per path.

    Parameters
    ----------
    path : str
        File to read.

    Returns
    -------
    list of list of str
        Fields of each line, in file order.
    """
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


def encode_file(exp, acc, ext):
    """Return the path of an ENCODE file under the configured directories.

    Parameters
    ----------
    exp : str
        Experiment accession (`ENCSR...`).
    acc : str
        File accession (`ENCFF...`).
    ext : str
        File extension, such as `bed.gz` or `bigWig`.

    Returns
    -------
    str
        The mirror path if it exists, else the path `download_encode` writes.
    """
    return encode_path(ENCODE_ROOT, DOWNLOAD_DIR, exp, acc, ext)


def hotspot_row(genome, row):
    """Return one row of a genome's hotspot list.

    Parameters
    ----------
    genome : str
        Genome key.
    row : int or str
        1-based row number.

    Returns
    -------
    list of str
        Experiment, peak file, signal file, then the other columns.
    """
    return read_rows(frozen_list(genome, "Hotspot-List.txt"))[int(row) - 1]


def hotspot_experiment(genome, peak):
    """Return the experiment of a hotspot peak file.

    Parameters
    ----------
    genome : str
        Genome key.
    peak : str
        Peak file accession, column 2 of the hotspot list.

    Returns
    -------
    str
        Experiment accession.
    """
    rows = read_rows(frozen_list(genome, "Hotspot-List.txt"))
    (exp,) = [r[0] for r in rows if r[1] == peak]
    return exp
