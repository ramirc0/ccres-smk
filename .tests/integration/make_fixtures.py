"""Build the integration test inputs: real data cut to one window per genome.

Run once from the repository root, after the frozen inputs exist:

    pixi run --manifest-path workflow/envs/pixi.toml python .tests/integration/make_fixtures.py

The lists are synthetic. Their rows point to real ENCODE files, and their
biosample names are chosen so that every branch of step 7 runs: toy_tf (DNase
and TF), toy_dnase (DNase, partial marks) and toy_nodnase (no DNase, an NA
H3K4me3 file). As in the real config, mm10 freezes its step 7 master list. The
ENCODE mirror at /data/projects/encode/data must exist.
"""

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import blake3

HERE = Path(__file__).parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
import master_cell_list
import tf_list

MIRROR = Path("/data/projects/encode/data")
PREVIOUS = Path(
    "/data/zusers/moorej3/moorej.ghpcc.project/ENCODE/Encyclopedia/Previous-Versions"
)
WINDOWS = {
    "hg38": ("chr21", 33_000_000, 33_500_000),
    "mm10": ("chr19", 5_000_000, 5_500_000),
}
TF_EXPERIMENTS = 6


def rows(path):
    """Read a tab-separated file.

    Parameters
    ----------
    path : pathlib.Path
        File to read.

    Returns
    -------
    list of list of str
        Fields per line.
    """
    with open(path) as f:
        return [line.rstrip("\n").split("\t") for line in f]


def write(path, table):
    """Write rows as tab-separated lines, creating the directory.

    Parameters
    ----------
    path : pathlib.Path
        File to write.
    table : list of list of str
        Rows.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.writelines("\t".join(r) + "\n" for r in table)


def in_window(fields, window):
    """Tell whether a BED row lies inside the window.

    Parameters
    ----------
    fields : list of str
        BED fields.
    window : tuple of (str, int, int)
        Chromosome, start and end.

    Returns
    -------
    bool
        True if the row starts and ends within the window.
    """
    chrom, start, end = window
    return fields[0] == chrom and int(fields[1]) >= start and int(fields[2]) <= end


def cut_bed(source, target, window, compress):
    """Copy the rows of a BED file that lie in the window.

    Parameters
    ----------
    source : pathlib.Path
        BED file, gzipped if its name ends in .gz.
    target : pathlib.Path
        BED file to write.
    window : tuple of (str, int, int)
        Chromosome, start and end.
    compress : bool
        Gzip the output.
    """
    opener = gzip.open if source.suffix == ".gz" else open
    with opener(source, "rt") as f:
        text = "".join(line for line in f if in_window(line.split("\t"), window))
    target.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 keeps the gzip bytes the same from run to run.
    data = text.encode()
    target.write_bytes(gzip.compress(data, mtime=0) if compress else data)


def cut_bigwig(source, target, window, sizes):
    """Cut a bigWig to the window with the UCSC tools.

    Parameters
    ----------
    source, target : pathlib.Path
        Input and output bigWig.
    window : tuple of (str, int, int)
        Chromosome, start and end.
    sizes : pathlib.Path
        Chromosome sizes.
    """
    chrom, start, end = window
    target.parent.mkdir(parents=True, exist_ok=True)
    graph = target.with_suffix(".bedGraph")
    subprocess.run(
        [
            "bigWigToBedGraph",
            f"-chrom={chrom}",
            f"-start={start}",
            f"-end={end}",
            source,
            graph,
        ],
        check=True,
    )
    subprocess.run(["bedGraphToBigWig", graph, sizes, target], check=True)
    graph.unlink()


def build(genome):
    """Write the fixture lists and files of one genome.

    Parameters
    ----------
    genome : str
        hg38 or mm10.

    Returns
    -------
    dict of str to str
        Frozen input name to fixture source path(s), relative to the test dir.
    """
    window = WINDOWS[genome]
    frozen = ROOT / "data" / "inputs" / genome
    res = HERE / "resources" / genome
    encode = HERE / "resources" / "encode"
    sizes = res / "chromInfo.txt"
    sizes.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(frozen / "chromInfo.txt", sizes)

    def encode_copy(exp, acc, ext):
        source, target = MIRROR / exp / f"{acc}.{ext}", encode / exp / f"{acc}.{ext}"
        if ext == "bigWig":
            cut_bigwig(source, target, window, sizes)
        else:
            cut_bed(source, target, window, compress=True)

    hotspots = rows(frozen / "Hotspot-List.txt")[:3]
    for r in hotspots:
        encode_copy(r[0], r[1], "bed.gz")
        encode_copy(r[0], r[2], "bigWig")
    write(res / "Hotspot-List.txt", hotspots)

    marks = {}
    marks["DNase"] = [
        [r[0], r[2], name] for r, name in zip(hotspots, ("toy_tf", "toy_dnase"))
    ]
    for mark, names in (
        ("H3K4me3", ("toy_tf",)),
        ("H3K27ac", ("toy_tf", "toy_dnase")),
        ("CTCF", ("toy_tf", "toy_nodnase")),
        ("ATAC", ("toy_atac",)),
    ):
        real = rows(frozen / f"{mark}-List.txt")[: len(names)]
        marks[mark] = [[r[0], r[1], name] for r, name in zip(real, names)]
    marks["H3K4me3"].append(["ENCSR999NAX", "NA", "toy_nodnase"])
    for mark, table in marks.items():
        for exp, acc, _ in table:
            if acc != "NA":
                encode_copy(exp, acc, "bigWig")
        write(res / f"{mark}-List.txt", table)

    tf_all = rows(frozen / "TF-List.All.txt")
    passing = [r for r in tf_all if tf_list.passes(r)]
    counts = []
    for r in passing[:200]:
        with gzip.open(MIRROR / r[0] / f"{r[1]}.bed.gz", "rt") as f:
            counts.append((sum(in_window(line.split("\t"), window) for line in f), r))
    chosen = [r for _, r in sorted(counts, key=lambda c: -c[0])[:TF_EXPERIMENTS]]
    for r in chosen:
        encode_copy(r[0], r[1], "bed.gz")
    write(res / "TF-List.All.txt", chosen + [r for r in tf_all if r[1] == "NA"][:1])
    tf_rows = [[*r, "ENCDO000TOY", "toy_tf"] for r in chosen[:2]]
    if genome == "mm10":
        # As in the real config, mm10's step 7 freezes its master list, and the
        # TF list comes as two column files.
        tf_files = ["TF-Filtered.txt", "TF-Biosamples.txt"]
        write(res / tf_files[0], [r[:-1] for r in tf_rows])
        write(res / tf_files[1], [r[-1:] for r in tf_rows])
        mark_lists = [res / f"{m}-List.txt" for m in master_cell_list.MARKS]
        write(
            res / "Master-Cell-List.txt",
            master_cell_list.master_rows(mark_lists, res / tf_files[1]),
        )
    else:
        tf_files = ["TF-Filtered-Mod.txt"]
        write(res / tf_files[0], tf_rows)
        write(res / "Master-Cell-List.txt", [["---"] * 8 + ["toy_tf", "Group1", "yes"]])

    # cDHSs cover only the first half, so that some rDHSs fall in the none class.
    chrom, start, end = window
    cut_bed(
        frozen / "cDHSs.bed",
        res / "cDHSs.bed",
        (chrom, start, (start + end) // 2),
        False,
    )
    for name in ("TSS.Basic.bed", "TSS.Basic.4K.bed"):
        cut_bed(frozen / name, res / name, window, compress=False)
    cut_bed(
        PREVIOUS / f"rDHS.{genome}.bed.gz", res / "rDHS-Previous.bed.gz", window, True
    )
    cut_bed(PREVIOUS / f"cCREs.{genome}.bed.gz", res / "cCREs.bed.gz", window, True)
    shutil.copy(frozen / "cres.as", res / "cres.as")

    rel = f"resources/{genome}"
    sources = {
        name: f"{rel}/{name}"
        for name in (
            "Hotspot-List.txt",
            *(f"{m}-List.txt" for m in marks),
            "TF-List.All.txt",
            "cDHSs.bed",
            "rDHS-Previous.bed.gz",
            "TSS.Basic.bed",
            "TSS.Basic.4K.bed",
            "chromInfo.txt",
            "cres.as",
        )
    }
    sources["cCREs-Previous.bed"] = f"{rel}/cCREs.bed.gz"
    sources["TF-List.Cell-Type.txt"] = [
        f"{rel}/{name}" for name in (*tf_files, "Master-Cell-List.txt")
    ]
    if genome == "mm10":
        sources["Master-Cell-List.txt"] = f"{rel}/Master-Cell-List.txt"
    return sources


def main():
    """Build both genomes' fixtures, then the test config and checksums."""
    sources = {genome: build(genome) for genome in WINDOWS}
    config = (ROOT / "config" / "config.yaml.template").read_text()
    head = config[: config.index("sources:")]
    head = head.replace(
        "encode_root: /data/projects/encode/data", "encode_root: resources/encode"
    )
    lines = [head, "sources:\n"]
    for genome, names in sources.items():
        lines.append(f"  {genome}:\n")
        for name, src in names.items():
            if isinstance(src, list):
                lines.append(f"    {name}:\n" + "".join(f"      - {s}\n" for s in src))
            else:
                lines.append(f"    {name}: {src}\n")
    (HERE / "config").mkdir(exist_ok=True)
    (HERE / "config" / "config.yaml").write_text("".join(lines))

    sums = []
    for genome, names in sources.items():
        for name, src in names.items():
            frozen = HERE / "data" / "inputs" / genome / name
            sums.append(f"{_digest(src, name, frozen)}  {genome}/{name}\n")
    (HERE / "config" / "inputs.b3").write_text("".join(sums))


def _digest(src, name, scratch):
    """Return the BLAKE3 digest a source freezes to.

    Parameters
    ----------
    src : str or list of str
        Source path(s), relative to the test dir.
    name : str
        Frozen input name.
    scratch : pathlib.Path
        Where to derive the input when it is not a plain copy.

    Returns
    -------
    str
        Hex digest.
    """
    import freeze_input

    scratch.parent.mkdir(parents=True, exist_ok=True)
    if name == "cCREs-Previous.bed":
        freeze_input.previous_ccres(HERE / src, scratch)
    elif name == "TF-List.Cell-Type.txt":
        freeze_input.cell_type_tf([HERE / s for s in src[:-1]], HERE / src[-1], scratch)
    else:
        shutil.copy(HERE / src, scratch)
    digest = blake3.blake3(scratch.read_bytes()).hexdigest()
    scratch.unlink()
    return digest


if __name__ == "__main__":
    main()
