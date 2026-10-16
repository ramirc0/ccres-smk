"""Step 8: assign promoter-like cCREs (PLS) to genes.

Replaces 8_Assign-PLS-Gene.sh with calculate-center-distance.py (assignment).
A PLS gets the gene of every TSS that step 6's TSS test accepts: a TSS base
inside it, or a TSS start within 200 bp of its center, bounds included. Pairs
are distinct and sorted.

The original took the TSSs a PLS overlaps. Only for a PLS that overlapped
none did it take the closest TSS, if its start lay within 200 bp of the
center. It missed the other TSSs near the center.
"""

import argparse

import _bed
import _cres
import polars as pl

CCRE = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "rdhs": pl.String,
    "ccre": pl.String,
    "class": pl.String,
}
MAX_DISTANCE = 200


def build_parser():
    """Build the pls_genes.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the pls_genes.py flags.
    """
    p = argparse.ArgumentParser(description="Assign PLS cCREs to genes.")
    p.add_argument("--ccres", required=True, help="cCREs-Unfiltered BED.")
    p.add_argument("--tss", required=True, help="TSS BED, zero-length.")
    p.add_argument("-o", "--output", required=True, help="cCRE-gene pairs to write.")
    return p


def main(argv=None):
    """Write the sorted, distinct cCRE-gene pairs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    pls = _bed.read_bed(args.ccres, CCRE).filter(pl.col("class") == "PLS")
    tss = _cres.read_tss(args.tss)
    center = (pl.col("start") + pl.col("end")) / 2
    # A TSS start s is near when center - 200 <= s <= center + 200.
    window = pls.select(
        "chrom",
        start=pl.min_horizontal("start", (center - MAX_DISTANCE).ceil().cast(pl.Int64)),
        end=pl.max_horizontal(
            "end", (center + MAX_DISTANCE).floor().cast(pl.Int64) + 1
        ),
    )
    pairs = _bed.overlap_pairs(window, tss)
    genes = pl.DataFrame(
        {"ccre": pls["ccre"][pairs["a"]], "gene": tss["gene"][pairs["b"]]}
    )
    _bed.write_bed(genes.unique().sort(genes.columns), args.output)


if __name__ == "__main__":
    main()
