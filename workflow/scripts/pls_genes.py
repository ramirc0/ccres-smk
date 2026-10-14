"""Step 8: assign promoter-like cCREs (PLS) to genes.

Replaces 8_Assign-PLS-Gene.sh with calculate-center-distance.py (assignment).
A PLS gets the gene of every TSS that step 6's TSS test accepts: a TSS inside
it, bounds included as bedtools reads a zero-length TSS, or within 200 bp of
its center, bounds included. Pairs go through `sort -u`.

The original took the TSSs a PLS overlaps. Only for a PLS that overlapped
none did it take the closest TSS, if its start lay within 200 bp of the
center. It missed the other TSSs near the center.
"""

import argparse

import _bed
import _cres
import numpy as np
import polars as pl

CCRE = ["chrom", "start", "end", "rdhs", "ccre", "class"]
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
    """Write the sorted, deduplicated cCRE-gene pairs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    pls = _bed.read_tsv(args.ccres, CCRE).filter(pl.col("class") == "PLS")
    tss = _bed.read_tsv(args.tss, _cres.TSS)
    pls_codes, tss_codes = _bed.shared_codes(pls["chrom"], tss["chrom"])
    start = pls["start"].cast(pl.Int64).to_numpy()
    end = pls["end"].cast(pl.Int64).to_numpy()
    center = (start + end) / 2
    # The center lies in [start, end], so this window is the union of the PLS,
    # end included, and the positions within 200 bp of the center.
    window = (
        pls_codes,
        np.minimum(start, np.ceil(center - MAX_DISTANCE)).astype(np.int64),
        np.maximum(end, np.floor(center + MAX_DISTANCE)).astype(np.int64) + 1,
    )
    position = tss["start"].cast(pl.Int64).to_numpy()
    ia, ib, _ = _bed.overlap_pairs(window, (tss_codes, position, position + 1))
    ccre, gene = pls["ccre"].to_list(), tss["gene"].to_list()
    pairs = [f"{ccre[i]}\t{gene[j]}" for i, j in zip(ia.tolist(), ib.tolist())]
    with open(args.output, "w") as f:
        f.writelines(f"{line}\n" for line in _bed.sort_unique(pairs))


if __name__ == "__main__":
    main()
