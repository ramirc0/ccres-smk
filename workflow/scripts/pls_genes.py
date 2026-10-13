"""Step 8: assign promoter-like cCREs (PLS) to genes.

Port of 8_Assign-PLS-Gene.sh with calculate-center-distance.py (assignment).
A PLS gets the gene of every TSS it overlaps. A PLS overlapping none but
within a TSS +-2 kb window gets the gene of each closest TSS (all ties) whose
start lies at most 200 bp from the PLS center. Pairs go through `sort -u`.

The published script tests `abs(diff) < 200`. The release's list was made with
`<= 200`: its 93 hg38 pairs at exactly 200 bp are kept here.
"""

import argparse

import _bed
import _cres
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
    p.add_argument("--tss-4k", required=True, help="TSS +-2 kb BED.")
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
    window = _bed.read_tsv(args.tss_4k, _cres.TSS)
    codes = _bed.shared_codes(pls["chrom"], tss["chrom"], window["chrom"])
    pls_iv, tss_iv, window_iv = (
        (c, t["start"].cast(pl.Int64).to_numpy(), t["end"].cast(pl.Int64).to_numpy())
        for c, t in zip(codes, (pls, tss, window))
    )
    ccre, gene = pls["ccre"].to_list(), tss["gene"].to_list()

    ia, ib, _ = _bed.overlap_pairs(pls_iv, tss_iv)
    pairs = [f"{ccre[i]}\t{gene[j]}" for i, j in zip(ia, ib)]

    on_tss = _bed.overlaps_any(pls_iv, tss_iv)
    near = (~on_tss) & _bed.overlaps_any(pls_iv, window_iv)
    rows = near.nonzero()[0]
    near_iv = tuple(a[rows] for a in pls_iv)
    ia, ib = _bed.closest_pairs(near_iv, tss_iv)
    # The strand only sets the sign of the distance, and the test is |distance| < 200.
    distance = abs(tss_iv[1][ib] - (near_iv[1][ia] + near_iv[2][ia]) / 2)
    pairs += [
        f"{ccre[rows[i]]}\t{gene[j]}"
        for i, j, d in zip(ia, ib, distance)
        if d <= MAX_DISTANCE
    ]
    with open(args.output, "w") as f:
        f.writelines(f"{line}\n" for line in _bed.sort_unique(pairs))


if __name__ == "__main__":
    main()
