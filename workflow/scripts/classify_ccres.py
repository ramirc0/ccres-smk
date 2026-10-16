"""Step 6: classify anchors into cell type-agnostic cCREs and accession them.

Port of 6_Classify-cCREs.sh. Anchors with DNase max z > 1.64 split into TSS,
proximal and distal groups (see `_cres.tss_groups`). Each group then splits by
max z > 1.64 of H3K4me3, H3K27ac and CTCF into PLS, pELS, dELS, CA-H3K4me3,
CA-CTCF and CA. CA anchors overlapping a TF cluster become CA-TF. Anchors with
DNase max z <= 1.64 that overlap a TF cluster become TF, whatever their
H3K4me3, H3K27ac and CTCF z-scores. The paper describes TF cCREs as low in
histone marks too. The cCREs are sorted by position and accessioned, reusing
the accessions of identical V1-V3 cCREs.
"""

import argparse

import _bed
import _cres
import _registry
import polars as pl
from _cres import high

PREVIOUS = {
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "accession": pl.String,
    "class": pl.String,
    "version": pl.String,
}
MARKS = ["dnase", "h3k4me3", "h3k27ac", "ctcf"]


def build_parser():
    """Build the classify_ccres.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the classify_ccres.py flags.
    """
    p = argparse.ArgumentParser(description="Classify cell type-agnostic cCREs.")
    p.add_argument("--anchors", required=True, help="Anchor Parquet.")
    for mark in ("DNase", "H3K4me3", "H3K27ac", "CTCF"):
        p.add_argument(f"--{mark.lower()}", required=True, help=f"{mark} maxZ table.")
    p.add_argument("--tss", required=True, help="TSS BED, zero-length.")
    p.add_argument("--tss-4k", required=True, help="TSS +-2 kb BED.")
    p.add_argument("--tf-clusters", required=True, help="TF cluster Parquet.")
    p.add_argument("--previous", required=True, help="V1-V3 cCRE BED.")
    p.add_argument("--genome", required=True, help="hg38 or mm10.")
    p.add_argument("-o", "--output", required=True, help="cCRE BED to write.")
    return p


def main(argv=None):
    """Write the accessioned cCREs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    anchors = _cres.with_z(
        pl.read_parquet(args.anchors),
        "accession",
        {mark: getattr(args, mark) for mark in MARKS},
    )
    anchors = anchors.with_columns(
        on_tf=_bed.overlaps_any(anchors, pl.read_parquet(args.tf_clusters))
    )
    dnase = anchors.filter(high("dnase"))
    groups = _cres.tss_groups(
        dnase, _cres.read_tss(args.tss), _bed.read_bed(args.tss_4k, _cres.TSS)
    )
    ca = pl.when("on_tf").then(pl.lit("CA-TF")).otherwise(pl.lit("CA"))
    classified = dnase.with_columns(groups).select(
        *_bed.COORDS, "accession", label=_cres.label(ca)
    )
    tf = anchors.filter("on_tf", ~high("dnase")).select(
        *_bed.COORDS, "accession", label=pl.lit("TF")
    )
    ccres = pl.concat([classified, tf]).sort(*_bed.COORDS, "accession")
    previous = _bed.read_bed(args.previous, PREVIOUS)
    ccre = _registry.accessions(ccres, args.genome, "cCRE", previous)
    _bed.write_bed(
        ccres.select(*_bed.COORDS, "accession", ccre.alias("ccre"), "label"),
        args.output,
    )


if __name__ == "__main__":
    main()
