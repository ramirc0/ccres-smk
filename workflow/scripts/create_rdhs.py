"""Step 2b: pick representative DHSs, accession them and keep the cDHS-supported ones.

Port of the rest of 2_Create-rDHSs.sh. The filtered DHSs are sorted with
`sort -k1,1 -k2,2n` and reduced to representative DHSs (rDHSs) by the
merge, pick and intersect loop. The sorted rDHSs get accessions, reusing
those of identical earlier rDHSs. Each rDHS then falls in one cDHS class:

- complete: it covers a whole cDHS (`intersect -F 1 -u`).
- none: it overlaps no cDHS.
- 135: it overlaps a cDHS by at least 135 bp.
- less: the rest. It is kept if its signal is at least the 90th percentile of
  the none class's signal.

The summary is the kept less rDHSs and the 135 rDHSs, each run through
`sort -u`, followed by the complete rDHSs.
"""

import argparse

import _bed
import _registry
import numpy as np
import polars as pl

COLUMNS = ["chrom", "start", "end", "name", "fdr", "signal"]
CDHS_COLUMNS = [f"c{i}" for i in range(1, 12)]
MIN_OVERLAP = 135
PERCENTILE = 90


def build_parser():
    """Build the create_rdhs.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the create_rdhs.py flags.
    """
    p = argparse.ArgumentParser(description="Create rDHSs from filtered DHSs.")
    p.add_argument("--filtered", required=True, help="Filtered DHS BED.")
    p.add_argument("--previous", required=True, help="Previous rDHSs, gzipped BED.")
    p.add_argument("--cdhs", required=True, help="Altius consensus DHS BED.")
    p.add_argument("--genome", required=True, help="hg38 or mm10.")
    p.add_argument("--prestam", required=True, help="Accessioned rDHSs to write.")
    p.add_argument("--complete", required=True, help="complete class to write.")
    p.add_argument("--none", required=True, help="none class to write.")
    p.add_argument("--summary", required=True, help="Kept rDHSs to write.")
    p.add_argument("--rdhs", required=True, help="rDHS BED to write.")
    return p


def gnu_sorted(df):
    """Sort a BED table like `sort -k1,1 -k2,2n` under en_US.UTF-8.

    Parameters
    ----------
    df : polars.DataFrame
        String columns, chrom and start among them.

    Returns
    -------
    polars.DataFrame
        Sorted rows.
    """
    return df[_bed.gnu_sort(df, text_keys=["chrom"], numeric_keys=["start"])]


def coordinates(df, codes):
    """Return the interval arrays of a table.

    Parameters
    ----------
    df : polars.DataFrame
        Table with start and end columns.
    codes : numpy.ndarray
        Chromosome code of each row.

    Returns
    -------
    tuple of numpy.ndarray
        Codes, starts and ends.
    """
    return (
        codes,
        df["start"].cast(pl.Int64).to_numpy(),
        df["end"].cast(pl.Int64).to_numpy(),
    )


def representatives(dhs):
    """Reduce sorted DHSs to representative DHSs.

    Parameters
    ----------
    dhs : polars.DataFrame
        DHSs in `sort -k1,1 -k2,2n` order.

    Returns
    -------
    polars.DataFrame
        Representative DHSs in the order the loop appends them.
    """
    picks = _registry.representative_peaks(
        *coordinates(dhs, _bed.chrom_codes(dhs["chrom"])),
        _bed.python_floats(dhs["signal"]),
        dhs["name"],
    )
    return dhs[picks]


def accessioned(rdhs, previous_path, genome):
    """Append rDHS accessions, as make-region-accession.py does.

    Parameters
    ----------
    rdhs : polars.DataFrame
        Sorted rDHSs.
    previous_path : str
        Previous rDHSs, gzipped BED of chrom, start, end, accession, version.
    genome : str
        Genome, hg38 or mm10.

    Returns
    -------
    polars.DataFrame
        The rDHSs with an accession column.
    """
    previous = _bed.read_tsv(previous_path, ["chrom", "start", "end", "acc", "version"])
    known = _registry.previous_accessions(previous.iter_rows(), "rDHS")
    regions = rdhs.select("chrom", "start", "end").iter_rows()
    return rdhs.with_columns(
        accession=pl.Series(_registry.accessions(regions, genome, "rDHS", known))
    )


def main(argv=None):
    """Write the accessioned rDHSs, their cDHS classes and the kept rDHSs.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    dhs = gnu_sorted(_bed.read_tsv(args.filtered, COLUMNS))
    zz = accessioned(gnu_sorted(representatives(dhs)), args.previous, args.genome)
    _bed.write_tsv(zz, args.prestam)

    cdhs = _bed.read_tsv(args.cdhs, CDHS_COLUMNS)
    zz_codes, cdhs_codes = _bed.shared_codes(zz["chrom"], cdhs["c1"])
    cdhs_iv = coordinates(cdhs.rename({"c2": "start", "c3": "end"}), cdhs_codes)
    complete = _bed.contains_any(coordinates(zz, zz_codes), cdhs_iv)
    _bed.write_tsv(zz.filter(complete), args.complete)

    rest, rest_codes = zz.filter(~complete), zz_codes[~complete]
    rest_iv = coordinates(rest, rest_codes)
    none = rest.filter(~_bed.overlaps_any(rest_iv, cdhs_iv))
    _bed.write_tsv(none, args.none)

    pairs, _, overlap = _bed.overlap_pairs(rest_iv, cdhs_iv)
    rest_lines = _bed.lines(rest).to_list()
    wide = np.unique(pairs[overlap >= MIN_OVERLAP])
    coverage_135 = _bed.sort_unique([rest_lines[i] for i in wide])
    names_135 = set(rest["name"].gather(wide).to_list())
    names = rest["name"].to_list()
    signal = _bed.python_floats(rest["signal"])
    # numpy prints the percentile with repr, which awk parses back exactly.
    cutoff = float(np.percentile(_bed.python_floats(none["signal"]), PERCENTILE))
    less = [
        rest_lines[i]
        for i in pairs.tolist()
        if names[i] not in names_135 and signal[i] >= cutoff
    ]
    summary = (
        _bed.sort_unique(less)
        + coverage_135
        + _bed.lines(zz.filter(complete)).to_list()
    )
    with open(args.summary, "w") as f:
        f.writelines(f"{line}\n" for line in summary)
    with open(args.rdhs, "w") as f:
        f.writelines(
            "\t".join([*fields[:3], fields[-1]]) + "\n"
            for fields in (line.split("\t") for line in summary)
        )


if __name__ == "__main__":
    main()
