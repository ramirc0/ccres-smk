"""Download one ENCODE file from the public S3 bucket and check its MD5.

The S3 key and MD5 come from a snapshot of the ENCODE file report, so the
download needs no portal API call.
"""

import argparse
import hashlib
import os
import urllib.request

import polars as pl

BUCKET = "s3://encode-public/"
HTTPS = "https://encode-public.s3.amazonaws.com/"


def build_parser():
    """Build the download_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the download_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Download one ENCODE file from S3.")
    p.add_argument("--accession", required=True, help="File accession, ENCFF...")
    p.add_argument(
        "--metadata", required=True, help="ENCODE file report snapshot (Parquet)."
    )
    p.add_argument("-o", "--output", required=True, help="File to write.")
    return p


def main(argv=None):
    """Download the file, then fail unless its MD5 matches the snapshot.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    record = (
        pl.scan_parquet(args.metadata)
        .filter(pl.col("accession") == args.accession)
        .select("s3_uri", "md5sum", "status", "restricted_file")
        .collect()
    )
    if record.height != 1:
        raise SystemExit(
            f"{args.accession}: {record.height} records in {args.metadata}."
        )
    s3_uri, md5sum, status, restricted = record.row(0)
    if restricted or not s3_uri or not s3_uri.startswith(BUCKET):
        raise SystemExit(
            f"{args.accession} ({status}) is not in the public bucket. Ask first."
        )

    md5 = hashlib.md5()
    part = f"{args.output}.part"
    url = HTTPS + s3_uri.removeprefix(BUCKET)
    with urllib.request.urlopen(url) as r, open(part, "wb") as o:
        while chunk := r.read(1 << 20):
            md5.update(chunk)
            o.write(chunk)
    if md5.hexdigest() != md5sum:
        raise SystemExit(f"{args.accession}: MD5 {md5.hexdigest()}, expected {md5sum}.")
    os.replace(part, args.output)


if __name__ == "__main__":
    main()
