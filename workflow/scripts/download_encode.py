"""Download one ENCODE file from the public portal and check its MD5."""

import argparse
import hashlib
import json
import os
import urllib.error
import urllib.request

PORTAL = "https://www.encodeproject.org"


def build_parser():
    """Build the download_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the download_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Download one ENCODE file.")
    p.add_argument("--accession", required=True, help="File accession, ENCFF...")
    p.add_argument("--extension", required=True, help="Extension, e.g. bigWig.")
    p.add_argument("-o", "--output", required=True, help="File to write.")
    return p


def main(argv=None):
    """Download the file, then fail unless its MD5 matches the portal record.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    acc = args.accession

    url = f"{PORTAL}/files/{acc}/@@download/{acc}.{args.extension}"
    md5 = hashlib.md5()
    part = f"{args.output}.part"
    try:
        with urllib.request.urlopen(f"{PORTAL}/files/{acc}/?format=json") as r:
            expected = json.load(r)["md5sum"]
        with urllib.request.urlopen(url) as r, open(part, "wb") as o:
            while chunk := r.read(1 << 20):
                md5.update(chunk)
                o.write(chunk)
    except urllib.error.HTTPError as e:
        if e.code in (401, 403):
            raise SystemExit(f"{acc} is restricted ({e.code}). Ask before using it.")
        raise
    if md5.hexdigest() != expected:
        raise SystemExit(f"{acc}: MD5 {md5.hexdigest()}, portal says {expected}.")
    os.replace(part, args.output)


if __name__ == "__main__":
    main()
