"""Freeze one pipeline input: copy or derive it, then check its BLAKE3 sum.

The sources sit in another user's tree and keep changing. A frozen copy whose
checksum is not recorded in the checksum file fails the job.
"""

import argparse
import gzip
import shutil
from contextlib import ExitStack

import blake3


def build_parser():
    """Build the freeze_input.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the freeze_input.py flags.
    """
    p = argparse.ArgumentParser(description="Freeze one pipeline input.")
    p.add_argument("--sources", nargs="+", required=True, help="Source file(s).")
    p.add_argument(
        "--kind",
        default="copy",
        choices=["copy", "previous-ccres", "cell-type-tf"],
        help="How to build the frozen file from its sources.",
    )
    p.add_argument(
        "--checksums", required=True, help="b3sum-format file of expected sums."
    )
    p.add_argument(
        "--key", required=True, help="Name of the input in the checksum file."
    )
    p.add_argument("-o", "--output", required=True, help="Frozen file to write.")
    return p


def previous_ccres(src, out):
    """Write the V1-V3 rows of the previous-versions cCRE file.

    The source now also holds the V4 release, built by the run this pipeline
    reproduces. Its accessions must not seed the V4 accessioning.

    Parameters
    ----------
    src : str
        Gzipped cCRE BED. Column 6 is the version that first reported the cCRE.
    out : str
        BED to write, uncompressed.
    """
    with gzip.open(src, "rt") as f, open(out, "w") as o:
        o.writelines(line for line in f if line.rstrip("\n").split("\t")[5] != "V4")


def cell_type_tf(tf_columns, master, out):
    """Write the TF experiments of the biosamples flagged TF=yes in May 2022.

    The step 7 TF list kept growing after the release. Its rows for the
    biosamples the release flagged keep the release's TF calls reproducible.

    Parameters
    ----------
    tf_columns : list of str
        TF experiment list, as files whose lines are joined side by side with
        tabs. The last column is the biosample.
    master : str
        Master cell list of the release. Column 9 is the biosample, the last
        column `yes` or `no`.
    out : str
        Rows of the TF list whose biosample is flagged `yes`, in file order.
    """
    with open(master) as f:
        flagged = {
            fields[8]
            for fields in (line.rstrip("\n").split("\t") for line in f)
            if fields[-1] == "yes"
        }
    with ExitStack() as stack, open(out, "w") as o:
        files = [stack.enter_context(open(path)) for path in tf_columns]
        for parts in zip(*files, strict=True):
            line = "\t".join(part.rstrip("\n") for part in parts)
            if line.split("\t")[-1] in flagged:
                o.write(f"{line}\n")


def read_checksums(path):
    """Read a b3sum-format checksum file.

    Parameters
    ----------
    path : str
        Lines of `<hex digest>  <key>`.

    Returns
    -------
    dict of str to str
        Digest per key.
    """
    with open(path) as f:
        return {key: digest for digest, key in (line.split() for line in f)}


def main(argv=None):
    """Build the frozen file and fail unless its checksum is the recorded one.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)

    if args.kind == "copy":
        (src,) = args.sources
        shutil.copyfile(src, args.output)
    elif args.kind == "previous-ccres":
        (src,) = args.sources
        previous_ccres(src, args.output)
    else:
        *tf_columns, master = args.sources
        cell_type_tf(tf_columns, master, args.output)

    digest = blake3.blake3(max_threads=blake3.blake3.AUTO).update_mmap(args.output)
    digest = digest.hexdigest()
    expected = read_checksums(args.checksums).get(args.key)
    if expected is None:
        raise SystemExit(f"{args.key}: no checksum in {args.checksums}. Got {digest}.")
    if digest != expected:
        raise SystemExit(f"{args.key}: checksum {digest}, expected {expected}.")
    print(f"{args.key}: {digest} OK")


if __name__ == "__main__":
    main()
