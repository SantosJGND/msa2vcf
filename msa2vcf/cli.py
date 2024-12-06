import argparse
from .msa2vcf import go
from .object_classes import Args


def parse_args() -> Args:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--keep_n", action="store_true", required=False, help="Retain N ALTs in vcf"
    )
    parser.add_argument("msa", help="Path to MSA")
    parser.add_argument("refname", help="Name of reference sequence")
    args = parser.parse_args()

    return Args(msa=args.msa, refname=args.refname, keep_n=args.keep_n, odir=".")


def main():
    go(parse_args())
