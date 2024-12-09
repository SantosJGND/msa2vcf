import argparse
from msa2vcf.object_classes import Args


def parse_args() -> Args:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--keep_n", action="store_true", required=False, help="Retain N ALTs in vcf"
    )
    parser.add_argument("--msa", help="Path to MSA")
    parser.add_argument("--refname", help="Name of reference sequence", required=True)
    parser.add_argument("--odir", default=".", help="Path to output directory")
    args = parser.parse_args()

    return Args(msa=args.msa, refname=args.refname, keep_n=args.keep_n, odir=args.odir)


if __name__ == "__main__":
    # main()

    import sys
    from msa2vcf.object_classes import Args
    from msa2vcf.msa2vcf_tools import MsaParser, Variant

    msa = "../data/Alignment_nt_All.fasta"
    odir = "work"
    refname = "SARS_CoV_2_Wuhan_Hu_1_MN908947"
    keep_n = False

    args = Args(msa=msa, odir=odir, refname=refname, keep_n=keep_n)

    msa_parser = MsaParser(args)
    refseq = msa_parser.get_ref_seq()

    if not refseq:
        print(args.refname + " not found\n")
        sys.exit(1)

    msa_parser.parse_msa()
    print(len(msa_parser.variants))

    variants = msa_parser.variants_clean()

    print("OIN")
    print(len(variants))

    total = 0

    for variant, list in sorted(variants.items(), key=lambda x: x[0].ref_pos):
        # print(variant, len(list))
        total += len(list)

    msa_parser.write_vcf()

    print("Total variants: ", total)
