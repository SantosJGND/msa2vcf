#!/usr/bin/env python3


import os
import re
import sys
from .readfq import readfq
from operator import itemgetter
from itertools import groupby
from .object_classes import Args
from typing import List
from typing import Dict


class Variant:

    def __init__(
        self, var_type, ref_base, ref_pos, alt_base, alt_pos, alt_len, alt_freq
    ):
        self.var_type = var_type
        self.ref_base = ref_base
        self.ref_pos = ref_pos
        self.alt_base = alt_base
        self.alt_pos = alt_pos
        self.alt_len = alt_len
        self.alt_freq = alt_freq

    def __str__(self):
        return f"{self.var_type} {self.ref_base} {self.ref_pos} {self.alt_base} {self.alt_pos} {self.alt_len} {self.alt_freq}"

    def set_seqname(self, seqname):
        self.seqname = seqname

    def __eq__(self, other):
        return (
            self.var_type == other.var_type
            and self.ref_base == other.ref_base
            and self.ref_pos == other.ref_pos
            and self.alt_base == other.alt_base
            and self.alt_pos == other.alt_pos
            and self.alt_len == other.alt_len
            and self.alt_freq == other.alt_freq
        )

    def __hash__(self):
        return hash(tuple(self.__dict__.values()))


def group_del_positions(dels):
    """
    Group deletions that are adjacent to each other."""

    grouped_dels = []

    for k, g in groupby(enumerate(dels), lambda x: x[0] - x[1]):
        group = map(itemgetter(1), g)
        group = list(map(int, group))
        grouped_dels.append((group[0] - 1, len(group) + 1))

    return grouped_dels


def make_vcf(snps: List[Variant], qname, rname, keep_n):

    vcflines = []

    # header
    vcflines.append("##fileformat=VCFv4.2")
    vcflines.append("##source=" + os.path.basename(sys.argv[0]))
    vcflines.append("##contig=<ID=" + rname + ">")
    vcflines.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">')
    vcflines.append('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
    vcflines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t" + qname)

    # variants
    for line in snps:
        vcf_line = None

        alt_no_n = re.sub(r"(N|-)+", "", line.alt_base)

        alt_no_gap = re.sub(r"-+", "", line.alt_base)

        if alt_no_n:
            if not line.ref_base == alt_no_n:
                vcf_line = "\t".join(
                    [
                        rname,
                        str(line.ref_pos + 1),
                        ".",
                        line.ref_base,
                        alt_no_gap,
                        ".",
                        "PASS",
                        "DP=1;AF=" + str(line.alt_freq),
                    ]
                )

        else:
            if keep_n:
                vcf_line = "\t".join(
                    [
                        rname,
                        str(line.ref_pos + 1),
                        ".",
                        line.ref_base,
                        alt_no_gap,
                        ".",
                        "PASS",
                        "DP=1;AF=" + str(line.alt_freq),
                    ]
                )

        if vcf_line:
            vcflines.append(vcf_line)

    return vcflines


def write_vcf(vcflines, qname, odir="."):

    output_file = os.path.join(odir, qname + ".vcf")
    with open(output_file, "w", encoding="utf-8") as f:
        for line in vcflines:
            f.write(line + "\n")


#########################################################
#########################################################


def iupac_to_base(base):

    lookup_iupac = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
    }

    if lookup_iupac.get(base):
        return lookup_iupac.get(base)

    else:
        return base


def remove_terminal_gapns(seq):
    """
    Remove terminal Ns and gaps from a sequence.
    """
    return re.sub(r"(N|-)*$", "", seq)


#########################################################
#########################################################


def group_ins_positions(ins):
    """
    Group insertions that are adjacent to each other."""
    res = [(el - 1, ins.count(el) + 1) for el in ins]

    return set(res)


def fix_complex_vars(all_vars: List[Variant], rseq, qseq):
    """
    Fix complex variants that are not properly detected by update_snps."""

    positions = [i.ref_pos for i in all_vars]

    duplicate_positions = set([x for x in positions if positions.count(x) > 1])

    for pos in duplicate_positions:
        duplicate_variants = [x for x in all_vars if x.ref_pos == pos]
        for variant in duplicate_variants:
            if variant.var_type == "snp":
                all_vars.remove(variant)

            elif variant.var_type == "ins":
                variant.alt_len += 1
                variant.ref_pos -= 1
                variant.alt_pos -= 1
                variant.ref_base = rseq[variant.alt_pos]
                variant.alt_base = qseq[
                    variant.alt_pos : variant.alt_pos + variant.alt_len
                ]

            elif variant.var_type == "del":

                variant.alt_len += 1
                variant.ref_pos -= 1
                variant.alt_pos -= 1
                variant.ref_base = rseq[variant.alt_pos]
                variant.alt_base = rseq[
                    variant.alt_pos : variant.alt_pos + variant.alt_len
                ]

    return all_vars


def sequence_compare_positions(qseq, rseq):
    """
    Compare two sequences and return the positions of differences.
    returns position ints for deletions and insertions, and tuples for snps.
    """

    dels = []
    ins = []
    snps = []

    for qpos, qbase in enumerate(qseq):
        ins_aware_pos = qpos - len(ins)

        if qbase != rseq[qpos]:
            if qbase == "-":
                # deletion
                dels.append(ins_aware_pos)

            elif rseq[qpos] == "-":
                # insertion
                ins.append(ins_aware_pos)

            else:
                # snp
                snps.append((ins_aware_pos, qbase))

        elif qbase == "-" and rseq[qpos] == "-":
            # insertion but not in this sequence
            ins.append(ins_aware_pos)

    return dels, ins, snps


def update_snps(qseq, rseq) -> List[Variant]:

    dels, ins, snps = sequence_compare_positions(qseq, rseq)

    all_vars: List[Variant] = []

    if dels:
        for start, length in group_del_positions(dels):
            ins_correction = len([i for i in ins if i < start])

            var = Variant(
                "del",
                rseq[start + ins_correction : start + ins_correction + length],
                start,
                rseq[start + ins_correction],
                start + ins_correction,
                length,
                1,
            )

            all_vars.append(var)

    if snps:
        for position, base in snps:
            ins_correction = len([i for i in ins if i < position])

            var = Variant(
                "snp",
                rseq[position + ins_correction],
                position,
                base,
                position + ins_correction,
                1,
                1,
            )
            all_vars.append(var)

    if ins:
        for start, length in group_ins_positions(ins):
            ins_correction = len([i for i in ins if i < start])
            if (
                rseq[start + ins_correction : start + ins_correction + length]
                == qseq[start : start + length]
            ):
                # insertion but not in this sequence
                continue

            var = Variant(
                "ins",
                rseq[start + ins_correction],
                start,
                qseq[start : start + length],
                start + ins_correction,
                length,
                1,
            )
            all_vars.append(var)

    if len(all_vars) >= 1:
        fixed_vars = fix_complex_vars(all_vars, rseq, qseq)

        for var in fixed_vars:
            deconvolute_IUPAC(var)
            # print(var)

        fixed_vars_s = sorted(fixed_vars, key=lambda x: x.ref_pos)

        return fixed_vars_s

    else:
        return None


def deconvolute_IUPAC(var: Variant):
    """
    Deconvolute IUPAC codes into multiple variants. Keep track of the allele frequency.
    """

    num_alts = 1

    for base in var.alt_base:
        no_iupac = iupac_to_base(base)
        if isinstance(no_iupac, list):
            num_alts = len(no_iupac)
            no_iupac.remove(var.ref_base)
            var.alt_base = ",".join(no_iupac)

    var.alt_freq = round(1 / num_alts, 2)

    return var


class MsaParser:

    def __init__(self, msa_parse_args: Args):
        self.msa = msa_parse_args.msa
        self.refname = msa_parse_args.refname
        self.keep_n = msa_parse_args.keep_n
        self.odir = msa_parse_args.odir

        self.refseq = self.get_ref_seq()
        if not self.refseq:
            print(self.refname + " not found\n")
            sys.exit(1)

        self.variants: Dict[int, List[Variant]] = {}
        self.seqs = []
        self.variant_dict = {}

        self.create_output_dir()

    def create_output_dir(self):
        if not os.path.exists(self.odir):
            os.makedirs(self.odir)

    def get_ref_seq(self):
        """
        Get the reference sequence from an MSA."""

        try:
            with open(self.msa, "r", encoding="utf-8") as f:
                for name, seq, qual in readfq(f):
                    if name == self.refname:
                        return seq.upper()
        except FileNotFoundError:
            print("MSA file not found, exiting")
            sys.exit(1)

        self.seqs.append(seq)

        return None

    def variant_add(self, variant: Variant, seqname=None):

        if variant in self.variants:
            self.variants[variant].append(seqname)
        else:
            self.variants[variant] = [seqname]

        # if variant.ref_pos in self.variants:
        #    self.variants[variant.ref_pos].append(variant)
        # else:
        #    self.variants[variant.ref_pos] = [variant]

    def register_variants(self, variants: List[Variant], seqname=None):
        for variant in variants:
            variant.set_seqname(seqname)
            self.variant_add(variant)

    def parse_msa(self):

        with open(self.msa) as f:
            for name, seq, _ in readfq(f):
                if name == self.refname:
                    continue
                else:
                    self.seqs.append(seq)
                    processed_qseq = remove_terminal_gapns(seq)

                    variants = update_snps(processed_qseq, self.refseq)
                    self.register_variants(variants)

    # def compile_vcf(self):


#
#    for pos, variants in self.variants.items():
#
#        grouped_variants = group(variants)
#        print(grouped_variants)


def go(args: Args):

    msa_parser = MsaParser(args)

    refseq = msa_parser.get_ref_seq()

    if not refseq:
        print(args.refname + " not found\n")
        sys.exit(1)

    msa_parser.parse_msa()
    print(len(msa_parser.variants))

    with open(args.msa) as f:
        for name, seq, _ in readfq(f):
            if name == args.refname:
                continue
            else:
                processed_qseq = remove_terminal_gapns(seq)
                variants = update_snps(processed_qseq, refseq)

                if variants:
                    vcflines = make_vcf(variants, name, args.refname, args.keep_n)
                    write_vcf(vcflines, name, odir=args.odir)
