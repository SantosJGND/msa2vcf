# msa2vcf

Turn a fasta-format MSA into a single VCF.

This script is a development on the original msa2vcf.py script by the Connor lab, [here](github.com/connor-lab/msa2vcf.git). WIP.

## Notes

VCF position is based on reference position.

IUPAC bases:

- Ref base is stripped from list of potentials
- Variant AF is given as 1/IUPAC redundancy
