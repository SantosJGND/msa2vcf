from dataclasses import dataclass


@dataclass
class Args:
    msa: str
    refname: str
    keep_n: bool
    odir: str
