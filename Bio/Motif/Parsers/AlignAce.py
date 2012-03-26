# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parsing AlignACE files: AlignAceParser
"""

from Bio.Motif import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


class Record(object):
    def __init__(self):
        self.motifs=[]
        self.current_motif=None
        self.param_dict = None


def read(handle):
    """read(handle)"""
    record = Record()
    record.ver = handle.next()
    record.cmd_line = handle.next()
    for line in handle:
        if line.strip() == "":
            pass
        elif line[:4]=="Para":
            record.param_dict={}
        elif line[0]=="#":
            seq_name = line.split("\t")[1]
            record.seq_dict.append(seq_name)
        elif "=" in line:
            par_name = line.split("=")[0].strip()
            par_value = line.split("=")[1].strip()
            record.param_dict[par_name]=par_value
        elif line[:5]=="Input":
            record.seq_dict=[]
        elif line[:5]=="Motif":
            record.current_motif = Motif()
            record.motifs.append(record.current_motif)
            record.current_motif.alphabet=IUPAC.unambiguous_dna
        elif line[:3]=="MAP":
            record.current_motif.score = float(line.split()[-1])
        elif len(line.split("\t"))==4:
            seq = Seq(line.split("\t")[0],IUPAC.unambiguous_dna)
            record.current_motif.add_instance(seq)
        elif "*" in line:
            record.current_motif.set_mask(line.strip("\n\c"))
        else:
            raise ValueError(line)
    return record


