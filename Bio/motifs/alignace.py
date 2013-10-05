# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parsing AlignACE output files
"""

from Bio.motifs import Motif, Instances
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


class Record(list):
    def __init__(self):
        self.parameters = None


def read(handle):
    """read(handle)"""
    record = Record()
    line = next(handle)
    record.version = line.strip()
    line = next(handle)
    record.command = line.strip()
    for line in handle:
        line = line.strip()
        if line=="":
            pass
        elif line[:4]=="Para":
            record.parameters={}
        elif line[0]=="#":
            seq_name = line.split("\t")[1]
            record.sequences.append(seq_name)
        elif "=" in line:
            par_name, par_value = line.split("=")
            par_name = par_name.strip()
            par_value = par_value.strip()
            record.parameters[par_name]=par_value
        elif line[:5]=="Input":
            record.sequences=[]
        elif line[:5]=="Motif":
            words = line.split()
            assert words[0]=="Motif"
            number = int(words[1])
            instances = []
        elif line[:3]=="MAP":
            alphabet = IUPAC.unambiguous_dna
            instances = Instances(instances, alphabet)
            motif = Motif(alphabet, instances)
            motif.score = float(line.split()[-1])
            motif.number = number
            motif.mask = mask
            record.append(motif)
        elif len(line.split("\t"))==4:
            seq = Seq(line.split("\t")[0], IUPAC.unambiguous_dna)
            instances.append(seq)
        elif "*" in line:
            mask = line.strip("\r\n")
        else:
            raise ValueError(line)
    return record
