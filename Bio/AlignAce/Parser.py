# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Classes for parsing AlignAce and CompareACE files
"""

from Bio.ParserSupport import *
from Scanner import AlignAceScanner,CompareAceScanner
from Motif import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


class AlignAceConsumer:
    """
    The general purpose consumer for the AlignAceScanner.

    Should be passed as the consumer to the feed method of the AlignAceScanner. After 'consuming' the file, it has the list of motifs in the motifs property.
    """
    def __init__(self):
        self.motifs=[]
        self.current_motif=None
        self.param_dict = None
    
    def parameters(self,line):
        self.param_dict={}

    def parameter(self,line):
        par_name = line.split("=")[0].strip()
        par_value = line.split("=")[1].strip()
        self.param_dict[par_name]=par_value
        
    def sequences(self,line):
        self.seq_dict=[]
        
    def sequence(self,line):
        seq_name = line.split("\t")[1]
        self.seq_dict.append(seq_name)
        
    def motif(self,line):
        self.current_motif = Motif()
        self.motifs.append(self.current_motif)
        self.current_motif.alphabet=IUPAC.unambiguous_dna
        
    def motif_hit(self,line):
        seq = Seq(line.split("\t")[0],IUPAC.unambiguous_dna)
        self.current_motif.add_instance(seq)
        
    def motif_score(self,line):
        self.current_motif.score = float(line.split()[-1])
        
    def motif_mask(self,line):
        self.current_motif.set_mask(line.strip("\n\c"))

    def noevent(self,line):
        pass
        
    def version(self,line):
        self.ver = line
        
    def command_line(self,line):
        self.cmd_line = line
    
class AlignAceParser(AbstractParser):
    """Parses AlignAce data into a sequence of Motifs.
    """
    def __init__(self):
        """__init__(self)"""
        self._scanner = AlignAceScanner()
        self._consumer = AlignAceConsumer()

    def parse(self, handle):
        """parse(self, handle)"""
        self._scanner.feed(handle, self._consumer)
        return self._consumer.motifs


class CompareAceConsumer:
    """
    The general purpose consumer for the CompareAceScanner.

    Should be passed as the consumer to the feed method of the CompareAceScanner. After 'consuming' the file, it has the list of motifs in the motifs property.
    """
    def __init__(self):
        pass
    def motif_score(self,line):
        self.data = float(line.split()[-1])
    
class CompareAceParser(AbstractParser):
    """Parses CompareAce output to usable form

    ### so far only in a very limited way
    """
    def __init__(self):
        """__init__(self)"""
        self._scanner = CompareAceScanner()
        self._consumer = CompareAceConsumer()

    def parse(self, handle):
        """parse(self, handle)"""
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data
