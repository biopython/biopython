# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""

This module provides code to work with the standalone version of AlignACE, 
for motif search in DNA sequences.

AlignACE homepage:

http://atlas.med.harvard.edu/

AlignACE Citations:

Computational identification of cis-regulatory elements associated with 
groups of functionally related genes in Saccharomyces cerevisiae, 
Hughes, JD, Estep, PW, Tavazoie S, & GM Church, Journal of Molecular 
Biology 2000 Mar 10;296(5):1205-14.

Finding DNA Regulatory Motifs within Unaligned Non-Coding Sequences 
Clustered by Whole-Genome mRNA Quantitation, 
Roth, FR, Hughes, JD, Estep, PE & GM Church, Nature Biotechnology 
1998 Oct;16(10):939-45. 

functions:
AlignAce - runs the AlignACE standalone prgram and returns the 
ApplicationResult object
CompareAce - runs the AlignACE standalone prgram and returns the ApplicationResult object

Clases for pparsing AlignAce and CompareACE files: AlignAceParser,CompareAceParser

"""
#changed string.atof to float, for compatibility with python 2.6 and 3k, BW

from Bio import File
from Bio.ParserSupport import *
from Bio.Motif import Motif
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import Application
from Bio.Application import _Option,_Argument
import os

class AlignAceCommandline(Application.AbstractCommandline):
    """Create a commandline for the AlignAce program.

    XXX This could use more checking for valid paramters to the program.
    """
    def __init__(self, cmd = "AlignACE"):

        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
          [
            _Option(["-i","input","Sequence File"],["input"],lambda x : x.__class__== str,1,
                    "Input Sequence file in FASTA format."),
            
            _Option(["-numcols","numcols","number of columns to align"],["input"],lambda x : x.__class__== int,0,
                    "Number of columns to align"),

            _Option(["-expect","expect","number of sites expected in model "],["input"],lambda x : x.__class__== int,0,
                    "number of sites expected in model "),
            
            _Option(["-gcback","gcback","background fractional GC content of input sequence"],["input"],lambda x : x.__class__== float,0,
                    "background fractional GC content of input sequence"),
            
            _Option(["-minpass","minpass","minimum number of non-improved passes in phase 1"],["input"],lambda x : x.__class__== int,0,
                    "minimum number of non-improved passes in phase 1"),
            
            _Option(["-seed","seed","set seed for random number generator (time)"],["input"],lambda x : x.__class__== int,0,
                    "set seed for random number generator (time)"),
            
            _Option(["-undersample","undersample","possible sites / (expect * numcols * seedings)"],["input"],lambda x : x.__class__== int,0,
                    "possible sites / (expect * numcols * seedings)"),

            _Option(["-oversample","oversample","1/undersample"],["input"],lambda x : x.__class__== int,0,
                    "1/undersample"),
          ]

    def run(self):
        return Application.generic_run(self)



class CompareAceCommandline(Application.AbstractCommandline):
    """Create a commandline for the CompareAce program.

    XXX This could use more checking for valid paramters to the program.
    """
    def __init__(self, cmd = "CompareACE"):

        import os.path
        Application.AbstractCommandline.__init__(self)
        self.program_name = cmd

        self.parameters = \
          [
            _Argument(["motif1"],["input","file"], os.path.exists,1,"name of file containing motif 1"),
            _Argument(["motif2"],["input","file"], os.path.exists,1,"name of file containing motif 2"),
          ]

    def run(self):
        return Application.generic_run(self)


def AlignAce(infile, cmd="AlignACE", **keywds):
    """Runs AlignACE and returns data.

    cmd == AlignACE executable
    infile == sequence file to process
    
    You may pass more parameters to **keywds to change the behavior of
    the search.  Otherwise, optional values will be chosen by blastall.

    numcols    	number of columns to align (10)
    expect     	number of sites expected in model (10)
    gcback     	background fractional GC content of input sequence (0.38)
    minpass    	minimum number of non-improved passes in phase 1 (200)
    seed       	set seed for random number generator (time)
    undersample	possible sites / (expect * numcols * seedings) (1)
    oversample	        1/undersample (1)
    """

    if not os.path.exists(cmd):
        raise IOError("Executable does not exist at %s" % cmd)

    if not os.path.exists(infile):
        raise IOError("Input file does not exist at %s" % infile)
    
    AlignCmd = AlignAceCommandline(cmd)

    AlignCmd.set_parameter("input",infile)
    
    for (par,val) in keywds.iteritems():
        AlignCmd.set_parameter(par,val)

    return AlignCmd.run()


def CompareAce( cmd="CompareACE", **keywds):
    """Runs CompareACE and returns data.

    motif1, motif2 == files containing AlignACE motifs
    """

    if not os.path.exists(cmd):
        raise IOError("Executable does not exist at %s" % cmd)
    
    CompareCmd = CompareAceCommandline(cmd)

    for (par,val) in keywds.iteritems():
        CompareCmd.set_parameter(par,val)

    return CompareCmd.run()


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
        return self._consumer

class AlignAceScanner:
    """Scannner for AlignACE output

    Methods:
    feed     Feed data into the scanner.

    The scanner generates (and calls the consumer) the following types of events:

    noevent - blank line

    version - AlignACE version number
    command_line - AlignACE command line string
    parameters - the begining of the parameters
    parameter - the line containing a parameter
    sequences - the begining of the sequences list
    sequence - line containing the name of the input sequence (and a respective number)
    motif - the begining of the motif (contains the number)
    motif_hit - one hit for a motif
    motif_mask - mask of the motif (space - gap, asterisk - significant position)
    motif_score - MAP score of the motif - approx. N * log R, where R == (num. of actual occur.) / (num. of occur. expected by random.)
    
    """
    def feed(self, handle, consumer):
        """S.feed(handle, consumer)

        Feed in a AlignACE report for scanning.  handle is a file-like
        object that contains the AlignACE report.  consumer is a Consumer
        object that will receive events as the report is scanned.
        """
        consumer.version(handle.readline())
        consumer.command_line(handle.readline())
        for line in handle:
            if line.strip() == "":
                consumer.noevent(line)
            elif line[:4]=="Para":
                consumer.parameters(line)
            elif line[0]=="#":
                consumer.sequence(line)
            elif "=" in line:
                consumer.parameter(line)
            elif line[:5]=="Input":
                consumer.sequences(line)
            elif line[:5]=="Motif":
                consumer.motif(line)
            elif line[:3]=="MAP":
                consumer.motif_score(line)
            elif len(line.split("\t"))==4:
                consumer.motif_hit(line)
            elif "*" in line:
                consumer.motif_mask(line)
            else:
                raise ValueError(line)

class CompareAceScanner:
    """Scannner for CompareACE output

    Methods:
    feed     Feed data into the scanner.

    The scanner generates (and calls the consumer) the following types of events:

    motif_score - CompareACE score of motifs

    ###### TO DO #############3
    extend the scanner to include other, more complex outputs.
    """
    def feed(self, handle, consumer):
        """S.feed(handle, consumer)

        Feed in a CompareACE report for scanning.  handle is a file-like
        object that contains the CompareACE report.  consumer is a Consumer
        object that will receive events as the report is scanned.
        """
        consumer.motif_score(handle.readline())


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
