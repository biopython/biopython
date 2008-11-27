# Copyright 2004 by James Casbon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""
Code to deal with COMPASS output, a program for profile/profile comparison.

Compass is described in:

Sadreyev R, Grishin N. COMPASS: a tool for comparison of multiple protein
alignments with assessment of statistical significance. J Mol Biol. 2003 Feb
7;326(1):317-36.

Tested with COMPASS 1.24.

Classes:
Record        One result of a compass file
_Scanner      Scan compass results
_Consumer     Consume scanner events
RecordParser  Parse one compass record
Iterator      Iterate through a number of compass records
"""
from Bio import File
from Bio.ParserSupport import *
import re

class Record:
    """
    Hold information from one compass hit.
    Ali1 one is the query, Ali2 the hit.
    """

    def __init__(self):
        self.query=''
        self.hit=''
        self.gap_threshold=0
        self.query_length=0
        self.query_filtered_length=0
        self.query_nseqs=0
        self.query_neffseqs=0
        self.hit_length=0
        self.hit_filtered_length=0
        self.hit_nseqs=0
        self.hit_neffseqs=0
        self.sw_score=0
        self.evalue=-1
        self.query_start=-1
        self.hit_start=-1
        self.query_aln=''
        self.hit_aln=''
        self.positives=''


    def query_coverage(self):
        """Return the length of the query covered in alignment"""
        s = self.query_aln.replace("=", "")
        return len(s)

    def hit_coverage(self):
        """Return the length of the hit covered in the alignment"""
        s = self.hit_aln.replace("=", "")
        return len(s)

class _Scanner:
    """Reads compass output and generate events"""

    def feed(self, handle, consumer):
        """Feed in COMPASS ouput"""

        if isinstance(handle, File.UndoHandle):
            pass
        else:
            handle = File.UndoHandle(handle)
                                        

        assert isinstance(handle, File.UndoHandle), \
               "handle must be an UndoHandle"
        if handle.peekline():
            self._scan_record(handle, consumer)
                        

    def _scan_record(self,handle,consumer):
        self._scan_names(handle, consumer)
        self._scan_threshold(handle, consumer)
        self._scan_lengths(handle,consumer)
        self._scan_profilewidth(handle, consumer)
        self._scan_scores(handle,consumer)
        self._scan_alignment(handle,consumer)

    def _scan_names(self,handle,consumer):
        """
        Ali1: 60456.blo.gz.aln  Ali2: allscop//14984.blo.gz.aln
        """
        read_and_call(handle, consumer.names, contains="Ali1:")

    def _scan_threshold(self,handle, consumer):
        """
        Threshold of effective gap content in columns: 0.5
        """
        read_and_call(handle, consumer.threshold, start="Threshold")

    def _scan_lengths(self,handle, consumer):
        """
        length1=388     filtered_length1=386    length2=145     filtered_length2=137
        """
        read_and_call(handle, consumer.lengths, start="length1=")

    def _scan_profilewidth(self,handle,consumer):
        """
        Nseqs1=399      Neff1=12.972    Nseqs2=1        Neff2=6.099
        """
        read_and_call(handle, consumer.profilewidth, contains="Nseqs1")
            
    def _scan_scores(self,handle, consumer):
        """
        Smith-Waterman score = 37        Evalue = 5.75e+02
        """
        read_and_call(handle, consumer.scores, start="Smith-Waterman")

    def _scan_alignment(self,handle, consumer):
        """
        QUERY   2      LSDRLELVSASEIRKLFDIAAGMKDVISLGIGEPDFDTPQHIKEYAKEALDKGLTHYGPN
                       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        QUERY   2      LSDRLELVSASEIRKLFDIAAGMKDVISLGIGEPDFDTPQHIKEYAKEALDKGLTHYGPN

        
        QUERY          IGLLELREAIAEKLKKQNGIEADPKTEIMVLLGANQAFLMGLSAFLKDGEEVLIPTPAFV
                       ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        QUERY          IGLLELREAIAEKLKKQNGIEADPKTEIMVLLGANQAFLMGLSAFLKDGEEVLIPTPAFV
                                      
        """
        while 1:
            line = handle.readline()
            if not line:
                break
            if is_blank_line(line):
                continue
            else:
                consumer.query_alignment(line)
                read_and_call(handle, consumer.positive_alignment)
                read_and_call(handle, consumer.hit_alignment)

class _Consumer:
    # all regular expressions used -- compile only once
    _re_names = re.compile("Ali1:\s+(\S+)\s+Ali2:\s+(\S+)\s+")
    _re_threshold = \
      re.compile("Threshold of effective gap content in columns: (\S+)")
    _re_lengths = \
      re.compile("length1=(\S+)\s+filtered_length1=(\S+)\s+length2=(\S+)"
        + "\s+filtered_length2=(\S+)")
    _re_profilewidth = \
      re.compile("Nseqs1=(\S+)\s+Neff1=(\S+)\s+Nseqs2=(\S+)\s+Neff2=(\S+)")
    _re_scores = re.compile("Smith-Waterman score = (\S+)\s+Evalue = (\S+)")
    _re_start = re.compile("(\d+)")
    _re_align = re.compile("^.{15}(\S+)")
    _re_positive_alignment = re.compile("^.{15}(.+)")

    def __init__(self):
        self.data = None

    def names(self, line):
        """
        Ali1: 60456.blo.gz.aln  Ali2: allscop//14984.blo.gz.aln
              ------query-----        -------hit-------------
        """        
        self.data = Record()
        m = self.__class__._re_names.search(line)
        self.data.query = m.group(1)
        self.data.hit = m.group(2)

    def threshold(self,line):
        m = self.__class__._re_threshold.search(line)
        self.data.gap_threshold = float(m.group(1))
                                        
    def lengths(self,line):
        m = self.__class__._re_lengths.search(line)
        self.data.query_length = int(m.group(1))
        self.data.query_filtered_length = float(m.group(2))
        self.data.hit_length = int(m.group(3))
        self.data.hit_filtered_length = float(m.group(4))

    def profilewidth(self,line):
        m = self.__class__._re_profilewidth.search(line)
        self.data.query_nseqs = int(m.group(1))
        self.data.query_neffseqs = float(m.group(2))
        self.data.hit_nseqs = int(m.group(3))
        self.data.hit_neffseqs = float(m.group(4))

    def scores(self, line):
        m = self.__class__._re_scores.search(line)
        if m:
            self.data.sw_score = int(m.group(1))
            self.data.evalue = float(m.group(2))
        else:
            self.data.sw_score = 0
            self.data.evalue = -1.0
                                    
    def query_alignment(self, line):
        m = self.__class__._re_start.search(line)
        if m:
            self.data.query_start = int(m.group(1))
        m = self.__class__._re_align.match(line)
        assert m!=None, "invalid match"
        self.data.query_aln = self.data.query_aln + m.group(1)
               
    def positive_alignment(self,line):
        m = self.__class__._re_positive_alignment.match(line)
        assert m!=None, "invalid match"
        self.data.positives = self.data.positives + m.group(1)

    def hit_alignment(self,line):
        m = self.__class__._re_start.search(line)
        if m:
            self.data.hit_start = int(m.group(1))
        m = self.__class__._re_align.match(line)
        assert m!=None, "invalid match"
        self.data.hit_aln = self.data.hit_aln + m.group(1)
    
class RecordParser(AbstractParser):
        """Parses compass results into a Record object.
        """
        def __init__(self):
            self._scanner = _Scanner()
            self._consumer = _Consumer()

            
        def parse(self, handle):
            if isinstance(handle, File.UndoHandle):
                uhandle = handle
            else:
                uhandle = File.UndoHandle(handle)
            self._scanner.feed(uhandle, self._consumer)
            return self._consumer.data
                
class Iterator:
    """Iterate through a file of compass results"""
    def __init__(self, handle):
        self._uhandle = File.UndoHandle(handle)
        self._parser = RecordParser()

    def next(self):
        lines = []
        while 1:
            line = self._uhandle.readline()
            if not line:
                break
            if line[0:4] == "Ali1" and lines:
                self._uhandle.saveline(line)
                break
                              
            lines.append(line)

     
        if not lines:
            return None

        data = ''.join(lines)
        return self._parser.parse(File.StringHandle(data))

