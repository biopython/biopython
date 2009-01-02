# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Patched by Brad Chapman.
# Chris Wroe added modifications for work in myGrid

"""
This module provides code to work with the WWW version of BLAST
provided by the NCBI.
http://blast.ncbi.nlm.nih.gov/

Functions:
qblast        Do a BLAST search using the QBLAST API.

Deprecated classes:
BlastParser   Parses output from WWW blast.
_Scanner      Scans output from NCBI's BLAST WWW server.


"""
import re

try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

from Bio.ParserSupport import *

class BlastParser(AbstractParser):
    """Parses WWW BLAST data into a Record.Blast object (DEPRECATED).

    This is a parser for the NCBI's HTML (web page) BLAST output.
    """
    def __init__(self):
        """Create a BlastParser object (DEPRECATED)."""
        import warnings      
        warnings.warn("Bio.Blast.NCBIWWW.BlastParser is deprecated." \
                      + " We recommend you use the XML output with" \
                      + " the parser in Bio.Blast.NCBIXML instead.",
                      DeprecationWarning)
                       
        import NCBIStandalone
        self._scanner = _Scanner()
        self._consumer = SGMLStrippingConsumer(NCBIStandalone._BlastConsumer())

    def parse(self, handle):
        """parse(self, handle)"""
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data
    
class _Scanner:
    """Scanner for the HTML BLAST parser (PRIVATE, DEPRECATED).
    
    Scan BLAST output from NCBI's web server at:
    http://www.ncbi.nlm.nih.gov/BLAST/
    
    Tested with BLAST v2.0.10

    Methods:
    feed     Feed data into the scanner.
    """
    def feed(self, handle, consumer):
        """S.feed(handle, consumer)

        Feed in a BLAST report for scanning.  handle is a file-like
        object that contains the BLAST report.  consumer is a Consumer
        object that will receive events as the report is scanned.

        """
        from Bio import File

        # This stuff appears in 2.0.12.
        # <p><!--
        # QBlastInfoBegin
        #         Status=READY
        # QBlastInfoEnd
        # --><p>

        # <HTML>
        # <HEAD>
        # <TITLE>BLAST Search Results </TITLE>
        # </HEAD>
        # <BODY BGCOLOR="#FFFFFF" LINK="#0000FF" VLINK="#660099" ALINK="#660099
        # <A HREF="http://www.ncbi.nlm.nih.gov/BLAST/blast_form.map"> <IMG SRC=
        # <BR><BR><PRE>

        # BLAST Formatted information
        
        # 
        # </BODY>
        # </HTML>
        # </BODY>
        # </HTML>
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)
        # Read HTML formatting up to the "BLAST" version line.
        read_and_call_until(uhandle, consumer.noevent,
                            has_re=re.compile(r'<b>.?BLAST'))

        self._scan_header(uhandle, consumer)
        self._scan_rounds(uhandle, consumer)
        self._scan_database_report(uhandle, consumer)
        self._scan_parameters(uhandle, consumer)

        # Read HTML footer information.
        while uhandle.peekline():
            read_and_call(uhandle, consumer.noevent)

    def _scan_header(self, uhandle, consumer):
        # <b>BLASTP 2.0.10 [Aug-26-1999]</b>
        # 
        # 
        # <b><a href="http://www.ncbi.nlm.nih.gov/htbin-
        # post/Entrez/query?uid=9254694&form=6&db=m&Dopt=r">Reference</a>:</b>
        # Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch&auml;ffer, 
        # Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), 
        # "Gapped BLAST and PSI-BLAST: a new generation of protein database sea
        # programs",  Nucleic Acids Res. 25:3389-3402.
        # <p>
        # <b>Query=</b> gi|120291|sp|P21297|FLBT_CAUCR FLBT PROTEIN.
        #          (141 letters)
        # 
        # <b>Database:</b> Non-redundant SwissProt sequences
        #            82,258 sequences; 29,652,561 total letters
        # 
        # <p> <p>If you have any problems or questions with the results of this

        # If there are hits, and Graphical Overview was selected:
        # <FORM NAME="BLASTFORM">
        # </PRE>
        # <CENTER>
        # <H3><a href="/BLAST/newoptions.html#graphical-overview"> Distribution
        # <input name=defline size=80 value="Mouse-over to show defline and sco
        # </CENTER>
        # <map name=img_map>
        # <area shape=rect coords=69,101,476,106 href="#120291" ONMOUSEOVER='do
        # <area shape=rect coords=156,108,305,113 href="#3024946" ONMOUSEOVER='
        # </map>
        # <CENTER>
        # <IMG WIDTH=529 HEIGHT=115 USEMAP=#img_map BORDER=1 SRC="nph-getgif.cg
        # <HR>
        # <PRE>  XXX
        consumer.start_header()

        # Read the "BLAST" version line and the following blanks.
        read_and_call(uhandle, consumer.version, contains='BLAST')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        # Read the reference lines and the '<p>' line.
        # TBLASTN 2.2.6 has a blank line instead of a "<p>".
        while 1:
            line = uhandle.readline()
            if line[:3] == '<p>' or not line.strip():
                consumer.noevent(line)
                break
            consumer.reference(line)

        # Read the RID line, for version 2.0.12 (2.0.11?) and above.
        attempt_read_and_call(uhandle, consumer.noevent, start='RID')
        # Brad Chapman noticed a '<p>' line in BLASTN 2.1.1; this line
        # seems to have disappeared again.
        # attempt_read_and_call(uhandle, consumer.noevent, start='<p>')
        attempt_read_and_call(uhandle, consumer.noevent)

        # Apparently, there's some discrepancy between whether the
        # Query or database comes first.  Usually the Query does, but
        # Brad noticed a case where the database came first.
        if uhandle.peekline().find("Query=") >= 0:
            self._scan_query_info(uhandle, consumer)
            self._scan_database_info(uhandle, consumer)
        else:
            self._scan_database_info(uhandle, consumer)
            self._scan_query_info(uhandle, consumer)
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        consumer.end_header()

    def _scan_blastform(self, uhandle, consumer):
        if attempt_read_and_call(uhandle, consumer.noevent,
                                 contains="BLASTFORM"):
            while 1:
                line = uhandle.peekline()
                if is_blank_line(line):
                    break
                elif "Query=" in line:
                    break
                consumer.noevent(uhandle.readline())

    def _scan_database_info(self, uhandle, consumer):
        attempt_read_and_call(uhandle, consumer.noevent, start='<p>')
        read_and_call(uhandle, consumer.database_info, contains='Database')
        # Sagar Damle reported that databases can consist of multiple lines.
        # But, trickily enough, sometimes the second line can also have the
        # word sequences in it. Try to use 'sequences;' (with a semicolon)
        read_and_call_until(uhandle, consumer.database_info,
                            contains='sequences;')
        read_and_call(uhandle, consumer.database_info, contains='sequences;')
        read_and_call(uhandle, consumer.noevent, blank=1)
        attempt_read_and_call(uhandle, consumer.noevent,
                              contains='problems or questions')
        self._scan_blastform(uhandle, consumer)
        
        attempt_read_and_call(uhandle, consumer.noevent, blank=1)
        if attempt_read_and_call(uhandle, consumer.noevent,
                                 start="<table border=0 width=600"):
            read_and_call_until(uhandle, consumer.noevent,
                            contains="</table>")
            consumer.noevent(uhandle.readline())
            read_and_call(uhandle, consumer.noevent, blank=1)

        attempt_read_and_call(uhandle, consumer.noevent, start="<p>")
        
        if attempt_read_and_call(uhandle, consumer.noevent,
                                 contains="Taxonomy reports"):
            read_and_call(uhandle, consumer.noevent, start="<BR>")
        attempt_read_and_call(uhandle, consumer.noevent, start="<PRE>")

        # </PRE>
        # <!-- Progress msg from the server 500 7-->
        # <!-- Progress msg from the server 1000 15-->
        # <!-- Progress msg from the server 1500 21-->
        # ...
        # <PRE><HR><BR><b>Query=</b> test
        #          (60 letters)
        if attempt_read_and_call(uhandle, consumer.noevent, start="</PRE>"):
            read_and_call_until(uhandle, consumer.noevent, start="<PRE>")
            while 1:
                line = uhandle.peekline()
                if not line[:5] == "<PRE>" or line.find("Query=") >= 0:
                    break
                read_and_call(uhandle, consumer.noevent, start="<PRE>")
            
        read_and_call_while(uhandle, consumer.noevent, blank=1)

    def _scan_query_info(self, uhandle, consumer):
        # Read the Query lines and the following blank line.
        read_and_call(uhandle, consumer.query_info, contains='Query=')
        read_and_call_until(uhandle, consumer.query_info, blank=1)
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        if attempt_read_and_call(uhandle, consumer.noevent, start="<PRE>"):
            read_and_call_while(uhandle, consumer.noevent, blank=1)
        self._scan_blastform(uhandle, consumer)
        
    def _scan_rounds(self, uhandle, consumer):
        self._scan_descriptions(uhandle, consumer)
        self._scan_alignments(uhandle, consumer)

    def _scan_descriptions(self, uhandle, consumer):
        consumer.start_descriptions()

        # Three things can happen here:
        # 1.  line contains 'Score     E'
        # 2.  line contains "No significant similarity"
        # 3.  no descriptions
        if not attempt_read_and_call(
            uhandle, consumer.description_header,
            has_re=re.compile(r"Score {4,5}E")):
            # Either case 2 or 3.  Look for "No hits found".
            attempt_read_and_call(uhandle, consumer.no_hits,
                                  contains='No significant similarity')
            read_and_call_while(uhandle, consumer.noevent, blank=1)
            consumer.end_descriptions()
            # Stop processing.
            return
        # Sequences producing significant alignments:                      
        # 
        # <a href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Ret
        # <a href="http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Ret
        # 
        # Read the score header lines and a blank line.
        read_and_call(uhandle, consumer.description_header,
                      start='Sequences producing')
        read_and_call(uhandle, consumer.noevent, blank=1)

        # Read the descriptions
        # The description contains at least an <a href> into the alignments.
        # What is no alignments are chosen?
        read_and_call_while(uhandle, consumer.description,
                            blank=0, contains='<a')

        # two choices here, either blank lines or a </PRE>
        if not attempt_read_and_call(uhandle, consumer.noevent,
                                     contains='</PRE>'):
            read_and_call_while(uhandle, consumer.noevent, blank=1)

        consumer.end_descriptions()

    def _scan_alignments(self, uhandle, consumer):
        # Check to see whether I'm at an alignment or database report.
        # Possibilities:
        # 1) BLASTP 2.0.14, pairwise alignment
        #   <CENTER><b><FONT color="green">Alignments</FONT></b></CENTER>
        #   <PRE>
        #   ><a name = 121837></a><a href="http://www.ncbi.nlm.nih.gov:80/entre
        # 2) BLASTP 2.0.10, pairwise alignment
        #   <PRE>
        #   <a name = 120291> </a><a href="http://www.ncbi.nlm.nih.gov:80/entre
        # 3) BLASTP 2.0.10, master-slave
        #   <PRE>
        #   blast_tmp 1    MFQQIGAVQAKSGTDEPAHPCEKFPPERKCEAVFWKPLPRHEAREILLAARK
        # 4) BLASTP 2.0.10, 2.0.14, database
        #   <PRE>
        #     Database: Non-redundant SwissProt sequences
        # 5) BLASTX 2.2.4, pairwise alignment
        #   <CENTER><b><FONT color="green">Alignments</FONT></b></CENTER>
        #   </form>
        #   <script src="blastResult.js"></script><table border="0"><tr><td><FO
        #   <PRE>
        # 6) Qblast 2.2.10, database (no 'Database' line)
        #  <PRE>
        #  Lambda     K      H

        # Get the first two lines and examine them.
        line1 = safe_readline(uhandle)
        line2 = safe_readline(uhandle)
        uhandle.saveline(line2)
        uhandle.saveline(line1)

        is_pairwise = is_masterslave = 0
        if 'Alignments' in line2:
            is_pairwise = 1
        elif line2.startswith('  Database'):
            pass
        elif line2.startswith('Lambda     K      H'):
            pass
        elif line2.startswith('blast_tmp'):
            is_masterslave = 1
        elif line1.startswith('<PRE>'):
            is_pairwise = 1
        else:
            raise ValueError("Cannot resolve location at lines:\n%s\n%s" \
                             % (line1, line2))

        if is_pairwise:
            self._scan_pairwise_alignments(uhandle, consumer)
        elif is_masterslave:
            self._scan_masterslave_alignment(uhandle, consumer)

    def _scan_pairwise_alignments(self, uhandle, consumer):
        while 1:
            read_and_call_until(uhandle, consumer.noevent, start='<PRE>')
        
            # The first line is <PRE>.  Check the second line to see if
            # I'm still at an alignment.
            line1 = safe_readline(uhandle)
            line2 = safe_readline(uhandle)
            uhandle.saveline(line2)
            uhandle.saveline(line1)
            # Lambda is for Q-blast results, which do not have a Database line
            if line1.find('Database') >= 0 or line2.find("Database") >= 0 \
                or line2.find('Lambda     K      H') >= 0:
                break

            # Occasionally, there's a bug where the alignment_header and
            # hsp_header are skipped, leaving only the hsp_alignment.
            # Detect this and handle it accordingly.
            if line2[:6] == 'Query:':
                self._scan_abbreviated_pairwise_alignment(uhandle, consumer)
            else:
                self._scan_one_pairwise_alignment(uhandle, consumer)

    def _scan_abbreviated_pairwise_alignment(self, uhandle, consumer):
        # Sometimes all header information is skipped, leaving
        # only the raw alignments.  I believe this is a bug because
        # without the header information, you lose vital information such
        # as score, target sequence id, etc.
        # Format:
        # <PRE>
        # hsp_alignment

        consumer.start_alignment()
        consumer.start_hsp()
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        self._scan_hsp_alignment(uhandle, consumer)
        consumer.end_hsp()
        consumer.end_alignment()
        
    def _scan_one_pairwise_alignment(self, uhandle, consumer):
        # Alignment format:
        # <CENTER><b><FONT color="green">Alignments</FONT></b></CENTER>
        #       (BLAST 2.0.14)
        # <PRE>
        # alignment_header
        #   hsp_header
        #   hsp_alignment
        #   [...]
        # The hsp_header and hsp_alignment blocks can be repeated.

        consumer.start_alignment()
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        self._scan_alignment_header(uhandle, consumer)

        # Scan a bunch of score/alignment's.
        while 1:
            # An HSP header starts with ' Score'.
            # However, if the HSP header is not the first one in the
            # alignment, there will be a '<PRE>' line first.  Therefore,
            # I will need to check either of the first two lines to
            # see if I'm at an HSP header.
            line1 = safe_readline(uhandle)
            line2 = safe_readline(uhandle)
            line3 = safe_readline(uhandle)
            uhandle.saveline(line3)
            uhandle.saveline(line2)
            uhandle.saveline(line1)
            # There can be <a> links in front of 'Score'
            rea = re.compile(r"</?a[^>]*>")
            line1 = rea.sub("", line1)
            line2 = rea.sub("", line2)
            line3 = rea.sub("", line3)
            if line1[:6] != ' Score' and line2[:6] != ' Score' and \
               line3[:6] != ' Score':
                break
            self._scan_hsp(uhandle, consumer)
                
        consumer.end_alignment()

    def _scan_alignment_header(self, uhandle, consumer):
        # <a name = 120291> </a><a href="http://www.ncbi.nlm.nih.gov:80/entrez/
        #            Length = 141
        #
        while 1:
            line = safe_readline(uhandle)
            if line.lstrip().startswith('Length ='):
                consumer.length(line)
                break
            elif is_blank_line(line):
                # Check to make sure I haven't missed the Length line
                raise ValueError("I missed the Length in an alignment header")
            consumer.title(line)

        if not attempt_read_and_call(uhandle, consumer.noevent,
                                     start='          '):
            read_and_call(uhandle, consumer.noevent, blank=1)

    def _scan_hsp(self, uhandle, consumer):
        consumer.start_hsp()
        self._scan_hsp_header(uhandle, consumer)
        self._scan_hsp_alignment(uhandle, consumer)
        consumer.end_hsp()

    def _scan_hsp_header(self, uhandle, consumer):
        # If the HSP is not the first one within an alignment, includes:
        # <PRE>
        
        #  Score = 22.7 bits (47), Expect = 2.5
        #  Identities = 10/36 (27%), Positives = 18/36 (49%)
        #  Strand = Plus / Plus
        #  Frame = +3
        #

        attempt_read_and_call(uhandle, consumer.noevent, start='<PRE>')
        attempt_read_and_call(uhandle, consumer.noevent, blank=1)
        read_and_call(uhandle, consumer.score,
                      has_re=re.compile(r'^ (<a[^>]*></a>)*Score'))
        read_and_call(uhandle, consumer.identities, start=' Identities')
        # BLASTN
        attempt_read_and_call(uhandle, consumer.strand, start = ' Strand')
        # BLASTX, TBLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.frame, start = ' Frame')
        read_and_call(uhandle, consumer.noevent, blank=1)

    def _scan_hsp_alignment(self, uhandle, consumer):
        # Query: 11 GRGVSACA-------TCDGFFYRNQKVAVIGGGNTAVEEALYLSNIASEVHLIHRRDGF
        #           GRGVS+         TC    Y  + + V GGG+ + EE   L     +   I R+
        # Sbjct: 12 GRGVSSVVRRCIHKPTCKE--YAVKIIDVTGGGSFSAEEVQELREATLKEVDILRKVSG
        # 
        # Query: 64 AEKILIKR 71
        #              I +K 
        # Sbjct: 70 PNIIQLKD 77
        # </PRE>
        #
        # 

        while 1:
            # Blastn adds an extra line filled with spaces before Query
            attempt_read_and_call(uhandle, consumer.noevent, start='     ')
            read_and_call(uhandle, consumer.query, start='Query')
            read_and_call(uhandle, consumer.align, start='     ')
            read_and_call(uhandle, consumer.sbjct, start='Sbjct')
            if not attempt_read_and_call(uhandle, consumer.noevent, blank=1):
                break
        read_and_call(uhandle, consumer.noevent, start='</PRE>')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

    def _scan_masterslave_alignment(self, uhandle, consumer):
        consumer.start_alignment()
        read_and_call(uhandle, consumer.noevent, start='<PRE>')
        while 1:
            line = safe_readline(uhandle)
            if is_blank_line(line):
                consumer.noevent(line)
            elif line[:6] == '</PRE>':
                consumer.noevent(line)
                break
            else:
                consumer.multalign(line)
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        consumer.end_alignment()

    def _scan_database_report(self, uhandle, consumer):
        # <PRE>
        #   Database: Non-redundant SwissProt sequences
        #     Posted date:  Dec 18, 1999  8:26 PM
        #   Number of letters in database: 29,652,561
        #   Number of sequences in database:  82,258
        #   
        # Lambda     K      H
        #    0.317    0.133    0.395 
        # 
        # Gapped
        # Lambda     K      H
        #    0.270   0.0470    0.230 
        # 

        # qblast (BLASTN 2.2.10) does not give the Database: bits before the Lambda
        # information, so that needs to be skipped

        consumer.start_database_report()

        # TBALSTN 2.2.6
        # <PRE>  Database: /tmp/affyA.fasta
        line = uhandle.peekline()
        # only look for database information if we aren't already at the
        # Lambda bits
        if line.find("Database") < 0:
            read_and_call(uhandle, consumer.noevent, start='<PRE>')
        line2 = uhandle.peekline()
        if line2.find("Lambda     K      H") < 0:
            read_and_call(uhandle, consumer.database, contains='  Database')
            read_and_call_until(uhandle, consumer.database, contains="Posted")
            read_and_call(uhandle, consumer.posted_date, start='    Posted')
            read_and_call(uhandle, consumer.num_letters_in_database,
                          start='  Number of letters')
            read_and_call(uhandle, consumer.num_sequences_in_database,
                          start='  Number of sequences')
            read_and_call(uhandle, consumer.noevent, start='  ')

        read_and_call(uhandle, consumer.noevent, start='Lambda')
        read_and_call(uhandle, consumer.ka_params)
        read_and_call(uhandle, consumer.noevent, blank=1)

        # not BLASTP
        attempt_read_and_call(uhandle, consumer.gapped, start='Gapped')
        # not TBLASTX
        if attempt_read_and_call(uhandle, consumer.noevent, start='Lambda'):
            read_and_call(uhandle, consumer.ka_params_gap)
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        consumer.end_database_report()

    def _scan_parameters(self, uhandle, consumer):
        # Matrix: BLOSUM62
        # Number of Hits to DB: 1st pass: 41542626, 2nd pass: 9765
        # Number of Sequences: 1st pass: 89405, 2nd pass: 84
        # Number of extensions: 1st pass: 500847, 2nd pass: 6747
        # Number of successful extensions: 1st pass: 14, 2nd pass: 49
        # Number of sequences better than 10.0: 20
        # length of query: 205
        # length of database: 10,955,950
        # effective HSP length: 46
        # effective length of query: 158
        # effective length of database: 6,843,320
        # effective search space: 1081244560
        # effective search space used: 1081244560
        # frameshift window, decay const: 50,  0.5
        # T: 13
        # A: 40
        # X1: 16 ( 7.3 bits)
        # X2: 0 ( 0.0 bits)
        # S1: 41 (21.7 bits)
        # S2: 52 (26.7 bits)
        # 
        # </PRE>
        
        # 6/3/2001, </PRE> is gone, replaced by </form>
        

        consumer.start_parameters()

        # qblast doesn't have Matrix line
        attempt_read_and_call(uhandle, consumer.matrix, start='Matrix')
        # not TBLASTX
        attempt_read_and_call(uhandle, consumer.gap_penalties, start='Gap')
        
        # in qblast the Number of Hits and Number of Sequences lines are
        # reversed
        if attempt_read_and_call(uhandle, consumer.num_hits,
                start='Number of Hits'):
            read_and_call(uhandle, consumer.num_sequences,
                          start='Number of Sequences')
        else:
            read_and_call(uhandle, consumer.num_sequences,
                          start='Number of Sequences')
            read_and_call(uhandle, consumer.num_hits,
                          start='Number of Hits')

        read_and_call(uhandle, consumer.num_extends,
                      start='Number of extensions')
        read_and_call(uhandle, consumer.num_good_extends,
                      start='Number of successful')

        read_and_call(uhandle, consumer.num_seqs_better_e,
                      start='Number of sequences')

        # not BLASTN, TBLASTX
        if attempt_read_and_call(uhandle, consumer.hsps_no_gap,
                                 start="Number of HSP's better"):
            # for qblast order of HSP info is changed
            if attempt_read_and_call(uhandle, consumer.hsps_prelim_gapped,
                    start="Number of HSP's successfully"):
                read_and_call(uhandle, consumer.hsps_prelim_gap_attempted,
                              start="Number of HSP's that")
                read_and_call(uhandle, consumer.hsps_gapped,
                              start="Number of HSP's gapped")
            else:
                read_and_call(uhandle, consumer.no_event,
                              start="Number of HSP's gapped")
                read_and_call(uhandle, consumer.no_event,
                              start="Number of HSP's successfully")
                read_and_call(uhandle, consumer.no_event,
                              start="Number of extra gapped")
        
        # QBlast has different capitalization on the Length info:
        if attempt_read_and_call(uhandle, consumer.query_length,
                start='Length of query'):
            read_and_call(uhandle, consumer.database_length,
                start='Length of database')
            read_and_call(uhandle, consumer.no_event,
                          start='Length adjustment')
            attempt_read_and_call(uhandle, consumer.effective_query_length,
                                  start='Effective length of query')
            read_and_call(uhandle, consumer.effective_database_length,
                          start='Effective length of database')
            attempt_read_and_call(uhandle, consumer.effective_search_space,
                                  start='Effective search space:')
            attempt_read_and_call(uhandle, consumer.effective_search_space_used,
                                  start='Effective search space used')

        else:
            attempt_read_and_call(uhandle, consumer.query_length,
                                  start='length of query')
            read_and_call(uhandle, consumer.database_length,
                start='length of database')
            read_and_call(uhandle, consumer.effective_hsp_length,
                          start='effective HSP')
            attempt_read_and_call(uhandle, consumer.effective_query_length,
                                  start='effective length of query')
            read_and_call(uhandle, consumer.effective_database_length,
                          start='effective length of database')
            attempt_read_and_call(uhandle, consumer.effective_search_space,
                                  start='effective search space:')
            attempt_read_and_call(uhandle, consumer.effective_search_space_used,
                                  start='effective search space used')

        # BLASTX, TBLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.frameshift, start='frameshift')
        attempt_read_and_call(uhandle, consumer.threshold, start='T')
        read_and_call(uhandle, consumer.window_size, start='A')
        read_and_call(uhandle, consumer.dropoff_1st_pass, start='X1')
        read_and_call(uhandle, consumer.gap_x_dropoff, start='X2')
        # not BLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.gap_x_dropoff_final,
                              start='X3')
        read_and_call(uhandle, consumer.gap_trigger, start='S1')
        attempt_read_and_call(uhandle, consumer.blast_cutoff, start='S2')

        attempt_read_and_call(uhandle, consumer.noevent, blank=1)
        attempt_read_and_call(uhandle, consumer.noevent, start="</PRE>")
        attempt_read_and_call(uhandle, consumer.noevent, start="</form>")

        consumer.end_parameters()

def qblast(program, database, sequence,
           auto_format=None,composition_based_statistics=None,
           db_genetic_code=None,endpoints=None,entrez_query='(none)',
           expect=10.0,filter=None,gapcosts=None,genetic_code=None,
           hitlist_size=50,i_thresh=None,layout=None,lcase_mask=None,
           matrix_name=None,nucl_penalty=None,nucl_reward=None,
           other_advanced=None,perc_ident=None,phi_pattern=None,
           query_file=None,query_believe_defline=None,query_from=None,
           query_to=None,searchsp_eff=None,service=None,threshold=None,
           ungapped_alignment=None,word_size=None,
           alignments=500,alignment_view=None,descriptions=500,
           entrez_links_new_window=None,expect_low=None,expect_high=None,
           format_entrez_query=None,format_object=None,format_type='XML',
           ncbi_gi=None,results_file=None,show_overview=None
           ):
    """Do a BLAST search using the QBLAST server at NCBI.

    Supports all parameters of the qblast API for Put and Get.
    Some useful parameters:
    program        blastn, blastp, blastx, tblastn, or tblastx (lower case)
    database       Which database to search against (e.g. "nr").
    sequence       The sequence to search.
    ncbi_gi        TRUE/FALSE whether to give 'gi' identifier.
    descriptions   Number of descriptions to show.  Def 500.
    alignments     Number of alignments to show.  Def 500.
    expect         An expect value cutoff.  Def 10.0.
    matrix_name    Specify an alt. matrix (PAM30, PAM70, BLOSUM80, BLOSUM45).
    filter         "none" turns off filtering.  Default no filtering
    format_type    "HTML", "Text", "ASN.1", or "XML".  Def. "XML".
    entrez_query   Entrez query to limit Blast search
    hitlist_size   Number of hits to return. Default 50

    This function does no checking of the validity of the parameters
    and passes the values to the server as is.  More help is available at:
    http://www.ncbi.nlm.nih.gov/BLAST/blast_overview.html

    """
    import urllib, urllib2
    import time

    assert program in ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

    # Format the "Put" command, which sends search requests to qblast.
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node5.html on 9 July 2007
    parameters = [
        ('AUTO_FORMAT',auto_format),
        ('COMPOSITION_BASED_STATISTICS',composition_based_statistics),
        ('DATABASE',database),
        ('DB_GENETIC_CODE',db_genetic_code),
        ('ENDPOINTS',endpoints),
        ('ENTREZ_QUERY',entrez_query),
        ('EXPECT',expect),
        ('FILTER',filter),
        ('GAPCOSTS',gapcosts),
        ('GENETIC_CODE',genetic_code),
        ('HITLIST_SIZE',hitlist_size),
        ('I_THRESH',i_thresh),
        ('LAYOUT',layout),
        ('LCASE_MASK',lcase_mask),
        ('MATRIX_NAME',matrix_name),
        ('NUCL_PENALTY',nucl_penalty),
        ('NUCL_REWARD',nucl_reward),
        ('OTHER_ADVANCED',other_advanced),
        ('PERC_IDENT',perc_ident),
        ('PHI_PATTERN',phi_pattern),
        ('PROGRAM',program),
        ('QUERY',sequence),
        ('QUERY_FILE',query_file),
        ('QUERY_BELIEVE_DEFLINE',query_believe_defline),
        ('QUERY_FROM',query_from),
        ('QUERY_TO',query_to),
        ('SEARCHSP_EFF',searchsp_eff),
        ('SERVICE',service),
        ('THRESHOLD',threshold),
        ('UNGAPPED_ALIGNMENT',ungapped_alignment),
        ('WORD_SIZE',word_size),
        ('CMD', 'Put'),
        ]
    query = [x for x in parameters if x[1] is not None]
    message = urllib.urlencode(query)

    # Send off the initial query to qblast.
    # Note the NCBI do not currently impose a rate limit here, other
    # than the request not to make say 50 queries at once using multiple
    # threads.
    request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi",
                              message,
                              {"User-Agent":"BiopythonClient"})
    handle = urllib2.urlopen(request)

    # Format the "Get" command, which gets the formatted results from qblast
    # Parameters taken from http://www.ncbi.nlm.nih.gov/BLAST/Doc/node6.html on 9 July 2007    
    rid, rtoe = _parse_qblast_ref_page(handle)
    parameters = [
        ('ALIGNMENTS',alignments),
        ('ALIGNMENT_VIEW',alignment_view),
        ('DESCRIPTIONS',descriptions),
        ('ENTREZ_LINKS_NEW_WINDOW',entrez_links_new_window),
        ('EXPECT_LOW',expect_low),
        ('EXPECT_HIGH',expect_high),
        ('FORMAT_ENTREZ_QUERY',format_entrez_query),
        ('FORMAT_OBJECT',format_object),
        ('FORMAT_TYPE',format_type),
        ('NCBI_GI',ncbi_gi),
        ('RID',rid),
        ('RESULTS_FILE',results_file),
        ('SERVICE',service),
        ('SHOW_OVERVIEW',show_overview),
        ('CMD', 'Get'),
        ]
    query = [x for x in parameters if x[1] is not None]
    message = urllib.urlencode(query)

    # Poll NCBI until the results are ready.  Use a 3 second wait
    delay = 3.0
    previous = time.time()
    while True:
        current = time.time()
        wait = previous + delay - current
        if wait > 0:
            time.sleep(wait)
            previous = current + wait
        else:
            previous = current

        request = urllib2.Request("http://blast.ncbi.nlm.nih.gov/Blast.cgi",
                                  message,
                                  {"User-Agent":"BiopythonClient"})
        handle = urllib2.urlopen(request)
        results = handle.read()
        # XML results don't have the Status tag when finished
        if results.find("Status=") < 0:
            break
        i = results.index("Status=")
        j = results.index("\n", i)
        status = results[i+len("Status="):j].strip()
        if status.upper() == "READY":
            break

    return StringIO.StringIO(results)

def _parse_qblast_ref_page(handle):
    """Extract a tuple of RID, RTOE from the 'please wait' page (PRIVATE).

    The NCBI FAQ pages use TOE for 'Time of Execution', so RTOE is proably
    'Request Time of Execution' and RID would be 'Request Identifier'.
    """
    s = handle.read()
    i = s.find("RID =")
    if i == -1 :
        rid = None
    else :
        j = s.find("\n", i)
        rid = s[i+len("RID ="):j].strip()

    i = s.find("RTOE =")
    if i == -1 :
        rtoe = None
    else :
        j = s.find("\n", i)
        rtoe = s[i+len("RTOE ="):j].strip()

    if not rid and not rtoe :
        #Can we reliably extract the error message from the HTML page?
        #e.g.  "Message ID#24 Error: Failed to read the Blast query:
        #       Nucleotide FASTA provided for protein sequence"
        #This occurs inside a <div class="error msInf"> entry so it might
        #be possible to grab this...
        raise ValueError("No RID and no RTOE found in the 'please wait' page."
                         " (there was probably a problem with your request)")
    elif not rid :
        #Can this happen?
        raise ValueError("No RID found in the 'please wait' page."
                         " (although RTOE = %s)" % repr(rtoe))
    elif not rtoe :
        #Can this happen?
        raise ValueError("No RTOE found in the 'please wait' page."
                         " (although RID = %s)" % repr(rid))

    try :
        return rid, int(rtoe)
    except ValueError :
        raise ValueError("A non-integer RTOE found in " \
                         +"the 'please wait' page, %s" % repr(rtoe))
