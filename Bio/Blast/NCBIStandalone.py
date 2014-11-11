# Copyright 1999-2000 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
# Patches by Mike Poidinger to support multiple databases.
# Updated by Peter Cock in 2007 to do a better job on BLAST 2.2.15

"""Code for calling standalone BLAST and parsing plain text output (DEPRECATED).

Rather than parsing the human readable plain text BLAST output (which seems to
change with every update to BLAST), we and the NBCI recommend you parse the
XML output instead. The plain text parser in this module still works at the
time of writing, but is considered obsolete and updating it to cope with the
latest versions of BLAST is not a priority for us.

This module also provides code to work with the "legacy" standalone version of
NCBI BLAST, tools blastall, rpsblast and blastpgp via three helper functions of
the same name. These functions are very limited for dealing with the output as
files rather than handles, for which the wrappers in Bio.Blast.Applications are
preferred. Furthermore, the NCBI themselves regard these command line tools as
"legacy", and encourage using the new BLAST+ tools instead. Biopython has
wrappers for these under Bio.Blast.Applications (see the tutorial).
"""

from __future__ import print_function

from Bio import BiopythonDeprecationWarning
import warnings
warnings.warn("This module has been deprecated. Consider Bio.SearchIO for "
              "parsing BLAST output instead.", BiopythonDeprecationWarning)

import os
import re
from Bio._py3k import StringIO

from Bio import File
from Bio.ParserSupport import *
from Bio.Blast import Record
from Bio.Application import _escape_filename

__docformat__ = "restructuredtext en"


class LowQualityBlastError(Exception):
    """Error caused by running a low quality sequence through BLAST.

    When low quality sequences (like GenBank entries containing only
    stretches of a single nucleotide) are BLASTed, they will result in
    BLAST generating an error and not being able to perform the BLAST.
    search. This error should be raised for the BLAST reports produced
    in this case.
    """
    pass


class ShortQueryBlastError(Exception):
    """Error caused by running a short query sequence through BLAST.

    If the query sequence is too short, BLAST outputs warnings and errors::

        Searching[blastall] WARNING:  [000.000]  AT1G08320: SetUpBlastSearch failed.
        [blastall] ERROR:  [000.000]  AT1G08320: Blast:
        [blastall] ERROR:  [000.000]  AT1G08320: Blast: Query must be at least wordsize
        done

    This exception is raised when that condition is detected.
    """
    pass


class _Scanner(object):
    """Scan BLAST output from blastall or blastpgp.

    Tested with blastall and blastpgp v2.0.10, v2.0.11

    Methods:
     - feed     Feed data into the scanner.
    """
    def feed(self, handle, consumer):
        """S.feed(handle, consumer)

        Feed in a BLAST report for scanning.  handle is a file-like
        object that contains the BLAST report.  consumer is a Consumer
        object that will receive events as the report is scanned.
        """
        if isinstance(handle, File.UndoHandle):
            uhandle = handle
        else:
            uhandle = File.UndoHandle(handle)

        # Try to fast-forward to the beginning of the blast report.
        read_and_call_until(uhandle, consumer.noevent, contains='BLAST')
        # Now scan the BLAST report.
        self._scan_header(uhandle, consumer)
        self._scan_rounds(uhandle, consumer)
        self._scan_database_report(uhandle, consumer)
        self._scan_parameters(uhandle, consumer)

    def _scan_header(self, uhandle, consumer):
        # BLASTP 2.0.10 [Aug-26-1999]
        #
        #
        # Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaf
        # Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997),
        # "Gapped BLAST and PSI-BLAST: a new generation of protein database sea
        # programs",  Nucleic Acids Res. 25:3389-3402.
        #
        # Query= test
        #          (140 letters)
        #
        # Database: sdqib40-1.35.seg.fa
        #            1323 sequences; 223,339 total letters
        #
        # ========================================================
        # This next example is from the online version of Blast,
        # note there are TWO references, an RID line, and also
        # the database is BEFORE the query line.
        # Note there possibleuse of non-ASCII in the author names.
        # ========================================================
        #
        # BLASTP 2.2.15 [Oct-15-2006]
        # Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Sch??ffer,
        # Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman
        # (1997), "Gapped BLAST and PSI-BLAST: a new generation of
        # protein database search programs", Nucleic Acids Res. 25:3389-3402.
        #
        # Reference: Sch??ffer, Alejandro A., L. Aravind, Thomas L. Madden, Sergei
        # Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and
        # Stephen F. Altschul (2001), "Improving the accuracy of PSI-BLAST
        # protein database searches with composition-based statistics
        # and other refinements", Nucleic Acids Res. 29:2994-3005.
        #
        # RID: 1166022616-19998-65316425856.BLASTQ1
        #
        #
        # Database: All non-redundant GenBank CDS
        # translations+PDB+SwissProt+PIR+PRF excluding environmental samples
        #            4,254,166 sequences; 1,462,033,012 total letters
        # Query=  gi:16127998
        # Length=428
        #

        consumer.start_header()

        read_and_call(uhandle, consumer.version, contains='BLAST')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        # There might be a <pre> line, for qblast output.
        attempt_read_and_call(uhandle, consumer.noevent, start="<pre>")

        # Read the reference(s)
        while attempt_read_and_call(uhandle,
                                consumer.reference, start='Reference'):
            # References are normally multiline terminated by a blank line
            # (or, based on the old code, the RID line)
            while True:
                line = uhandle.readline()
                if is_blank_line(line):
                    consumer.noevent(line)
                    break
                elif line.startswith("RID"):
                    break
                else:
                    # More of the reference
                    consumer.reference(line)

        # Deal with the optional RID: ...
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        attempt_read_and_call(uhandle, consumer.reference, start="RID:")
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        # blastpgp may have a reference for compositional score matrix
        # adjustment (see Bug 2502):
        if attempt_read_and_call(
           uhandle, consumer.reference, start="Reference"):
            read_and_call_until(uhandle, consumer.reference, blank=1)
            read_and_call_while(uhandle, consumer.noevent, blank=1)

        # blastpgp has a Reference for composition-based statistics.
        if attempt_read_and_call(
           uhandle, consumer.reference, start="Reference"):
            read_and_call_until(uhandle, consumer.reference, blank=1)
            read_and_call_while(uhandle, consumer.noevent, blank=1)

        line = uhandle.peekline()
        assert line.strip() != ""
        assert not line.startswith("RID:")
        if line.startswith("Query="):
            # This is an old style query then database...

            # Read the Query lines and the following blank line.
            read_and_call(uhandle, consumer.query_info, start='Query=')
            read_and_call_until(uhandle, consumer.query_info, blank=1)
            read_and_call_while(uhandle, consumer.noevent, blank=1)

            # Read the database lines and the following blank line.
            read_and_call_until(uhandle, consumer.database_info, end='total letters')
            read_and_call(uhandle, consumer.database_info, contains='sequences')
            read_and_call_while(uhandle, consumer.noevent, blank=1)
        elif line.startswith("Database:"):
            # This is a new style database then query...
            read_and_call_until(uhandle, consumer.database_info, end='total letters')
            read_and_call(uhandle, consumer.database_info, contains='sequences')
            read_and_call_while(uhandle, consumer.noevent, blank=1)

            # Read the Query lines and the following blank line.
            # Or, on BLAST 2.2.22+ there is no blank link - need to spot
            # the "... Score     E" line instead.
            read_and_call(uhandle, consumer.query_info, start='Query=')
            # BLAST 2.2.25+ has a blank line before Length=
            read_and_call_until(uhandle, consumer.query_info, start='Length=')
            while True:
                line = uhandle.peekline()
                if not line.strip() or "Score     E" in line:
                    break
                # It is more of the query (and its length)
                read_and_call(uhandle, consumer.query_info)
            read_and_call_while(uhandle, consumer.noevent, blank=1)
        else:
            raise ValueError("Invalid header?")

        consumer.end_header()

    def _scan_rounds(self, uhandle, consumer):
        # Scan a bunch of rounds.
        # Each round begins with either a "Searching......" line
        # or a 'Score     E' line followed by descriptions and alignments.
        # The email server doesn't give the "Searching....." line.
        # If there is no 'Searching.....' line then you'll first see a
        # 'Results from round' line

        while not self._eof(uhandle):
            line = safe_peekline(uhandle)
            if not line.startswith('Searching') and \
               not line.startswith('Results from round') and \
               re.search(r"Score +E", line) is None and \
               'No hits found' not in line:
                break
            self._scan_descriptions(uhandle, consumer)
            self._scan_alignments(uhandle, consumer)

    def _scan_descriptions(self, uhandle, consumer):
        # Searching..................................................done
        # Results from round 2
        #
        #
        #                                                                    Sc
        # Sequences producing significant alignments:                        (b
        # Sequences used in model and found again:
        #
        # d1tde_2 3.4.1.4.4 (119-244) Thioredoxin reductase [Escherichia ...
        # d1tcob_ 1.31.1.5.16 Calcineurin regulatory subunit (B-chain) [B...
        # d1symb_ 1.31.1.2.2 Calcyclin (S100) [RAT (RATTUS NORVEGICUS)]
        #
        # Sequences not found previously or not previously below threshold:
        #
        # d1osa__ 1.31.1.5.11 Calmodulin [Paramecium tetraurelia]
        # d1aoza3 2.5.1.3.3 (339-552) Ascorbate oxidase [zucchini (Cucurb...
        #

        # If PSI-BLAST, may also have:
        #
        # CONVERGED!

        consumer.start_descriptions()

        # Read 'Searching'
        # This line seems to be missing in BLASTN 2.1.2 (others?)
        attempt_read_and_call(uhandle, consumer.noevent, start='Searching')

        # blastpgp 2.0.10 from NCBI 9/19/99 for Solaris sometimes crashes here.
        # If this happens, the handle will yield no more information.
        if not uhandle.peekline():
            raise ValueError("Unexpected end of blast report.  " +
                  "Looks suspiciously like a PSI-BLAST crash.")

        # BLASTN 2.2.3 sometimes spews a bunch of warnings and errors here:
        # Searching[blastall] WARNING:  [000.000]  AT1G08320: SetUpBlastSearch
        # [blastall] ERROR:  [000.000]  AT1G08320: Blast:
        # [blastall] ERROR:  [000.000]  AT1G08320: Blast: Query must be at leas
        # done
        # Reported by David Weisman.
        # Check for these error lines and ignore them for now.  Let
        # the BlastErrorParser deal with them.
        line = uhandle.peekline()
        if "ERROR:" in line or line.startswith("done"):
            read_and_call_while(uhandle, consumer.noevent, contains="ERROR:")
            read_and_call(uhandle, consumer.noevent, start="done")

        # Check to see if this is PSI-BLAST.
        # If it is, the 'Searching' line will be followed by:
        # (version 2.0.10)
        #     Searching.............................
        #     Results from round 2
        # or (version 2.0.11)
        #     Searching.............................
        #
        #
        #     Results from round 2

        # Skip a bunch of blank lines.
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        # Check for the results line if it's there.
        if attempt_read_and_call(uhandle, consumer.round, start='Results'):
            read_and_call_while(uhandle, consumer.noevent, blank=1)

        # Three things can happen here:
        # 1.  line contains 'Score     E'
        # 2.  line contains "No hits found"
        # 3.  no descriptions
        # The first one begins a bunch of descriptions.  The last two
        # indicates that no descriptions follow, and we should go straight
        # to the alignments.
        if not attempt_read_and_call(
           uhandle, consumer.description_header,
           has_re=re.compile(r'Score +E')):
            # Either case 2 or 3.  Look for "No hits found".
            attempt_read_and_call(uhandle, consumer.no_hits,
                                  contains='No hits found')
            try:
                read_and_call_while(uhandle, consumer.noevent, blank=1)
            except ValueError as err:
                if str(err) != "Unexpected end of stream.":
                    raise err

            consumer.end_descriptions()
            # Stop processing.
            return

        # Read the score header lines
        read_and_call(uhandle, consumer.description_header,
                      start='Sequences producing')

        # If PSI-BLAST, read the 'Sequences used in model' line.
        attempt_read_and_call(uhandle, consumer.model_sequences,
                              start='Sequences used in model')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        # In BLAT, rather than a "No hits found" line, we just
        # get no descriptions (and no alignments). This can be
        # spotted because the next line is the database block:
        if safe_peekline(uhandle).startswith("  Database:"):
            consumer.end_descriptions()
            # Stop processing.
            return

        # Read the descriptions and the following blank lines, making
        # sure that there are descriptions.
        if not uhandle.peekline().startswith('Sequences not found'):
            read_and_call_until(uhandle, consumer.description, blank=1)
            read_and_call_while(uhandle, consumer.noevent, blank=1)

        # If PSI-BLAST, read the 'Sequences not found' line followed
        # by more descriptions.  However, I need to watch out for the
        # case where there were no sequences not found previously, in
        # which case there will be no more descriptions.
        if attempt_read_and_call(uhandle, consumer.nonmodel_sequences,
                                 start='Sequences not found'):
            # Read the descriptions and the following blank lines.
            read_and_call_while(uhandle, consumer.noevent, blank=1)
            l = safe_peekline(uhandle)
            # Brad -- added check for QUERY. On some PSI-BLAST outputs
            # there will be a 'Sequences not found' line followed by no
            # descriptions. Check for this case since the first thing you'll
            # get is a blank line and then 'QUERY'
            if not l.startswith('CONVERGED') and l[0] != '>' \
                    and not l.startswith('QUERY'):
                read_and_call_until(uhandle, consumer.description, blank=1)
                read_and_call_while(uhandle, consumer.noevent, blank=1)

        attempt_read_and_call(uhandle, consumer.converged, start='CONVERGED')
        read_and_call_while(uhandle, consumer.noevent, blank=1)

        consumer.end_descriptions()

    def _scan_alignments(self, uhandle, consumer):
        if self._eof(uhandle):
            return

        # qblast inserts a helpful line here.
        attempt_read_and_call(uhandle, consumer.noevent, start="ALIGNMENTS")

        # First, check to see if I'm at the database report.
        line = safe_peekline(uhandle)
        if not line:
            # EOF
            return
        elif line.startswith('  Database') or line.startswith("Lambda"):
            return
        elif line[0] == '>':
            # XXX make a better check here between pairwise and masterslave
            self._scan_pairwise_alignments(uhandle, consumer)
        elif line.startswith('Effective'):
            return
        else:
            # XXX put in a check to make sure I'm in a masterslave alignment
            self._scan_masterslave_alignment(uhandle, consumer)

    def _scan_pairwise_alignments(self, uhandle, consumer):
        while not self._eof(uhandle):
            line = safe_peekline(uhandle)
            if line[0] != '>':
                break
            self._scan_one_pairwise_alignment(uhandle, consumer)

    def _scan_one_pairwise_alignment(self, uhandle, consumer):
        if self._eof(uhandle):
            return
        consumer.start_alignment()

        self._scan_alignment_header(uhandle, consumer)

        # Scan a bunch of score/alignment pairs.
        while True:
            if self._eof(uhandle):
                # Shouldn't have issued that _scan_alignment_header event...
                break
            line = safe_peekline(uhandle)
            if not line.startswith(' Score'):
                break
            self._scan_hsp(uhandle, consumer)
        consumer.end_alignment()

    def _scan_alignment_header(self, uhandle, consumer):
        # >d1rip__ 2.24.7.1.1 Ribosomal S17 protein [Bacillus
        #           stearothermophilus]
        #           Length = 81
        #
        # Or, more recently with different white space:
        #
        # >gi|15799684|ref|NP_285696.1| threonine synthase ...
        #  gi|15829258|ref|NP_308031.1| threonine synthase
        #  ...
        # Length=428
        read_and_call(uhandle, consumer.title, start='>')
        while True:
            line = safe_readline(uhandle)
            if line.lstrip().startswith('Length =') \
            or line.lstrip().startswith('Length='):
                consumer.length(line)
                break
            elif is_blank_line(line):
                # Check to make sure I haven't missed the Length line
                raise ValueError("I missed the Length in an alignment header")
            consumer.title(line)

        # Older versions of BLAST will have a line with some spaces.
        # Version 2.0.14 (maybe 2.0.13?) and above print a true blank line.
        if not attempt_read_and_call(uhandle, consumer.noevent,
                                     start='          '):
            read_and_call(uhandle, consumer.noevent, blank=1)

    def _scan_hsp(self, uhandle, consumer):
        consumer.start_hsp()
        self._scan_hsp_header(uhandle, consumer)
        self._scan_hsp_alignment(uhandle, consumer)
        consumer.end_hsp()

    def _scan_hsp_header(self, uhandle, consumer):
        #  Score = 22.7 bits (47), Expect = 2.5
        #  Identities = 10/36 (27%), Positives = 18/36 (49%)
        #  Strand = Plus / Plus
        #  Frame = +3
        #

        read_and_call(uhandle, consumer.score, start=' Score')
        read_and_call(uhandle, consumer.identities, start=' Identities')
        # BLASTN
        attempt_read_and_call(uhandle, consumer.strand, start=' Strand')
        # BLASTX, TBLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.frame, start=' Frame')
        read_and_call(uhandle, consumer.noevent, blank=1)

    def _scan_hsp_alignment(self, uhandle, consumer):
        # Query: 11 GRGVSACA-------TCDGFFYRNQKVAVIGGGNTAVEEALYLSNIASEVHLIHRRDGF
        #           GRGVS+         TC    Y  + + V GGG+ + EE   L     +   I R+
        # Sbjct: 12 GRGVSSVVRRCIHKPTCKE--YAVKIIDVTGGGSFSAEEVQELREATLKEVDILRKVSG
        #
        # Query: 64 AEKILIKR 71
        #              I +K
        # Sbjct: 70 PNIIQLKD 77
        #

        while True:
            # Blastn adds an extra line filled with spaces before Query
            attempt_read_and_call(uhandle, consumer.noevent, start='     ')
            read_and_call(uhandle, consumer.query, start='Query')
            read_and_call(uhandle, consumer.align, start='     ')
            read_and_call(uhandle, consumer.sbjct, start='Sbjct')
            try:
                read_and_call_while(uhandle, consumer.noevent, blank=1)
            except ValueError as err:
                if str(err) != "Unexpected end of stream.":
                    raise err
                # End of File (well, it looks like it with recent versions
                # of BLAST for multiple queries after the Iterator class
                # has broken up the whole file into chunks).
                break
            line = safe_peekline(uhandle)
            # Alignment continues if I see a 'Query' or the spaces for Blastn.
            if not (line.startswith('Query') or line.startswith('     ')):
                break

    def _scan_masterslave_alignment(self, uhandle, consumer):
        consumer.start_alignment()
        while True:
            line = safe_readline(uhandle)
            # Check to see whether I'm finished reading the alignment.
            # This is indicated by 1) database section, 2) next psi-blast
            # round, which can also be a 'Results from round' if no
            # searching line is present
            # patch by chapmanb
            if line.startswith('Searching') or \
                    line.startswith('Results from round'):
                uhandle.saveline(line)
                break
            elif line.startswith('  Database'):
                uhandle.saveline(line)
                break
            elif is_blank_line(line):
                consumer.noevent(line)
            else:
                consumer.multalign(line)
        read_and_call_while(uhandle, consumer.noevent, blank=1)
        consumer.end_alignment()

    def _eof(self, uhandle):
        try:
            line = safe_peekline(uhandle)
        except ValueError as err:
            if str(err) != "Unexpected end of stream.":
                raise err
            line = ""
        return not line

    def _scan_database_report(self, uhandle, consumer):
        #   Database: sdqib40-1.35.seg.fa
        #     Posted date:  Nov 1, 1999  4:25 PM
        #   Number of letters in database: 223,339
        #   Number of sequences in database:  1323
        #
        # Lambda     K      H
        #    0.322    0.133    0.369
        #
        # Gapped
        # Lambda     K      H
        #    0.270   0.0470    0.230
        #
        ##########################################
        # Or, more recently Blast 2.2.15 gives less blank lines
        ##########################################
        #   Database: All non-redundant GenBank CDS translations+PDB+SwissProt+PIR+PRF excluding
        # environmental samples
        #     Posted date:  Dec 12, 2006  5:51 PM
        #   Number of letters in database: 667,088,753
        #   Number of sequences in database:  2,094,974
        # Lambda     K      H
        #    0.319    0.136    0.395
        # Gapped
        # Lambda     K      H
        #    0.267   0.0410    0.140

        if self._eof(uhandle):
            return

        consumer.start_database_report()

        # Subset of the database(s) listed below
        #    Number of letters searched: 562,618,960
        #    Number of sequences searched:  228,924
        if attempt_read_and_call(uhandle, consumer.noevent, start="  Subset"):
            read_and_call(uhandle, consumer.noevent, contains="letters")
            read_and_call(uhandle, consumer.noevent, contains="sequences")
            read_and_call(uhandle, consumer.noevent, start="  ")

        # Sameet Mehta reported seeing output from BLASTN 2.2.9 that
        # was missing the "Database" stanza completely.
        while attempt_read_and_call(uhandle, consumer.database,
                start='  Database'):
            # BLAT output ends abruptly here, without any of the other
            # information.  Check to see if this is the case.  If so,
            # then end the database report here gracefully.
            if not uhandle.peekline().strip() \
            or uhandle.peekline().startswith("BLAST"):
                consumer.end_database_report()
                return

            # Database can span multiple lines.
            read_and_call_until(uhandle, consumer.database, start='    Posted')
            read_and_call(uhandle, consumer.posted_date, start='    Posted')
            read_and_call(uhandle, consumer.num_letters_in_database,
                       start='  Number of letters')
            read_and_call(uhandle, consumer.num_sequences_in_database,
                       start='  Number of sequences')
            # There may not be a line starting with spaces...
            attempt_read_and_call(uhandle, consumer.noevent, start='  ')

            line = safe_readline(uhandle)
            uhandle.saveline(line)
            if 'Lambda' in line:
                break

        try:
            read_and_call(uhandle, consumer.noevent, start='Lambda')
            read_and_call(uhandle, consumer.ka_params)
        except:
            pass

        # This blank line is optional:
        attempt_read_and_call(uhandle, consumer.noevent, blank=1)

        # not BLASTP
        attempt_read_and_call(uhandle, consumer.gapped, start='Gapped')
        # not TBLASTX
        if attempt_read_and_call(uhandle, consumer.noevent, start='Lambda'):
            read_and_call(uhandle, consumer.ka_params_gap)

        # Blast 2.2.4 can sometimes skip the whole parameter section.
        # Thus, I need to be careful not to read past the end of the
        # file.
        try:
            read_and_call_while(uhandle, consumer.noevent, blank=1)
        except ValueError as x:
            if str(x) != "Unexpected end of stream.":
                raise
        consumer.end_database_report()

    def _scan_parameters(self, uhandle, consumer):
        # Matrix: BLOSUM62
        # Gap Penalties: Existence: 11, Extension: 1
        # Number of Hits to DB: 50604
        # Number of Sequences: 1323
        # Number of extensions: 1526
        # Number of successful extensions: 6
        # Number of sequences better than 10.0: 5
        # Number of HSP's better than 10.0 without gapping: 5
        # Number of HSP's successfully gapped in prelim test: 0
        # Number of HSP's that attempted gapping in prelim test: 1
        # Number of HSP's gapped (non-prelim): 5
        # length of query: 140
        # length of database: 223,339
        # effective HSP length: 39
        # effective length of query: 101
        # effective length of database: 171,742
        # effective search space: 17345942
        # effective search space used: 17345942
        # T: 11
        # A: 40
        # X1: 16 ( 7.4 bits)
        # X2: 38 (14.8 bits)
        # X3: 64 (24.9 bits)
        # S1: 41 (21.9 bits)
        # S2: 42 (20.8 bits)
        ##########################################
        # Or, more recently Blast(x) 2.2.15 gives
        ##########################################
        # Matrix: BLOSUM62
        # Gap Penalties: Existence: 11, Extension: 1
        # Number of Sequences: 4535438
        # Number of Hits to DB: 2,588,844,100
        # Number of extensions: 60427286
        # Number of successful extensions: 126433
        # Number of sequences better than  2.0: 30
        # Number of HSP's gapped: 126387
        # Number of HSP's successfully gapped: 35
        # Length of query: 291
        # Length of database: 1,573,298,872
        # Length adjustment: 130
        # Effective length of query: 161
        # Effective length of database: 983,691,932
        # Effective search space: 158374401052
        # Effective search space used: 158374401052
        # Neighboring words threshold: 12
        # Window for multiple hits: 40
        # X1: 16 ( 7.3 bits)
        # X2: 38 (14.6 bits)
        # X3: 64 (24.7 bits)
        # S1: 41 (21.7 bits)
        # S2: 32 (16.9 bits)

        # Blast 2.2.4 can sometimes skip the whole parameter section.
        # BLAT also skips the whole parameter section.
        # Thus, check to make sure that the parameter section really
        # exists.
        if not uhandle.peekline().strip():
            return

        # BLASTN 2.2.9 looks like it reverses the "Number of Hits" and
        # "Number of Sequences" lines.
        consumer.start_parameters()

        # Matrix line may be missing in BLASTN 2.2.9
        attempt_read_and_call(uhandle, consumer.matrix, start='Matrix')
        # not TBLASTX
        attempt_read_and_call(uhandle, consumer.gap_penalties, start='Gap')

        attempt_read_and_call(uhandle, consumer.num_sequences,
                              start='Number of Sequences')
        attempt_read_and_call(uhandle, consumer.num_hits,
                      start='Number of Hits')
        attempt_read_and_call(uhandle, consumer.num_sequences,
                              start='Number of Sequences')
        attempt_read_and_call(uhandle, consumer.num_extends,
                      start='Number of extensions')
        attempt_read_and_call(uhandle, consumer.num_good_extends,
                      start='Number of successful')

        attempt_read_and_call(uhandle, consumer.num_seqs_better_e,
                      start='Number of sequences')

        # not BLASTN, TBLASTX
        if attempt_read_and_call(uhandle, consumer.hsps_no_gap,
                                 start="Number of HSP's better"):
            # BLASTN 2.2.9
            if attempt_read_and_call(uhandle, consumer.noevent,
                                     start="Number of HSP's gapped:"):
                read_and_call(uhandle, consumer.noevent,
                              start="Number of HSP's successfully")
                # This is omitted in 2.2.15
                attempt_read_and_call(uhandle, consumer.noevent,
                              start="Number of extra gapped extensions")
            else:
                read_and_call(uhandle, consumer.hsps_prelim_gapped,
                              start="Number of HSP's successfully")
                read_and_call(uhandle, consumer.hsps_prelim_gap_attempted,
                              start="Number of HSP's that")
                read_and_call(uhandle, consumer.hsps_gapped,
                              start="Number of HSP's gapped")
        # e.g. BLASTX 2.2.15 where the "better" line is missing
        elif attempt_read_and_call(uhandle, consumer.noevent,
                                     start="Number of HSP's gapped"):
            read_and_call(uhandle, consumer.noevent,
                          start="Number of HSP's successfully")

        # not in blastx 2.2.1
        attempt_read_and_call(uhandle, consumer.query_length,
                              has_re=re.compile(r"[Ll]ength of query"))
        # Not in BLASTX 2.2.22+
        attempt_read_and_call(uhandle, consumer.database_length,
                          has_re=re.compile(r"[Ll]ength of \s*[Dd]atabase"))

        # BLASTN 2.2.9
        attempt_read_and_call(uhandle, consumer.noevent,
                              start="Length adjustment")
        attempt_read_and_call(uhandle, consumer.effective_hsp_length,
                              start='effective HSP')
        # Not in blastx 2.2.1
        attempt_read_and_call(
            uhandle, consumer.effective_query_length,
            has_re=re.compile(r'[Ee]ffective length of query'))

        # This is not in BLASTP 2.2.15
        attempt_read_and_call(
            uhandle, consumer.effective_database_length,
            has_re=re.compile(r'[Ee]ffective length of \s*[Dd]atabase'))
        # Not in blastx 2.2.1, added a ':' to distinguish between
        # this and the 'effective search space used' line
        attempt_read_and_call(
            uhandle, consumer.effective_search_space,
            has_re=re.compile(r'[Ee]ffective search space:'))
        # Does not appear in BLASTP 2.0.5
        attempt_read_and_call(
            uhandle, consumer.effective_search_space_used,
            has_re=re.compile(r'[Ee]ffective search space used'))

        # BLASTX, TBLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.frameshift, start='frameshift')

        # not in BLASTN 2.2.9
        attempt_read_and_call(uhandle, consumer.threshold, start='T')
        # In BLASTX 2.2.15 replaced by: "Neighboring words threshold: 12"
        attempt_read_and_call(uhandle, consumer.threshold, start='Neighboring words threshold')

        # not in BLASTX 2.2.15
        attempt_read_and_call(uhandle, consumer.window_size, start='A')
        # get this instead: "Window for multiple hits: 40"
        attempt_read_and_call(uhandle, consumer.window_size, start='Window for multiple hits')

        # not in BLASTX 2.2.22+
        attempt_read_and_call(uhandle, consumer.dropoff_1st_pass, start='X1')
        # not TBLASTN
        attempt_read_and_call(uhandle, consumer.gap_x_dropoff, start='X2')

        # not BLASTN, TBLASTX
        attempt_read_and_call(uhandle, consumer.gap_x_dropoff_final,
                              start='X3')

        # not TBLASTN
        attempt_read_and_call(uhandle, consumer.gap_trigger, start='S1')
        # not in blastx 2.2.1
        # first we make sure we have additional lines to work with, if
        # not then the file is done and we don't have a final S2
        if not is_blank_line(uhandle.peekline(), allow_spaces=1):
            read_and_call(uhandle, consumer.blast_cutoff, start='S2')

        consumer.end_parameters()


class BlastParser(AbstractParser):
    """Parses BLAST data into a Record.Blast object.

    """
    def __init__(self):
        """__init__(self)"""
        self._scanner = _Scanner()
        self._consumer = _BlastConsumer()

    def parse(self, handle):
        """parse(self, handle)"""
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data


class PSIBlastParser(AbstractParser):
    """Parses BLAST data into a Record.PSIBlast object.

    """
    def __init__(self):
        """__init__(self)"""
        self._scanner = _Scanner()
        self._consumer = _PSIBlastConsumer()

    def parse(self, handle):
        """parse(self, handle)"""
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data


class _HeaderConsumer(object):
    def start_header(self):
        self._header = Record.Header()

    def version(self, line):
        c = line.split()
        self._header.application = c[0]
        self._header.version = c[1]
        if len(c) > 2:
            # The date is missing in the new C++ output from blastx 2.2.22+
            # Just get "BLASTX 2.2.22+\n" and that's all.
            self._header.date = c[2][1:-1]

    def reference(self, line):
        if line.startswith('Reference: '):
            self._header.reference = line[11:]
        else:
            self._header.reference = self._header.reference + line

    def query_info(self, line):
        if line.startswith('Query= '):
            self._header.query = line[7:].lstrip()
        elif line.startswith('Length='):
            # New style way to give the query length in BLAST 2.2.22+ (the C++ code)
            self._header.query_letters = _safe_int(line[7:].strip())
        elif not line.startswith('       '):  # continuation of query_info
            self._header.query = "%s%s" % (self._header.query, line)
        else:
            # Hope it is the old style way to give the query length:
            letters, = _re_search(
                r"([0-9,]+) letters", line,
                "I could not find the number of letters in line\n%s" % line)
            self._header.query_letters = _safe_int(letters)

    def database_info(self, line):
        line = line.rstrip()
        if line.startswith('Database: '):
            self._header.database = line[10:]
        elif not line.endswith('total letters'):
            if self._header.database:
                # Need to include a space when merging multi line datase descr
                self._header.database = self._header.database + " " + line.strip()
            else:
                self._header.database = line.strip()
        else:
            sequences, letters = _re_search(
                r"([0-9,]+) sequences; ([0-9,-]+) total letters", line,
                "I could not find the sequences and letters in line\n%s" % line)
            self._header.database_sequences = _safe_int(sequences)
            self._header.database_letters = _safe_int(letters)

    def end_header(self):
        # Get rid of the trailing newlines
        self._header.reference = self._header.reference.rstrip()
        self._header.query = self._header.query.rstrip()


class _DescriptionConsumer(object):
    def start_descriptions(self):
        self._descriptions = []
        self._model_sequences = []
        self._nonmodel_sequences = []
        self._converged = 0
        self._type = None
        self._roundnum = None

        self.__has_n = 0   # Does the description line contain an N value?

    def description_header(self, line):
        if line.startswith('Sequences producing'):
            cols = line.split()
            if cols[-1] == 'N':
                self.__has_n = 1

    def description(self, line):
        dh = self._parse(line)
        if self._type == 'model':
            self._model_sequences.append(dh)
        elif self._type == 'nonmodel':
            self._nonmodel_sequences.append(dh)
        else:
            self._descriptions.append(dh)

    def model_sequences(self, line):
        self._type = 'model'

    def nonmodel_sequences(self, line):
        self._type = 'nonmodel'

    def converged(self, line):
        self._converged = 1

    def no_hits(self, line):
        pass

    def round(self, line):
        if not line.startswith('Results from round'):
            raise ValueError("I didn't understand the round line\n%s" % line)
        self._roundnum = _safe_int(line[18:].strip())

    def end_descriptions(self):
        pass

    def _parse(self, description_line):
        line = description_line  # for convenience
        dh = Record.Description()

        # I need to separate the score and p-value from the title.
        # sp|P21297|FLBT_CAUCR FLBT PROTEIN     [snip]         284  7e-77
        # sp|P21297|FLBT_CAUCR FLBT PROTEIN     [snip]         284  7e-77  1
        # special cases to handle:
        #   - title must be preserved exactly (including whitespaces)
        #   - score could be equal to e-value (not likely, but what if??)
        #   - sometimes there's an "N" score of '1'.
        cols = line.split()
        if len(cols) < 3:
            raise ValueError(
                  "Line does not appear to contain description:\n%s" % line)
        if self.__has_n:
            i = line.rfind(cols[-1])        # find start of N
            i = line.rfind(cols[-2], 0, i)  # find start of p-value
            i = line.rfind(cols[-3], 0, i)  # find start of score
        else:
            i = line.rfind(cols[-1])        # find start of p-value
            i = line.rfind(cols[-2], 0, i)  # find start of score
        if self.__has_n:
            dh.title, dh.score, dh.e, dh.num_alignments = \
                      line[:i].rstrip(), cols[-3], cols[-2], cols[-1]
        else:
            dh.title, dh.score, dh.e, dh.num_alignments = \
                      line[:i].rstrip(), cols[-2], cols[-1], 1
        dh.num_alignments = _safe_int(dh.num_alignments)
        dh.score = _safe_int(dh.score)
        dh.e = _safe_float(dh.e)
        return dh


class _AlignmentConsumer(object):
    # This is a little bit tricky.  An alignment can either be a
    # pairwise alignment or a multiple alignment.  Since it's difficult
    # to know a-priori which one the blast record will contain, I'm going
    # to make one class that can parse both of them.
    def start_alignment(self):
        self._alignment = Record.Alignment()
        self._multiple_alignment = Record.MultipleAlignment()

    def title(self, line):
        if self._alignment.title:
            self._alignment.title += " "
        self._alignment.title += line.strip()

    def length(self, line):
        # e.g. "Length = 81" or more recently, "Length=428"
        parts = line.replace(" ", "").split("=")
        assert len(parts) == 2, "Unrecognised format length line"
        self._alignment.length = parts[1]
        self._alignment.length = _safe_int(self._alignment.length)

    def multalign(self, line):
        # Standalone version uses 'QUERY', while WWW version uses blast_tmp.
        if line.startswith('QUERY') or line.startswith('blast_tmp'):
            # If this is the first line of the multiple alignment,
            # then I need to figure out how the line is formatted.

            # Format of line is:
            # QUERY 1   acttg...gccagaggtggtttattcagtctccataagagaggggacaaacg 60
            try:
                name, start, seq, end = line.split()
            except ValueError:
                raise ValueError("I do not understand the line\n%s" % line)
            self._start_index = line.index(start, len(name))
            self._seq_index = line.index(seq,
                                         self._start_index + len(start))
            # subtract 1 for the space
            self._name_length = self._start_index - 1
            self._start_length = self._seq_index - self._start_index - 1
            self._seq_length = line.rfind(end) - self._seq_index - 1

            # self._seq_index = line.index(seq)
            # # subtract 1 for the space
            # self._seq_length = line.rfind(end) - self._seq_index - 1
            # self._start_index = line.index(start)
            # self._start_length = self._seq_index - self._start_index - 1
            # self._name_length = self._start_index

        # Extract the information from the line
        name = line[:self._name_length]
        name = name.rstrip()
        start = line[self._start_index:self._start_index + self._start_length]
        start = start.rstrip()
        if start:
            start = _safe_int(start)
        end = line[self._seq_index + self._seq_length:].rstrip()
        if end:
            end = _safe_int(end)
        seq = line[self._seq_index:self._seq_index + self._seq_length].rstrip()
        # right pad the sequence with spaces if necessary
        if len(seq) < self._seq_length:
            seq += ' ' * (self._seq_length - len(seq))

        # I need to make sure the sequence is aligned correctly with the query.
        # First, I will find the length of the query.  Then, if necessary,
        # I will pad my current sequence with spaces so that they will line
        # up correctly.

        # Two possible things can happen:
        # QUERY
        # 504
        #
        # QUERY
        # 403
        #
        # Sequence 504 will need padding at the end.  Since I won't know
        # this until the end of the alignment, this will be handled in
        # end_alignment.
        # Sequence 403 will need padding before being added to the alignment.

        align = self._multiple_alignment.alignment  # for convenience
        align.append((name, start, seq, end))

        # This is old code that tried to line up all the sequences
        # in a multiple alignment by using the sequence title's as
        # identifiers.  The problem with this is that BLAST assigns
        # different HSP's from the same sequence the same id.  Thus,
        # in one alignment block, there may be multiple sequences with
        # the same id.  I'm not sure how to handle this, so I'm not
        # going to.

        # # If the sequence is the query, then just add it.
        # if name == 'QUERY':
        #     if len(align) == 0:
        #         align.append((name, start, seq))
        #     else:
        #         aname, astart, aseq = align[0]
        #         if name != aname:
        #             raise ValueError, "Query is not the first sequence"
        #         aseq = aseq + seq
        #         align[0] = aname, astart, aseq
        # else:
        #     if len(align) == 0:
        #         raise ValueError, "I could not find the query sequence"
        #     qname, qstart, qseq = align[0]
        #
        #     # Now find my sequence in the multiple alignment.
        #     for i in range(1, len(align)):
        #         aname, astart, aseq = align[i]
        #         if name == aname:
        #             index = i
        #             break
        #     else:
        #         # If I couldn't find it, then add a new one.
        #         align.append((None, None, None))
        #         index = len(align)-1
        #         # Make sure to left-pad it.
        #         aname, astart, aseq = name, start, ' '*(len(qseq)-len(seq))
        #
        #     if len(qseq) != len(aseq) + len(seq):
        #         # If my sequences are shorter than the query sequence,
        #         # then I will need to pad some spaces to make them line up.
        #         # Since I've already right padded seq, that means aseq
        #         # must be too short.
        #         aseq = aseq + ' '*(len(qseq)-len(aseq)-len(seq))
        #     aseq = aseq + seq
        #     if astart is None:
        #         astart = start
        #     align[index] = aname, astart, aseq

    def end_alignment(self):
        # Remove trailing newlines
        if self._alignment:
            self._alignment.title = self._alignment.title.rstrip()

        # This code is also obsolete.  See note above.
        # If there's a multiple alignment, I will need to make sure
        # all the sequences are aligned.  That is, I may need to
        # right-pad the sequences.
        # if self._multiple_alignment is not None:
        #     align = self._multiple_alignment.alignment
        #     seqlen = None
        #     for i in range(len(align)):
        #         name, start, seq = align[i]
        #         if seqlen is None:
        #             seqlen = len(seq)
        #         else:
        #             if len(seq) < seqlen:
        #                 seq = seq + ' '*(seqlen - len(seq))
        #                 align[i] = name, start, seq
        #             elif len(seq) > seqlen:
        #                 raise ValueError, \
        #                       "Sequence %s is longer than the query" % name

        # Clean up some variables, if they exist.
        try:
            del self._seq_index
            del self._seq_length
            del self._start_index
            del self._start_length
            del self._name_length
        except AttributeError:
            pass


class _HSPConsumer(object):
    def start_hsp(self):
        self._hsp = Record.HSP()

    def score(self, line):
        self._hsp.bits, self._hsp.score = _re_search(
            r"Score =\s*([0-9.e+]+) bits \(([0-9]+)\)", line,
            "I could not find the score in line\n%s" % line)
        self._hsp.score = _safe_float(self._hsp.score)
        self._hsp.bits = _safe_float(self._hsp.bits)

        x, y = _re_search(
            r"Expect\(?(\d*)\)? = +([0-9.e\-|\+]+)", line,
            "I could not find the expect in line\n%s" % line)
        if x:
            self._hsp.num_alignments = _safe_int(x)
        else:
            self._hsp.num_alignments = 1
        self._hsp.expect = _safe_float(y)

    def identities(self, line):
        x, y = _re_search(
            r"Identities = (\d+)\/(\d+)", line,
            "I could not find the identities in line\n%s" % line)
        self._hsp.identities = _safe_int(x), _safe_int(y)
        self._hsp.align_length = _safe_int(y)

        if 'Positives' in line:
            x, y = _re_search(
                r"Positives = (\d+)\/(\d+)", line,
                "I could not find the positives in line\n%s" % line)
            self._hsp.positives = _safe_int(x), _safe_int(y)
            assert self._hsp.align_length == _safe_int(y)

        if 'Gaps' in line:
            x, y = _re_search(
                r"Gaps = (\d+)\/(\d+)", line,
                "I could not find the gaps in line\n%s" % line)
            self._hsp.gaps = _safe_int(x), _safe_int(y)
            assert self._hsp.align_length == _safe_int(y)

    def strand(self, line):
        self._hsp.strand = _re_search(
            r"Strand\s?=\s?(\w+)\s?/\s?(\w+)", line,
            "I could not find the strand in line\n%s" % line)

    def frame(self, line):
        # Frame can be in formats:
        # Frame = +1
        # Frame = +2 / +2
        if '/' in line:
            self._hsp.frame = _re_search(
                r"Frame\s?=\s?([-+][123])\s?/\s?([-+][123])", line,
                "I could not find the frame in line\n%s" % line)
        else:
            self._hsp.frame = _re_search(
                r"Frame = ([-+][123])", line,
                "I could not find the frame in line\n%s" % line)

    # Match a space, if one is available.  Masahir Ishikawa found a
    # case where there's no space between the start and the sequence:
    # Query: 100tt 101
    # line below modified by Yair Benita, Sep 2004
    # Note that the colon is not always present. 2006
    _query_re = re.compile(r"Query(:?) \s*(\d+)\s*(.+) (\d+)")

    def query(self, line):
        m = self._query_re.search(line)
        if m is None:
            raise ValueError("I could not find the query in line\n%s" % line)

        # line below modified by Yair Benita, Sep 2004.
        # added the end attribute for the query
        colon, start, seq, end = m.groups()
        self._hsp.query = self._hsp.query + seq
        if self._hsp.query_start is None:
            self._hsp.query_start = _safe_int(start)

        # line below added by Yair Benita, Sep 2004.
        # added the end attribute for the query
        self._hsp.query_end = _safe_int(end)

        # Get index for sequence start (regular expression element 3)
        self._query_start_index = m.start(3)
        self._query_len = len(seq)

    def align(self, line):
        seq = line[self._query_start_index:].rstrip()
        if len(seq) < self._query_len:
            # Make sure the alignment is the same length as the query
            seq += ' ' * (self._query_len - len(seq))
        elif len(seq) < self._query_len:
            raise ValueError("Match is longer than the query in line\n%s"
                             % line)
        self._hsp.match = self._hsp.match + seq

    # To match how we do the query, cache the regular expression.
    # Note that the colon is not always present.
    _sbjct_re = re.compile(r"Sbjct(:?) \s*(\d+)\s*(.+) (\d+)")

    def sbjct(self, line):
        m = self._sbjct_re.search(line)
        if m is None:
            raise ValueError("I could not find the sbjct in line\n%s" % line)
        colon, start, seq, end = m.groups()
        # mikep 26/9/00
        # On occasion, there is a blast hit with no subject match
        # so far, it only occurs with 1-line short "matches"
        # I have decided to let these pass as they appear
        if not seq.strip():
            seq = ' ' * self._query_len
        self._hsp.sbjct = self._hsp.sbjct + seq
        if self._hsp.sbjct_start is None:
            self._hsp.sbjct_start = _safe_int(start)

        self._hsp.sbjct_end = _safe_int(end)
        if len(seq) != self._query_len:
            raise ValueError(
                  "QUERY and SBJCT sequence lengths don't match in line\n%s"
                  % line)

        del self._query_start_index   # clean up unused variables
        del self._query_len

    def end_hsp(self):
        pass


class _DatabaseReportConsumer(object):

    def start_database_report(self):
        self._dr = Record.DatabaseReport()

    def database(self, line):
        m = re.search(r"Database: (.+)$", line)
        if m:
            self._dr.database_name.append(m.group(1))
        elif self._dr.database_name:
            # This must be a continuation of the previous name.
            self._dr.database_name[-1] = "%s%s" % (self._dr.database_name[-1],
                                                   line.strip())

    def posted_date(self, line):
        self._dr.posted_date.append(_re_search(
            r"Posted date:\s*(.+)$", line,
            "I could not find the posted date in line\n%s" % line))

    def num_letters_in_database(self, line):
        letters, = _get_cols(
            line, (-1,), ncols=6, expected={2: "letters", 4: "database:"})
        self._dr.num_letters_in_database.append(_safe_int(letters))

    def num_sequences_in_database(self, line):
        sequences, = _get_cols(
            line, (-1,), ncols=6, expected={2: "sequences", 4: "database:"})
        self._dr.num_sequences_in_database.append(_safe_int(sequences))

    def ka_params(self, line):
        self._dr.ka_params = [_safe_float(x) for x in line.split()]

    def gapped(self, line):
        self._dr.gapped = 1

    def ka_params_gap(self, line):
        self._dr.ka_params_gap = [_safe_float(x) for x in line.split()]

    def end_database_report(self):
        pass


class _ParametersConsumer(object):
    def start_parameters(self):
        self._params = Record.Parameters()

    def matrix(self, line):
        self._params.matrix = line[8:].rstrip()

    def gap_penalties(self, line):
        self._params.gap_penalties = [_safe_float(x) for x in _get_cols(
            line, (3, 5), ncols=6, expected={2: "Existence:", 4: "Extension:"})]

    def num_hits(self, line):
        if '1st pass' in line:
            x, = _get_cols(line, (-4,), ncols=11, expected={2: "Hits"})
            self._params.num_hits = _safe_int(x)
        else:
            x, = _get_cols(line, (-1,), ncols=6, expected={2: "Hits"})
            self._params.num_hits = _safe_int(x)

    def num_sequences(self, line):
        if '1st pass' in line:
            x, = _get_cols(line, (-4,), ncols=9, expected={2: "Sequences:"})
            self._params.num_sequences = _safe_int(x)
        else:
            x, = _get_cols(line, (-1,), ncols=4, expected={2: "Sequences:"})
            self._params.num_sequences = _safe_int(x)

    def num_extends(self, line):
        if '1st pass' in line:
            x, = _get_cols(line, (-4,), ncols=9, expected={2: "extensions:"})
            self._params.num_extends = _safe_int(x)
        else:
            x, = _get_cols(line, (-1,), ncols=4, expected={2: "extensions:"})
            self._params.num_extends = _safe_int(x)

    def num_good_extends(self, line):
        if '1st pass' in line:
            x, = _get_cols(line, (-4,), ncols=10, expected={3: "extensions:"})
            self._params.num_good_extends = _safe_int(x)
        else:
            x, = _get_cols(line, (-1,), ncols=5, expected={3: "extensions:"})
            self._params.num_good_extends = _safe_int(x)

    def num_seqs_better_e(self, line):
        self._params.num_seqs_better_e, = _get_cols(
            line, (-1,), ncols=7, expected={2: "sequences"})
        self._params.num_seqs_better_e = _safe_int(
            self._params.num_seqs_better_e)

    def hsps_no_gap(self, line):
        self._params.hsps_no_gap, = _get_cols(
            line, (-1,), ncols=9, expected={3: "better", 7: "gapping:"})
        self._params.hsps_no_gap = _safe_int(self._params.hsps_no_gap)

    def hsps_prelim_gapped(self, line):
        self._params.hsps_prelim_gapped, = _get_cols(
            line, (-1,), ncols=9, expected={4: "gapped", 6: "prelim"})
        self._params.hsps_prelim_gapped = _safe_int(
            self._params.hsps_prelim_gapped)

    def hsps_prelim_gapped_attempted(self, line):
        self._params.hsps_prelim_gapped_attempted, = _get_cols(
            line, (-1,), ncols=10, expected={4: "attempted", 7: "prelim"})
        self._params.hsps_prelim_gapped_attempted = _safe_int(
            self._params.hsps_prelim_gapped_attempted)

    def hsps_gapped(self, line):
        self._params.hsps_gapped, = _get_cols(
            line, (-1,), ncols=6, expected={3: "gapped"})
        self._params.hsps_gapped = _safe_int(self._params.hsps_gapped)

    def query_length(self, line):
        self._params.query_length, = _get_cols(
            line.lower(), (-1,), ncols=4, expected={0: "length", 2: "query:"})
        self._params.query_length = _safe_int(self._params.query_length)

    def database_length(self, line):
        self._params.database_length, = _get_cols(
            line.lower(), (-1,), ncols=4, expected={0: "length", 2: "database:"})
        self._params.database_length = _safe_int(self._params.database_length)

    def effective_hsp_length(self, line):
        self._params.effective_hsp_length, = _get_cols(
            line, (-1,), ncols=4, expected={1: "HSP", 2: "length:"})
        self._params.effective_hsp_length = _safe_int(
            self._params.effective_hsp_length)

    def effective_query_length(self, line):
        self._params.effective_query_length, = _get_cols(
            line, (-1,), ncols=5, expected={1: "length", 3: "query:"})
        self._params.effective_query_length = _safe_int(
            self._params.effective_query_length)

    def effective_database_length(self, line):
        self._params.effective_database_length, = _get_cols(
            line.lower(), (-1,), ncols=5, expected={1: "length", 3: "database:"})
        self._params.effective_database_length = _safe_int(
            self._params.effective_database_length)

    def effective_search_space(self, line):
        self._params.effective_search_space, = _get_cols(
            line, (-1,), ncols=4, expected={1: "search"})
        self._params.effective_search_space = _safe_int(
            self._params.effective_search_space)

    def effective_search_space_used(self, line):
        self._params.effective_search_space_used, = _get_cols(
            line, (-1,), ncols=5, expected={1: "search", 3: "used:"})
        self._params.effective_search_space_used = _safe_int(
            self._params.effective_search_space_used)

    def frameshift(self, line):
        self._params.frameshift = _get_cols(
           line, (4, 5), ncols=6, expected={0: "frameshift", 2: "decay"})

    def threshold(self, line):
        if line[:2] == "T:":
            # Assume its an old stlye line like "T: 123"
            self._params.threshold, = _get_cols(
                line, (1,), ncols=2, expected={0: "T:"})
        elif line[:28] == "Neighboring words threshold:":
            self._params.threshold, = _get_cols(
                line, (3,), ncols=4, expected={0: "Neighboring", 1: "words", 2: "threshold:"})
        else:
            raise ValueError("Unrecognised threshold line:\n%s" % line)
        self._params.threshold = _safe_int(self._params.threshold)

    def window_size(self, line):
        if line[:2] == "A:":
            self._params.window_size, = _get_cols(
                line, (1,), ncols=2, expected={0: "A:"})
        elif line[:25] == "Window for multiple hits:":
            self._params.window_size, = _get_cols(
                line, (4,), ncols=5, expected={0: "Window", 2: "multiple", 3: "hits:"})
        else:
            raise ValueError("Unrecognised window size line:\n%s" % line)
        self._params.window_size = _safe_int(self._params.window_size)

    def dropoff_1st_pass(self, line):
        score, bits = _re_search(
            r"X1: (\d+) \(\s*([0-9,.]+) bits\)", line,
            "I could not find the dropoff in line\n%s" % line)
        self._params.dropoff_1st_pass = _safe_int(score), _safe_float(bits)

    def gap_x_dropoff(self, line):
        score, bits = _re_search(
            r"X2: (\d+) \(\s*([0-9,.]+) bits\)", line,
            "I could not find the gap dropoff in line\n%s" % line)
        self._params.gap_x_dropoff = _safe_int(score), _safe_float(bits)

    def gap_x_dropoff_final(self, line):
        score, bits = _re_search(
            r"X3: (\d+) \(\s*([0-9,.]+) bits\)", line,
            "I could not find the gap dropoff final in line\n%s" % line)
        self._params.gap_x_dropoff_final = _safe_int(score), _safe_float(bits)

    def gap_trigger(self, line):
        score, bits = _re_search(
            r"S1: (\d+) \(\s*([0-9,.]+) bits\)", line,
            "I could not find the gap trigger in line\n%s" % line)
        self._params.gap_trigger = _safe_int(score), _safe_float(bits)

    def blast_cutoff(self, line):
        score, bits = _re_search(
            r"S2: (\d+) \(\s*([0-9,.]+) bits\)", line,
            "I could not find the blast cutoff in line\n%s" % line)
        self._params.blast_cutoff = _safe_int(score), _safe_float(bits)

    def end_parameters(self):
        pass


class _BlastConsumer(AbstractConsumer,
                     _HeaderConsumer,
                     _DescriptionConsumer,
                     _AlignmentConsumer,
                     _HSPConsumer,
                     _DatabaseReportConsumer,
                     _ParametersConsumer
                     ):
    # This Consumer is inherits from many other consumer classes that handle
    # the actual dirty work.  An alternate way to do it is to create objects
    # of those classes and then delegate the parsing tasks to them in a
    # decorator-type pattern.  The disadvantage of that is that the method
    # names will need to be resolved in this classes.  However, using
    # a decorator will retain more control in this class (which may or
    # may not be a bad thing).  In addition, having each sub-consumer as
    # its own object prevents this object's dictionary from being cluttered
    # with members and reduces the chance of member collisions.
    def __init__(self):
        self.data = None

    def round(self, line):
        # Make sure nobody's trying to pass me PSI-BLAST data!
        raise ValueError("This consumer doesn't handle PSI-BLAST data")

    def start_header(self):
        self.data = Record.Blast()
        _HeaderConsumer.start_header(self)

    def end_header(self):
        _HeaderConsumer.end_header(self)
        self.data.__dict__.update(self._header.__dict__)

    def end_descriptions(self):
        self.data.descriptions = self._descriptions

    def end_alignment(self):
        _AlignmentConsumer.end_alignment(self)
        if self._alignment.hsps:
            self.data.alignments.append(self._alignment)
        if self._multiple_alignment.alignment:
            self.data.multiple_alignment = self._multiple_alignment

    def end_hsp(self):
        _HSPConsumer.end_hsp(self)
        try:
            self._alignment.hsps.append(self._hsp)
        except AttributeError:
            raise ValueError("Found an HSP before an alignment")

    def end_database_report(self):
        _DatabaseReportConsumer.end_database_report(self)
        self.data.__dict__.update(self._dr.__dict__)

    def end_parameters(self):
        _ParametersConsumer.end_parameters(self)
        self.data.__dict__.update(self._params.__dict__)


class _PSIBlastConsumer(AbstractConsumer,
                        _HeaderConsumer,
                        _DescriptionConsumer,
                        _AlignmentConsumer,
                        _HSPConsumer,
                        _DatabaseReportConsumer,
                        _ParametersConsumer
                        ):
    def __init__(self):
        self.data = None

    def start_header(self):
        self.data = Record.PSIBlast()
        _HeaderConsumer.start_header(self)

    def end_header(self):
        _HeaderConsumer.end_header(self)
        self.data.__dict__.update(self._header.__dict__)

    def start_descriptions(self):
        self._round = Record.Round()
        self.data.rounds.append(self._round)
        _DescriptionConsumer.start_descriptions(self)

    def end_descriptions(self):
        _DescriptionConsumer.end_descriptions(self)
        self._round.number = self._roundnum
        if self._descriptions:
            self._round.new_seqs.extend(self._descriptions)
        self._round.reused_seqs.extend(self._model_sequences)
        self._round.new_seqs.extend(self._nonmodel_sequences)
        if self._converged:
            self.data.converged = 1

    def end_alignment(self):
        _AlignmentConsumer.end_alignment(self)
        if self._alignment.hsps:
            self._round.alignments.append(self._alignment)
        if self._multiple_alignment:
            self._round.multiple_alignment = self._multiple_alignment

    def end_hsp(self):
        _HSPConsumer.end_hsp(self)
        try:
            self._alignment.hsps.append(self._hsp)
        except AttributeError:
            raise ValueError("Found an HSP before an alignment")

    def end_database_report(self):
        _DatabaseReportConsumer.end_database_report(self)
        self.data.__dict__.update(self._dr.__dict__)

    def end_parameters(self):
        _ParametersConsumer.end_parameters(self)
        self.data.__dict__.update(self._params.__dict__)


class Iterator(object):
    """Iterates over a file of multiple BLAST results.

    Methods:
    next   Return the next record from the stream, or None.

    """
    def __init__(self, handle, parser=None):
        """__init__(self, handle, parser=None)

        Create a new iterator.  handle is a file-like object.  parser
        is an optional Parser object to change the results into another form.
        If set to None, then the raw contents of the file will be returned.

        """
        try:
            handle.readline
        except AttributeError:
            raise ValueError(
                "I expected a file handle or file-like object, got %s"
                % type(handle))
        self._uhandle = File.UndoHandle(handle)
        self._parser = parser
        self._header = []

    def __next__(self):
        """next(self) -> object

        Return the next Blast record from the file.  If no more records,
        return None.

        """
        lines = []
        query = False
        while True:
            line = self._uhandle.readline()
            if not line:
                break
            # If I've reached the next one, then put the line back and stop.
            if lines and (line.startswith('BLAST')
                          or line.startswith('BLAST', 1)
                          or line.startswith('<?xml ')):
                self._uhandle.saveline(line)
                break
            # New style files omit the BLAST line to mark a new query:
            if line.startswith("Query="):
                if not query:
                    if not self._header:
                        self._header = lines[:]
                    query = True
                else:
                    # Start of another record
                    self._uhandle.saveline(line)
                    break
            lines.append(line)

        if query and "BLAST" not in lines[0]:
            # Cheat and re-insert the header
            # print "-"*50
            # print "".join(self._header)
            # print "-"*50
            # print "".join(lines)
            # print "-"*50
            lines = self._header + lines

        if not lines:
            return None

        data = ''.join(lines)
        if self._parser is not None:
            return self._parser.parse(StringIO(data))
        return data

    if sys.version_info[0] < 3:
        def next(self):
            """Python 2 style alias for Python 3 style __next__ method."""
            return self.__next__()

    def __iter__(self):
        return iter(self.__next__, None)


def _re_search(regex, line, error_msg):
    m = re.search(regex, line)
    if not m:
        raise ValueError(error_msg)
    return m.groups()


def _get_cols(line, cols_to_get, ncols=None, expected={}):
    cols = line.split()

    # Check to make sure number of columns is correct
    if ncols is not None and len(cols) != ncols:
        raise ValueError("I expected %d columns (got %d) in line\n%s"
                         % (ncols, len(cols), line))

    # Check to make sure columns contain the correct data
    for k in expected:
        if cols[k] != expected[k]:
            raise ValueError("I expected '%s' in column %d in line\n%s"
                             % (expected[k], k, line))

    # Construct the answer tuple
    results = []
    for c in cols_to_get:
        results.append(cols[c])
    return tuple(results)


def _safe_int(str):
    try:
        return int(str)
    except ValueError:
        # Something went wrong.  Try to clean up the string.
        # Remove all commas from the string
        str = str.replace(',', '')
    # try again after removing commas.
    # Note int() will return a long rather than overflow
    try:
        return int(str)
    except ValueError:
        pass
    # Call float to handle things like "54.3", note could lose precision, e.g.
    # >>> int("5399354557888517312")
    # 5399354557888517312
    # >>> int(float("5399354557888517312"))
    # 5399354557888517120
    return int(float(str))


def _safe_float(str):
    # Thomas Rosleff Soerensen (rosleff@mpiz-koeln.mpg.de) noted that
    # float('e-172') does not produce an error on his platform.  Thus,
    # we need to check the string for this condition.

    # Sometimes BLAST leaves of the '1' in front of an exponent.
    if str and str[0] in ['E', 'e']:
        str = '1' + str
    try:
        return float(str)
    except ValueError:
        # Remove all commas from the string
        str = str.replace(',', '')
    # try again.
    return float(str)


class _BlastErrorConsumer(_BlastConsumer):
    def __init__(self):
        _BlastConsumer.__init__(self)

    def noevent(self, line):
        if 'Query must be at least wordsize' in line:
            raise ShortQueryBlastError("Query must be at least wordsize")
        # Now pass the line back up to the superclass.
        method = getattr(_BlastConsumer, 'noevent',
                         _BlastConsumer.__getattr__(self, 'noevent'))
        method(line)


class BlastErrorParser(AbstractParser):
    """Attempt to catch and diagnose BLAST errors while parsing.

    This utilizes the BlastParser module but adds an additional layer
    of complexity on top of it by attempting to diagnose ValueErrors
    that may actually indicate problems during BLAST parsing.

    Current BLAST problems this detects are:
    o LowQualityBlastError - When BLASTing really low quality sequences
    (ie. some GenBank entries which are just short streches of a single
    nucleotide), BLAST will report an error with the sequence and be
    unable to search with this. This will lead to a badly formatted
    BLAST report that the parsers choke on. The parser will convert the
    ValueError to a LowQualityBlastError and attempt to provide useful
    information.
    """
    def __init__(self, bad_report_handle=None):
        """Initialize a parser that tries to catch BlastErrors.

        Arguments:
        o bad_report_handle - An optional argument specifying a handle
        where bad reports should be sent. This would allow you to save
        all of the bad reports to a file, for instance. If no handle
        is specified, the bad reports will not be saved.
        """
        self._bad_report_handle = bad_report_handle

        # self._b_parser = BlastParser()
        self._scanner = _Scanner()
        self._consumer = _BlastErrorConsumer()

    def parse(self, handle):
        """Parse a handle, attempting to diagnose errors.
        """
        results = handle.read()

        try:
            self._scanner.feed(StringIO(results), self._consumer)
        except ValueError:
            # if we have a bad_report_file, save the info to it first
            if self._bad_report_handle:
                # send the info to the error handle
                self._bad_report_handle.write(results)

            # now we want to try and diagnose the error
            self._diagnose_error(
                StringIO(results), self._consumer.data)

            # if we got here we can't figure out the problem
            # so we should pass along the syntax error we got
            raise
        return self._consumer.data

    def _diagnose_error(self, handle, data_record):
        """Attempt to diagnose an error in the passed handle.

        Arguments:
        o handle - The handle potentially containing the error
        o data_record - The data record partially created by the consumer.
        """
        line = handle.readline()

        while line:
            # 'Searchingdone' instead of 'Searching......done' seems
            # to indicate a failure to perform the BLAST due to
            # low quality sequence
            if line.startswith('Searchingdone'):
                raise LowQualityBlastError("Blast failure occurred on query: ",
                                           data_record.query)
            line = handle.readline()
