# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ plain text output formats.

At the moment this is a wrapper around Biopython's NCBIStandalone text
parser.

"""

from Bio.Blast import NCBIStandalone
from Bio.SearchIO._objects import QueryResult, Hit, HSP


class BlastTextIterator(object):

    """Parser for the Blast tabular format."""

    def __init__(self, handle):
        self.handle = handle
        blast_parser = NCBIStandalone.BlastParser()
        self.blast_iter = NCBIStandalone.Iterator(handle, blast_parser)

    def __iter__(self):
        for rec in self.blast_iter:
            # set attributes to SearchIO's
            # get id and desc
            if rec.query.startswith('>'):
                rec.query = rec.query[1:]
            try:
                qid, qdesc = rec.query.split(' ', 1)
            except ValueError:
                qid, qdesc = rec.query, ''
            qdesc = qdesc.replace('\n', '').replace('\r', '')

            qresult = QueryResult(qid)
            qresult.program = rec.application.lower()
            qresult.target = rec.database
            qresult.seq_len = rec.query_letters
            qresult.version = rec.version

            # iterate over the 'alignments' (hits) and the hit table
            for idx, aln in enumerate(rec.alignments):
                # get id and desc
                if aln.title.startswith('> '):
                    aln.title = aln.title[2:]
                elif aln.title.startswith('>'):
                    aln.title = aln.title[1:]
                try:
                    hid, hdesc = aln.title.split(' ', 1)
                except ValueError:
                    hid, hdesc = aln.title, ''
                hdesc = hdesc.replace('\n', '').replace('\r', '')

                hit = Hit(hid, qid)
                hit.seq_len = aln.length
                # get hit-level e-value and score, if available
                try:
                    hit.evalue = rec.descriptions[idx].e
                    hit.score = rec.descriptions[idx].score
                except IndexError:
                    pass

                # iterate over the hsps and group them in the current hit
                for bhsp in aln.hsps:
                    hsp = HSP(hid, qid)
                    hsp.evalue = bhsp.expect
                    hsp.bitscore = bhsp.bits
                    hsp.bitscore_raw = bhsp.score
                    # set alignment length
                    hsp.aln_span = bhsp.identities[1]
                    # set frames
                    try:
                        hsp.query_frame = int(bhsp.frame[0])
                    except IndexError:
                        if qresult.program in ('blastp', 'tblastn'):
                            hsp.query_frame = 0
                        else:
                            hsp.query_frame = 1
                    try:
                        hsp.hit_frame = int(bhsp.frame[1])
                    except IndexError:
                        if qresult.program in ('blastp', 'tblastn'):
                            hsp.hit_frame = 0
                        else:
                            hsp.hit_frame = 1
                    # set gap
                    try:
                        hsp.gap_num = bhsp.gaps[0]
                    except IndexError:
                        hsp.gap_num = 0
                    # set identity
                    hsp.ident_num = bhsp.identities[0]
                    hsp.pos_num = bhsp.positives[0]
                    if hsp.pos_num is None:
                        hsp.pos_num = hsp.aln_span
                    # set query coordinates
                    startfunc = min if hsp.query_strand >= 0 else max
                    endfunc = max if hsp.query_strand >= 0 else min
                    hsp.query_start = startfunc(bhsp.query_start, \
                            bhsp.query_end) - 1
                    hsp.query_end = endfunc(bhsp.query_start, bhsp.query_end)
                    # set hit coordinates
                    startfunc = min if hsp.hit_strand >= 0 else max
                    endfunc = max if hsp.hit_strand >= 0 else min
                    hsp.hit_start = startfunc(bhsp.sbjct_start, \
                            bhsp.sbjct_end) - 1
                    hsp.hit_end = endfunc(bhsp.sbjct_start, bhsp.sbjct_end)

                    # set query, hit sequences and its annotation
                    qseq = ''
                    hseq = ''
                    midline = ''
                    for seqtrio in zip(bhsp.query, bhsp.sbjct, bhsp.match):
                        qchar, hchar, mchar = seqtrio
                        if qchar == ' ' or hchar == ' ':
                            assert all([' ' == x for x in seqtrio])
                        else:
                            qseq += qchar
                            hseq += hchar
                            midline += mchar
                    hsp.query, hsp.hit = qseq, hseq
                    hsp.alignment_annotation['homology'] = midline

                    hit.append(hsp)

                hit.desc = hdesc
                qresult.append(hit)

            qresult.desc = qdesc
            yield qresult


def _test():
    """Run the Bio.SearchIO.BlastIO module's doctests.

    This will try and locate the unit tests directory, and run the doctests
    from there in order that the relative paths used in the examples work.
    """
    import doctest
    import os

    test_dir = 'Tests'

    if os.path.isdir(os.path.join('..', '..', test_dir)):
        print "Runing doctests..."
        cur_dir = os.path.abspath(os.curdir)
        os.chdir(os.path.join('..', '..', test_dir))
        doctest.testmod()
        os.chdir(cur_dir)
        print "Done"


# if not used as a module, run the doctest
if __name__ == "__main__":
    _test()
