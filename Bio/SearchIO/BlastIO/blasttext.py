# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for BLAST+ plain text output formats.

At the moment this is a wrapper around Biopython's NCBIStandalone text
parser.

"""

from Bio.Blast import NCBIStandalone
from Bio.SearchIO._objects import QueryResult, Hit, HSP, HSPFragment


class BlastTextParser(object):

    """Parser for the BLAST text format."""

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

                # iterate over the hsps and group them in a list
                hit_list = []
                for bhsp in aln.hsps:
                    frag = HSPFragment(hid, qid)
                    # set alignment length
                    frag.aln_span = bhsp.identities[1]
                    # set frames
                    try:
                        frag.query_frame = int(bhsp.frame[0])
                    except IndexError:
                        if qresult.program in ('blastp', 'tblastn'):
                            frag.query_frame = 0
                        else:
                            frag.query_frame = 1
                    try:
                        frag.hit_frame = int(bhsp.frame[1])
                    except IndexError:
                        if qresult.program in ('blastp', 'tblastn'):
                            frag.hit_frame = 0
                        else:
                            frag.hit_frame = 1
                    # set query coordinates
                    startfunc = min if frag.query_strand >= 0 else max
                    endfunc = max if frag.query_strand >= 0 else min
                    frag.query_start = startfunc(bhsp.query_start, \
                            bhsp.query_end) - 1
                    frag.query_end = endfunc(bhsp.query_start, bhsp.query_end)
                    # set hit coordinates
                    startfunc = min if frag.hit_strand >= 0 else max
                    endfunc = max if frag.hit_strand >= 0 else min
                    frag.hit_start = startfunc(bhsp.sbjct_start, \
                            bhsp.sbjct_end) - 1
                    frag.hit_end = endfunc(bhsp.sbjct_start, bhsp.sbjct_end)
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
                    frag.query, frag.hit = qseq, hseq
                    frag.alignment_annotation['homology'] = midline

                    # create HSP object with the fragment
                    hsp = HSP([frag])
                    hsp.evalue = bhsp.expect
                    hsp.bitscore = bhsp.bits
                    hsp.bitscore_raw = bhsp.score
                    # set gap
                    try:
                        hsp.gap_num = bhsp.gaps[0]
                    except IndexError:
                        hsp.gap_num = 0
                    # set identity
                    hsp.ident_num = bhsp.identities[0]
                    hsp.pos_num = bhsp.positives[0]
                    if hsp.pos_num is None:
                        hsp.pos_num = hsp[0].aln_span

                    hit_list.append(hsp)

                hit = Hit(hid, qid, hit_list)
                hit.seq_len = aln.length
                # get hit-level e-value and score, if available
                try:
                    hit.evalue = rec.descriptions[idx].e
                    hit.score = rec.descriptions[idx].score
                except IndexError:
                    pass

                hit.description = hdesc
                qresult.append(hit)

            qresult.description = qdesc
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
