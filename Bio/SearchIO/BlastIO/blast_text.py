# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.SearchIO parser for BLAST+ plain text output formats.

At the moment this is a wrapper around Biopython's NCBIStandalone text
parser (which is now deprecated).

"""

from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.SearchIO._legacy import NCBIStandalone


__all__ = ("BlastTextParser",)


class BlastTextParser:
    """Parser for the BLAST text format."""

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        blast_parser = NCBIStandalone.BlastParser()
        self.blast_iter = NCBIStandalone.Iterator(handle, blast_parser)

    def __iter__(self):
        """Iterate over BlastTextParser, yields query results."""
        for rec in self.blast_iter:
            # set attributes to SearchIO's
            # get id and desc
            if rec.query.startswith(">"):
                rec.query = rec.query[1:]
            try:
                qid, qdesc = rec.query.split(" ", 1)
            except ValueError:
                qid, qdesc = rec.query, ""
            qdesc = qdesc.replace("\n", "").replace("\r", "")

            qresult = QueryResult(id=qid)
            qresult.program = rec.application.lower()
            qresult.target = rec.database
            qresult.seq_len = rec.query_letters
            qresult.version = rec.version

            # determine molecule_type based on program
            if qresult.program == "blastn":
                molecule_type = "DNA"
            elif qresult.program in ["blastp", "blastx", "tblastn", "tblastx"]:
                molecule_type = "protein"

            # iterate over the 'alignments' (hits) and the hit table
            for idx, aln in enumerate(rec.alignments):
                # get id and desc
                if aln.title.startswith("> "):
                    aln.title = aln.title[2:]
                elif aln.title.startswith(">"):
                    aln.title = aln.title[1:]
                try:
                    hid, hdesc = aln.title.split(" ", 1)
                except ValueError:
                    hid, hdesc = aln.title, ""
                hdesc = hdesc.replace("\n", "").replace("\r", "")

                # iterate over the hsps and group them in a list
                hsp_list = []
                for bhsp in aln.hsps:
                    frag = HSPFragment(hid, qid)
                    frag.molecule_type = molecule_type
                    # set alignment length
                    frag.aln_span = bhsp.identities[1]
                    # set frames
                    try:
                        frag.query_frame = int(bhsp.frame[0])
                    except IndexError:
                        if qresult.program in ("blastp", "tblastn"):
                            frag.query_frame = 0
                        else:
                            frag.query_frame = 1
                    try:
                        frag.hit_frame = int(bhsp.frame[1])
                    except IndexError:
                        if qresult.program in ("blastp", "tblastn"):
                            frag.hit_frame = 0
                        else:
                            frag.hit_frame = 1
                    # set query coordinates
                    frag.query_start = min(bhsp.query_start, bhsp.query_end) - 1
                    frag.query_end = max(bhsp.query_start, bhsp.query_end)
                    # set hit coordinates
                    frag.hit_start = min(bhsp.sbjct_start, bhsp.sbjct_end) - 1
                    frag.hit_end = max(bhsp.sbjct_start, bhsp.sbjct_end)
                    # set query, hit sequences and its annotation
                    qseq = ""
                    hseq = ""
                    midline = ""
                    for seqtrio in zip(bhsp.query, bhsp.sbjct, bhsp.match):
                        qchar, hchar, mchar = seqtrio
                        if qchar == " " or hchar == " ":
                            assert all(" " == x for x in seqtrio)
                        else:
                            qseq += qchar
                            hseq += hchar
                            midline += mchar
                    frag.query, frag.hit = qseq, hseq
                    frag.aln_annotation["similarity"] = midline

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

                    hsp_list.append(hsp)

                hit = Hit(hsp_list)
                hit.seq_len = aln.length
                hit.description = hdesc
                qresult.append(hit)

            qresult.description = qdesc
            yield qresult


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
