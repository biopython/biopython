# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO parser for BLAT output formats.

This module adds support for parsing BLAT outputs. BLAT (BLAST-Like Alignment
Tool) is a sequence similarity search program initially built for annotating
the human genome.

Bio.SearchIO.BlastIO was tested using standalone BLAT version 34, psLayout
version 3. It should be able to parse psLayout version 4 without problems.

More information on BLAT is available from these sites:

    - Publication: http://genome.cshlp.org/content/12/4/656
    - User guide: http://genome.ucsc.edu/goldenPath/help/blatSpec.html
    - Source download: http://www.soe.ucsc.edu/~kent/src
    - Executable download: http://hgdownload.cse.ucsc.edu/admin/exe/
    - Blat score calculation: http://genome.ucsc.edu/FAQ/FAQblat.html#blat4


Supported Formats
=================

BlatIO supports parsing, indexing, and writing for both PSL and PSLX output
formats, with or without header. To parse, index, or write PSLX files, use the
'pslx' keyword argument and set it to True.

    # blat-psl defaults to PSL files
    >>> from Bio import SearchIO
    >>> psl = 'Blat/psl_34_004.psl'
    >>> qresult = SearchIO.read(psl, 'blat-psl')
    >>> qresult
    QueryResult(id='hg19_dna', 10 hits)

    # set the pslx flag to parse PSLX files
    >>> pslx = 'Blat/pslx_34_004.pslx'
    >>> qresult = SearchIO.read(pslx, 'blat-psl', pslx=True)
    >>> qresult
    QueryResult(id='hg19_dna', 10 hits)

For parsing and indexing, you do not need to specify whether the file has a
header or not. For writing, if you want to write a header, you can set the
'header' keyword argument to True. This will write a 'psLayout version 3' header
to your output file.

    from Bio import SearchIO
    qresult = SearchIO.read(psl, 'blat-psl')
    SearchIO.write(qresult, 'header.psl', header=True)
    <stdout> (1, 10, 19, 23)

Note that the number of HSPFragments written may exceed the number of HSP
objects. This is because in PSL files, it is possible to have single matches
consisting of noncontiguous sequence fragments. This is where the HSPFragment
object comes into play. These fragments are grouped into a single HSP because
they share the same statistics (e.g. match numbers, BLAT score, etc.). However,
they do not share the same sequence attributes, such as the start and end
coordinates, making them distinct objects.

In addition to parsing PSL(X) files, BlatIO also computes the percent identities
and scores of your search results. This is done using the calculation formula
posted here: http://genome.ucsc.edu/FAQ/FAQblat.html#blat4. It mimics the score
and percent identity calculation done by UCSC's web BLAT service.

Since BlatIO parses the file in a single pass, it expects all results from
the same query to be in consecutive rows. If the results from one query are
spread in nonconsecutive rows, BlatIO will consider them to be separate
QueryResult objects.

In most cases, the PSL(X) format uses the same coordinate system as Python
(zero-based, half open). These coordinates are anchored on the plus strand.
However, if the query aligns on the minus strand, BLAT will anchor the qStarts
coordinates on the minus strand instead. BlatIO is aware of this, and will
re-anchor the qStarts coordinates to the plus strand whenever it sees a minus
strand query match. Conversely, when you write out to a PSL(X) file, BlatIO will
reanchor qStarts to the minus strand again.

BlatIO provides the following attribute-column mapping:

+----------------+-------------------------+-----------------------------------+
| Object         | Attribute               | Column Name, Value                |
+================+=========================+===================================+
| QueryResult    | id                      | Q name, query sequence ID         |
|                +-------------------------+-----------------------------------+
|                | seq_len                 | Q size, query sequence full       |
|                |                         | length                            |
+----------------+-------------------------+-----------------------------------+
| Hit            | id                      | T name, hit sequence ID           |
|                +-------------------------+-----------------------------------+
|                | seq_len                 | T size, hit sequence full length  |
+----------------+-------------------------+-----------------------------------+
| HSP            | hit_end                 | T end, end coordinate of the last |
|                |                         | hit fragment                      |
|                +-------------------------+-----------------------------------+
|                | hit_gap_num             | T gap bases, number of bases      |
|                |                         | inserted in hit                   |
|                +-------------------------+-----------------------------------+
|                | hit_gapopen_num         | T gap count, number of hit gap    |
|                |                         | inserts                           |
|                +-------------------------+-----------------------------------+
|                | hit_span_all            | blockSizes, sizes of each         |
|                |                         | fragment                          |
|                +-------------------------+-----------------------------------+
|                | hit_start               | T start, start coordinate of the  |
|                |                         | first hit fragment                |
|                +-------------------------+-----------------------------------+
|                | hit_start_all           | tStarts, start coordinate of each |
|                |                         | hit fragment                      |
|                +-------------------------+-----------------------------------+
|                | match_num               | match, number of non-repeat       |
|                |                         | matches                           |
|                +-------------------------+-----------------------------------+
|                | mismatch_num            | mismatch, number of mismatches    |
|                +-------------------------+-----------------------------------+
|                | match_rep_num           | rep. match, number of matches     |
|                |                         | that are part of repeats          |
|                +-------------------------+-----------------------------------+
|                | n_num                   | N's, number of N bases            |
|                +-------------------------+-----------------------------------+
|                | query_end               | Q end, end coordinate of the last |
|                +-------------------------+-----------------------------------+
|                |                         | query fragment                    |
|                | query_gap_num           | Q gap bases, number of bases      |
|                |                         | inserted in query                 |
|                +-------------------------+-----------------------------------+
|                | query_gapopen_num       | Q gap count, number of query gap  |
|                |                         | inserts                           |
|                +-------------------------+-----------------------------------+
|                | query_span_all          | blockSizes, sizes of each         |
|                |                         | fragment                          |
|                +-------------------------+-----------------------------------+
|                | query_start             | Q start, start coordinate of the  |
|                |                         | first query block                 |
|                +-------------------------+-----------------------------------+
|                | query_start_all         | qStarts, start coordinate of each |
|                |                         | query fragment                    |
|                +-------------------------+-----------------------------------+
|                | len [*]_                | block count, the number of blocks |
|                |                         | in the alignment                  |
+----------------+-------------------------+-----------------------------------+
| HSPFragment    | hit                     | hit sequence, if present          |
|                +-------------------------+-----------------------------------+
|                | hit_strand              | strand, hit sequence strand       |
|                +-------------------------+-----------------------------------+
|                | query                   | query sequence, if present        |
|                +-------------------------+-----------------------------------+
|                | query_strand            | strand, query sequence strand     |
+----------------+-------------------------+-----------------------------------+

In addition to the column mappings above, BlatIO also provides the following
object attributes:

+----------------+-------------------------+-----------------------------------+
| Object         | Attribute               | Value                             |
+================+=========================+===================================+
| HSP            | gapopen_num             | Q gap count + T gap count, total  |
|                |                         |  number of gap openings           |
|                +-------------------------+-----------------------------------+
|                | ident_num               | matches + repmatches, total       |
|                |                         | number of identical residues      |
|                +-------------------------+-----------------------------------+
|                | ident_pct               | percent identity, calculated      |
|                |                         | using UCSC's formula              |
|                +-------------------------+-----------------------------------+
|                | query_is_protein        | boolean, whether the query        |
|                |                         | sequence is a protein             |
|                +-------------------------+-----------------------------------+
|                | score                   | HSP score, calculated using       |
|                |                         | UCSC's formula                    |
+----------------+-------------------------+-----------------------------------+

Finally, the default HSP and HSPFragment properties are also provided. See the
HSP and HSPFragment documentation for more details on these properties.


.. [*] You can obtain the number of blocks / fragments in the HSP by invoking
   ``len`` on the HSP

"""
import re
from math import log

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment


__all__ = ("BlatPslParser", "BlatPslIndexer", "BlatPslWriter")


# precompile regex patterns
_PTR_ROW_CHECK = r"^\d+\s+\d+\s+\d+\s+\d+"
_RE_ROW_CHECK = re.compile(_PTR_ROW_CHECK)
_RE_ROW_CHECK_IDX = re.compile(_PTR_ROW_CHECK.encode())


def _list_from_csv(csv_string, caster=None):
    """Transform the given comma-separated string into a list (PRIVATE).

    :param csv_string: comma-separated input string
    :type csv_string: string
    :param caster: function used to cast each item in the input string
                   to its intended type
    :type caster: callable, accepts string, returns object

    """
    if caster is None:
        return [x for x in csv_string.split(",") if x]
    else:
        return [caster(x) for x in csv_string.split(",") if x]


def _reorient_starts(starts, blksizes, seqlen, strand):
    """Reorients block starts into the opposite strand's coordinates (PRIVATE).

    :param starts: start coordinates
    :type starts: list [int]
    :param blksizes: block sizes
    :type blksizes: list [int]
    :param seqlen: sequence length
    :type seqlen: int
    :param strand: sequence strand
    :type strand: int, choice of -1, 0, or 1

    """
    if len(starts) != len(blksizes):
        raise RuntimeError(
            "Unequal start coordinates and block sizes list (%r vs %r)"
            % (len(starts), len(blksizes))
        )
    # see: http://genome.ucsc.edu/goldenPath/help/blatSpec.html
    # no need to reorient if it's already the positive strand
    if strand >= 0:
        return starts
    else:
        # the plus-oriented coordinate is calculated by this:
        # plus_coord = length - minus_coord - block_size
        return [seqlen - start - blksize for start, blksize in zip(starts, blksizes)]


def _is_protein(psl):
    """Validate if psl is protein (PRIVATE)."""
    # check if query is protein or not
    # adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4
    if len(psl["strand"]) == 2:
        if psl["strand"][1] == "+":
            return psl["tend"] == psl["tstarts"][-1] + 3 * psl["blocksizes"][-1]
        elif psl["strand"][1] == "-":
            return psl["tstart"] == psl["tsize"] - (
                psl["tstarts"][-1] + 3 * psl["blocksizes"][-1]
            )

    return False


def _calc_millibad(psl, is_protein):
    """Calculate millibad (PRIVATE)."""
    # adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4
    size_mul = 3 if is_protein else 1
    millibad = 0

    qali_size = size_mul * (psl["qend"] - psl["qstart"])
    tali_size = psl["tend"] - psl["tstart"]
    ali_size = min(qali_size, tali_size)
    if ali_size <= 0:
        return 0

    size_dif = qali_size - tali_size
    size_dif = 0 if size_dif < 0 else size_dif

    total = size_mul * (psl["matches"] + psl["repmatches"] + psl["mismatches"])
    if total != 0:
        millibad = (
            1000
            * (
                psl["mismatches"] * size_mul
                + psl["qnuminsert"]
                + round(3 * log(1 + size_dif))
            )
        ) / total

    return millibad


def _calc_score(psl, is_protein):
    """Calculate score (PRIVATE)."""
    # adapted from http://genome.ucsc.edu/FAQ/FAQblat.html#blat4
    size_mul = 3 if is_protein else 1
    return (
        size_mul * (psl["matches"] + (psl["repmatches"] >> 1))
        - size_mul * psl["mismatches"]
        - psl["qnuminsert"]
        - psl["tnuminsert"]
    )


def _create_hsp(hid, qid, psl):
    """Create high scoring pair object (PRIVATE)."""
    # protein flag
    is_protein = _is_protein(psl)
    # strand
    # if query is protein, strand is 0
    if is_protein:
        qstrand = 0
    else:
        qstrand = 1 if psl["strand"][0] == "+" else -1
    # try to get hit strand, if it exists
    try:
        hstrand = 1 if psl["strand"][1] == "+" else -1
    except IndexError:
        hstrand = 1  # hit strand defaults to plus

    blocksize_multiplier = 3 if is_protein else 1
    # query block starts
    qstarts = _reorient_starts(psl["qstarts"], psl["blocksizes"], psl["qsize"], qstrand)
    # hit block starts
    if len(psl["strand"]) == 2:
        hstarts = _reorient_starts(
            psl["tstarts"],
            [blocksize_multiplier * i for i in psl["blocksizes"]],
            psl["tsize"],
            hstrand,
        )
    else:
        hstarts = psl["tstarts"]
    # set query and hit coords
    # this assumes each block has no gaps (which seems to be the case)
    assert len(qstarts) == len(hstarts) == len(psl["blocksizes"])
    query_range_all = list(
        zip(qstarts, [x + y for x, y in zip(qstarts, psl["blocksizes"])])
    )
    hit_range_all = list(
        zip(
            hstarts,
            [x + y * blocksize_multiplier for x, y in zip(hstarts, psl["blocksizes"])],
        )
    )
    # check length of sequences and coordinates, all must match
    if "tseqs" in psl and "qseqs" in psl:
        assert (
            len(psl["tseqs"])
            == len(psl["qseqs"])
            == len(query_range_all)
            == len(hit_range_all)
        )
    else:
        assert len(query_range_all) == len(hit_range_all)

    frags = []
    # iterating over query_range_all, but hit_range_all works just as well
    for idx, qcoords in enumerate(query_range_all):
        hseqlist = psl.get("tseqs")
        hseq = "" if not hseqlist else hseqlist[idx]
        qseqlist = psl.get("qseqs")
        qseq = "" if not qseqlist else qseqlist[idx]
        frag = HSPFragment(hid, qid, hit=hseq, query=qseq)
        # set molecule type
        frag.molecule_type = "DNA"
        # set coordinates
        frag.query_start = qcoords[0]
        frag.query_end = qcoords[1]
        frag.hit_start = hit_range_all[idx][0]
        frag.hit_end = hit_range_all[idx][1]
        # and strands
        frag.query_strand = qstrand
        frag.hit_strand = hstrand
        frags.append(frag)

    # create hsp object
    hsp = HSP(frags)
    # check if start and end are set correctly
    assert hsp.query_start == psl["qstart"]
    assert hsp.query_end == psl["qend"]
    assert hsp.hit_start == psl["tstart"]
    assert hsp.hit_end == psl["tend"]
    # and check block spans as well
    hit_spans = [span / blocksize_multiplier for span in hsp.hit_span_all]
    assert hit_spans == hsp.query_span_all == psl["blocksizes"]
    # set its attributes
    hsp.match_num = psl["matches"]
    hsp.mismatch_num = psl["mismatches"]
    hsp.match_rep_num = psl["repmatches"]
    hsp.n_num = psl["ncount"]
    hsp.query_gapopen_num = psl["qnuminsert"]
    hsp.query_gap_num = psl["qbaseinsert"]
    hsp.hit_gapopen_num = psl["tnuminsert"]
    hsp.hit_gap_num = psl["tbaseinsert"]

    hsp.ident_num = psl["matches"] + psl["repmatches"]
    hsp.gapopen_num = psl["qnuminsert"] + psl["tnuminsert"]
    hsp.gap_num = psl["qbaseinsert"] + psl["tbaseinsert"]
    hsp.query_is_protein = is_protein
    hsp.ident_pct = 100.0 - _calc_millibad(psl, is_protein) * 0.1
    hsp.score = _calc_score(psl, is_protein)
    # helper flag, for writing
    hsp._has_hit_strand = len(psl["strand"]) == 2

    return hsp


class BlatPslParser:
    """Parser for the BLAT PSL format."""

    def __init__(self, handle, pslx=False):
        """Initialize the class."""
        self.handle = handle
        self.line = self.handle.readline()
        self.pslx = pslx

    def __iter__(self):
        """Iterate over BlatPslParser, yields query results."""
        # break out if it's an empty file
        if not self.line:
            return

        # read through header
        # this assumes that the result row match the regex
        while not re.search(_RE_ROW_CHECK, self.line.strip()):
            self.line = self.handle.readline()
            if not self.line:
                return

        # parse into query results
        for qresult in self._parse_qresult():
            qresult.program = "blat"
            yield qresult

    def _parse_row(self):
        """Return a dictionary of parsed column values (PRIVATE)."""
        assert self.line
        cols = [x for x in self.line.strip().split("\t") if x]
        self._validate_cols(cols)

        psl = {}
        psl["qname"] = cols[9]  # qName
        psl["qsize"] = int(cols[10])  # qSize
        psl["tname"] = cols[13]  # tName
        psl["tsize"] = int(cols[14])  # tSize
        psl["matches"] = int(cols[0])  # matches
        psl["mismatches"] = int(cols[1])  # misMatches
        psl["repmatches"] = int(cols[2])  # repMatches
        psl["ncount"] = int(cols[3])  # nCount
        psl["qnuminsert"] = int(cols[4])  # qNumInsert
        psl["qbaseinsert"] = int(cols[5])  # qBaseInsert
        psl["tnuminsert"] = int(cols[6])  # tNumInsert
        psl["tbaseinsert"] = int(cols[7])  # tBaseInsert
        psl["strand"] = cols[8]  # strand
        psl["qstart"] = int(cols[11])  # qStart
        psl["qend"] = int(cols[12])  # qEnd
        psl["tstart"] = int(cols[15])  # tStart
        psl["tend"] = int(cols[16])  # tEnd
        psl["blockcount"] = int(cols[17])  # blockCount
        psl["blocksizes"] = _list_from_csv(cols[18], int)  # blockSizes
        psl["qstarts"] = _list_from_csv(cols[19], int)  # qStarts
        psl["tstarts"] = _list_from_csv(cols[20], int)  # tStarts
        if self.pslx:
            psl["qseqs"] = _list_from_csv(cols[21])  # query sequence
            psl["tseqs"] = _list_from_csv(cols[22])  # hit sequence

        return psl

    def _validate_cols(self, cols):
        """Validate column's length of PSL or PSLX (PRIVATE)."""
        if not self.pslx:
            if len(cols) != 21:
                raise ValueError(
                    "Invalid PSL line: %r. Expected 21 tab-separated columns, found %i"
                    % (self.line, len(cols))
                )
        else:
            if len(cols) != 23:
                raise ValueError(
                    "Invalid PSLX line: %r. Expected 23 tab-separated columns, found %i"
                    % (self.line, len(cols))
                )

    def _parse_qresult(self):
        """Yield QueryResult objects (PRIVATE)."""
        # state values, determines what to do for each line
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        state_HIT_NEW = 2
        state_HIT_SAME = 4
        # initial dummy values
        qres_state = None
        file_state = None
        cur_qid, cur_hid = None, None
        prev_qid, prev_hid = None, None
        cur, prev = None, None
        hit_list, hsp_list = [], []

        while True:
            # store previous line's parsed values for all lines after the first
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            # only parse the result row if it's not EOF
            if self.line:
                cur = self._parse_row()
                cur_qid = cur["qname"]
                cur_hid = cur["tname"]
            else:
                file_state = state_EOF
                # mock values, since we have nothing to parse
                cur_qid, cur_hid = None, None

            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME
            # new hits are hits with different ids or hits in a new qresult
            if prev_hid != cur_hid or qres_state == state_QRES_NEW:
                hit_state = state_HIT_NEW
            else:
                hit_state = state_HIT_SAME

            if prev is not None:
                # create fragment and HSP and set their attributes
                hsp = _create_hsp(prev_hid, prev_qid, prev)
                hsp_list.append(hsp)

                if hit_state == state_HIT_NEW:
                    # create Hit and set its attributes
                    hit = Hit(hsp_list)
                    hit.seq_len = prev["tsize"]
                    hit_list.append(hit)
                    hsp_list = []

                # create qresult and yield if we're at a new qresult or at EOF
                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    qresult = QueryResult(id=prev_qid)
                    for hit in hit_list:
                        qresult.absorb(hit)
                    qresult.seq_len = prev["qsize"]
                    yield qresult
                    # if we're at EOF, break
                    if file_state == state_EOF:
                        break
                    hit_list = []

            self.line = self.handle.readline()


class BlatPslIndexer(SearchIndexer):
    """Indexer class for BLAT PSL output."""

    _parser = BlatPslParser

    def __init__(self, filename, pslx=False):
        """Initialize the class."""
        SearchIndexer.__init__(self, filename, pslx=pslx)

    def __iter__(self):
        """Iterate over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        # denotes column location for query identifier
        query_id_idx = 9
        qresult_key = None
        tab_char = b"\t"

        start_offset = handle.tell()
        line = handle.readline()
        # read through header
        # this assumes that the result row match the regex
        while not re.search(_RE_ROW_CHECK_IDX, line.strip()):
            start_offset = handle.tell()
            line = handle.readline()
            if not line:
                return

        # and index the qresults
        while True:
            end_offset = handle.tell()

            cols = [x for x in line.strip().split(tab_char) if x]
            if qresult_key is None:
                qresult_key = cols[query_id_idx]
            else:
                curr_key = cols[query_id_idx]

                if curr_key != qresult_key:
                    yield qresult_key.decode(), start_offset, end_offset - start_offset
                    qresult_key = curr_key
                    start_offset = end_offset - len(line)

            line = handle.readline()
            if not line:
                yield qresult_key.decode(), start_offset, end_offset - start_offset
                break

    def get_raw(self, offset):
        """Return raw bytes string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        query_id_idx = 9
        qresult_key = None
        qresult_raw = b""
        tab_char = b"\t"

        while True:
            line = handle.readline()
            if not line:
                break
            cols = [x for x in line.strip().split(tab_char) if x]
            if qresult_key is None:
                qresult_key = cols[query_id_idx]
            else:
                curr_key = cols[query_id_idx]
                if curr_key != qresult_key:
                    break
            qresult_raw += line

        return qresult_raw


class BlatPslWriter:
    """Writer for the blat-psl format."""

    def __init__(self, handle, header=False, pslx=False):
        """Initialize the class."""
        self.handle = handle
        # flag for writing header or not
        self.header = header
        self.pslx = pslx

    def write_file(self, qresults):
        """Write query results to file."""
        handle = self.handle
        qresult_counter, hit_counter, hsp_counter, frag_counter = 0, 0, 0, 0

        if self.header:
            handle.write(self._build_header())

        for qresult in qresults:
            if qresult:
                handle.write(self._build_row(qresult))
                qresult_counter += 1
                hit_counter += len(qresult)
                hsp_counter += sum(len(hit) for hit in qresult)
                frag_counter += sum(len(hit.fragments) for hit in qresult)

        return qresult_counter, hit_counter, hsp_counter, frag_counter

    def _build_header(self):
        """Build header, tab-separated string (PRIVATE)."""
        # for now, always use the psLayout version 3
        header = "psLayout version 3\n"

        # adapted from BLAT's source: lib/psl.c#L496
        header += (
            "\nmatch\tmis- \trep. \tN's\tQ gap\tQ gap\tT gap\tT "
            "gap\tstrand\tQ        \tQ   \tQ    \tQ  \tT        \tT   "
            "\tT    \tT  \tblock\tblockSizes \tqStarts\t tStarts"
            "\n     \tmatch\tmatch\t   \tcount\tbases\tcount\tbases"
            "\t      \tname     \tsize\tstart\tend\tname     \tsize"
            "\tstart\tend\tcount\n%s\n" % ("-" * 159)
        )

        return header

    def _build_row(self, qresult):
        """Return a string or one row or more of the QueryResult object (PRIVATE)."""
        # For now, our writer writes the row according to the order in
        # the QueryResult and Hit objects.
        # This is different from BLAT's native output, where the rows are
        # grouped by strand.
        # Should we tweak the behavior to better mimic the native output?
        qresult_lines = []

        for hit in qresult:
            for hsp in hit.hsps:
                query_is_protein = getattr(hsp, "query_is_protein", False)
                blocksize_multiplier = 3 if query_is_protein else 1

                line = []
                line.append(hsp.match_num)
                line.append(hsp.mismatch_num)
                line.append(hsp.match_rep_num)
                line.append(hsp.n_num)
                line.append(hsp.query_gapopen_num)
                line.append(hsp.query_gap_num)
                line.append(hsp.hit_gapopen_num)
                line.append(hsp.hit_gap_num)

                # check spans
                eff_query_spans = [blocksize_multiplier * s for s in hsp.query_span_all]
                if hsp.hit_span_all != eff_query_spans:
                    raise ValueError("HSP hit span and query span values do not match.")
                block_sizes = hsp.query_span_all

                # set strand and starts
                if hsp[0].query_strand >= 0:  # since it may be a protein seq
                    strand = "+"
                else:
                    strand = "-"
                qstarts = _reorient_starts(
                    [x[0] for x in hsp.query_range_all],
                    hsp.query_span_all,
                    qresult.seq_len,
                    hsp[0].query_strand,
                )

                if hsp[0].hit_strand == 1:
                    hstrand = 1
                    # only write hit strand if it was present in the source file
                    if hsp._has_hit_strand:
                        strand += "+"
                else:
                    hstrand = -1
                    strand += "-"
                hstarts = _reorient_starts(
                    [x[0] for x in hsp.hit_range_all],
                    hsp.hit_span_all,
                    hit.seq_len,
                    hstrand,
                )

                line.append(strand)
                line.append(qresult.id)
                line.append(qresult.seq_len)
                line.append(hsp.query_start)
                line.append(hsp.query_end)
                line.append(hit.id)
                line.append(hit.seq_len)
                line.append(hsp.hit_start)
                line.append(hsp.hit_end)
                line.append(len(hsp))
                line.append(",".join(str(x) for x in block_sizes) + ",")
                line.append(",".join(str(x) for x in qstarts) + ",")
                line.append(",".join(str(x) for x in hstarts) + ",")

                if self.pslx:
                    line.append(",".join(str(x.seq) for x in hsp.query_all) + ",")
                    line.append(",".join(str(x.seq) for x in hsp.hit_all) + ",")

                qresult_lines.append("\t".join(str(x) for x in line))

        return "\n".join(qresult_lines) + "\n"


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
