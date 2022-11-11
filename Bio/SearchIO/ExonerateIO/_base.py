# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO abstract base parser for Exonerate standard output format."""

import re
from functools import reduce
from abc import ABC, abstractmethod

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.SeqUtils import seq1


# strand char-value mapping
_STRAND_MAP = {"+": 1, "-": -1, ".": 0}

_RE_SHIFTS = re.compile(r"(#+)")
# regex for checking whether a vulgar line has protein/translated components
_RE_TRANS = re.compile(r"[53ISCF]")


def _set_frame(frag):
    """Set the HSPFragment frames (PRIVATE)."""
    frag.hit_frame = (frag.hit_start % 3 + 1) * frag.hit_strand
    frag.query_frame = (frag.query_start % 3 + 1) * frag.query_strand


def _make_triplets(seq, phase=0):
    """Select a valid amino acid sequence given a 3-letter code input (PRIVATE).

    This function takes a single three-letter amino acid sequence and the phase
    of the sequence to return the longest intact amino acid sequence possible.
    Parts of the input sequence before and after the selected sequence are also
    returned.

    This is an internal private function and is meant for parsing Exonerate's
    three-letter amino acid output.

    >>> from Bio.SearchIO.ExonerateIO._base import _make_triplets
    >>> _make_triplets('GlyThrSerAlaPro')
    ('', ['Gly', 'Thr', 'Ser', 'Ala', 'Pro'], '')
    >>> _make_triplets('yThrSerAla', phase=1)
    ('y', ['Thr', 'Ser', 'Ala'], '')
    >>> _make_triplets('yThrSerAlaPr', phase=1)
    ('y', ['Thr', 'Ser', 'Ala'], 'Pr')

    """
    pre = seq[:phase]
    np_seq = seq[phase:]
    non_triplets = len(np_seq) % 3
    post = "" if not non_triplets else np_seq[-1 * non_triplets :]
    intacts = [np_seq[3 * i : 3 * (i + 1)] for i in range(len(np_seq) // 3)]
    return pre, intacts, post


def _get_fragments_coord(frags):
    """Return the letter coordinate of the given list of fragments (PRIVATE).

    This function takes a list of three-letter amino acid sequences and
    returns a list of coordinates for each fragment had all the input
    sequences been flattened.

    This is an internal private function and is meant for parsing Exonerate's
    three-letter amino acid output.

    >>> from Bio.SearchIO.ExonerateIO._base import _get_fragments_coord
    >>> _get_fragments_coord(['Thr', 'Ser', 'Ala'])
    [0, 3, 6]
    >>> _get_fragments_coord(['Thr', 'SerAlaPro', 'GlyLeu'])
    [0, 3, 12]
    >>> _get_fragments_coord(['Thr', 'SerAlaPro', 'GlyLeu', 'Cys'])
    [0, 3, 12, 18]

    """
    if not frags:
        return []
    # first fragment always starts from position 0
    init = [0]
    return reduce(lambda acc, frag: acc + [acc[-1] + len(frag)], frags[:-1], init)


def _get_fragments_phase(frags):
    """Return the phases of the given list of 3-letter amino acid fragments (PRIVATE).

    This is an internal private function and is meant for parsing Exonerate's
    three-letter amino acid output.

    >>> from Bio.SearchIO.ExonerateIO._base import _get_fragments_phase
    >>> _get_fragments_phase(['Thr', 'Ser', 'Ala'])
    [0, 0, 0]
    >>> _get_fragments_phase(['ThrSe', 'rAla'])
    [0, 1]
    >>> _get_fragments_phase(['ThrSe', 'rAlaLeu', 'ProCys'])
    [0, 1, 0]
    >>> _get_fragments_phase(['ThrSe', 'rAlaLeuP', 'roCys'])
    [0, 1, 2]
    >>> _get_fragments_phase(['ThrSe', 'rAlaLeuPr', 'oCys'])
    [0, 1, 1]

    """
    return [(3 - (x % 3)) % 3 for x in _get_fragments_coord(frags)]


def _adjust_aa_seq(fraglist):
    """Transform 3-letter AA codes of input fragments to one-letter codes (PRIVATE).

    Argument fraglist should be a list of HSPFragments objects.
    """
    custom_map = {"***": "*", "<->": "-"}
    hsp_hstart = fraglist[0].hit_start
    hsp_qstart = fraglist[0].query_start
    frag_phases = _get_fragments_phase(fraglist)
    for frag, phase in zip(fraglist, frag_phases):
        assert frag.query_strand == 0 or frag.hit_strand == 0
        # hit step may be -1 as we're aligning to DNA
        hstep = 1 if frag.hit_strand >= 0 else -1

        # set fragment phase
        frag.phase = phase

        # fragment should have a length that is a multiple of 3
        # assert len(frag) % 3 == 0
        qseq = str(frag.query.seq)
        q_triplets_pre, q_triplets, q_triplets_post = _make_triplets(qseq, phase)

        hseq = str(frag.hit.seq)
        h_triplets_pre, h_triplets, h_triplets_post = _make_triplets(hseq, phase)

        # get one letter codes
        # and replace gap codon markers and termination characters
        hseq1_pre = "X" if h_triplets_pre else ""
        hseq1_post = "X" if h_triplets_post else ""
        hseq1 = seq1("".join(h_triplets), custom_map=custom_map)
        hstart = hsp_hstart + (len(hseq1_pre) * hstep)
        hend = hstart + len(hseq1.replace("-", "")) * hstep

        qseq1_pre = "X" if q_triplets_pre else ""
        qseq1_post = "X" if q_triplets_post else ""
        qseq1 = seq1("".join(q_triplets), custom_map=custom_map)
        qstart = hsp_qstart + len(qseq1_pre)
        qend = qstart + len(qseq1.replace("-", ""))

        # replace the old frag sequences with the new ones
        frag.hit = None
        frag.query = None
        frag.hit = hseq1_pre + hseq1 + hseq1_post
        frag.query = qseq1_pre + qseq1 + qseq1_post

        # set coordinates for the protein sequence
        if frag.query_strand == 0:
            frag.query_start, frag.query_end = qstart, qend
        elif frag.hit_strand == 0:
            frag.hit_start, frag.hit_end = hstart, hend

        # update alignment annotation
        # by turning them into list of triplets
        for annot, annotseq in frag.aln_annotation.items():
            pre, intact, post = _make_triplets(annotseq, phase)
            frag.aln_annotation[annot] = (
                list(filter(None, [pre])) + intact + list(filter(None, [post]))
            )

        # update values for next iteration
        hsp_hstart, hsp_qstart = hend, qend

    return fraglist


def _split_fragment(frag):
    """Split one HSPFragment containing frame-shifted alignment into two (PRIVATE)."""
    # given an HSPFragment object with frameshift(s), this method splits it
    # into fragments without frameshifts by sequentially chopping it off
    # starting from the beginning
    simil = frag.aln_annotation["similarity"]
    # we should have at least 1 frame shift for splitting
    assert simil.count("#") > 0

    split_frags = []
    qstep = 1 if frag.query_strand >= 0 else -1
    hstep = 1 if frag.hit_strand >= 0 else -1
    qpos = min(frag.query_range) if qstep >= 0 else max(frag.query_range)
    hpos = min(frag.hit_range) if hstep >= 0 else max(frag.hit_range)
    abs_pos = 0
    # split according to hit, then query
    while simil:

        try:
            shifts = re.search(_RE_SHIFTS, simil).group(1)
            s_start = simil.find(shifts)
            s_stop = s_start + len(shifts)
            split = frag[abs_pos : abs_pos + s_start]
        except AttributeError:  # no '#' in simil, i.e. last frag
            shifts = ""
            s_start = 0
            s_stop = len(simil)
            split = frag[abs_pos:]

        # coordinates for the split strand
        qstart, hstart = qpos, hpos
        qpos += (
            len(split) - sum(split.query.seq.count(x) for x in ("-", "<", ">"))
        ) * qstep
        hpos += (
            len(split) - sum(split.hit.seq.count(x) for x in ("-", "<", ">"))
        ) * hstep

        split.hit_start = min(hstart, hpos)
        split.query_start = min(qstart, qpos)
        split.hit_end = max(hstart, hpos)
        split.query_end = max(qstart, qpos)

        # account for frameshift length
        abs_slice = slice(abs_pos + s_start, abs_pos + s_stop)
        if len(frag.aln_annotation) == 2:
            seqs = (frag[abs_slice].query.seq, frag[abs_slice].hit.seq)
        elif len(frag.aln_annotation) == 3:
            seqs = (
                frag[abs_slice].aln_annotation["query_annotation"],
                frag[abs_slice].aln_annotation["hit_annotation"],
            )
        if "#" in seqs[0]:
            qpos += len(shifts) * qstep
        elif "#" in seqs[1]:
            hpos += len(shifts) * hstep

        # set frame
        _set_frame(split)
        split_frags.append(split)
        # set similarity string and absolute position for the next loop
        simil = simil[s_stop:]
        abs_pos += s_stop

    return split_frags


def _create_hsp(hid, qid, hspd):
    """Return a list of HSP objects from the given parsed HSP values (PRIVATE)."""
    frags = []
    # we are iterating over query_ranges, but hit_ranges works just as well
    for idx, qcoords in enumerate(hspd["query_ranges"]):
        # get sequences, create object
        hseqlist = hspd.get("hit")
        hseq = "" if hseqlist is None else hseqlist[idx]
        qseqlist = hspd.get("query")
        qseq = "" if qseqlist is None else qseqlist[idx]
        frag = HSPFragment(hid, qid, hit=hseq, query=qseq)
        # coordinates
        frag.query_start = qcoords[0]
        frag.query_end = qcoords[1]
        frag.hit_start = hspd["hit_ranges"][idx][0]
        frag.hit_end = hspd["hit_ranges"][idx][1]
        # alignment annotation
        try:
            aln_annot = hspd.get("aln_annotation", {})
            for key, value in aln_annot.items():
                frag.aln_annotation[key] = value[idx]
        except IndexError:
            pass
        # strands
        frag.query_strand = hspd["query_strand"]
        frag.hit_strand = hspd["hit_strand"]
        # and append the hsp object to the list
        if frag.aln_annotation.get("similarity") is not None:
            if "#" in frag.aln_annotation["similarity"]:
                frags.extend(_split_fragment(frag))
                continue
        # try to set frame if there are translation in the alignment
        if (
            len(frag.aln_annotation) > 1
            or frag.query_strand == 0
            or ("vulgar_comp" in hspd and re.search(_RE_TRANS, hspd["vulgar_comp"]))
        ):
            _set_frame(frag)

        frags.append(frag)

    # if the query is protein, we need to change the hit and query sequences
    # from three-letter amino acid codes to one letter, and adjust their
    # coordinates accordingly
    if len(frags[0].aln_annotation) == 2:  # 2 annotations == protein query
        frags = _adjust_aa_seq(frags)

    hsp = HSP(frags)
    # set hsp-specific attributes
    for attr in (
        "score",
        "hit_split_codons",
        "query_split_codons",
        "model",
        "vulgar_comp",
        "cigar_comp",
        "molecule_type",
    ):
        if attr in hspd:
            setattr(hsp, attr, hspd[attr])

    return hsp


def _parse_hit_or_query_line(line):
    """Parse the 'Query:' line of exonerate alignment outputs (PRIVATE)."""
    try:
        mark, id, desc = line.split(" ", 2)
    except ValueError:  # no desc
        mark, id = line.split(" ", 1)
        desc = ""

    return id, desc


def _get_strand_from_desc(desc, is_protein, modify_desc=True):
    """Determine the strand from the description (PRIVATE).

    Exonerate appends ``:[revcomp]`` (versions <= 2.2) or ``[revcomp]``
    (versions > 2.2) to the query and/or hit description string. This function
    outputs '-' if the description has such modifications or '+' if not. If the
    query and/or hit is a protein sequence, a '.' is output instead.

    Aside from the strand, the input description value is also returned. It is
    returned unmodified if ``modify_desc`` is ``False``. Otherwise, the appended
    ``:[revcomp]`` or ``[revcomp]`` is removed.

    """
    if is_protein:
        return ".", desc

    suffix = ""
    if desc.endswith("[revcomp]"):
        suffix = ":[revcomp]" if desc.endswith(":[revcomp]") else "[revcomp]"

    if not suffix:
        return "+", desc

    if modify_desc:
        return "-", desc[: -len(suffix)]

    return "-", desc


class _BaseExonerateParser(ABC):
    """Abstract base class iterator for exonerate format."""

    _ALN_MARK = None

    def __init__(self, handle):
        self.handle = handle
        self.has_c4_alignment = False

    def __iter__(self):
        # read line until the first alignment block or cigar/vulgar lines
        while True:
            self.line = self.handle.readline()
            # flag for human-readable alignment block
            if self.line.startswith("C4 Alignment:") and not self.has_c4_alignment:
                self.has_c4_alignment = True
            if (
                self.line.startswith("C4 Alignment:")
                or self.line.startswith("vulgar:")
                or self.line.startswith("cigar:")
            ):
                break
            elif not self.line or self.line.startswith("-- completed "):
                return

        for qresult in self._parse_qresult():
            qresult.program = "exonerate"
            # HACK: so that all descriptions are set
            qresult.description = qresult.description
            for hit in qresult:
                hit.description = hit.description
            yield qresult

    def read_until(self, bool_func):
        """Read the file handle until the given bool function returns True."""
        while True:
            if not self.line or bool_func(self.line):
                return
            else:
                self.line = self.handle.readline()

    @abstractmethod
    def parse_alignment_block(self, header):
        raise NotImplementedError

    def _parse_alignment_header(self):
        # read all header lines and store them
        aln_header = []
        # header is everything before the first empty line
        while self.line.strip():
            aln_header.append(self.line.strip())
            self.line = self.handle.readline()
        # then parse them
        qresult, hit, hsp = {}, {}, {}
        for line in aln_header:
            # query line
            if line.startswith("Query:"):
                qresult["id"], qresult["description"] = _parse_hit_or_query_line(line)
            # target line
            elif line.startswith("Target:"):
                hit["id"], hit["description"] = _parse_hit_or_query_line(line)
            # model line
            elif line.startswith("Model:"):
                qresult["model"] = line.split(" ", 1)[1]
            # score line
            elif line.startswith("Raw score:"):
                hsp["score"] = line.split(" ", 2)[2]
            # query range line
            elif line.startswith("Query range:"):
                # line is always 'Query range: \d+ -> \d+', so we can pluck
                # the numbers directly
                hsp["query_start"], hsp["query_end"] = line.split(" ", 4)[2:5:2]
            # hit range line
            elif line.startswith("Target range:"):
                # same logic with query range
                hsp["hit_start"], hsp["hit_end"] = line.split(" ", 4)[2:5:2]

        # determine strand
        qresult_strand, qresult_desc = _get_strand_from_desc(
            desc=qresult["description"],
            is_protein="protein2" in qresult["model"],
            modify_desc=True,
        )
        hsp["query_strand"] = qresult_strand
        qresult["description"] = qresult_desc

        hit_strand, hit_desc = _get_strand_from_desc(
            desc=hit["description"],
            is_protein="2protein" in qresult["model"],
            modify_desc=True,
        )
        hsp["hit_strand"] = hit_strand
        hit["description"] = hit_desc

        # NOTE: we haven't processed the coordinates types
        # and the strands are not yet Biopython's standard (1 / -1 / 0)
        # since it's easier if we do the conversion later

        return {"qresult": qresult, "hit": hit, "hsp": hsp}

    def _parse_qresult(self):
        # state values
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        state_HIT_NEW = 2
        state_HIT_SAME = 4
        # initial dummies
        qres_state, hit_state = None, None
        file_state = None
        cur_qid, cur_hid = None, None
        prev_qid, prev_hid = None, None
        cur, prev = None, None
        hit_list, hsp_list = [], []
        # if the file has c4 alignments, use that as the alignment mark
        if self.has_c4_alignment:
            self._ALN_MARK = "C4 Alignment:"

        while True:
            self.read_until(lambda line: line.startswith(self._ALN_MARK))
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            # only parse the result row if it's not EOF
            if self.line:
                assert self.line.startswith(self._ALN_MARK), self.line
                # create temp dicts for storing parsed values
                header = {"qresult": {}, "hit": {}, "hsp": {}}
                # if the file has c4 alignments, try to parse the header
                if self.has_c4_alignment:
                    self.read_until(lambda line: line.strip().startswith("Query:"))
                    header = self._parse_alignment_header()
                # parse the block contents
                cur = self.parse_alignment_block(header)
                cur_qid = cur["qresult"]["id"]
                cur_hid = cur["hit"]["id"]
            elif not self.line or self.line.startswith("-- completed "):
                file_state = state_EOF
                cur_qid, cur_hid = None, None

            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME
            # new hits are hits with different ids or hits in a new query
            if prev_hid != cur_hid or qres_state == state_QRES_NEW:
                hit_state = state_HIT_NEW
            else:
                hit_state = state_HIT_SAME

            if prev is not None:
                hsp = _create_hsp(prev_hid, prev_qid, prev["hsp"])
                hsp_list.append(hsp)

                if hit_state == state_HIT_NEW:
                    hit = Hit(hsp_list)
                    for attr, value in prev["hit"].items():
                        setattr(hit, attr, value)
                    hit_list.append(hit)
                    hsp_list = []

                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    qresult = QueryResult(id=prev_qid)
                    for hit in hit_list:
                        # not using append since Exonerate may separate the
                        # same hit if it has different strands
                        qresult.absorb(hit)
                    for attr, value in prev["qresult"].items():
                        setattr(qresult, attr, value)
                    yield qresult
                    if file_state == state_EOF:
                        break
                    hit_list = []

            # only readline() here if we're not parsing C4 alignments
            # C4 alignments readline() is handled by its parse_alignment_block
            # function
            if not self.has_c4_alignment:
                self.line = self.handle.readline()


class _BaseExonerateIndexer(SearchIndexer):
    """Indexer class for Exonerate plain text."""

    _parser = None  # should be defined by subclass
    _query_mark = None  # this one too

    def get_qresult_id(self, pos):
        raise NotImplementedError("Should be defined by subclass")

    def __iter__(self):
        """Iterate over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        qresult_key = None

        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if line.startswith(self._query_mark):
                if qresult_key is None:
                    qresult_key = self.get_qresult_id(start_offset)
                    qresult_offset = start_offset
                else:
                    curr_key = self.get_qresult_id(start_offset)
                    if curr_key != qresult_key:
                        yield qresult_key, qresult_offset, start_offset - qresult_offset
                        qresult_key = curr_key
                        qresult_offset = start_offset
                        handle.seek(qresult_offset)
            elif not line:
                yield qresult_key, qresult_offset, start_offset - qresult_offset
                break


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
