# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.


"""Bio.SearchIO parser for Infernal plain text output format."""

import operator
import re

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult
from Bio.SearchIO._utils import read_forward

from ._base import _BaseInfernalParser

__all__ = ("InfernalTextParser", "InfernalTextIndexer")

# precompile regex patterns for faster processing
# program name
_RE_PROGRAM = re.compile(r"^# .*?(\w?cm\w+) :: .*$")
# version string
_RE_VERSION = re.compile(r"# INFERNAL+ ([\w+\.]{2,}) .*$")
# option string
_RE_OPT = re.compile(r"^# (.+):\s+(.+)$")
# numbers only
_RE_NUMERIC = re.compile(r"\d+")
_RE_NOT_NUMERIC = re.compile(r"\D")
# letters only
_RE_LETTERS = re.compile(r"[^A-Za-z]")
# missing segments in the alignment block
_RE_SPLIT_ALN = re.compile(r"(\*\[[ 0-9]+\]\*)+")

# divider for Infernal plain text output, these divider are
# always taken from the beginning of the line to be usable in
# self.line.startswith() method
# header options
_DIV_HEADER_OPT = "# - - -"
# query
_DIV_QUERY_END = "//"
_DIV_QUERY_START = "Query:"
# hits
_DIV_HITS_END_CM = "Internal CM pipeline statistics summary:"
_DIV_HITS_END_HMM = "Internal HMM-only pipeline statistics summary:"
_DIV_HIT_SCORE_TABLE = "Hit scores:"
_DIV_HIT_ALIGNMENT = "Hit alignments:"
_DIV_NO_HIT = "   [No hits detected that satisfy reporting thresholds]"
_DIV_TABLE_START = " ----   --------- ------"
_DIV_INC_THRESHOLD = " ------ inclusion threshold ------"
# hit alignment
_DIV_ALIGNMENT_START = ">> "


class InfernalTextParser(_BaseInfernalParser):
    """Parser for the Infernal text output."""

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        self.line = read_forward(self.handle)
        self._meta = self._parse_header()

    def __iter__(self):
        """Iterate over query results."""
        yield from self._parse_qresult()

    def _read_until(self, bool_func):
        """Read the file handle until the given function returns True (PRIVATE)."""
        while True:
            if not self.line or bool_func(self.line):
                return
            self.line = read_forward(self.handle)

    def _parse_header(self):
        """Parse Infernal header (PRIVATE)."""
        meta = {}
        # set the default value for the presence of alignment
        # as this information is important for the hit section parsing
        meta["show alignments in output"] = "yes"
        # bool flag for storing state ~ whether we are parsing the option
        # lines or not
        has_opts = False
        while True:
            # no pound sign means we've left the preamble
            if not self.line.startswith("#"):
                break
            # dashes could either mean we are entering or leaving the options
            # section, so it's a switch for the has_opts flag
            elif self.line.startswith(_DIV_HEADER_OPT):
                if not has_opts:
                    # if flag is false, that means we're entering opts
                    # so switch the flag accordingly
                    has_opts = True
                else:
                    # if flag is true, that means we've reached the end of opts
                    # so we can break out of the function
                    break
            elif not has_opts:
                # try parsing program
                regx = re.search(_RE_PROGRAM, self.line)
                if regx:
                    meta["program"] = regx.group(1)
                # try parsing version
                regx = re.search(_RE_VERSION, self.line)
                if regx:
                    meta["version"] = regx.group(1)
            elif has_opts:
                regx = re.search(_RE_OPT, self.line)
                # if target in regx.group(1), then we store the key as target
                if "target" in regx.group(1):
                    meta["target"] = regx.group(2).strip()
                else:
                    meta[regx.group(1)] = regx.group(2)

            self.line = read_forward(self.handle)

        return meta

    def _parse_qresult(self):
        """Parse a Infernal query block (PRIVATE)."""
        self._read_until(lambda line: line.startswith(_DIV_QUERY_START))

        while self.line:
            # create qresult for query id
            if self.line.startswith(_DIV_QUERY_START):
                qid = self.line.strip().split()[1]
                qlen = int(re.sub(_RE_NOT_NUMERIC, "", self.line.strip().split()[-1]))

                # store qresult attributes
                qresult_attrs = {
                    "id": qid,
                    "seq_len": qlen,
                    "program": self._meta.get("program"),
                    "version": self._meta.get("version"),
                    "target": self._meta.get("target"),
                }
            else:
                self.line = read_forward(self.handle)

            # get description and accession, if they exist
            qdesc = "<unknown description>"  # placeholder
            while not self.line.startswith(_DIV_HIT_SCORE_TABLE):
                self.line = read_forward(self.handle)

                if self.line.startswith("Accession:"):
                    acc = self.line.strip().split(" ", 1)[1]
                    qresult_attrs["accession"] = acc.strip()
                elif self.line.startswith("Description:"):
                    qdesc = self.line.strip().split(" ", 1)[1].strip()
                    qresult_attrs["description"] = qdesc

            # parse the query hits
            # initializing hit_list directly to handle empty hits
            hit_list = []
            while self.line and not self.line.startswith(_DIV_QUERY_END):
                hit_list = self._parse_hit(qid, qdesc)
                # read through the statistics summary
                if self.line.startswith(_DIV_HITS_END_CM) or self.line.startswith(
                    _DIV_HITS_END_HMM
                ):
                    while self.line and not self.line.startswith(_DIV_QUERY_END):
                        self.line = read_forward(self.handle)
            # create qresult, set its attributes and yield
            qresult = QueryResult(id=qid, hits=hit_list)
            for attr, value in qresult_attrs.items():
                setattr(qresult, attr, value)
            yield qresult
            self.line = read_forward(self.handle)

            # Infernal outputs '[ok]' at the end of all results file,
            # which means we can break the main loop when we see the line
            if "[ok]" in self.line:
                break

    def _parse_hit(self, qid, qdesc):
        """Parse an Infernal hit (PRIVATE)."""
        # state values, determines what to do for each block
        hit_end = False
        in_score_table = False
        parsing_hits = False
        prev_line = None
        # hits are not guaranteed to be in order, and duplicate hit
        # IDs will raise a ValueError. to avoid this, hits are stored in
        # a temporary dictionary and the hit list required to create
        # the QueryResult is generated at the end of the query block
        hit_dict = {}

        # set the divider based on the output type (with or without alignment)
        if self._meta["show alignments in output"] == "no":
            div_hit_start = _DIV_HIT_SCORE_TABLE
        else:
            div_hit_start = _DIV_ALIGNMENT_START

        while True:
            if not self.line:
                raise ValueError("Unexpected end of file")
            # if there are no hits, forward-read to the end of the query
            elif self.line.startswith(_DIV_NO_HIT):
                while True:
                    self.line = read_forward(self.handle)
                    if self.line.startswith(_DIV_HITS_END_CM) or self.line.startswith(
                        _DIV_HITS_END_HMM
                    ):
                        hit_end = True
                        return []
            # entering hit alignment block
            elif self.line.startswith(div_hit_start):
                # for --noali output, move to the beginning of the hit score table
                if self._meta["show alignments in output"] == "no":
                    assert not in_score_table
                    self._read_until(lambda line: line.startswith(_DIV_TABLE_START))
                    self.line = read_forward(self.handle)
                    parsing_hits = in_score_table = True
                else:
                    self._parse_hit_from_alignment(qid, hit_dict)
                    parsing_hits = True
            # we've reached the end of the hit section
            elif self.line.startswith(_DIV_HITS_END_CM) or self.line.startswith(
                _DIV_HITS_END_HMM
            ):
                hit_end = True
                parsing_hits = False
            # Read through the scores table. For regular output this information
            # is not needed, so we can skip it. For --noali output, this table
            # contains the HSP information
            else:
                # for --noali output, keep the line
                if in_score_table:
                    if self.line.strip():
                        prev_line = self.line
                    else:
                        parsing_hits = in_score_table = False

                # read one line at the time
                self.line = self.handle.readline()

            # for --noali output, parse the scores table row
            if in_score_table and prev_line is not None:
                if not prev_line.startswith(_DIV_INC_THRESHOLD):
                    self._parse_scores_table_row(prev_line, qid, hit_dict)

            # we've reached the end of the query block hits
            # creating the Hit objects for this query
            if hit_end:
                return self._hit_to_list(hit_dict)

    def _parse_hit_from_alignment(self, qid, hit_dict):
        """Parse an Infernal hit alignment (PRIVATE)."""
        hid, hdesc = self.line[len(_DIV_ALIGNMENT_START) :].split("  ", 1)
        hdesc = hdesc.strip()

        # read through the hit table header and move one more line
        self._read_until(lambda line: line.startswith(_DIV_TABLE_START))
        self.line = read_forward(self.handle)

        # parse the hit table
        row = [x for x in self.line.strip().split() if x]
        assert len(row) == 16

        # create hit and append to hit container
        hit_attrs = {"id": hid, "query_id": qid, "description": hdesc}
        hsp_attrs = {
            "evalue": float(row[2]),
            "bitscore": float(row[3]),
            "bias": float(row[4]),
            "model": row[5],
            "truncated": row[14],
            "gc": float(row[15]),
            "avg_acc": float(row[13]),
            "query_endtype": row[8],
            "hit_endtype": row[12],
            "is_included": True if row[1] == "!" else False,
        }
        query_start = int(row[6])
        query_end = int(row[7])
        hit_start = int(row[9]) if row[11] == "+" else int(row[10])
        hit_end = int(row[10]) if row[11] == "+" else int(row[9])
        hit_strand = 0 if row[11] == "+" else -1

        # move to the HSP alignment block
        self.line = read_forward(self.handle)

        # create the hsp
        frag_list = self._parse_aln_block(
            hit_attrs["id"],
            hit_attrs["query_id"],
            hsp_attrs["model"],
            query_start,
            query_end,
            hit_start,
            hit_end,
            hit_strand,
        )
        hsp = HSP(frag_list)
        for attr, value in hsp_attrs.items():
            setattr(hsp, attr, value)

        # add the hit to the container
        self._add_hit_to_dict(hit_attrs, hsp, hit_dict)

    def _parse_scores_table_row(self, row, qid, hit_dict):
        """Parse an Infernal hit scores table (when used with --noali) (PRIVATE)."""

        # parse the columns into a list
        row = [x for x in row.strip().split(" ") if x]
        # join the description words if it's >1 word
        if len(row) > 12:
            row[12] = " ".join(row[12:])
        # if there's no description, set it to an empty string
        elif len(row) < 12:
            row.append("")
            assert len(row) == 12

        # parse the attributes
        hit_attrs = {"id": row[5], "query_id": qid, "description": row[12]}
        hsp_attrs = {
            "evalue": float(row[2]),
            "bitscore": float(row[3]),
            "bias": float(row[4]),
            "model": row[9],
            "truncated": row[10],
            "gc": float(row[11]),
            "is_included": True if row[1] == "!" else False,
        }
        hsp_frag_attrs = {
            "hit_start": int(row[6]) if row[8] == "+" else int(row[7]),
            "hit_end": int(row[7]) if row[8] == "+" else int(row[6]),
            "hit_strand": 0 if row[8] == "+" else -1,
        }

        # create the hsp fragment and set it's attributes
        hsp_frag = HSPFragment(row[5], qid)
        for attr, value in hsp_frag_attrs.items():
            setattr(hsp_frag, attr, value)

        # create the hsp and set it's attributes
        hsp = HSP([hsp_frag])
        for attr, value in hsp_attrs.items():
            setattr(hsp, attr, value)

        # add the hit to the container
        self._add_hit_to_dict(hit_attrs, hsp, hit_dict)

    def _parse_aln_block(
        self, hid, qid, model, query_start, query_end, hit_start, hit_end, hit_strand
    ):
        """Parse a Infernal HSP alignment block (PRIVATE)."""
        frag_list = []
        model_seq = ""
        hit_seq = ""
        if model == "cm":
            annot = {"NC": "", "CS": "", "similarity": "", "PP": ""}
        else:
            annot = {"CS": "", "similarity": "", "PP": ""}
        while True:
            # we've reached the end of the alignment section
            if (
                self.line.startswith(_DIV_ALIGNMENT_START)
                or self.line.startswith(_DIV_HITS_END_CM)
                or self.line.startswith(_DIV_HITS_END_HMM)
            ):
                # Process local end in infernal hit alignment. Local end are
                # large insertion or deletion indicated by *[NN]* where N is
                # the number of model positions are deleted or the number of
                # residues are inserted in the sequence. We split the
                # alignment block in HSPs on these local ends.

                # get local ends string (*[NN]*) indexes in the model sequence
                # there can be more than one local ends back-to-back
                local_aln_idx = [
                    (0, 0)
                ]  # there is always at least one local alignment starting at 0
                local_aln_idx += [
                    (m.start(0), m.end(0))
                    for m in re.finditer(_RE_SPLIT_ALN, model_seq)
                ]

                prev_hit_start = hit_start if hit_strand == 0 else hit_end
                prev_model_start = query_start
                hsps = []

                for i in range(len(local_aln_idx)):
                    local_start = local_aln_idx[i][1]
                    local_end = (
                        local_aln_idx[i + 1][0] if i + 1 < len(local_aln_idx) else None
                    )

                    # hsp hit position. forward strand moves up, reverse strand moves down
                    op = operator.add if hit_strand == 0 else operator.sub
                    cur_hit_seq = hit_seq[local_start:local_end]
                    cur_hit_gap_size = self._local_aln_gap_size(
                        local_aln_idx[i], hit_seq
                    )
                    cur_hit_start = op(prev_hit_start, cur_hit_gap_size)
                    cur_hit_end = op(
                        cur_hit_start, len(re.sub(_RE_LETTERS, "", cur_hit_seq))
                    )
                    # adjust start and end elements positon
                    if hit_strand == 0 and i == 0:
                        cur_hit_end -= 1
                    if hit_strand == -1 and i == 0:
                        cur_hit_end += 1
                    prev_hit_start = cur_hit_end
                    # hsp model position
                    cur_model_seq = model_seq[local_start:local_end].replace(".", "-")
                    cur_model_gap_size = self._local_aln_gap_size(
                        local_aln_idx[i], model_seq
                    )
                    cur_model_start = prev_model_start + cur_model_gap_size
                    cur_model_end = cur_model_start + len(
                        re.sub(_RE_LETTERS, "", cur_model_seq)
                    )
                    if i == 0:
                        cur_model_end -= 1
                    prev_model_start = cur_model_end
                    # annotations
                    cur_annot = {k: v[local_start:local_end] for k, v in annot.items()}

                    # create the hsp fragment and add to container
                    frag = HSPFragment(hid, qid)
                    frag.query = cur_model_seq
                    frag.hit = cur_hit_seq
                    frag.query_start = cur_model_start
                    frag.query_end = cur_model_end
                    frag.hit_start = cur_hit_start if hit_strand == 0 else cur_hit_end
                    frag.hit_end = cur_hit_end if hit_strand == 0 else cur_hit_start
                    frag.hit_strand = hit_strand
                    frag.aln_annotation = cur_annot

                    frag_list.append(frag)

                return frag_list

            # parse the alignment blocks in the hsp
            # each block have 4 (hmmonly) or 5 (cm) lines followed by an empty line
            block_size = 6 if model == "cm" else 5
            offset = 1 if model == "cm" else 0  # offset for the annotation line indexes
            lines = [None] * block_size
            for i in range(block_size):
                lines[i] = self.line
                self.line = read_forward(self.handle)

            # get the position of the alignment in the string using the PP line
            blklen = len(lines[4 + offset].strip().split()[0])
            blkstart = len(lines[4 + offset]) - blklen - 4
            blkend = len(lines[4 + offset]) - 4
            model_seq += lines[1 + offset][blkstart:blkend]
            hit_seq += lines[3 + offset][blkstart:blkend]
            # NC line is specific to cm model searches
            if model == "cm":
                annot["NC"] += lines[0][blkstart:blkend]
            annot["CS"] += lines[0 + offset][blkstart:blkend]
            annot["similarity"] += lines[2 + offset][blkstart:blkend]
            annot["PP"] += lines[4 + offset][blkstart:blkend]

    def _local_aln_gap_size(self, cur_aln_idx, seq):
        """Calculate the gap size between the local alignments (PRIVATE)."""
        gap_len = 0
        if cur_aln_idx[1] > 0:
            gap_len = sum(
                [
                    int(n)
                    for n in re.findall(
                        _RE_NUMERIC, seq[cur_aln_idx[0] : cur_aln_idx[1]]
                    )
                ]
            )
            assert gap_len > 0
        return gap_len


class InfernalTextIndexer(SearchIndexer):
    """Indexer class for Infernal plain text output."""

    _parser = InfernalTextParser

    def __init__(self, *args, **kwargs):
        """Initialize the class."""
        super().__init__(*args, **kwargs)
        self._preamble = b""

    def __iter__(self):
        """Iterate over InfernalTextIndexer; yields query results' key, offsets, 0."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()

        while True:
            line = read_forward(handle)
            end_offset = handle.tell()

            if line.startswith(_DIV_QUERY_START.encode()):
                qresult_key = line.strip().split()[1]
                # qresult start offset is the offset of this line
                # (starts with the start mark)
                start_offset = end_offset - len(line)
            elif line.startswith(_DIV_QUERY_END.encode()):
                yield qresult_key.decode(), start_offset, 0
                start_offset = end_offset
            elif not line:
                break

    def get_raw(self, offset):
        """Return the raw record from the file as a bytes string."""
        handle = self._handle
        qresult_raw = b""

        # read header
        if not self._preamble:
            handle.seek(0)
            while True:
                line = handle.readline()
                if line.startswith(_DIV_QUERY_START.encode()):
                    break
                self._preamble += line

        qresult_raw += self._preamble

        # read the qresult raw string
        handle.seek(offset)
        while True:
            # preserve whitespace, don't use read_forward
            line = handle.readline()
            qresult_raw += line

            # break when we've reached qresult end
            if line.startswith(_DIV_QUERY_END.encode()) or not line:
                break

        return qresult_raw


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
