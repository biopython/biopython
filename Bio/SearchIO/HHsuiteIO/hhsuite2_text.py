# Copyright 2019 by Jens Thomas.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO parser for HHSUITE version 2 and 3 plain text output format."""

import re
import warnings

from Bio.SearchIO._utils import read_forward
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment

__all__ = ("Hhsuite2TextParser",)

# precompile regex patterns for faster processing
# regex for query name capture
_RE_QUERY = re.compile(r"^Query\s+(.+)\s?$")

# regex for version string capture
_RE_HIT_BLOCK_START = re.compile(r"^No +(\d+)\s+$")

# id and full description
_RE_HIT_BLOCK_DESC = re.compile(r">(\S+)\s+(.*)$")

# sequence alignment line
# Q sp|Q9BSU1|CP07  229 DAKMRVFERSVYFGDSCQDVLSMLGSPHKV  258 (422)
_RE_MATCH_BLOCK_QUERY_SEQ = re.compile(r"^Q\s+(.+) +(\d+) +([A-Z-]+) +(\d+) +\(\d+\)$")
_RE_MATCH_BLOCK_HIT_SEQ = re.compile(r"^T\s+(.+) +(\d+) +([A-Z-]+) +(\d+) +\(\d+\)$")

_END_OF_FILE_MARKER = "Done!"

_PROGRAM = "HHSUITE"

# Maximum number of lines to read before expecting a hit block
# This determines the maximum number of hits that would be allowed in
# the initial hit table.
MAX_READ_UNTIL = 5000


class Hhsuite2TextParser:
    """Parser for the HHSUITE version 2 and 3 text output."""

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        self.line = read_forward(self.handle)
        self.done = False
        self.query_id = None
        self.seq_len = None

    def __iter__(self):
        """Iterate over query results - there will only ever be one."""
        yield from self._parse_qresult()

    def _read_until(self, bool_func, stop_on_blank=True, max_read_until=MAX_READ_UNTIL):
        """Read the file handle until the given function returns True (PRIVATE)."""
        count = 0
        while True:
            if stop_on_blank and not self.line:
                return
            if bool_func(self.line):
                return
            else:
                self.line = read_forward(self.handle)
            count += 1
            if count >= max_read_until:
                raise RuntimeError("Exceeded max_read_until in _read_until")

    def _parse_qresult(self):
        """Parse HHSUITE output file (PRIVATE)."""
        hit_block_data = []
        self._parse_preamble()
        self._read_until(
            lambda line: re.search(_RE_HIT_BLOCK_START, line), stop_on_blank=False
        )
        while not self.done:
            hit_dict = self._parse_hit_block()
            hit_block_data.append(hit_dict)
        return self._create_qresult(hit_block_data)

    def _parse_preamble(self):
        """Parse metadata about query (PRIVATE)."""
        meta = {}
        while self.line:
            regx = re.search(_RE_QUERY, self.line)
            if regx:
                self.query_id = regx.group(1)
            if self.line.startswith("Match_columns"):
                self.seq_len = int(self.line.strip().split()[1])
            self.line = self.handle.readline().strip()
        return meta

    def _parse_hit_block(self):
        """Parse a hit block (PRIVATE)."""
        self.line = read_forward(self.handle)
        match = re.search(_RE_HIT_BLOCK_DESC, self.line)
        if not match:
            raise RuntimeError(
                f"Unexpected content in HIT_BLOCK_DESC line'{self.line}'"
            )
        hit_data = {
            "hit_id": match.group(1),
            "description": match.group(2).lstrip(" ;"),
            "evalue": None,
            "hit_start": None,
            "hit_end": None,
            "hit_seq": "",
            "prob": None,
            "query_start": None,
            "query_end": None,
            "query_seq": "",
            "score": None,
        }
        self.line = self.handle.readline()
        self._process_score_line(self.line, hit_data)
        while True:
            self.line = read_forward(self.handle)
            if not self.line.strip() or self.line.startswith(_END_OF_FILE_MARKER):
                # _END_OF_FILE_MARKER isn't always present
                self.done = True
                return hit_data
            elif re.search(_RE_HIT_BLOCK_START, self.line):
                return hit_data
            else:
                self._parse_hit_match_block(hit_data)

    @staticmethod
    def _process_score_line(line, hit_data):
        """Parse the scores from the line and populate hit_data dict (PRIVATE).

        Lines are of the form:
        Probab=99.95  E-value=3.7e-34  Score=210.31  Aligned_cols=171  Identities=100%  Similarity=2.050  Sum_probs=166.9

        E-value could be in decimal or scientific notation, so split the string rather then use regexp - this
        also means we should be tolerant of additional fields being added/removed
        """
        score_map = {"E-value": "evalue", "Score": "score", "Probab": "prob"}
        for score_pair in line.strip().split():
            key, value = score_pair.split("=")
            if key in score_map:
                try:
                    hit_data[score_map[key]] = float(value)
                except KeyError:
                    # We trigger warnings here as it's not a big enough problem to crash, but indicates something unexpected.
                    warnings.warn(
                        f"HHsuite parser: unable to extract {key} from line: {line}"
                    )

    def _parse_hit_match_block(self, hit_match_data):
        """Parse a single block of hit sequence data (PRIVATE).

        Parses block such as ::

            Q ss_pred             ceecchHHHHHHHHHHHHHHHHHHHhhhhhcCCCCccc
            Q 4P79:A|PDBID|C  160 YELGPALYLGWSASLLSILGGICVFSTAAASSKEEPAT  197 (198)
            Q Consensus       160 ~~~g~sf~l~~~~~~l~~~~~~l~~~~~~~~~~~~~~~  197 (198)
                                  .++|||||++|++.++.+++++++++..+..++++..+
            T Consensus       327 ~~~GwS~~l~~~s~~l~lia~~l~~~~~~~~~~~~~~~  364 (364)
            T 5B2G_A          327 REMGASLYVGWAASGLLLLGGGLLCCSGPSSGENLYFQ  364 (364)
            T ss_dssp             EEECTHHHHHHHHHHHHHHHHHHHHCC-----------
            T ss_pred             cccchHHHHHHHHHHHHHHHHHHHHhcCCCCCCccccC

        """

        def match_is_valid(match):
            """Return True if match is not a Consensus column (PRIVATE).

            It's not possible to distinguish a sequence line from a Consensus line with
            a regexp, so need to check the ID column.
            """
            return match.group(1).strip() != "Consensus"

        while True:
            if not self.line.strip():  # blank lines indicate the end of a hit block
                return
            match = re.match(_RE_MATCH_BLOCK_QUERY_SEQ, self.line)
            if match and match_is_valid(match):
                hit_match_data["query_seq"] += match.group(3).strip()
                if hit_match_data["query_start"] is None:
                    hit_match_data["query_start"] = int(match.group(2))
                hit_match_data["query_end"] = int(match.group(4))
            else:
                match = re.match(_RE_MATCH_BLOCK_HIT_SEQ, self.line)
                if match and match_is_valid(match):
                    hit_match_data["hit_seq"] += match.group(3).strip()
                    if hit_match_data["hit_start"] is None:
                        hit_match_data["hit_start"] = int(match.group(2))
                    hit_match_data["hit_end"] = int(match.group(4))
            self.line = self.handle.readline()

    def _create_qresult(self, hit_blocks):
        """Create the Biopython data structures from the parsed data (PRIVATE)."""
        query_id = self.query_id
        hit_dict = {}

        for output_index, block in enumerate(hit_blocks):
            hit_id = block["hit_id"]

            frag = HSPFragment(hit_id, query_id)
            frag.molecule_type = "protein"
            frag.query_start = block["query_start"] - 1
            frag.query_end = block["query_end"]
            frag.hit_start = block["hit_start"] - 1
            frag.hit_end = block["hit_end"]
            frag.hit = block["hit_seq"]
            frag.query = block["query_seq"]

            hsp = HSP([frag])
            hsp.hit_id = hit_id
            hsp.output_index = output_index
            hsp.query_id = query_id
            hsp.hit_description = block["description"]
            is_included = True  # Should everything should be included?
            hsp.is_included = is_included
            hsp.evalue = block["evalue"]
            hsp.score = block["score"]
            hsp.prob = block["prob"]

            if hit_id not in hit_dict:
                hit = Hit([hsp], hit_id)
                hit.description = block["description"]
                hit.is_included = is_included
                hit.evalue = block["evalue"]
                hit.score = block["score"]
                hit_dict[hit_id] = hit
            else:
                hit_dict[hit_id].append(hsp)

        qresult = QueryResult(hit_dict.values(), query_id)
        qresult.program = _PROGRAM
        qresult.seq_len = self.seq_len
        return [qresult]
