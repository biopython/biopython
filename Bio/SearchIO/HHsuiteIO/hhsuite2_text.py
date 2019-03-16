# Copyright 2019 by Jens Thomas.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO parser for HHSUITE version 2 plain text output format."""

import re

from Bio._utils import read_forward
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment
from Bio.Alphabet import generic_protein

__all__ = ('Hhsuite2TextParser')

# precompile regex patterns for faster processing
# regex for query name capture
_RE_QUERY = re.compile(r'^Query\s+(.+)\s?$')

# regex for version string capture
_RE_HIT_BLOCK_START = re.compile(r'^No +(\d+)\s+$')

# id and full description
_RE_HIT_BLOCK_DESC = re.compile(r'>(\S+)\s+(.*)$')

# sequence alignment line
# Q sp|Q9BSU1|CP07  229 DAKMRVFERSVYFGDSCQDVLSMLGSPHKV  258 (422)
_RE_MATCH_BLOCK_QUERY_SEQ = re.compile(r'^Q\s+(.+) +(\d+) +([A-Z-]+) +(\d+) +\(\d+\)$')
_RE_MATCH_BLOCK_HIT_SEQ = re.compile(r'^T\s+(.+) +(\d+) +([A-Z-]+) +(\d+) +\(\d+\)$')

_END_OF_FILE_MARKER = 'Done!'

_PROGRAM = 'HHSUITE'

# Maximum number of lines to read before expecting a hit block
# This determines the maximum numnber of hits that would be allowed in
# the initial hit table.
MAX_READ_UNTIL = 5000


class Hhsuite2TextParser(object):
    """Parser for the HMMER 3.0 text output."""

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        self.line = read_forward(self.handle)
        self.done = False
        self.query_id = None
        self.seq_len = None

    def __iter__(self):
        """Iterate over query results - there will only ever be one."""
        for qresult in self._parse_qresult():
            yield qresult

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
        self._read_until(lambda line: re.search(_RE_HIT_BLOCK_START, line), stop_on_blank=False)
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
            if self.line.startswith('Match_columns'):
                self.seq_len = int(self.line.strip().split()[1])
            self.line = self.handle.readline().strip()
        return meta

    def _parse_hit_block(self):
        """Parse a hit block (PRIVATE)."""
        self.line = read_forward(self.handle)
        match = re.search(_RE_HIT_BLOCK_DESC, self.line)
        if not match:
            raise RuntimeError("Unexpected content in HIT_BLOCK_DESC line'{}'".format(self.line))
        hit_data = {
            'hit_id': match.group(1),
            'description': match.group(2).lstrip(' ;'),
            'query_start': None,
            'query_end': None,
            'query_seq': '',
            'hit_start': None,
            'hit_end': None,
            'hit_seq': ''
        }
        self.line = self.handle.readline()
        # E-value could be in decimal or scientific notation, so split the string rather then use regexp - this
        # also means we should be tolerant of additional fields being added/removed
        # Probab=99.95  E-value=3.7e-34  Score=210.31  Aligned_cols=171  Identities=100%  Similarity=2.050  Sum_probs=166.9
        scores = {}
        for score_pair in self.line.strip().split():
            key, value = score_pair.split('=')
            scores[key] = value
        hit_data['evalue'] = float(scores['E-value'])
        hit_data['score'] = float(scores['Score'])
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
            """Return True if match is not a Consensus column.

            It's not possible to distinguish a sequence line from a Consensus line with
            a regexp, so need to check the ID column.
            """
            return match.group(1).strip() != 'Consensus'

        while True:
            if not self.line.strip():  # blank lines indicate the end of a hit block
                return
            match = re.match(_RE_MATCH_BLOCK_QUERY_SEQ, self.line)
            if match and match_is_valid(match):
                hit_match_data['query_seq'] += match.group(3).strip()
                if hit_match_data['query_start'] is None:
                    hit_match_data['query_start'] = int(match.group(2))
                hit_match_data['query_end'] = int(match.group(4))
            else:
                match = re.match(_RE_MATCH_BLOCK_HIT_SEQ, self.line)
                if match and match_is_valid(match):
                    hit_match_data['hit_seq'] += match.group(3).strip()
                    if hit_match_data['hit_start'] is None:
                        hit_match_data['hit_start'] = int(match.group(2))
                    hit_match_data['hit_end'] = int(match.group(4))
            self.line = self.handle.readline()

    def _create_qresult(self, hit_blocks):
        """Create the Biopython data structures from the parsed data (PRIVATE)."""
        query_id = self.query_id
        hit_list = []
        hit_ids = set()

        count = 0
        for block in hit_blocks:
            count += 1
            hit_id = self._unique_hit_id(block['hit_id'], hit_ids)

            frag = HSPFragment(hit_id, query_id)
            frag.alphabet = generic_protein
            frag.query_start = block['query_start']
            frag.query_end = block['query_end']
            frag.hit_start = block['hit_start']
            frag.hit_end = block['hit_end']
            frag.hit = block['hit_seq']
            frag.query = block['query_seq']

            hsp = HSP([frag])
            hsp.hit_id = hit_id
            hsp.query_id = query_id
            hsp.is_included = True  # Should everything should be included?
            hsp.evalue = block['evalue']
            hsp.score = block['score']

            hit = Hit([hsp], hit_id)
            hit.description = block['description']
            hit.is_included = True  # Should everything should be included?
            hit.evalue = block['evalue']
            hit.score = block['score']
            hit_list.append(hit)

        qresult = QueryResult(hit_list, query_id)
        qresult.program = _PROGRAM
        qresult.seq_len = self.seq_len
        return [qresult]

    def _unique_hit_id(self, hit_id, existing_ids, separator='_'):
        """Return a unique hit id (PRIVATE).

        Always append a numeric id to each hit as there may be multiple with the same id.
        """
        i = 1
        new_id = "{}{}{}".format(hit_id, separator, i)
        while True:
            if new_id not in existing_ids:
                existing_ids.add(new_id)
                return new_id
            fields = new_id.split(separator)
            hit_id = separator.join(fields[0:-1])
            idx = int(fields[-1])
            new_id = "{}{}{}".format(hit_id, separator, idx + i)
            i += 1
