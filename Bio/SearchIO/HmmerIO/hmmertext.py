# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER plain text output format."""

import re

from Bio._py3k import _as_bytes, _bytes_to_string
from Bio.SearchIO._objects import QueryResult, Hit, HSP
from Bio.SearchIO._index import SearchIndexer


# precompile regex patterns for faster processing
# regex for program name capture
_RE_PROGRAM = re.compile(r'^# (\w*hmm\w+) :: .*$')
# regex for version string capture
_RE_VERSION = re.compile(r'# \w+ ([\w+\.]+) .*; http.*$')
# regex for option string capture
_RE_OPT = re.compile(r'^# (.+):\s+(.+)$')
# regex for parsing query id and length
_QRE_ID_LEN = re.compile(r'^Query:\s*(.*)\s+\[\w=(\d+)\]')
# regex for hsp validation
_HRE_VALIDATE = re.compile(r'score:\s(-?\d+\.?\d+)\sbits.*value:\s(.*)')
# regexes for parsing hsp alignment blocks
_HRE_ANNOT_LINE = re.compile(r'^(\s+)(.+)\s(\w+)')
_HRE_ID_LINE = re.compile(r'^(\s+\S+\s+\d+ )(.+) (\d+)')


def read_forward(handle):
    """Reads through whitespaces, returns the first non-whitespace line."""
    while True:
        line = handle.readline()
        # if line has characters and stripping does not remove them,
        # return the line
        if line and line.strip():
            return line
        # if line ends, return None
        elif not line:
            return line


class HmmerTextIterator(object):

    """Iterator for the HMMER 3.0 text output."""

    def __init__(self, handle):
        self.handle = handle
        self.line = read_forward(self.handle)
        self._meta = self.parse_preamble()

    def __iter__(self):
        for qresult in self.parse_qresult():
            qresult.program = self._meta.get('program')
            qresult.target = self._meta.get('target')
            yield qresult

    def read_until(self, bool_func):
        """Reads the file handle until the given function returns True."""
        while True:
            if not self.line or bool_func(self.line):
                return
            else:
                self.line = read_forward(self.handle)

    def parse_preamble(self):
        """Parses HMMER preamble (lines beginning with '#')."""
        meta = {}
        # bool flag for storing state ~ whether we are parsing the option self.lines
        # or not
        has_opts = False
        # store self.line as a local var
        while True:
            # no pound sign means we've left the preamble
            if not self.line.startswith('#'):
                break
            # dashes could either mean we are entering or leaving the options
            # section ~ so it's a switch for the has_opts flag
            elif '- - -' in self.line:
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
                    meta['program'] = regx.group(1)
                # try parsing version
                regx = re.search(_RE_VERSION, self.line)
                if regx:
                    meta['version'] = regx.group(1)
            elif has_opts:
                regx = re.search(_RE_OPT, self.line)
                # if target in regx.group(1), then we store the key as target
                if 'target' in regx.group(1):
                    meta['target'] = regx.group(2)
                else:
                    meta[regx.group(1)] = regx.group(2)

            self.line = read_forward(self.handle)

        return meta

    def parse_qresult(self):
        """Parses a HMMER3 query block."""

        self.read_until(lambda line: line.startswith('Query:'))

        while self.line:

            # get query id and length
            regx = re.search(_QRE_ID_LEN, self.line)
            id = regx.group(1).strip()
            seq_len = regx.group(2)
            # create qresult object
            self.qresult = QueryResult(id)
            self.qresult.seq_len = seq_len
            self.qresult.program = self._meta.get('program')
            self.qresult.version = self._meta.get('version')
            self.qresult.target = self._meta.get('target')
            
            # get description and accession, if they exist
            desc = '' # placeholder
            while not self.line.startswith('Scores for '):
                self.line = read_forward(self.handle)

                if self.line.startswith('Accession:'):
                    acc = self.line.strip().split(' ', 1)[1]
                    self.qresult.acc = acc.strip()
                elif self.line.startswith('Description:'):
                    desc = self.line.strip().split(' ', 1)[1]

            # parse the query hits
            while self.line and '//' not in self.line:
                self.parse_hit()
                # read through the statistics summary
                # TODO: parse and store this information?
                if self.line.startswith('Internal pipeline'):
                    while self.line and '//' not in self.line:
                        self.line = read_forward(self.handle)

            # append desc here, so hsp attributes are also set
            self.qresult.desc = desc.strip()
            yield self.qresult
            self.line = read_forward(self.handle)

    def parse_hit(self):
        """Parses a HMMER3 hit block, beginning with the hit table."""
        # get to the end of the hit table delimiter and read one more line
        self.read_until(lambda line: \
                line.startswith('    ------- ------ -----'))
        self.line = read_forward(self.handle)

        # assume every hit is in inclusion threshold until the inclusion
        # threshold line is encountered
        is_included = True

        # parse the hit table
        while True:
            if not self.line:
                break
            elif self.line.startswith('  ------ inclusion'):
                is_included = False
                self.line = read_forward(self.handle)
            # if there are no hits, then there are no hsps
            # so we forward-read until 'Internal pipeline..'
            elif self.line.startswith('   [No hits detected that satisfy reporting'):
                while True:
                    self.line = read_forward(self.handle)
                    if self.line.startswith('Internal pipeline'):
                        return
            elif self.line.startswith('Domain annotation for each '):
                self.parse_hsp()
                break
            # entering hit results row
            # parse the columns into a list
            row = filter(None, self.line.strip().split(' '))
            # join the description words if it's >1 word
            if len(row) > 10:
                row[9] = ' '.join(row[9:])
            # if there's no description, set it to an empty string
            elif len(row) < 10:
                row.append('')
                assert len(row) == 10
            # create the hit object
            hit_id = row[8]
            hit = Hit(hit_id, self.qresult.id)
            # store the parsed results appropriately
            hit.evalue = row[0]
            hit.bitscore = row[1]
            hit.bias = row[2]
            # row[3:6] is not parsed, since the info is available the the HSP level
            hit.domain_exp_num = row[6]
            hit.domain_obs_num = row[7]
            hit.desc = row[9]
            # don't forget to attach the boolean is_included
            hit.is_included = is_included

            self.qresult.append(hit)

            self.line = read_forward(self.handle)

    def parse_hsp(self):
        """Parses a HMMER3 hsp block, beginning with the hsp table."""
        # read through until the beginning of the hsp block
        self.read_until(lambda line: line.startswith('Internal pipeline') \
                or line.startswith('>>'))

        # start parsing the hsp block
        while True:
            if self.line.startswith('Internal pipeline'):
                break
            assert self.line.startswith('>>')
            hid, hdesc = self.line[len('>> '):].split('  ', 1)
            # ensure that the parsed hit ID is in the qresult object
            assert hid in self.qresult, "Unexpected hit ID: %s" % hid

            # read through the hsp table header and move one more line
            self.read_until(lambda line: \
                    line.startswith(' ---   ------ ----- --------'))
            self.line = read_forward(self.handle)

            # parse the hsp table for the current hit
            while True:
                # break out of hsp parsing if there are no hits, it's the last hsp
                # or it's the start of a new hit
                if self.line.startswith('   [No targets detected that satisfy') or \
                        self.line.startswith('Internal pipeline statistics summary:') or \
                        self.line.startswith('  Alignments for each domain:') or \
                        self.line.startswith('>>'):
                    break

                parsed = filter(None, self.line.strip().split(' '))
                assert len(parsed) == 16
                # parsed column order:
                # index, is_included, bitscore, bias, evalue_cond, evalue
                # hmmfrom, hmmto, query_ends, hit_ends, alifrom, alito,
                # envfrom, envto, acc_avg
                hsp = HSP(hid, self.qresult.id)
                hsp.domain_index = parsed[0]
                hsp.is_included = parsed[1] == '!'
                hsp.bitscore = parsed[2]
                hsp.bias = parsed[3]
                hsp.evalue_cond = parsed[4]
                hsp.evalue = parsed[5]
                # depending on whether the program is hmmsearch, hmmscan, or phmmer
                # {hmm,ali}{from,to} can either be hit_{from,to} or query_{from,to}
                # for hmmscan, hit is the hmm profile, query is the sequence
                if self._meta.get('program') == 'hmmscan':
                    # adjust 'from' and 'to' coordinates to 0-based ones
                    hsp.hit_from = int(parsed[6]) - 1
                    hsp.hit_to = int(parsed[7]) - 1
                    hsp.query_from = int(parsed[9]) - 1
                    hsp.query_to = int(parsed[10]) - 1
                    hsp.hit_endtype = parsed[8]
                    hsp.query_endtype = parsed[11]
                elif self._meta.get('program') in ['hmmsearch', 'phmmer']:
                    # adjust 'from' and 'to' coordinates to 0-based ones
                    hsp.hit_from = int(parsed[9]) - 1
                    hsp.hit_to = int(parsed[10]) - 1
                    hsp.query_from = int(parsed[6]) - 1
                    hsp.query_to = int(parsed[7]) - 1
                    hsp.hit_endtype = parsed[11]
                    hsp.query_endtype = parsed[8]
                # adjust 'from' and 'to' coordinates to 0-based ones
                hsp.env_from = int(parsed[12]) - 1
                hsp.env_to = int(parsed[13]) - 1
                hsp.env_endtype = parsed[14]
                hsp.acc_avg = parsed[15]

                self.qresult[hid].append(hsp)
                self.line = read_forward(self.handle)

            # parse the hsp alignments
            if self.line.startswith('  Alignments for each domain:'):
                self.parse_hsp_alignment(hid)

    def parse_hsp_alignment(self, hid):
        """Parses a HMMER3 HSP alignment block."""
        self.line = read_forward(self.handle)
        dom_counter = 0
        while True:
            if self.line.startswith('>>') or \
                    self.line.startswith('Internal pipeline'):
                break
            assert self.line.startswith('  == domain %i' % (dom_counter + 1))
            # alias hsp to local var
            # but note that we're still changing the attrs of the actual
            # hsp inside the qresult as we're not creating a copy
            hsp = self.qresult[hid][dom_counter]
            # XXX: should we validate again here? regex is expensive..
            #regx = re.search(_HRE_VALIDATE, self.line)
            #assert hsp.bitscore == float(regx.group(1))
            #assert hsp.evalue_cond == float(regx.group(2))
            hmmseq = ''
            aliseq = ''
            annot = {}
            self.line = self.handle.readline()

            # parse all the alignment blocks in the hsp
            while True:

                regx = None

                # check for hit or query line
                # we don't check for the hit or query id specifically
                # to anticipate special cases where query id == hit id
                regx = re.search(_HRE_ID_LINE, self.line)
                if regx:
                    # the first hit/query self.line we encounter is the hmmseq
                    if len(hmmseq) == len(aliseq):
                        hmmseq += regx.group(2)
                    # and for subsequent self.lines, len(hmmseq) is either > or == len(aliseq)
                    elif len(hmmseq) > len(aliseq):
                        aliseq += regx.group(2)
                    assert len(hmmseq) >= len(aliseq)
                # check for start of new domain
                elif self.line.startswith('  == domain') or \
                        self.line.startswith('>>') or \
                        self.line.startswith('Internal pipeline'):
                    hsp.alignment_annotation = annot
                    if self._meta.get('program') == 'hmmscan':
                        hsp.hit = hmmseq
                        hsp.query = aliseq
                    elif self._meta.get('program') in ['hmmsearch', 'phmmer']:
                        hsp.hit = aliseq
                        hsp.query = hmmseq
                    hsp.hit.description = self.qresult[hid].desc
                    hsp.query.description = self.qresult.desc
                    dom_counter += 1
                    hmmseq = ''
                    aliseq = ''
                    annot = {}
                    break
                # otherwise check if it's an annotation line and parse it
                # len(hmmseq) is only != len(aliseq) when the cursor is parsing
                # the homology character. Since we're not parsing that, we
                # check for when the condition is False (i.e. when it's ==)
                elif len(hmmseq) == len(aliseq):
                    regx = re.search(_HRE_ANNOT_LINE, self.line)
                    if regx:
                        annot_name = regx.group(3)
                        if annot_name in annot:
                            annot[annot_name] += regx.group(2)
                        else:
                            annot[annot_name] = regx.group(2)

                self.line = self.handle.readline()


class HmmerTextIndexer(SearchIndexer):

    """Indexer class for HMMER plain text output."""

    _parser = HmmerTextIterator
    qresult_start = _as_bytes('Query: ')
    qresult_end = _as_bytes('//')

    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()

        while True:
            line = read_forward(handle)
            end_offset = handle.tell()

            if line.startswith(self.qresult_start):
                regx = re.search(_QRE_ID_LEN, line)
                qresult_key = regx.group(1).strip()
                # qresult start offset is the offset of this line
                # (starts with the start mark)
                start_offset = end_offset - len(line)
            elif line.startswith(self.qresult_end):
                yield _bytes_to_string(qresult_key), start_offset, \
                        end_offset - start_offset
                start_offset = end_offset
            elif not line:
                break

    def get_raw(self, offset):
        handle = self._handle
        qresult_raw = ''

        # read header first
        handle.seek(0)
        while True:
            line = handle.readline()
            if line.startswith(self.qresult_start):
                break
            qresult_raw += line

        # and read the qresult raw string
        handle.seek(offset)
        while True:
            # preserve whitespace, don't use read_forward
            line = handle.readline()
            qresult_raw += line

            # break when we've reached qresult end
            if line.startswith(self.qresult_end):
                break

        return qresult_raw


def _test():
    """Run the Bio.SearchIO.HmmerIO.hmmertext module's doctests.

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
