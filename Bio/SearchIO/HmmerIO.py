# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SearchIO parser for HMMER output formats.

This module adds support for parsing HMMER outputs, from version 3.0 onwards.
HMMER is a suite of programs implementing the profile hidden Markov models
to find homology across protein sequences.

Specifically, this module supports the following HMMER output formats:

  - Plain text - 'hmmer-text'

And the following HMMER programs: hmmersearch, hmmerscan

More information are available through these links:
  - Web page: http://hmmer.janelia.org/
  - User guide: ftp://selab.janelia.org/pub/software/hmmer3/3.0/Userguide.pdf

"""

import re

from Bio.SearchIO._objects import QueryResult, Hit, HSP, SearchIndexer


# precompile regex patterns for faster processing
# regex for program name capture
re_program = re.compile(r'^# (\w*hmm\w+) :: .*$')
# regex for version string capture
re_version = re.compile(r'# \w+ ([\w+\.]+) .*; http.*$')
# regex for option string capture
re_opt = re.compile(r'^# (.+):\s+(.+)$')
# regex for parsing query id and length
qre_id_len = re.compile(r'^Query:\s*(.*)\s+\[\w=(\d+)\]')
# regex for hsp validation
hre_validate = re.compile(r'score:\s(-?\d+\.?\d+)\sbits.*value:\s(.*)')
# regexes for parsing hsp alignment blocks
hre_annot_line = re.compile(r'^(\s+)(.+)\s(\w+)')
hre_id_line = re.compile(r'^(\s+\S+\s+\d+ )(.+) (\d+)')


def hmmer_text_iterator(handle):
    """Generator function to parse HMMER plain text output as QueryResult objects.

    handle -- Handle to the file, or the filename as string.

    """
    iterator = HmmerTextIterator(handle)
    for qresult in iterator:
        yield qresult


class HmmerTextIterator(object):

    """Iterator for the HMMER 3.0 text output."""

    def __init__(self, handle):
        self.handle = handle
        self.meta = {}
        self.line = self.read_forward()
        self.parse_preamble()

    def __iter__(self):
        for qresult in self.parse_qresult():
            qresult.meta = self.meta
            qresult.program = self.meta['program']
            qresult.target = self.meta['target']
            yield qresult

    def read_forward(self):
        """Reads through whitespaces at the beginning of handle, returns the first
        non-whitespace line."""
        while True:
            line = self.handle.readline()
            # if line has characters and stripping does not remove them,
            # return the line
            if line and line.strip():
                return line
            elif line and not line.strip():
                continue
            # if line ends, return None
            elif not line:
                return

    def parse_preamble(self):
        """Parses HMMER preamble (lines beginning with '#')."""
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
                # XXX: is there a better way to do the regex check?
                # try parsing program
                regx = re.search(re_program, self.line)
                if regx:
                    self.meta['program'] = regx.group(1)
                # try parsing version
                regx = re.search(re_version, self.line)
                if regx:
                    self.meta['program_version'] = regx.group(1)
            elif has_opts:
                regx = re.search(re_opt, self.line)
                # if target in regx.group(1), then we store the key as target
                if 'target' in regx.group(1):
                    self.meta['target'] = regx.group(2)
                else:
                    self.meta[regx.group(1)] = regx.group(2)

            self.line = self.read_forward()

    def parse_qresult(self):
        """Parses a HMMER3 query block."""
        while True:
            if self.line.startswith('Query:'):
                break
            else:
                self.line = self.read_forward()

        while True:

            assert self.line.startswith('Query:')

            # get query id and length
            regx = re.search(qre_id_len, self.line)
            id = regx.group(1).strip()
            seq_len = regx.group(2)
            # create qresult object
            self.qresult = QueryResult(id)
            self.qresult.seq_len = seq_len
            self.qresult.meta = self.meta
            self.qresult.program = self.meta['program']
            self.qresult.target = self.meta['target']
            # get description, if it exists
            self.line = self.read_forward()
            if self.line.startswith('Description:'):
                desc = self.line.split(' ', 1)[1]
                self.qresult.desc = desc

            while True:
                self.parse_hit()
                if not self.line:
                    break
                if self.line.startswith('Internal pipeline'):
                    while True:
                        self.line = self.read_forward()
                        if '//' in self.line or not self.line:
                            break
                if '//' in self.line or not self.line:
                    break

            yield self.qresult
            self.line = self.read_forward()
            if not self.line:
                break

    def parse_hit(self):
        """Parses a HMMER3 hit block, beginning with the hit table."""
        # get to the end of the hit table delimiter
        while True:
            if not self.line:
                break
            elif self.line.startswith('    ------- ------ -----'):
                self.line = self.read_forward()
                break
            self.line = self.read_forward()

        # assume every hit is in inclusion threshold until the inclusion
        # threshold line is encountered
        is_in_inclusion = True

        # parse the hit table
        while True:
            if not self.line:
                break
            elif self.line.startswith('  ------ inclusion'):
                is_in_inclusion = False
                self.line = self.read_forward()
                continue
            # if there are no hits, then there are no hsps
            # so we forward-read until 'Internal pipeline..'
            elif self.line.startswith('   [No hits detected that satisfy reporting'):
                while True:
                    self.line = self.read_forward()
                    if self.line.startswith('Internal pipeline'):
                        return
            elif self.line.startswith('Domain annotation for each '):
                self.parse_hsp()
                break
            # entering hit results row
            # parse the columns into a list
            row = filter(None, self.line.strip().split(' '))
            try:
                assert len(row) == 10
            except AssertionError:
                # If length is not 10, then there are >1 words in the
                # description column ~ which should be concatenated with
                # the first description word
                extra_desc = ' '.join(row[10:])
                row = row[:11]
                row[9] = '%s %s' % (row[9], extra_desc)
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
            # don't forget to attach the boolean is_in_inclusion
            hit.is_in_inclusion = is_in_inclusion

            self.qresult.append(hit)

            self.line = self.read_forward()

    def parse_hsp(self):
        """Parses a HMMER3 hsp block, beginning with the hsp table."""
        # read through until the beginning of the hsp block
        while True:
            if not self.line or self.line.startswith('Internal pipeline'):
                break
            elif self.line.startswith('Domain annotation for each '):
                self.line = self.read_forward()
                break
            self.line = self.read_forward()

        if not self.line or self.line.startswith('   [No targets detected'):
            return
        elif self.line.startswith('Internal pipeline'):
            return

        # start parsing the hsp block
        while True:
            if self.line.startswith('Internal pipeline'):
                break
            assert self.line.startswith('>>')
            hid, hdesc = self.line[len('>> '):].split('  ', 1)
            # ensure that the parsed hit ID is in the qresult object
            assert hid in self.qresult, "Unexpected hit ID: %s" % hid

            # read through the hsp table header
            while True:
                if self.line.startswith(' ---   ------ ----- --------'):
                    self.line = self.read_forward()
                    break
                self.line = self.read_forward()

            # parse the hsp table for the current hit
            while True:
                if self.line.startswith('   [No targets detected that satisfy') or \
                        self.line.startswith('Internal pipeline statistics summary:') or \
                        self.line.startswith('  Alignments for each domain:'):
                    break

                parsed = filter(None, self.line.strip().split(' '))
                assert len(parsed) == 16
                # parsed column order:
                # index, is_in_inclusion, bitscore, bias, evalue_cond, evalue
                # hmmfrom, hmmto, query_ends, hit_ends, alifrom, alito,
                # envfrom, envto, acc_avg
                hsp = HSP(hid, self.qresult.id)
                hsp.index = parsed[0]
                hsp.is_in_inclusion = parsed[1] == '!'
                hsp.bitscore = parsed[2]
                hsp.bias = parsed[3]
                hsp.evalue_cond = parsed[4]
                hsp.evalue = parsed[5]
                # depending on whether the program is hmmsearch or hmmscan,
                # {hmm,ali}{from,to} can either be hit_{from,to} or query_{from,to}
                # for hmmscan, hit is the hmm profile, query is the sequence
                if self.meta['program'] == 'hmmscan':
                    hsp.hit_from = parsed[6]
                    hsp.hit_to = parsed[7]
                    hsp.hit_endtype = parsed[8]
                    hsp.query_from = parsed[9]
                    hsp.query_to = parsed[10]
                    hsp.query_endtype = parsed[11]
                else:
                    hsp.hit_from = parsed[9]
                    hsp.hit_to = parsed[10]
                    hsp.hit_endtype = parsed[11]
                    hsp.query_from = parsed[6]
                    hsp.query_to = parsed[7]
                    hsp.query_endtype = parsed[8]
                hsp.env_from = parsed[12]
                hsp.env_to = parsed[13]
                hsp.env_endtype = parsed[14]
                hsp.acc_avg = parsed[15]

                self.qresult[hid].append(hsp)
                self.line = self.read_forward()

            # parse the hsp alignments
            if self.line.startswith('  Alignments for each domain:'):
                self.parse_hsp_alignment(hid)

    def parse_hsp_alignment(self, hid):
        """Parses a HMMER3 HSP alignment block."""
        self.line = self.read_forward()
        dom_counter = 0
        while True:
            if self.line.startswith('>>') or \
                    self.line.startswith('Internal pipeline'):
                break
            assert self.line.startswith('  == domain %i' % (dom_counter + 1))
            # alias hsp to local var
            # but not that we're still changing the attrs of the actual
            # hsp inside the qresult as we're not creating a copy
            hsp = self.qresult[hid][dom_counter]
            # XXX: should we validate again here? regex is expensive..
            #regx = re.search(hre_validate, self.line)
            #assert hsp.bitscore == float(regx.group(1))
            #assert hsp.evalue_cond == float(regx.group(2))
            spc_len = None
            ali_len = None
            hmmseq = ''
            aliseq = ''
            annot = {}
            self.line = self.handle.readline()

            # parse all the alignment blocks in the hsp
            while True:

                regx = None

                if spc_len is None:
                    # assume that self.line is a query or hit line
                    try:
                        regx = re.search(hre_id_line, self.line)
                        spc_len = len(regx.group(1))
                    # if not, then the self.line is an annotation line
                    except AttributeError:
                        regx = re.search(hre_annot_line, self.line)
                        spc_len = len(regx.group(1))

                # check for hit or query line
                # we don't check for the hit or query id specifically yet
                # to anticipate special cases where query id == hit id

                # check for query id
                regx = re.search(hre_id_line, self.line)
                if regx:
                    # the first hit/query self.line we encounter is the hmmseq
                    if len(hmmseq) == len(aliseq):
                        hmmseq += regx.group(2)
                    # and for subsequent self.lines, len(hmmseq) is either > or == len(aliseq)
                    elif len(hmmseq) > len(aliseq):
                        aliseq += regx.group(2)
                    assert len(hmmseq) >= len(aliseq)
                # homology self.line is only encountered when len(hmmseq) > len(aliseq)
                # and its length is len(hmmseq) - len(aliseq)
                elif len(hmmseq) > len(aliseq):
                    # note that ali_len changes depending on the difference of hmmseq
                    # and aliseq
                    ali_len = len(hmmseq) - len(aliseq)
                    if 'homology' in annot:
                        annot['homology'] += self.line[spc_len:spc_len+ali_len]
                    else:
                        annot['homology'] = self.line[spc_len:spc_len+ali_len]
                # check for blank self.line
                elif not self.line.strip():
                    pass
                # check for start of new domain
                elif self.line.startswith('  == domain') or \
                        self.line.startswith('>>') or \
                        self.line.startswith('Internal pipeline'):
                    hsp.alignment_annotation = annot
                    if self.meta['program'] == 'hmmscan':
                        hsp.add_alignment(hmmseq, aliseq)
                        hsp.hit.description = 'HMM sequence'
                        hsp.query.description = 'protein sequence'
                    else:
                        hsp.add_alignment(aliseq, hmmseq)
                        hsp.hit.description = 'protein sequence'
                        hsp.query.description = 'HMM sequence'
                    dom_counter += 1
                    hmmseq = ''
                    aliseq = ''
                    spc_len = None
                    ali_len = None
                    annot = {}
                    break
                # otherwise it's an annotation self.line
                else:
                    regx = re.search(hre_annot_line, self.line)
                    annot_name = regx.group(3)
                    if annot_name in annot:
                        annot[annot_name] += regx.group(2)
                    else:
                        annot[annot_name] = regx.group(2)

                self.line = self.handle.readline()


class HmmerTextIndexer(SearchIndexer):

    """Indexer class for HMMER plain text output."""

    def __init__(self, handle):
        pass



def _test():
    """Run the Bio.SearchIO.HmmerIO module's doctests.

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
