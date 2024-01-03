# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO parser for HMMER plain text output format."""

import re

from Bio.SearchIO._utils import read_forward, removesuffix
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment

from ._base import _BaseHmmerTextIndexer

__all__ = ("Hmmer3TextParser", "Hmmer3TextIndexer")


# precompile regex patterns for faster processing
# regex for program name capture
_RE_PROGRAM = re.compile(r"^# .*?(\w?hmm\w+) :: .*$")
# regex for version string capture
_RE_VERSION = re.compile(r"# \w+ ([\w+\.]+) .*; http.*$")
# regex for option string capture
_RE_OPT = re.compile(r"^# (.+):\s+(.+)$")
# regex for parsing query id and length, for parsing
_QRE_ID_LEN_PTN = r"^Query:\s*(.*)\s+\[\w=(\d+)\]"
_QRE_ID_LEN = re.compile(_QRE_ID_LEN_PTN)
# regex for hsp validation
_HRE_VALIDATE = re.compile(r"score:\s(-?\d+\.?\d+)\sbits.*value:\s(.*)")
# regexes for parsing hsp alignment blocks
_HRE_ANNOT_LINE = re.compile(r"^(\s+)(.+)\s(\w+)")
_HRE_ID_LINE = re.compile(r"^(\s+\S+\s+[0-9-]+ )(.+?)(\s+[0-9-]+)")


class Hmmer3TextParser:
    """Parser for the HMMER 3.0 text output."""

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        self.line = read_forward(self.handle)
        self._meta = self._parse_preamble()

    def __iter__(self):
        """Iterate over query results."""
        yield from self._parse_qresult()

    def _read_until(self, bool_func):
        """Read the file handle until the given function returns True (PRIVATE)."""
        while True:
            if not self.line or bool_func(self.line):
                return
            else:
                self.line = read_forward(self.handle)

    def _parse_preamble(self):
        """Parse HMMER preamble (lines beginning with '#') (PRIVATE)."""
        meta = {}
        # bool flag for storing state ~ whether we are parsing the option
        # lines or not
        has_opts = False
        while True:
            # no pound sign means we've left the preamble
            if not self.line.startswith("#"):
                break
            # dashes could either mean we are entering or leaving the options
            # section ~ so it's a switch for the has_opts flag
            elif "- - -" in self.line:
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
        """Parse a HMMER3 query block (PRIVATE)."""
        self._read_until(lambda line: line.startswith("Query:"))

        while self.line:
            regx = re.search(_QRE_ID_LEN, self.line)

            while not regx:
                self.line = read_forward(self.handle)
                regx = re.search(_QRE_ID_LEN, self.line)

            # get query id and length
            qid = regx.group(1).strip()
            # store qresult attributes
            qresult_attrs = {
                "seq_len": int(regx.group(2)),
                "program": self._meta.get("program"),
                "version": self._meta.get("version"),
                "target": self._meta.get("target"),
            }

            # get description and accession, if they exist
            qdesc = "<unknown description>"  # placeholder
            while not self.line.startswith("Scores for "):
                self.line = read_forward(self.handle)

                if self.line.startswith("Accession:"):
                    acc = self.line.strip().split(" ", 1)[1]
                    qresult_attrs["accession"] = acc.strip()
                elif self.line.startswith("Description:"):
                    qdesc = self.line.strip().split(" ", 1)[1].strip()
                    qresult_attrs["description"] = qdesc

            # parse the query hits
            while self.line and "//" not in self.line:
                hit_list = self._parse_hit(qid, qdesc)
                # read through the statistics summary
                # TODO: parse and store this information?
                if self.line.startswith("Internal pipeline"):
                    while self.line and "//" not in self.line:
                        self.line = read_forward(self.handle)

            # create qresult, set its attributes and yield
            # not initializing hit_list directly to handle empty hits
            # (i.e. need to set its query description manually)
            qresult = QueryResult(id=qid, hits=hit_list)
            for attr, value in qresult_attrs.items():
                setattr(qresult, attr, value)
            yield qresult
            self.line = read_forward(self.handle)

            # Skip line beginning with '# Alignment of', which are output
            # when running phmmer with the '-A' flag.
            if self.line.startswith("#"):
                self.line = self.handle.readline()

            # HMMER >= 3.1 outputs '[ok]' at the end of all results file,
            # which means we can break the main loop when we see the line
            if "[ok]" in self.line:
                break

    def _parse_hit(self, qid, qdesc):
        """Parse a HMMER3 hit block, beginning with the hit table (PRIVATE)."""
        # get to the end of the hit table delimiter and read one more line
        self._read_until(lambda line: line.startswith("    ------- ------ -----"))
        self.line = read_forward(self.handle)

        # assume every hit is in inclusion threshold until the inclusion
        # threshold line is encountered
        is_included = True

        # parse the hit table
        hit_attr_list = []
        while True:
            if not self.line:
                return []
            elif self.line.startswith("  ------ inclusion"):
                is_included = False
                self.line = read_forward(self.handle)
            # if there are no hits, then there are no hsps
            # so we forward-read until 'Internal pipeline..'
            elif self.line.startswith("   [No hits detected that satisfy reporting"):
                while True:
                    self.line = read_forward(self.handle)
                    if self.line.startswith("Internal pipeline"):
                        assert len(hit_attr_list) == 0
                        return []
            elif self.line.startswith("Domain annotation for each "):
                hit_list = self._create_hits(hit_attr_list, qid, qdesc)
                return hit_list
            # entering hit results row
            # parse the columns into a list
            row = [x for x in self.line.strip().split(" ") if x]
            # join the description words if it's >1 word
            if len(row) > 10:
                row[9] = " ".join(row[9:])
            # if there's no description, set it to an empty string
            elif len(row) < 10:
                row.append("")
                assert len(row) == 10
            # create the hit object
            hit_attrs = {
                "id": row[8],
                "query_id": qid,
                "evalue": float(row[0]),
                "bitscore": float(row[1]),
                "bias": float(row[2]),
                # row[3:6] is not parsed, since the info is available
                # at the HSP level
                "domain_exp_num": float(row[6]),
                "domain_obs_num": int(row[7]),
                "description": row[9],
                "is_included": is_included,
            }
            hit_attr_list.append(hit_attrs)

            self.line = read_forward(self.handle)

    def _create_hits(self, hit_attrs, qid, qdesc):
        """Parse a HMMER3 hsp block, beginning with the hsp table (PRIVATE)."""
        # read through until the beginning of the hsp block
        self._read_until(lambda line: line.startswith(("Internal pipeline", ">>")))

        # start parsing the hsp block
        hit_list = []
        while True:
            if self.line.startswith("Internal pipeline"):
                # by this time we should've emptied the hit attr list
                assert len(hit_attrs) == 0
                return hit_list
            assert self.line.startswith(">>")
            hid, hdesc = self.line[len(">> ") :].split("  ", 1)
            hdesc = hdesc.strip()

            # read through the hsp table header and move one more line
            self._read_until(
                lambda line: line.startswith(
                    (" ---   ------ ----- --------", "   [No individual domains")
                )
            )
            self.line = read_forward(self.handle)

            # parse the hsp table for the current hit
            hsp_list = []
            while True:
                # break out of hsp parsing if there are no hits, it's the last hsp
                # or it's the start of a new hit
                if (
                    self.line.startswith("   [No targets detected that satisfy")
                    or self.line.startswith("   [No individual domains")
                    or self.line.startswith("Internal pipeline statistics summary:")
                    or self.line.startswith("  Alignments for each domain:")
                    or self.line.startswith(">>")
                ):
                    hit_attr = hit_attrs.pop(0)
                    hit = Hit(hsp_list)
                    for attr, value in hit_attr.items():
                        if attr == "description":
                            cur_val = getattr(hit, attr)
                            if cur_val and value and cur_val.startswith(value):
                                continue
                        setattr(hit, attr, value)
                    if not hit:
                        hit.query_description = qdesc
                    hit_list.append(hit)
                    break

                parsed = [x for x in self.line.strip().split(" ") if x]
                assert len(parsed) == 16
                # parsed column order:
                # index, is_included, bitscore, bias, evalue_cond, evalue
                # hmmfrom, hmmto, query_ends, hit_ends, alifrom, alito,
                # envfrom, envto, acc_avg
                frag = HSPFragment(hid, qid)
                # set query and hit descriptions if they are defined / nonempty string
                if qdesc:
                    frag.query_description = qdesc
                if hdesc:
                    frag.hit_description = hdesc
                # HMMER3 results are always protein
                frag.molecule_type = "protein"
                # depending on whether the program is hmmsearch, hmmscan, or phmmer
                # {hmm,ali}{from,to} can either be hit_{from,to} or query_{from,to}
                # for hmmscan, hit is the hmm profile, query is the sequence
                if self._meta.get("program") == "hmmscan":
                    # adjust 'from' and 'to' coordinates to 0-based ones
                    frag.hit_start = int(parsed[6]) - 1
                    frag.hit_end = int(parsed[7])
                    frag.query_start = int(parsed[9]) - 1
                    frag.query_end = int(parsed[10])
                elif self._meta.get("program") in ["hmmsearch", "phmmer"]:
                    # adjust 'from' and 'to' coordinates to 0-based ones
                    frag.hit_start = int(parsed[9]) - 1
                    frag.hit_end = int(parsed[10])
                    frag.query_start = int(parsed[6]) - 1
                    frag.query_end = int(parsed[7])
                # strand is always 0, since HMMER now only handles protein
                frag.hit_strand = frag.query_strand = 0

                hsp = HSP([frag])
                hsp.domain_index = int(parsed[0])
                hsp.is_included = parsed[1] == "!"
                hsp.bitscore = float(parsed[2])
                hsp.bias = float(parsed[3])
                hsp.evalue_cond = float(parsed[4])
                hsp.evalue = float(parsed[5])
                if self._meta.get("program") == "hmmscan":
                    # adjust 'from' and 'to' coordinates to 0-based ones
                    hsp.hit_endtype = parsed[8]
                    hsp.query_endtype = parsed[11]
                elif self._meta.get("program") in ["hmmsearch", "phmmer"]:
                    # adjust 'from' and 'to' coordinates to 0-based ones
                    hsp.hit_endtype = parsed[11]
                    hsp.query_endtype = parsed[8]
                # adjust 'from' and 'to' coordinates to 0-based ones
                hsp.env_start = int(parsed[12]) - 1
                hsp.env_end = int(parsed[13])
                hsp.env_endtype = parsed[14]
                hsp.acc_avg = float(parsed[15])

                hsp_list.append(hsp)
                self.line = read_forward(self.handle)

            # parse the hsp alignments
            if self.line.startswith("  Alignments for each domain:"):
                self._parse_aln_block(hid, hit.hsps)

    def _parse_aln_block(self, hid, hsp_list):
        """Parse a HMMER3 HSP alignment block (PRIVATE)."""
        self.line = read_forward(self.handle)
        dom_counter = 0
        while True:
            if self.line.startswith(">>") or self.line.startswith("Internal pipeline"):
                return hsp_list
            assert self.line.startswith("  == domain %i" % (dom_counter + 1))
            # alias hsp to local var
            # but note that we're still changing the attrs of the actual
            # hsp inside the qresult as we're not creating a copy
            frag = hsp_list[dom_counter][0]
            # XXX: should we validate again here? regex is expensive..
            # regx = re.search(_HRE_VALIDATE, self.line)
            # assert hsp.bitscore == float(regx.group(1))
            # assert hsp.evalue_cond == float(regx.group(2))
            hmmseq = ""
            aliseq = ""
            annot = {}
            aln_prefix_len = None
            self.line = self.handle.readline()

            # parse all the alignment blocks in the hsp
            while True:
                regx = None

                # check for hit or query line
                # we don't check for the hit or query id specifically
                # to anticipate special cases where query id == hit id
                regx = re.search(_HRE_ID_LINE, self.line)
                if regx:
                    # Try to capture the block prefix len, to parse the similarity
                    # string later.
                    if aln_prefix_len is None:
                        aln_prefix_len = len(regx.group(1))
                    else:
                        assert aln_prefix_len == len(regx.group(1))
                    # the first hit/query self.line we encounter is the hmmseq
                    if len(hmmseq) == len(aliseq):
                        hmmseq += regx.group(2)
                    # and for subsequent self.lines, len(hmmseq) is either
                    # > or == len(aliseq)
                    elif len(hmmseq) > len(aliseq):
                        aliseq += regx.group(2)
                    assert len(hmmseq) >= len(aliseq)
                # check for start of new domain
                elif (
                    self.line.startswith("  == domain")
                    or self.line.startswith(">>")
                    or self.line.startswith("Internal pipeline")
                ):
                    frag.aln_annotation = annot
                    if self._meta.get("program") == "hmmscan":
                        frag.hit = hmmseq
                        frag.query = aliseq
                    elif self._meta.get("program") in ["hmmsearch", "phmmer"]:
                        frag.hit = aliseq
                        frag.query = hmmseq
                    dom_counter += 1
                    hmmseq = ""
                    aliseq = ""
                    annot = {}
                    aln_prefix_len = None
                    break
                # check if it's an annotation line and parse it
                # len(hmmseq) is only != len(aliseq) when the cursor is parsing
                # the similarity character.
                elif len(hmmseq) == len(aliseq):
                    regx = re.search(_HRE_ANNOT_LINE, self.line)
                    if regx:
                        annot_name = regx.group(3)
                        if annot_name in annot:
                            annot[annot_name] += regx.group(2)
                        else:
                            annot[annot_name] = regx.group(2)
                # otherwise, assume we are seeing the similarity string (the string
                # between the two alignments). We do that by checking if we have
                # defined the aln_prefix_len variable, which is the length of string
                # from the beginning of the line to the start of the actual alignment
                # (sans the IDs). We need this string to see how many characters we
                # should strip from the beginning of the line. We can't call `strip()`
                # since the similarity string might start with empty characters.
                elif aln_prefix_len is not None:
                    similarity = removesuffix(
                        removesuffix(self.line[aln_prefix_len:], "\n"),
                        "\r",
                    )
                    if "similarity" not in annot:
                        annot["similarity"] = similarity
                    else:
                        annot["similarity"] += similarity

                self.line = self.handle.readline()


class Hmmer3TextIndexer(_BaseHmmerTextIndexer):
    """Indexer class for HMMER plain text output."""

    _parser = Hmmer3TextParser
    qresult_start = b"Query: "
    qresult_end = b"//"

    def __iter__(self):
        """Iterate over Hmmer3TextIndexer; yields query results' key, offsets, 0."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()
        regex_id = re.compile(_QRE_ID_LEN_PTN.encode())

        while True:
            line = read_forward(handle)
            end_offset = handle.tell()

            if line.startswith(self.qresult_start):
                regx = re.search(regex_id, line)
                qresult_key = regx.group(1).strip()
                # qresult start offset is the offset of this line
                # (starts with the start mark)
                start_offset = end_offset - len(line)
            elif line.startswith(self.qresult_end):
                yield qresult_key.decode(), start_offset, 0
                start_offset = end_offset
            elif not line:
                break


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
