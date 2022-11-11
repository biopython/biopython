# Copyright 2012 by Kai Blin.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO parser for HMMER 2 text output."""

import re

from Bio.SearchIO._utils import read_forward
from Bio.SearchIO._model import QueryResult, Hit, HSP, HSPFragment

from ._base import _BaseHmmerTextIndexer

__all__ = ("Hmmer2TextParser", "Hmmer2TextIndexer")


_HSP_ALIGN_LINE = re.compile(r"(\S+):\s+domain (\d+) of (\d+)")


class _HitPlaceholder:
    def createHit(self, hsp_list):
        hit = Hit(hsp_list)
        hit.id_ = self.id_
        hit.evalue = self.evalue
        hit.bitscore = self.bitscore
        if self.description:
            hit.description = self.description
        hit.domain_obs_num = self.domain_obs_num
        return hit


class Hmmer2TextParser:
    """Iterator for the HMMER 2.0 text output."""

    def __init__(self, handle):
        """Initialize the class."""
        self.handle = handle
        self.buf = []
        self._meta = self.parse_preamble()

    def __iter__(self):
        """Iterate over Hmmer2TextParser, yields query results."""
        for qresult in self.parse_qresult():
            qresult.program = self._meta.get("program")
            qresult.target = self._meta.get("target")
            qresult.version = self._meta.get("version")
            yield qresult

    def read_next(self, rstrip=True):
        """Return the next non-empty line, trailing whitespace removed."""
        if len(self.buf) > 0:
            return self.buf.pop()
        self.line = self.handle.readline()
        while self.line and rstrip and not self.line.strip():
            self.line = self.handle.readline()
        if self.line:
            if rstrip:
                self.line = self.line.rstrip()
        return self.line

    def push_back(self, line):
        """Un-read a line that should not be parsed yet."""
        self.buf.append(line)

    def parse_key_value(self):
        """Parse key-value pair separated by colon."""
        key, value = self.line.split(":", 1)
        return key.strip(), value.strip()

    def parse_preamble(self):
        """Parse HMMER2 preamble."""
        meta = {}
        state = "GENERIC"
        while self.read_next():
            if state == "GENERIC":
                if self.line.startswith("hmm"):
                    meta["program"] = self.line.split("-")[0].strip()
                elif self.line.startswith("HMMER is"):
                    continue
                elif self.line.startswith("HMMER"):
                    meta["version"] = self.line.split()[1]
                elif self.line.count("-") == 36:
                    state = "OPTIONS"
                continue

            assert state == "OPTIONS"
            assert "program" in meta

            if self.line.count("-") == 32:
                break

            key, value = self.parse_key_value()
            if meta["program"] == "hmmsearch":
                if key == "Sequence database":
                    meta["target"] = value
                    continue
            elif meta["program"] == "hmmpfam":
                if key == "HMM file":
                    meta["target"] = value
                    continue
            meta[key] = value

        return meta

    def parse_qresult(self):
        """Parse a HMMER2 query block."""
        while self.read_next():
            if not self.line.startswith("Query"):
                return
            _, id_ = self.parse_key_value()
            self.qresult = QueryResult(id=id_)

            description = None

            while self.read_next() and not self.line.startswith("Scores"):
                if self.line.startswith("Accession"):
                    self.qresult.accession = self.parse_key_value()[1]
                if self.line.startswith("Description"):
                    description = self.parse_key_value()[1]

            hit_placeholders = self.parse_hits()
            if len(hit_placeholders) > 0:
                self.parse_hsps(hit_placeholders)
                self.parse_hsp_alignments()

            while not self.line.startswith("Query"):
                self.read_next()
                if not self.line:
                    break
            self.buf.append(self.line)

            if description is not None:
                self.qresult.description = description
            yield self.qresult

    def parse_hits(self):
        """Parse a HMMER2 hit block, beginning with the hit table."""
        hit_placeholders = []
        while self.read_next():
            if self.line.startswith("Parsed"):
                break
            if self.line.find("no hits") > -1:
                break

            if (
                self.line.startswith("Sequence")
                or self.line.startswith("Model")
                or self.line.startswith("-------- ")
            ):
                continue

            fields = self.line.split()
            id_ = fields.pop(0)
            domain_obs_num = int(fields.pop())
            evalue = float(fields.pop())
            bitscore = float(fields.pop())
            description = " ".join(fields).strip()

            hit = _HitPlaceholder()
            hit.id_ = id_
            hit.evalue = evalue
            hit.bitscore = bitscore
            hit.description = description
            hit.domain_obs_num = domain_obs_num
            hit_placeholders.append(hit)

        return hit_placeholders

    def parse_hsps(self, hit_placeholders):
        """Parse a HMMER2 hsp block, beginning with the hsp table."""
        # HSPs may occur in different order than the hits
        # so store Hit objects separately first
        unordered_hits = {}
        while self.read_next():
            if (
                self.line.startswith("Alignments")
                or self.line.startswith("Histogram")
                or self.line == "//"
            ):
                break
            if (
                self.line.startswith("Model")
                or self.line.startswith("Sequence")
                or self.line.startswith("--------")
            ):
                continue

            (
                id_,
                domain,
                seq_f,
                seq_t,
                seq_compl,
                hmm_f,
                hmm_t,
                hmm_compl,
                score,
                evalue,
            ) = self.line.split()

            frag = HSPFragment(id_, self.qresult.id)
            frag.molecule_type = "protein"
            if self._meta["program"] == "hmmpfam":
                frag.hit_start = int(hmm_f) - 1
                frag.hit_end = int(hmm_t)
                frag.query_start = int(seq_f) - 1
                frag.query_end = int(seq_t)
            elif self._meta["program"] == "hmmsearch":
                frag.query_start = int(hmm_f) - 1
                frag.query_end = int(hmm_t)
                frag.hit_start = int(seq_f) - 1
                frag.hit_end = int(seq_t)

            hsp = HSP([frag])
            hsp.evalue = float(evalue)
            hsp.bitscore = float(score)
            hsp.domain_index = int(domain.split("/")[0])
            if self._meta["program"] == "hmmpfam":
                hsp.hit_endtype = hmm_compl
                hsp.query_endtype = seq_compl
            elif self._meta["program"] == "hmmsearch":
                hsp.query_endtype = hmm_compl
                hsp.hit_endtype = seq_compl

            if id_ not in unordered_hits:
                placeholder = [p for p in hit_placeholders if p.id_ == id_][0]
                hit = placeholder.createHit([hsp])
                unordered_hits[id_] = hit
            else:
                hit = unordered_hits[id_]
                hsp.hit_description = hit.description
                hit.append(hsp)

        # The placeholder list is in the correct order, so use that order for
        # the Hit objects in the qresult
        for p in hit_placeholders:
            self.qresult.append(unordered_hits[p.id_])

    def parse_hsp_alignments(self):
        """Parse a HMMER2 HSP alignment block."""
        if not self.line.startswith("Alignments"):
            return

        while self.read_next():
            if self.line == "//" or self.line.startswith("Histogram"):
                break

            match = re.search(_HSP_ALIGN_LINE, self.line)
            if match is None:
                continue

            id_ = match.group(1)
            idx = int(match.group(2))
            num = int(match.group(3))

            hit = self.qresult[id_]
            if hit.domain_obs_num != num:
                continue

            frag = hit[idx - 1][0]

            hmmseq = ""
            consensus = ""
            otherseq = ""
            structureseq = ""
            pad = 0
            while self.read_next() and self.line.startswith(" "):
                # if there's structure information, parse that
                if self.line[16:18] == "CS":
                    structureseq += self.line[19:].strip()

                    if not self.read_next():
                        break

                # skip the *-> start marker if it exists
                if self.line[19:22] == "*->":
                    seq = self.line[22:]
                    pad = 3
                else:
                    seq = self.line[19:]
                    pad = 0

                hmmseq += seq
                line_len = len(seq)
                if not self.read_next(rstrip=False):
                    break
                consensus += self.line[19 + pad : 19 + pad + line_len]
                # If there's no consensus sequence, hmmer2 doesn't
                # bother to put spaces here, so add extra padding
                extra_padding = len(hmmseq) - len(consensus)
                consensus += " " * extra_padding

                if not self.read_next():
                    break

                # if we have a line break in the end marker, we get a
                # whitespace-only otherseq line, making split()[0] return
                # the end coordinate. That'll be a -, which is a valid character
                # in the sequence, meaning we can't just strip it.
                parts = self.line[19:].split()
                if len(parts) == 2:
                    otherseq += self.line[19:].split()[0].strip()

            self.push_back(self.line)

            # get rid of the end marker
            if hmmseq.endswith("<-*"):
                hmmseq = hmmseq[:-3]
                consensus = consensus[:-3]

            # add similarity sequence to annotation
            frag.aln_annotation["similarity"] = consensus

            # if there's structure information, add it to the fragment
            if structureseq:
                frag.aln_annotation["CS"] = structureseq

            if self._meta["program"] == "hmmpfam":
                frag.hit = hmmseq
                frag.query = otherseq
            else:
                frag.hit = otherseq
                frag.query = hmmseq


class Hmmer2TextIndexer(_BaseHmmerTextIndexer):
    """Indexer for hmmer2-text format."""

    _parser = Hmmer2TextParser
    qresult_start = b"Query"
    # qresults_ends for hmmpfam and hmmsearch
    # need to anticipate both since hmmsearch have different query end mark
    qresult_end = b"//"

    def __iter__(self):
        """Iterate over Hmmer2TextIndexer; yields query results' key, offsets, 0."""
        handle = self._handle
        handle.seek(0)
        start_offset = handle.tell()
        regex_id = re.compile(rb"Query\s*(?:sequence|HMM)?:\s*(.*)")

        # determine flag for hmmsearch
        is_hmmsearch = False
        line = read_forward(handle)
        if line.startswith(b"hmmsearch"):
            is_hmmsearch = True

        while True:
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
                # HACK: since hmmsearch can only have one query result
                if is_hmmsearch:
                    yield qresult_key.decode(), start_offset, 0
                break

            line = read_forward(handle)


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
