# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO parser for Exonerate vulgar output format."""

import re

from ._base import _BaseExonerateParser, _BaseExonerateIndexer, _STRAND_MAP


__all__ = ("ExonerateVulgarParser", "ExonerateVulgarIndexer")


# precompile regex
_RE_VULGAR = re.compile(
    r"""^vulgar:\s+
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # query: ID, start, end, strand
        (\S+)\s+(\d+)\s+(\d+)\s+([\+-\.])\s+  # hit: ID, start, end, strand
        (\d+)(\s+.*)$                         # score, vulgar components
        """,
    re.VERBOSE,
)

_RE_VCOMP = re.compile(
    r"""
        \s+(\S+) # vulgar label (C/M: codon/match, G: gap, N: ner, 5/3: splice
                 #               site, I: intron, S: split codon, F: frameshift)
        \s+(\d+) # how many residues to advance in query sequence
        \s+(\d+) # how many residues to advance in hit sequence
        """,
    re.VERBOSE,
)


def parse_vulgar_comp(hsp, vulgar_comp):
    """Parse the vulgar components present in the hsp dictionary."""
    # containers for block coordinates
    qstarts = [hsp["query_start"]]
    qends = []
    hstarts = [hsp["hit_start"]]
    hends = []
    # containers for split codons
    hsp["query_split_codons"] = []
    hsp["hit_split_codons"] = []
    # containers for ner blocks
    hsp["query_ner_ranges"] = []
    hsp["hit_ner_ranges"] = []
    # sentinels for tracking query and hit positions
    qpos = hsp["query_start"]
    hpos = hsp["hit_start"]
    # multiplier for determining sentinel movement
    qmove = 1 if hsp["query_strand"] >= 0 else -1
    hmove = 1 if hsp["hit_strand"] >= 0 else -1

    vcomps = re.findall(_RE_VCOMP, vulgar_comp)
    for idx, match in enumerate(vcomps):
        label, qstep, hstep = match[0], int(match[1]), int(match[2])
        # check for label, must be recognized
        assert label in "MCGF53INS", "Unexpected vulgar label: %r" % label
        # match, codon, or gaps
        if label in "MCGS":
            # if the previous comp is not an MCGS block, it's the
            # start of a new block
            if vcomps[idx - 1][0] not in "MCGS":
                qstarts.append(qpos)
                hstarts.append(hpos)
        # other labels
        # store the values in the hsp dict as a tuple of (start, stop)
        # we're not doing anything if the label is in '53IN', as these
        # basically tell us what the inter-block coordinates are and
        # inter-block coordinates are automatically calculated by
        # and HSP property
        if label == "S":
            # get start and stop from parsed values
            qstart, hstart = qpos, hpos
            qend = qstart + qstep * qmove
            hend = hstart + hstep * hmove
            # adjust the start-stop ranges
            sqstart, sqend = min(qstart, qend), max(qstart, qend)
            shstart, shend = min(hstart, hend), max(hstart, hend)
            # split codons
            # XXX: is it possible to have a frameshift that introduces
            # a codon split? If so, this may need a different treatment..
            hsp["query_split_codons"].append((sqstart, sqend))
            hsp["hit_split_codons"].append((shstart, shend))

        # move sentinels accordingly
        qpos += qstep * qmove
        hpos += hstep * hmove

        # append to ends if the next comp is not an MCGS block or
        # if it's the last comp
        if idx == len(vcomps) - 1 or (
            label in "MCGS" and vcomps[idx + 1][0] not in "MCGS"
        ):
            qends.append(qpos)
            hends.append(hpos)

    # adjust coordinates
    for seq_type in ("query_", "hit_"):
        strand = hsp[seq_type + "strand"]
        # switch coordinates if strand is < 0
        if strand < 0:
            # switch the starts and ends
            hsp[seq_type + "start"], hsp[seq_type + "end"] = (
                hsp[seq_type + "end"],
                hsp[seq_type + "start"],
            )
            if seq_type == "query_":
                qstarts, qends = qends, qstarts
            else:
                hstarts, hends = hends, hstarts

    # set start and end ranges
    hsp["query_ranges"] = list(zip(qstarts, qends))
    hsp["hit_ranges"] = list(zip(hstarts, hends))
    return hsp


class ExonerateVulgarParser(_BaseExonerateParser):
    """Parser for Exonerate vulgar strings."""

    _ALN_MARK = "vulgar"

    def parse_alignment_block(self, header):
        """Parse alignment block for vulgar format, return query results, hits, hsps."""
        qresult = header["qresult"]
        hit = header["hit"]
        hsp = header["hsp"]
        self.read_until(lambda line: line.startswith("vulgar"))
        vulgars = re.search(_RE_VULGAR, self.line)
        # if the file has c4 alignments
        # check if vulgar values match our previously parsed header values
        if self.has_c4_alignment:
            assert qresult["id"] == vulgars.group(1)
            assert hsp["query_start"] == vulgars.group(2)
            assert hsp["query_end"] == vulgars.group(3)
            assert hsp["query_strand"] == vulgars.group(4)
            assert hit["id"] == vulgars.group(5)
            assert hsp["hit_start"] == vulgars.group(6)
            assert hsp["hit_end"] == vulgars.group(7)
            assert hsp["hit_strand"] == vulgars.group(8)
            assert hsp["score"] == vulgars.group(9)
        else:
            qresult["id"] = vulgars.group(1)
            hsp["query_start"] = vulgars.group(2)
            hsp["query_end"] = vulgars.group(3)
            hsp["query_strand"] = vulgars.group(4)
            hit["id"] = vulgars.group(5)
            hsp["hit_start"] = vulgars.group(6)
            hsp["hit_end"] = vulgars.group(7)
            hsp["hit_strand"] = vulgars.group(8)
            hsp["score"] = vulgars.group(9)

        # adjust strands
        hsp["hit_strand"] = _STRAND_MAP[hsp["hit_strand"]]
        hsp["query_strand"] = _STRAND_MAP[hsp["query_strand"]]
        # cast coords into ints
        hsp["query_start"] = int(hsp["query_start"])
        hsp["query_end"] = int(hsp["query_end"])
        hsp["hit_start"] = int(hsp["hit_start"])
        hsp["hit_end"] = int(hsp["hit_end"])
        # cast score into int
        hsp["score"] = int(hsp["score"])
        # store vulgar line and parse it
        # rstrip to remove line endings (otherwise gives errors in Windows)
        hsp["vulgar_comp"] = vulgars.group(10).rstrip()
        hsp = parse_vulgar_comp(hsp, hsp["vulgar_comp"])

        return {"qresult": qresult, "hit": hit, "hsp": hsp}


class ExonerateVulgarIndexer(_BaseExonerateIndexer):
    """Indexer class for exonerate vulgar lines."""

    _parser = ExonerateVulgarParser
    _query_mark = b"vulgar"

    def get_qresult_id(self, pos):
        """Return the query ID of the nearest vulgar line."""
        handle = self._handle
        handle.seek(pos)
        # get line, check if it's a vulgar line, and get query ID
        line = handle.readline()
        assert line.startswith(self._query_mark), line
        id = re.search(_RE_VULGAR, line.decode())
        return id.group(1)

    def get_raw(self, offset):
        """Return the raw bytes string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        qresult_key = None
        qresult_raw = b""

        while True:
            line = handle.readline()
            if not line:
                break
            elif line.startswith(self._query_mark):
                cur_pos = handle.tell() - len(line)
                if qresult_key is None:
                    qresult_key = self.get_qresult_id(cur_pos)
                else:
                    curr_key = self.get_qresult_id(cur_pos)
                    if curr_key != qresult_key:
                        break
            qresult_raw += line

        return qresult_raw


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
