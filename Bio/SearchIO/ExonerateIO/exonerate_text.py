# Copyright 2012 by Wibowo Arindrarto.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SearchIO parser for Exonerate plain text output format."""

import re
from itertools import chain


from ._base import (
    _BaseExonerateParser,
    _BaseExonerateIndexer,
    _STRAND_MAP,
    _parse_hit_or_query_line,
)
from .exonerate_vulgar import _RE_VULGAR


__all__ = ("ExonerateTextParser", "ExonerateTextIndexer")


# for capturing sequences in alignment blocks
# e.g. ' 529 : ATCCCTTATCTCTTTATCTTGTA :    472'
_RE_ALN_ROW = re.compile(r"\s*\d+\s+: (.*) :\s+\d+")
# for splitting the line based on intron annotations
# e.g. '  >>>> Target Intron 1 >>>>  ' or 'gt.........................ag'
_RE_EXON = re.compile(
    r"[atgc ]{2}?(?:(?:[<>]+ \w+ Intron \d+ [<>]+)|(?:\.+))[atgc ]{2}?"
)
# captures the intron length
# from e.g. '61 bp // 154295 bp' (joint intron lengths) or '177446 bp'
_RE_EXON_LEN = re.compile(r"(?:(\d+) bp // (\d+) bp)|(?:(\d+) bp)")
# for splitting lines in the NER model
_RE_NER = re.compile(r"--<\s+\d+\s+>--")
# for capturing NER gap lengths
_RE_NER_LEN = re.compile(r"--<\s+(\d+)\s+>--")
# regexes for capturing the letters inside curly braces
# no. of letters is either 1 or 2, since they are split codons
_RE_SCODON_START = re.compile(r"\{(\w{1,2})\}$")
_RE_SCODON_END = re.compile(r"^\{(\w{1,2})\}")


def _flip_codons(codon_seq, target_seq):
    """Flips the codon characters from one seq to another (PRIVATE)."""
    a, b = "", ""
    for char1, char2 in zip(codon_seq, target_seq):
        # no need to do anything if the codon seq line has nothing
        if char1 == " ":
            a += char1
            b += char2
        else:
            a += char2
            b += char1

    return a, b


def _get_block_coords(parsed_seq, row_dict, has_ner=False):
    """Return a list of start, end coordinates for each given block in the sequence (PRIVATE)."""
    start = 0
    coords = []
    if not has_ner:
        splitter = _RE_EXON
    else:
        splitter = _RE_NER

    # use the query line for reference
    seq = parsed_seq[row_dict["query"]]

    for block in re.split(splitter, seq):
        start += seq[start:].find(block)
        end = start + len(block)
        coords.append((start, end))

    return coords


def _get_inter_coords(coords, strand=1):
    """Return list of pairs covering intervening ranges (PRIVATE).

    From the given pairs of coordinates, returns a list of pairs
    covering the intervening ranges.
    """
    # adapted from Python's itertools guide
    # if strand is -1, adjust coords to the ends and starts are chained
    if strand == -1:
        sorted_coords = [(max(a, b), min(a, b)) for a, b in coords]
        inter_coords = list(chain(*sorted_coords))[1:-1]
        return list(zip(inter_coords[1::2], inter_coords[::2]))
    else:
        inter_coords = list(chain(*coords))[1:-1]
        return list(zip(inter_coords[::2], inter_coords[1::2]))


def _stitch_rows(raw_rows):
    """Stitches together the parsed alignment rows and returns them in a list (PRIVATE)."""
    # deal with possible codon surprise!
    # (i.e. alignments with codons using cdna2genome model)
    # by creating additional rows to contain the codons
    try:
        max_len = max(len(x) for x in raw_rows)
        for row in raw_rows:
            assert len(row) == max_len
    except AssertionError:
        for idx, row in enumerate(raw_rows):
            if len(row) != max_len:
                # codons must be present in the query and hit (so +2)
                assert len(row) + 2 == max_len
                # add additional empty lines to contain codons
                raw_rows[idx] = [" " * len(row[0])] + row + [" " * len(row[0])]

    cmbn_rows = []
    for idx, row in enumerate(raw_rows[0]):
        cmbn_row = "".join(aln_row[idx] for aln_row in raw_rows)
        cmbn_rows.append(cmbn_row)

    # the real aligned sequence is always the 'outer' one, so we want
    # to flip them with their 'inner' pairs
    if len(cmbn_rows) == 5:
        # flip query sequence
        cmbn_rows[0], cmbn_rows[1] = _flip_codons(cmbn_rows[0], cmbn_rows[1])
        # flip hit sequence
        cmbn_rows[4], cmbn_rows[3] = _flip_codons(cmbn_rows[4], cmbn_rows[3])

    return cmbn_rows


def _get_row_dict(row_len, model):
    """Return a dictionary of row indices for parsing alignment blocks (PRIVATE)."""
    idx = {}
    # 3 lines, usually in dna vs dna models
    if row_len == 3:
        idx["query"] = 0
        idx["midline"] = 1
        idx["hit"] = 2
        idx["qannot"], idx["hannot"] = None, None
    # 4 lines, in protein vs dna models or dna vs protein models
    # TODO: currently we check this from the model string; is there
    # a better way to do it?
    elif row_len == 4:
        if "protein2" in model:
            idx["query"] = 0
            idx["midline"] = 1
            idx["hit"] = 2
            idx["hannot"] = 3
            idx["qannot"] = None
        elif "2protein" in model:
            idx["query"] = 1
            idx["midline"] = 2
            idx["hit"] = 3
            idx["hannot"] = None
            idx["qannot"] = 0
        else:
            raise ValueError("Unexpected model: " + model)
    # 5 lines, translated dna vs translated dna
    elif row_len == 5:
        # set sequence indexes
        idx["qannot"] = 0
        idx["query"] = 1
        idx["midline"] = 2
        idx["hit"] = 3
        idx["hannot"] = 4
    else:
        raise ValueError("Unexpected row count in alignment block: %i" % row_len)
    return idx


def _get_blocks(rows, coords, idx):
    """Return a list of dictionaries of sequences split by the coordinates (PRIVATE)."""
    for idx_name in ("query", "hit", "midline", "qannot", "hannot"):
        assert idx_name in idx
    blocks = []
    for start, end in coords:
        block = {}
        # get seqs according to index
        block["query"] = rows[idx["query"]][start:end]
        block["hit"] = rows[idx["hit"]][start:end]
        block["similarity"] = rows[idx["midline"]][start:end]
        if idx["qannot"] is not None:
            block["query_annotation"] = rows[idx["qannot"]][start:end]
        if idx["hannot"] is not None:
            block["hit_annotation"] = rows[idx["hannot"]][start:end]
        blocks.append(block)

    return blocks


def _get_scodon_moves(tmp_seq_blocks):
    """Get a dictionary of split codon locations relative to each fragment end (PRIVATE)."""
    scodon_moves = {"query": [], "hit": []}
    for seq_type in scodon_moves:
        scoords = []
        for block in tmp_seq_blocks:
            # check both ends of the sequence for residues in curly braces
            m_start = re.search(_RE_SCODON_START, block[seq_type])
            m_end = re.search(_RE_SCODON_END, block[seq_type])
            if m_start:
                m_start = len(m_start.group(1))
                scoords.append((m_start, 0))
            else:
                scoords.append((0, 0))
            if m_end:
                m_end = len(m_end.group(1))
                scoords.append((0, m_end))
            else:
                scoords.append((0, 0))
        scodon_moves[seq_type] = scoords

    return scodon_moves


def _clean_blocks(tmp_seq_blocks):
    """Remove curly braces (split codon markers) from the given sequences (PRIVATE)."""
    seq_blocks = []
    for seq_block in tmp_seq_blocks:
        for line_name in seq_block:
            seq_block[line_name] = (
                seq_block[line_name].replace("{", "").replace("}", "")
            )
        seq_blocks.append(seq_block)

    return seq_blocks


def _comp_intron_lens(seq_type, inter_blocks, raw_inter_lens):
    """Return the length of introns between fragments (PRIVATE)."""
    # set opposite type, for setting introns
    opp_type = "hit" if seq_type == "query" else "query"
    # list of flags to denote if an intron follows a block
    # it reads e.g. this line:
    # "ATGTT{TT}  >>>> Target Intron 1 >>>>  {G}TGTGTGTACATT"
    # and sets the opposing sequence type's intron (since this
    # line is present on the opposite sequence type line)
    has_intron_after = ["Intron" in x[seq_type] for x in inter_blocks]
    assert len(has_intron_after) == len(raw_inter_lens)
    # create list containing coord adjustments incorporating
    # intron lengths
    inter_lens = []
    for flag, parsed_len in zip(has_intron_after, raw_inter_lens):
        if flag:
            # joint introns
            if all(parsed_len[:2]):
                # intron len is [0] if opp_type is query, otherwise it's [1]
                intron_len = (
                    int(parsed_len[0]) if opp_type == "query" else int(parsed_len[1])
                )
            # single hit/query introns
            elif parsed_len[2]:
                intron_len = int(parsed_len[2])
            else:
                raise ValueError("Unexpected intron parsing result: %r" % parsed_len)
        else:
            intron_len = 0

        inter_lens.append(intron_len)

    return inter_lens


def _comp_coords(hsp, seq_type, inter_lens):
    """Fill the block coordinates of the given hsp dictionary (PRIVATE)."""
    assert seq_type in ("hit", "query")
    # manually fill the first coord
    seq_step = 1 if hsp["%s_strand" % seq_type] >= 0 else -1
    fstart = hsp["%s_start" % seq_type]
    # fend is fstart + number of residues in the sequence, minus gaps
    fend = (
        fstart
        + len(hsp[seq_type][0].replace("-", "").replace(">", "").replace("<", ""))
        * seq_step
    )
    coords = [(fstart, fend)]
    # and start from the second block, after the first inter seq
    for idx, block in enumerate(hsp[seq_type][1:]):
        bstart = coords[-1][1] + inter_lens[idx] * seq_step
        bend = bstart + seq_step * len(block.replace("-", ""))
        coords.append((bstart, bend))

    # adjust the coords so the smallest is [0], if strand is -1
    # couldn't do this in the previous steps since we need the initial
    # block ordering
    if seq_step != 1:
        for idx, coord in enumerate(coords):
            coords[idx] = coords[idx][1], coords[idx][0]

    return coords


def _comp_split_codons(hsp, seq_type, scodon_moves):
    """Compute positions of split codons, store in given HSP dictionary (PRIVATE)."""
    scodons = []
    for idx in range(len(scodon_moves[seq_type])):
        pair = scodon_moves[seq_type][idx]
        if not any(pair):
            continue
        else:
            assert not all(pair)
        a, b = pair
        anchor_pair = hsp["%s_ranges" % seq_type][idx // 2]
        strand = 1 if hsp["%s_strand" % seq_type] >= 0 else -1

        if a:
            func = max if strand == 1 else min
            anchor = func(anchor_pair)
            start_c, end_c = anchor + a * strand * -1, anchor
        elif b:
            func = min if strand == 1 else max
            anchor = func(anchor_pair)
            start_c, end_c = anchor + b * strand, anchor
        scodons.append((min(start_c, end_c), max(start_c, end_c)))

    return scodons


class ExonerateTextParser(_BaseExonerateParser):
    """Parser for Exonerate plain text output."""

    _ALN_MARK = "C4 Alignment:"

    def parse_alignment_block(self, header):
        """Parse alignment block, return query result, hits, hsps."""
        qresult = header["qresult"]
        hit = header["hit"]
        hsp = header["hsp"]
        # check for values that must have been set by previous methods
        for val_name in (
            "query_start",
            "query_end",
            "hit_start",
            "hit_end",
            "query_strand",
            "hit_strand",
        ):
            assert val_name in hsp, hsp

        # get the alignment rows
        # and stitch them so we have the full sequences in single strings
        raw_aln_blocks, vulgar_comp = self._read_alignment()
        # cmbn_rows still has split codon markers (curly braces)
        cmbn_rows = _stitch_rows(raw_aln_blocks)
        row_dict = _get_row_dict(len(cmbn_rows), qresult["model"])
        # get the sequence blocks
        has_ner = "NER" in qresult["model"].upper()
        seq_coords = _get_block_coords(cmbn_rows, row_dict, has_ner)
        tmp_seq_blocks = _get_blocks(cmbn_rows, seq_coords, row_dict)
        # get split codon temp coords for later use
        # this result in pairs of base movement for both ends of each row
        scodon_moves = _get_scodon_moves(tmp_seq_blocks)
        # remove the split codon markers
        seq_blocks = _clean_blocks(tmp_seq_blocks)

        # adjust strands
        hsp["query_strand"] = _STRAND_MAP[hsp["query_strand"]]
        hsp["hit_strand"] = _STRAND_MAP[hsp["hit_strand"]]
        # cast coords into ints
        hsp["query_start"] = int(hsp["query_start"])
        hsp["query_end"] = int(hsp["query_end"])
        hsp["hit_start"] = int(hsp["hit_start"])
        hsp["hit_end"] = int(hsp["hit_end"])
        # cast score into ints
        hsp["score"] = int(hsp["score"])
        # set sequences
        hsp["query"] = [x["query"] for x in seq_blocks]
        hsp["hit"] = [x["hit"] for x in seq_blocks]
        hsp["aln_annotation"] = {}
        # set the molecule type
        # currently only limited to models with protein queries
        if (
            "protein2" in qresult["model"]
            or "coding2" in qresult["model"]
            or "2protein" in qresult["model"]
        ):
            hsp["molecule_type"] = "protein"
        # get the annotations if they exist
        for annot_type in ("similarity", "query_annotation", "hit_annotation"):
            try:
                hsp["aln_annotation"][annot_type] = [x[annot_type] for x in seq_blocks]
            except KeyError:
                pass

        # use vulgar coordinates if vulgar line is present and return
        # if vulgar_comp is not None:
        #    hsp = parse_vulgar_comp(hsp, vulgar_comp)

        #    return {'qresult': qresult, 'hit': hit, 'hsp': hsp}

        # otherwise we need to get the coordinates from the alignment
        # get the intervening blocks first, so we can use them
        # to adjust the coordinates
        if not has_ner:
            # get intervening coordinates and blocks, only if model is not ner
            # ner models have a much more simple coordinate calculation
            inter_coords = _get_inter_coords(seq_coords)
            inter_blocks = _get_blocks(cmbn_rows, inter_coords, row_dict)
            # returns a three-component tuple of intron lengths
            # first two component filled == intron in hit and query
            # last component filled == intron in hit or query
            raw_inter_lens = re.findall(_RE_EXON_LEN, cmbn_rows[row_dict["midline"]])

        # compute start and end coords for each block
        for seq_type in ("query", "hit"):

            # ner blocks and intron blocks require different adjustments
            if not has_ner:
                opp_type = "hit" if seq_type == "query" else "query"
                inter_lens = _comp_intron_lens(seq_type, inter_blocks, raw_inter_lens)
            else:
                # for NER blocks, the length of the inter-fragment gaps is
                # written on the same strand, so opp_type is seq_type
                opp_type = seq_type
                inter_lens = [
                    int(x)
                    for x in re.findall(_RE_NER_LEN, cmbn_rows[row_dict[seq_type]])
                ]

            # check that inter_lens's length is len opp_type block - 1
            if len(inter_lens) != len(hsp[opp_type]) - 1:
                raise ValueError(
                    "Length mismatch: %r vs %r"
                    % (len(inter_lens), len(hsp[opp_type]) - 1)
                )
            # fill the hsp query and hit coordinates
            hsp["%s_ranges" % opp_type] = _comp_coords(hsp, opp_type, inter_lens)
            # and fill the split codon coordinates, if model != ner
            # can't do this in the if-else clause above since we need to
            # compute the ranges first
            if not has_ner:
                hsp["%s_split_codons" % opp_type] = _comp_split_codons(
                    hsp, opp_type, scodon_moves
                )

        # now that we've finished parsing coords, we can set the hit and start
        # coord according to Biopython's convention (start <= end)
        for seq_type in ("query", "hit"):
            if hsp["%s_strand" % seq_type] == -1:
                n_start = "%s_start" % seq_type
                n_end = "%s_end" % seq_type
                hsp[n_start], hsp[n_end] = hsp[n_end], hsp[n_start]

        return {"qresult": qresult, "hit": hit, "hsp": hsp}

    def _read_alignment(self):
        """Read the raw alignment block strings, returns them in a list (PRIVATE)."""
        raw_aln_blocks = []
        # flag to check whether we're in an alignment row
        in_aln_row = False
        # flag for vulgar line, if present, we can parse coordinates from it
        vulgar_comp = None
        while True:

            match = re.search(_RE_ALN_ROW, self.line.strip())
            # if we have a match, set flags and values
            if match and not in_aln_row:
                start_idx = self.line.index(match.group(1))
                row_len = len(match.group(1))
                in_aln_row = True
                raw_aln_block = []
            # if we're in an alignment row, grab the sequence
            if in_aln_row:
                raw_aln_block.append(self.line[start_idx : start_idx + row_len])
            # reset flags and values if the line matches, we're in an alignment
            # row, and there are more than 1 line in rows
            if match and in_aln_row and len(raw_aln_block) > 1:
                raw_aln_blocks.append(raw_aln_block)
                start_idx = None
                row_len = None
                in_aln_row = False

            self.line = self.handle.readline()
            # try to parse vulgar line if present
            if self.line.startswith("vulgar"):
                vulgar = re.search(_RE_VULGAR, self.line)
                vulgar_comp = vulgar.group(10)
            if not self.line or self.line.startswith(self._ALN_MARK):
                # HACK: this is so that the parse_qresult method does not
                # yield the objects before appending the last HSP. We are doing
                # this to keep the parser compatible with outputs without
                # human-readable alignment outputs. This also relies on the
                # fact that repeated readline() always returns '' on EOF.
                if not self.line:
                    self.line = "mock"
                break

        return raw_aln_blocks, vulgar_comp


class ExonerateTextIndexer(_BaseExonerateIndexer):
    """Indexer class for Exonerate plain text."""

    _parser = ExonerateTextParser
    _query_mark = b"C4 Alignment"

    def get_qresult_id(self, pos):
        """Return the query ID from the nearest "Query:" line."""
        handle = self._handle
        handle.seek(pos)
        sentinel = b"Query:"

        while True:
            line = handle.readline().strip()
            if line.startswith(sentinel):
                break
            if not line:
                raise StopIteration
        qid, desc = _parse_hit_or_query_line(line.decode())

        return qid

    def get_raw(self, offset):
        """Return the raw string of a QueryResult object from the given offset."""
        handle = self._handle
        handle.seek(offset)
        qresult_key = None
        qresult_raw = b""

        while True:
            line = handle.readline()
            if not line:
                break
            elif line.startswith(self._query_mark):
                cur_pos = handle.tell()
                if qresult_key is None:
                    qresult_key = self.get_qresult_id(cur_pos)
                else:
                    curr_key = self.get_qresult_id(cur_pos)
                    if curr_key != qresult_key:
                        break
                handle.seek(cur_pos)
            qresult_raw += line

        return qresult_raw


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
