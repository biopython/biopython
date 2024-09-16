# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.SearchIO parser for Infernal plain text output format."""

import re
import operator

from Bio.SearchIO._model import Hit
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult
from Bio.SearchIO._utils import read_forward

__all__ = ("InfernalTextParser")

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
_DIV_HITS_END = "Internal CM pipeline statistics summary:"
_DIV_HIT_SCORE = "Hit scores:"
_DIV_HIT_ALIGNMENT = "Hit alignments:"
_DIV_NO_HIT = "   [No hits detected that satisfy reporting thresholds]"
_DIV_TABLE_START = " ----   --------- ------"
# hit alignment
_DIV_ALIGNMENT_START = ">> "

class InfernalTextParser:
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
            else:
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
            while not self.line.startswith(_DIV_HIT_SCORE):
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
                if self.line.startswith(_DIV_HITS_END):
                    while self.line and not self.line.startswith(_DIV_QUERY_END):
                        self.line = read_forward(self.handle)
            # create qresult, set its attributes and yield   
            #print(hit_list)
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
        """Parse an Infernal hit section (PRIVATE)."""
        # empty container
        hit_list, hsp_list = [], []
        # dummie for hit informations
        hit_attrs = None
        prev_hid = None
        cur_hid = None

        # the hit score and alignments tables are redundant, if the output 
        # does not include the alignment, we must parse the hit score table
        # otherwise we can skip it and parse the hit alignment only
        if self._meta["show alignments in output"] == "no":
            # parse the hit table
            while True:
                if not self.line:
                    return []
                # if there are no hits, forward-read to the end of the query
                elif self.line.startswith(_DIV_NO_HIT):
                    while True:
                        self.line = read_forward(self.handle)
                        if self.line.startswith(_DIV_HITS_END):
                            return []
                elif self.line.startswith(_DIV_HIT_SCORE):
                    # read through the header 
                    self._read_until(lambda line: line.startswith(_DIV_TABLE_START))
                    self.line = read_forward(self.handle)

                    # parse the hit score table
                    while True:
                        # we've reached the end of the hit score table
                        if self.line.startswith(_DIV_HITS_END):
                            # create the last hit 
                            hit = Hit(hsp_list)
                            for attr, value in hit_attrs.items():
                                setattr(hit, attr, value)
                            hit_list.append(hit)
                            return hit_list

                        # parse the columns into a list
                        row = [x for x in self.line.strip().split(" ") if x]
                        # join the description words if it's >1 word
                        if len(row) > 12:
                            row[12] = " ".join(row[12:])
                        # if there's no description, set it to an empty string
                        elif len(row) < 12:
                            row.append("")
                            assert len(row) == 12
                        
                        # create hit and append to hit container
                        cur_hid = row[5]
                        if prev_hid is not None and cur_hid != prev_hid:
                            hit = Hit(hsp_list)
                            for attr, value in hit_attrs.items():
                                setattr(hit, attr, value)
                            hit_list.append(hit)
                            hsp_list = []
                        
                        # parse the attributes
                        hit_attrs = {
                            "id": cur_hid,
                            "query_id": qid,
                            "description": row[12]
                        }
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
                            "hit_strand": 0 if row[8] == "+" else -1
                        }

                        # create the hsp fragment and set it's attributes
                        hsp_frag = HSPFragment(row[5], qid)
                        for attr, value in hsp_frag_attrs.items():
                            setattr(hsp_frag, attr, value)
                        
                        # create the hsp and set it's attributes
                        hsp = HSP([hsp_frag])
                        for attr, value in hsp_attrs.items():
                            setattr(hsp, attr, value)
                        hsp_list.append(hsp)

                        prev_hid = hit_attrs["id"]
                        self.line = read_forward(self.handle)
        else:
            # skip the hit score table 
            self._read_until(lambda line: line.startswith(_DIV_HIT_ALIGNMENT))
            self.line = read_forward(self.handle)

            while True:
                if not self.line:
                    return []
                # if there are no hits, forward-read to the end of the query
                elif self.line.startswith(_DIV_NO_HIT):
                    while True:
                        self.line = read_forward(self.handle)
                        if self.line.startswith(_DIV_HITS_END):
                            return []
                # we've reached the end of the hit section
                elif self.line.startswith(_DIV_HITS_END):
                    # create the last hit 
                    hit = Hit(hsp_list)
                    for attr, value in hit_attrs.items():
                        setattr(hit, attr, value)
                    hit_list.append(hit)
                    return hit_list
                # entering hit alignment table
                elif self.line.startswith(_DIV_ALIGNMENT_START):
                    hid, hdesc = self.line[len(_DIV_ALIGNMENT_START) :].split("  ", 1)
                    hdesc = hdesc.strip()

                    # read through the hit table header and move one more line
                    self._read_until(lambda line: line.startswith(_DIV_TABLE_START))
                    self.line = read_forward(self.handle)

                    # parse the hit table
                    row = [x for x in self.line.strip().split() if x]
                    assert len(row) == 16

                    # create hit and append to hit container
                    cur_hid = hid
                    if prev_hid is not None and cur_hid != prev_hid:
                        hit = Hit(hsp_list)
                        for attr, value in hit_attrs.items():
                            setattr(hit, attr, value)
                        hit_list.append(hit)
                        hsp_list = []

                    hit_attrs = {
                        "id": hid,
                        "query_id": qid,
                        "description": hdesc
                    }
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
                        "is_included": True if row[1] == "!" else False
                    }
                    query_start = int(row[6])
                    query_end = int(row[7])
                    hit_start = int(row[9]) if row[11] == "+" else int(row[10])
                    hit_end = int(row[10]) if row[11] == "+" else int(row[9])
                    hit_strand = 0 if row[11] == "+" else -1

                    # move to the HSP alignment block
                    self.line = read_forward(self.handle)

                    # create the hsp
                    frag_list = self._parse_aln_block(hit_attrs["id"], hit_attrs["query_id"], \
                        hsp_attrs["model"], query_start, query_end, hit_start, hit_end, hit_strand)
                    hsp = HSP(frag_list)
                    for attr, value in hsp_attrs.items():
                        setattr(hsp, attr, value)
                    hsp_list.append(hsp)

                    prev_hid = hit_attrs["id"]
                    # create hit and append to hit container
                    #hit = Hit([hsp])
                    #for attr, value in hit_attrs.items():
                    #    setattr(hit, attr, value)
                    #hit_list.append(hit)
        
        return hit_list


    def _parse_aln_block(self, hid, qid, model, query_start, query_end, hit_start, hit_end, hit_strand):
        """Parse a Infernal HSP alignment block (PRIVATE)."""
        frag_list = []
        model_seq = ""
        hit_seq = ""
        annot = {"NC": "", "CS": "", "similarity": "", "PP": ""}
        while True:
            # we've reached the end of the alignment section
            if self.line.startswith(_DIV_ALIGNMENT_START) or self.line.startswith(_DIV_HITS_END):
                # Process local end in infernal hit alignment. Local end are
                # large insertion or deletion indicated by *[NN]* where N is 
                # the number of model positions are deleted or the number of
                # residues are inserted in the sequence. We split the
                # alignment block in HSPs on these local ends. 

                # get local ends string (*[NN]*) indexes in the model sequence
                # there can be more than one local ends back-to-back
                local_aln_idx = [(0,0)] # there is always at least one local alignment starting at 0
                local_aln_idx += [(m.start(0), m.end(0)) for m in re.finditer(_RE_SPLIT_ALN, model_seq)]
                                
                prev_hit_start = hit_start if hit_strand == 0 else hit_end
                prev_model_start = query_start
                hsps = []
                
                for i in range(len(local_aln_idx)):
                    local_start = local_aln_idx[i][1]
                    local_end = local_aln_idx[i+1][0] if i+1 < len(local_aln_idx) else None

                    # hsp hit position. forward strand moves up, reverse strand moves down
                    op = operator.add if hit_strand == 0 else operator.sub
                    cur_hit_seq = hit_seq[local_start:local_end]
                    cur_hit_gap_size = self._local_aln_gap_size(local_aln_idx[i], hit_seq)
                    cur_hit_start = op(prev_hit_start, cur_hit_gap_size)
                    cur_hit_end = op(cur_hit_start, len(re.sub(_RE_LETTERS, "", cur_hit_seq)))
                    # adjust start and end elements positon 
                    if hit_strand == 0 and i == 0:
                        cur_hit_end -= 1
                    if hit_strand == -1 and i == 0:
                        cur_hit_end += 1
                    prev_hit_start = cur_hit_end
                    # hsp model position
                    cur_model_seq = model_seq[local_start:local_end].replace(".", "-")
                    cur_model_gap_size = self._local_aln_gap_size(local_aln_idx[i], model_seq)
                    cur_model_start = prev_model_start + cur_model_gap_size
                    cur_model_end = cur_model_start + len(re.sub(_RE_LETTERS, "", cur_model_seq))
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
                break
            
            # parse the alignment blocks in the hsp
            # each block have 4 (hmmonly) or 5 (cm) lines followed by an empty line
            block_size = 6 if model == "cm" else 5
            lines = [None] * block_size
            for i in range(block_size):
                lines[i] = self.line
                self.line = read_forward(self.handle)
                                    
            # get the position of the alignment in the string using the pp line
            blklen = len(lines[5].strip().split()[0])
            blkstart = len(lines[5]) - blklen - 4
            blkend = len(lines[5]) - 4
            model_seq += lines[2][blkstart:blkend]
            hit_seq += lines[4][blkstart:blkend]
            annot["NC"] += lines[0][blkstart:blkend]
            annot["CS"] += lines[1][blkstart:blkend]
            annot["similarity"] += lines[3][blkstart:blkend]
            annot["PP"] += lines[5][blkstart:blkend]

        return frag_list


    def _local_aln_gap_size(self, cur_aln_idx, seq):
        """Calculate the gap size between the local alignments (PRIVATE)."""
        gap_len = 0
        if cur_aln_idx[1] > 0:
            gap_len = sum([int(n) for n in re.findall(_RE_NUMERIC, seq[cur_aln_idx[0]:cur_aln_idx[1]])])
            assert gap_len > 0
        return gap_len

            

# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
