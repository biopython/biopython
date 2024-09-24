# Copyright 2024 by Samuel Prince. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Bio.SearchIO parser for Infernal tabular output format."""

from Bio.SearchIO._index import SearchIndexer
from Bio.SearchIO._model import HSP
from Bio.SearchIO._model import HSPFragment
from Bio.SearchIO._model import QueryResult

from ._base import _BaseInfernalParser

__all__ = ("InfernalTabParser", "InfernalTabIndexer")

# tabular format column names
_TAB_FORMAT = {
    1: (
        "target_name",
        "target_acc",
        "query_name",
        "query_acc",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "evalue",
        "inc",
        "description",
    ),
    2: (
        "idx",
        "target_name",
        "target_acc",
        "query_name",
        "query_acc",
        "clan",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "evalue",
        "inc",
        "olp",
        "anyidx",
        "afrct1",
        "afrct2",
        "winidx",
        "wfrct1",
        "wfrct2",
        "mdl_len",
        "seq_len",
        "description",
    ),
    3: (
        "target_name",
        "target_acc",
        "query_name",
        "query_acc",
        "mdl",
        "mdl_from",
        "mdl_to",
        "seq_from",
        "seq_to",
        "strand",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "evalue",
        "inc",
        "mdl_len",
        "seq_len",
        "description",
    ),
}
_DEFAULT_TAB_FORMAT = 1

# column to class attribute map
_COLUMN_QRESULT = {
    "query_name": ("id", str),
    "query_acc": ("accession", str),
    "seq_len": ("seq_len", int),
    "clan": ("clan", str),
    "mdl": ("model", str),
}
_COLUMN_HIT = {
    "target_name": ("id", str),
    "target_acc": ("accession", str),
    "query_name": ("query_id", str),
    "description": ("description", str),
    "mdl_len": ("seq_len", int),
}
_COLUMN_HSP = {
    "score": ("bitscore", float),
    "evalue": ("evalue", float),
    "bias": ("bias", float),
    "gc": ("gc", float),
    "trunc": ("truncated", str),
    "pass": ("pipeline_pass", int),
    "inc": ("is_included", str),
    "olp": ("olp", str),
    "anyidx": ("anyidx", str),
    "afrct1": ("afrct1", str),
    "afrct2": ("afrct2", str),
    "winidx": ("winidx", str),
    "wfrct1": ("wfrct1", str),
    "wfrct2": ("wfrct2", str),
}
_COLUMN_FRAG = {
    "mdl_from": ("query_start", int),
    "mdl_to": ("query_end", int),
    "seq_from": ("hit_start", int),
    "seq_to": ("hit_end", int),
    "strand": ("hit_strand", str),
}


class InfernalTabParser(_BaseInfernalParser):
    """Parser for the Infernal tabular format."""

    def __init__(self, handle, fmt=_DEFAULT_TAB_FORMAT):
        """Initialize the class."""
        self.handle = handle
        self.line = self.handle.readline().strip()
        if not isinstance(fmt, int):
            raise TypeError
        if not 1 <= fmt <= 3:
            raise ValueError("Invalid tabular format number, must be 1, 2 or 3.")
        self.fmt = fmt

    def __iter__(self):
        """Iterate over InfernalTabParser, yields query results."""
        # read through the header and footer
        while self.line.startswith("#"):
            self.line = self.handle.readline()
        # if we have result rows, parse it
        if self.line:
            yield from self._parse_qresult()

    def _parse_row(self):
        """Return a dictionary of parsed row values (PRIVATE)."""
        cols = [x for x in self.line.strip().split(" ") if x]
        if len(cols) < len(_TAB_FORMAT[self.fmt]):
            raise ValueError(
                f"Less columns than expected for format {self.fmt}, only {len(cols)}"
            )
        # combine extra description columns into one string
        cols[len(_TAB_FORMAT[self.fmt]) - 1] = " ".join(
            cols[len(_TAB_FORMAT[self.fmt]) - 1 :]
        )

        qresult, hit, hsp, frag = {}, {}, {}, {}
        for sname, value in zip(
            _TAB_FORMAT[self.fmt], cols[: len(_TAB_FORMAT[self.fmt])]
        ):
            # iterate over each dict, mapping pair to determine
            # attribute name and value of each column
            for parsed_dict, mapping in (
                (qresult, _COLUMN_QRESULT),
                (hit, _COLUMN_HIT),
                (hsp, _COLUMN_HSP),
                (frag, _COLUMN_FRAG),
            ):
                # process parsed value according to mapping
                if sname in mapping:
                    attr_name, caster = mapping[sname]
                    if caster is not str:
                        value = caster(value)
                    parsed_dict[attr_name] = value

        # adjust start and end coordinates according to strand
        self._adjust_coords(frag)
        # convert inclusion string to a bool
        self._inclusion_str_to_bool(hsp)

        return {"qresult": qresult, "hit": hit, "hsp": hsp, "frag": frag}

    def _adjust_coords(self, frag):
        """Adjust start and end coordinates according to strand (PRIVATE)."""
        strand = frag["hit_strand"]
        assert strand is not None
        # switch start <--> end coordinates if strand is -1 and the strand to an integer (0 or -1)
        if strand == "-":
            hit_start = frag["hit_start"]
            hit_end = frag["hit_end"]
            frag["hit_start"] = hit_end
            frag["hit_end"] = hit_start
            frag["hit_strand"] = -1
        else:
            frag["hit_strand"] = 0

    def _inclusion_str_to_bool(self, hsp):
        """Convert inclusion string to a bool (PRIVATE)."""
        is_included = hsp["is_included"]
        hsp["is_included"] = True if is_included == "!" else False

    def _parse_qresult(self):
        """Yield QueryResult objects (PRIVATE)."""
        # state values, determines what to do for each line
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        state_HIT_NEW = 2
        state_HIT_SAME = 4
        # dummies for initial states
        qres_state = None
        hit_state = None
        file_state = None
        cur_qid = None
        cur_hid = None
        # dummies for initial id caches
        prev_qid = None
        prev_hid = None
        # dummies for initial parsed value containers
        cur, prev = None, None
        hit_dict = {}

        while True:
            # store previous line's parsed values for all lines after the first
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            # only parse the result row if it's not EOF or a comment line
            if self.line and not self.line.startswith("#"):
                cur = self._parse_row()
                cur_qid = cur["qresult"]["id"]
                cur_hid = cur["hit"]["id"]
            else:
                file_state = state_EOF
                # mock values for cur_qid and cur_hid since the line is empty
                cur_qid, cur_hid = None, None

            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME
            # new hits are hits with different id or hits in a new qresult
            if prev_hid != cur_hid or qres_state == state_QRES_NEW:
                hit_state = state_HIT_NEW
            else:
                hit_state = state_HIT_SAME

            # creating objects for the previously parsed line(s), so nothing is done
            # in the first parsed line (prev == None)
            if prev is not None:
                # create fragment and HSP and set their attributes
                # and append to hit container
                frag = HSPFragment(prev_hid, prev_qid)
                for attr, value in prev["frag"].items():
                    setattr(frag, attr, value)
                hsp = HSP([frag])
                for attr, value in prev["hsp"].items():
                    setattr(hsp, attr, value)
                self._add_hit_to_dict(prev["hit"], hsp, hit_dict)

                # create qresult and yield if we're at a new qresult or at EOF
                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    # create the hit list from the temporary container
                    qresult = QueryResult(self._hit_to_list(hit_dict), prev_qid)
                    for attr, value in prev["qresult"].items():
                        setattr(qresult, attr, value)
                    yield qresult
                    # if we're at EOF, break
                    if file_state == state_EOF:
                        break
                    hit_dict = {}

            self.line = self.handle.readline()


class InfernalTabIndexer(SearchIndexer):
    """Indexer class for Infernal tabular output."""

    _parser = InfernalTabParser

    def __init__(self, filename, fmt=_DEFAULT_TAB_FORMAT):
        """Initialize the class."""
        SearchIndexer.__init__(self, filename, fmt=fmt)
        self._query_id_idx = 3 if fmt == 2 else 2

    def __iter__(self):
        """Iterate over the file handle; yields key, start offset, and length."""
        handle = self._handle
        handle.seek(0)
        qresult_key = None
        start_offset = handle.tell()

        # set line with initial mock value
        line = b"#"

        # read through header
        while line.startswith(b"#"):
            start_offset = handle.tell()
            line = handle.readline()

        # index the qresults
        while True:
            # get end offset here since we only know a qresult ends after
            # encountering the next one
            end_offset = handle.tell()

            if not line:
                break

            cols = [x for x in line.strip().split(b" ") if x]

            if qresult_key is None:
                qresult_key = cols[self._query_id_idx]
            else:
                curr_key = cols[self._query_id_idx]
                if curr_key != qresult_key:
                    adj_end = end_offset - len(line)
                    yield (qresult_key.decode(), start_offset, adj_end - start_offset)
                    qresult_key = curr_key
                    start_offset = adj_end

            line = handle.readline()
            if not line or line.startswith(b"#"):
                yield (qresult_key.decode(), start_offset, end_offset - start_offset)
                break

    def get_raw(self, offset):
        """Return the raw bytes string of a QueryResult object from the given offset."""
        handle = self._handle
        qresult_key = None
        qresult_raw = b""

        # read the qresult raw string
        handle.seek(offset)
        while True:
            line = handle.readline()
            if not line or line.startswith(b"#"):
                break
            cols = [x for x in line.strip().split(b" ") if x]
            if qresult_key is None:
                qresult_key = cols[self._query_id_idx]
            else:
                curr_key = cols[self._query_id_idx]
                if curr_key != qresult_key:
                    break
            qresult_raw += line

        return qresult_raw


# if not used as a module, run the doctest
if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
