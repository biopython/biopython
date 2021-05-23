# Copyright 2008-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for "fasta-m10" output from Bill Pearson's FASTA tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

This module contains a parser for the pairwise alignments produced by Bill
Pearson's FASTA tools, for use from the Bio.AlignIO interface where it is
referred to as the "fasta-m10" file format (as we only support the machine
readable output format selected with the -m 10 command line option).

This module does NOT cover the generic "fasta" file format originally
developed as an input format to the FASTA tools.  The Bio.AlignIO and
Bio.SeqIO both use the Bio.SeqIO.FastaIO module to deal with these files,
which can also be used to store a multiple sequence alignments.
"""
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _extract_alignment_region(alignment_seq_with_flanking, annotation):
    """Extract alignment region (PRIVATE).

    Helper function for the main parsing code.

    To get the actual pairwise alignment sequences, we must first
    translate the un-gapped sequence based coordinates into positions
    in the gapped sequence (which may have a flanking region shown
    using leading - characters).  To date, I have never seen any
    trailing flanking region shown in the m10 file, but the
    following code should also cope with that.

    Note that this code seems to work fine even when the "sq_offset"
    entries are present as a result of using the -X command line option.
    """
    align_stripped = alignment_seq_with_flanking.strip("-")
    display_start = int(annotation["al_display_start"])
    if int(annotation["al_start"]) <= int(annotation["al_stop"]):
        start = int(annotation["al_start"]) - display_start
        end = int(annotation["al_stop"]) - display_start + 1
    else:
        # FASTA has flipped this sequence...
        start = display_start - int(annotation["al_start"])
        end = display_start - int(annotation["al_stop"]) + 1

    end += align_stripped.count("-")
    if start < 0 or start >= end or end > len(align_stripped):
        raise ValueError(
            "Problem with sequence start/stop,\n%s[%i:%i]\n%s"
            % (alignment_seq_with_flanking, start, end, annotation)
        )
    return align_stripped[start:end]


def FastaM10Iterator(handle, seq_count=None):
    """Alignment iterator for the FASTA tool's pairwise alignment output.

    This is for reading the pairwise alignments output by Bill Pearson's
    FASTA program when called with the -m 10 command line option for machine
    readable output.  For more details about the FASTA tools, see the website
    http://fasta.bioch.virginia.edu/ and the paper:

         W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448

    This class is intended to be used via the Bio.AlignIO.parse() function
    by specifying the format as "fasta-m10" as shown in the following code::

        from Bio import AlignIO
        handle = ...
        for a in AlignIO.parse(handle, "fasta-m10"):
            assert len(a) == 2, "Should be pairwise!"
            print("Alignment length %i" % a.get_alignment_length())
            for record in a:
                print("%s %s %s" % (record.seq, record.name, record.id))

    Note that this is not a full blown parser for all the information
    in the FASTA output - for example, most of the header and all of the
    footer is ignored.  Also, the alignments are not batched according to
    the input queries.

    Also note that there can be up to about 30 letters of flanking region
    included in the raw FASTA output as contextual information.  This is NOT
    part of the alignment itself, and is not included in the resulting
    MultipleSeqAlignment objects returned.
    """
    state_PREAMBLE = -1
    state_NONE = 0
    state_QUERY_HEADER = 1
    state_ALIGN_HEADER = 2
    state_ALIGN_QUERY = 3
    state_ALIGN_MATCH = 4
    state_ALIGN_CONS = 5

    def build_hsp():
        if not query_tags and not match_tags:
            raise ValueError("No data for query %r, match %r" % (query_id, match_id))
        assert query_tags, query_tags
        assert match_tags, match_tags
        evalue = align_tags.get("fa_expect")
        tool = global_tags.get("tool", "").upper()

        q = _extract_alignment_region(query_seq, query_tags)
        if tool in ["TFASTX"] and len(match_seq) == len(q):
            m = match_seq
            # Quick hack until I can work out how -, * and / characters
            # and the apparent mix of aa and bp coordinates works.
        else:
            m = _extract_alignment_region(match_seq, match_tags)
        if len(q) != len(m):
            raise ValueError(
                f"""\
Darn... amino acids vs nucleotide coordinates?
tool: {tool}
query_seq: {query_seq}
query_tags: {query_tags}
{q} length: {len(q)}
match_seq: {match_seq}
match_tags: {match_tags}
{m} length: {len(m)}
handle.name: {handle.name}
"""
            )

        annotations = {}
        records = []

        # Want to record both the query header tags, and the alignment tags.
        annotations.update(header_tags)
        annotations.update(align_tags)

        # Query
        # =====
        record = SeqRecord(
            Seq(q),
            id=query_id,
            name="query",
            description=query_descr,
            annotations={"original_length": int(query_tags["sq_len"])},
        )
        # TODO - handle start/end coordinates properly. Short term hack for now:
        record._al_start = int(query_tags["al_start"])
        record._al_stop = int(query_tags["al_stop"])

        # TODO - Can FASTA output RNA?
        if "sq_type" in query_tags:
            if query_tags["sq_type"] == "D":
                record.annotations["molecule_type"] = "DNA"
            elif query_tags["sq_type"] == "p":
                record.annotations["molecule_type"] = "protein"

        records.append(record)

        # Match
        # =====
        record = SeqRecord(
            Seq(m),
            id=match_id,
            name="match",
            description=match_descr,
            annotations={"original_length": int(match_tags["sq_len"])},
        )
        # TODO - handle start/end coordinates properly. Short term hack for now:
        record._al_start = int(match_tags["al_start"])
        record._al_stop = int(match_tags["al_stop"])

        if "sq_type" in match_tags:
            if match_tags["sq_type"] == "D":
                record.annotations["molecule_type"] = "DNA"
            elif match_tags["sq_type"] == "p":
                record.annotations["molecule_type"] = "protein"

        records.append(record)

        return MultipleSeqAlignment(records, annotations=annotations)

    state = state_PREAMBLE
    query_id = None
    match_id = None
    query_descr = ""
    match_descr = ""
    global_tags = {}
    header_tags = {}
    align_tags = {}
    query_tags = {}
    match_tags = {}
    query_seq = ""
    match_seq = ""
    cons_seq = ""
    for line in handle:
        if ">>>" in line and not line.startswith(">>>"):
            if query_id and match_id:
                # This happens on old FASTA output which lacked an end of
                # query >>><<< marker line.
                yield build_hsp()
            state = state_NONE
            query_descr = line[line.find(">>>") + 3 :].strip()
            query_id = query_descr.split(None, 1)[0]
            match_id = None
            header_tags = {}
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
        elif line.startswith("!! No "):
            # e.g.
            # !! No library sequences with E() < 0.5
            # or on more recent versions,
            # No sequences with E() < 0.05
            assert state == state_NONE
            assert not header_tags
            assert not align_tags
            assert not match_tags
            assert not query_tags
            assert match_id is None
            assert not query_seq
            assert not match_seq
            assert not cons_seq
            query_id = None
        elif line.strip() in [">>><<<", ">>>///"]:
            # End of query, possible end of all queries
            if query_id and match_id:
                yield build_hsp()
            state = state_NONE
            query_id = None
            match_id = None
            header_tags = {}
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
        elif line.startswith(">>>"):
            # Should be start of a match!
            assert query_id is not None
            assert line[3:].split(", ", 1)[0] == query_id, line
            assert match_id is None
            assert not header_tags
            assert not align_tags
            assert not query_tags
            assert not match_tags
            assert not match_seq
            assert not query_seq
            assert not cons_seq
            state = state_QUERY_HEADER
        elif line.startswith(">>"):
            # Should now be at start of a match alignment!
            if query_id and match_id:
                yield build_hsp()
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
            match_descr = line[2:].strip()
            match_id = match_descr.split(None, 1)[0]
            state = state_ALIGN_HEADER
        elif line.startswith(">--"):
            # End of one HSP
            assert query_id and match_id, line
            yield build_hsp()
            # Clean up read for next HSP
            # but reuse header_tags
            align_tags = {}
            query_tags = {}
            match_tags = {}
            query_seq = ""
            match_seq = ""
            cons_seq = ""
            state = state_ALIGN_HEADER
        elif line.startswith(">"):
            if state == state_ALIGN_HEADER:
                # Should be start of query alignment seq...
                assert query_id is not None, line
                assert match_id is not None, line
                assert query_id.startswith(line[1:].split(None, 1)[0]), line
                state = state_ALIGN_QUERY
            elif state == state_ALIGN_QUERY:
                # Should be start of match alignment seq
                assert query_id is not None, line
                assert match_id is not None, line
                assert match_id.startswith(line[1:].split(None, 1)[0]), line
                state = state_ALIGN_MATCH
            elif state == state_NONE:
                # Can get > as the last line of a histogram
                pass
            else:
                raise RuntimeError("state %i got %r" % (state, line))
        elif line.startswith("; al_cons"):
            assert state == state_ALIGN_MATCH, line
            state = state_ALIGN_CONS
            # Next line(s) should be consensus seq...
        elif line.startswith("; "):
            if ": " in line:
                key, value = [s.strip() for s in line[2:].split(": ", 1)]
            else:
                import warnings
                from Bio import BiopythonParserWarning

                # Seen in lalign36, specifically version 36.3.4 Apr, 2011
                # Fixed in version 36.3.5b Oct, 2011(preload8)
                warnings.warn(
                    "Missing colon in line: %r" % line, BiopythonParserWarning
                )
                try:
                    key, value = [s.strip() for s in line[2:].split(" ", 1)]
                except ValueError:
                    raise ValueError("Bad line: %r" % line) from None
            if state == state_QUERY_HEADER:
                header_tags[key] = value
            elif state == state_ALIGN_HEADER:
                align_tags[key] = value
            elif state == state_ALIGN_QUERY:
                query_tags[key] = value
            elif state == state_ALIGN_MATCH:
                match_tags[key] = value
            else:
                raise RuntimeError("Unexpected state %r, %r" % (state, line))
        elif state == state_ALIGN_QUERY:
            query_seq += line.strip()
        elif state == state_ALIGN_MATCH:
            match_seq += line.strip()
        elif state == state_ALIGN_CONS:
            cons_seq += line.strip("\n")
        elif state == state_PREAMBLE:
            if line.startswith("#"):
                global_tags["command"] = line[1:].strip()
            elif line.startswith(" version "):
                global_tags["version"] = line[9:].strip()
            elif " compares a " in line:
                global_tags["tool"] = line[: line.find(" compares a ")].strip()
            elif " searches a " in line:
                global_tags["tool"] = line[: line.find(" searches a ")].strip()
        else:
            pass
