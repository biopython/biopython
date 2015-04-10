# Copyright 2009-2011 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Dictionary like indexing of sequence files (PRIVATE).

You are not expected to access this module, or any of its code, directly. This
is all handled internally by the Bio.SeqIO.index(...) function which is the
public interface for this functionality.

The basic idea is that we scan over a sequence file, looking for new record
markers. We then try and extract the string that Bio.SeqIO.parse/read would
use as the record id, ideally without actually parsing the full record. We
then use a subclassed Python dictionary to record the file offset for the
record start against the record id.

Note that this means full parsing is on demand, so any invalid or problem
record may not trigger an exception until it is accessed. This is by design.

This means our dictionary like objects have in memory ALL the keys (all the
record identifiers), which shouldn't be a problem even with second generation
sequencing. If this is an issue later on, storing the keys and offsets in a
temp lookup file might be one idea (e.g. using SQLite or an OBDA style index).
"""

from __future__ import print_function

import re
from Bio._py3k import StringIO
from Bio._py3k import _bytes_to_string, _as_bytes

from Bio import SeqIO
from Bio import Alphabet
from Bio.File import _IndexedSeqFileProxy, _open_for_random_access

__docformat__ = "restructuredtext en"


class SeqFileRandomAccess(_IndexedSeqFileProxy):
    def __init__(self, filename, format, alphabet):
        self._handle = _open_for_random_access(filename)
        self._alphabet = alphabet
        self._format = format
        # Load the parser class/function once an avoid the dict lookup in each
        # __getitem__ call:
        i = SeqIO._FormatToIterator[format]
        # The following alphabet code is a bit nasty... duplicates logic in
        # Bio.SeqIO.parse()
        if alphabet is None:
            def _parse(handle):
                """Dynamically generated parser function (PRIVATE)."""
                return next(i(handle))
        else:
            # TODO - Detect alphabet support ONCE at __init__
            def _parse(handle):
                """Dynamically generated parser function (PRIVATE)."""
                try:
                    return next(i(handle, alphabet=alphabet))
                except TypeError:
                    return next(SeqIO._force_alphabet(i(handle), alphabet))
        self._parse = _parse

    def get(self, offset):
        """Returns SeqRecord."""
        # Should be overridden for binary file formats etc:
        return self._parse(StringIO(_bytes_to_string(self.get_raw(offset))))


####################
# Special indexers #
####################
# Anything where the records cannot be read simply by parsing from
# the record start. For example, anything requiring information from
# a file header - e.g. SFF files where we would need to know the
# number of flows.
class SffRandomAccess(SeqFileRandomAccess):
    """Random access to a Standard Flowgram Format (SFF) file."""
    def __init__(self, filename, format, alphabet):
        SeqFileRandomAccess.__init__(self, filename, format, alphabet)
        header_length, index_offset, index_length, number_of_reads, \
            self._flows_per_read, self._flow_chars, self._key_sequence \
            = SeqIO.SffIO._sff_file_header(self._handle)

    def __iter__(self):
        """Load any index block in the file, or build it the slow way (PRIVATE)."""
        if self._alphabet is None:
            self._alphabet = Alphabet.generic_dna
        handle = self._handle
        handle.seek(0)
        # Alread did this in __init__ but need handle in right place
        header_length, index_offset, index_length, number_of_reads, \
            self._flows_per_read, self._flow_chars, self._key_sequence \
            = SeqIO.SffIO._sff_file_header(handle)
        if index_offset and index_length:
            # There is an index provided, try this the fast way:
            count = 0
            max_offset = 0
            try:
                for name, offset in SeqIO.SffIO._sff_read_roche_index(handle):
                    max_offset = max(max_offset, offset)
                    yield name, offset, 0
                    count += 1
                assert count == number_of_reads, \
                    "Indexed %i records, expected %i" \
                    % (count, number_of_reads)
                # If that worked, call _check_eof ...
            except ValueError as err:
                import warnings
                from Bio import BiopythonParserWarning
                warnings.warn("Could not parse the SFF index: %s" % err,
                              BiopythonParserWarning)
                assert count == 0, "Partially populated index"
                handle.seek(0)
                # Drop out to the slow way...
            else:
                # Fast way worked, check EOF
                if index_offset + index_length <= max_offset:
                    # Can have an index at start (or mid-file)
                    handle.seek(max_offset)
                    # Parse the final read,
                    SeqIO.SffIO._sff_read_raw_record(handle, self._flows_per_read)
                    # Should now be at the end of the file!
                SeqIO.SffIO._check_eof(handle, index_offset, index_length)
                return
        # We used to give a warning in this case, but Ion Torrent's
        # SFF files don't have an index so that would be annoying.
        # Fall back on the slow way!
        count = 0
        for name, offset in SeqIO.SffIO._sff_do_slow_index(handle):
            yield name, offset, 0
            count += 1
        assert count == number_of_reads, \
            "Indexed %i records, expected %i" % (count, number_of_reads)
        SeqIO.SffIO._check_eof(handle, index_offset, index_length)

    def get(self, offset):
        handle = self._handle
        handle.seek(offset)
        return SeqIO.SffIO._sff_read_seq_record(handle,
                                                self._flows_per_read,
                                                self._flow_chars,
                                                self._key_sequence,
                                                self._alphabet)

    def get_raw(self, offset):
        handle = self._handle
        handle.seek(offset)
        return SeqIO.SffIO._sff_read_raw_record(handle, self._flows_per_read)


class SffTrimedRandomAccess(SffRandomAccess):
    def get(self, offset):
        handle = self._handle
        handle.seek(offset)
        return SeqIO.SffIO._sff_read_seq_record(handle,
                                                self._flows_per_read,
                                                self._flow_chars,
                                                self._key_sequence,
                                                self._alphabet,
                                                trim=True)


###################
# Simple indexers #
###################

class SequentialSeqFileRandomAccess(SeqFileRandomAccess):
    def __init__(self, filename, format, alphabet):
        SeqFileRandomAccess.__init__(self, filename, format, alphabet)
        marker = {"ace": "CO ",
                  "embl": "ID ",
                  "fasta": ">",
                  "genbank": "LOCUS ",
                  "gb": "LOCUS ",
                  "imgt": "ID ",
                  "phd": "BEGIN_SEQUENCE",
                  "pir": ">..;",
                  "qual": ">",
                  "qual": ">",
                  "swiss": "ID ",
                  "uniprot-xml": "<entry ",
                  }[format]
        self._marker = marker
        self._marker_re = re.compile(_as_bytes("^%s" % marker))

    def __iter__(self):
        """Returns (id,offset) tuples."""
        marker_offset = len(self._marker)
        marker_re = self._marker_re
        handle = self._handle
        handle.seek(0)
        # Skip and header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        # Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            # Here we can assume the record.id is the first word after the
            # marker. This is generally fine... but not for GenBank, EMBL, Swiss
            id = line[marker_offset:].strip().split(None, 1)[0]
            length = len(line)
            while True:
                end_offset = handle.tell()
                line = handle.readline()
                if marker_re.match(line) or not line:
                    yield _bytes_to_string(id), start_offset, length
                    start_offset = end_offset
                    break
                else:
                    # Track this explicitly as can't do file offset difference on BGZF
                    length += len(line)
        assert not line, repr(line)

    def get_raw(self, offset):
        """Similar to the get method, but returns the record as a raw string."""
        # For non-trivial file formats this must be over-ridden in the subclass
        handle = self._handle
        marker_re = self._marker_re
        handle.seek(offset)
        lines = [handle.readline()]
        while True:
            line = handle.readline()
            if marker_re.match(line) or not line:
                # End of file, or start of next record => end of this record
                break
            lines.append(line)
        return _as_bytes("").join(lines)


#######################################
# Fiddly indexers: GenBank, EMBL, ... #
#######################################

class GenBankRandomAccess(SequentialSeqFileRandomAccess):
    """Indexed dictionary like access to a GenBank file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        dot_char = _as_bytes(".")
        accession_marker = _as_bytes("ACCESSION ")
        version_marker = _as_bytes("VERSION ")
        # Skip and header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        # Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            # We cannot assume the record.id is the first word after LOCUS,
            # normally the first entry on the VERSION or ACCESSION line is used.
            key = None
            length = len(line)
            while True:
                end_offset = handle.tell()
                line = handle.readline()
                if marker_re.match(line) or not line:
                    if not key:
                        raise ValueError(
                            "Did not find ACCESSION/VERSION lines")
                    yield _bytes_to_string(key), start_offset, length
                    start_offset = end_offset
                    break
                elif line.startswith(accession_marker):
                    key = line.rstrip().split()[1]
                elif line.startswith(version_marker):
                    version_id = line.rstrip().split()[1]
                    if version_id.count(dot_char) == 1 and version_id.split(dot_char)[1].isdigit():
                        # This should mimic the GenBank parser...
                        key = version_id
                length += len(line)
        assert not line, repr(line)


class EmblRandomAccess(SequentialSeqFileRandomAccess):
    """Indexed dictionary like access to an EMBL file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        semi_char = _as_bytes(";")
        dot_char = _as_bytes(".")
        sv_marker = _as_bytes("SV ")
        ac_marker = _as_bytes("AC ")
        # Skip any header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        # Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            # We cannot assume the record.id is the first word after ID,
            # normally the SV line is used.
            setbysv = False  # resets sv as false
            length = len(line)
            if line[2:].count(semi_char) == 6:
                # Looks like the semi colon separated style introduced in 2006
                parts = line[3:].rstrip().split(semi_char)
                if parts[1].strip().startswith(sv_marker):
                    # The SV bit gives the version
                    key = parts[0].strip() + dot_char + \
                        parts[1].strip().split()[1]
                    setbysv = True
                else:
                    key = parts[0].strip()
            elif line[2:].count(semi_char) == 3:
                # Looks like the pre 2006 style, take first word only
                key = line[3:].strip().split(None, 1)[0]
                if key.endswith(semi_char):
                    key = key[:-1]
            else:
                raise ValueError(
                    'Did not recognise the ID line layout:\n' + line)
            while True:
                end_offset = handle.tell()
                line = handle.readline()
                if marker_re.match(line) or not line:
                    end_offset = handle.tell() - len(line)
                    yield _bytes_to_string(key), start_offset, length
                    start_offset = end_offset
                    break
                elif line.startswith(ac_marker) and not setbysv:
                    key = line.rstrip().split()[1]
                    if key.endswith(semi_char):
                        key = key[:-1]
                elif line.startswith(sv_marker):
                    key = line.rstrip().split()[1]
                    setbysv = True
                length += len(line)
        assert not line, repr(line)


class SwissRandomAccess(SequentialSeqFileRandomAccess):
    """Random access to a SwissProt file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        semi_char = _as_bytes(";")
        # Skip any header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        # Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            length = len(line)
            # We cannot assume the record.id is the first word after ID,
            # normally the following AC line is used.
            line = handle.readline()
            length += len(line)
            assert line.startswith(_as_bytes("AC "))
            key = line[3:].strip().split(semi_char)[0].strip()
            while True:
                end_offset = handle.tell()
                line = handle.readline()
                if marker_re.match(line) or not line:
                    yield _bytes_to_string(key), start_offset, length
                    start_offset = end_offset
                    break
                length += len(line)
        assert not line, repr(line)


class UniprotRandomAccess(SequentialSeqFileRandomAccess):
    """Random access to a UniProt XML file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        start_acc_marker = _as_bytes("<accession>")
        end_acc_marker = _as_bytes("</accession>")
        end_entry_marker = _as_bytes("</entry>")
        less_than = _as_bytes("<")
        # Skip any header before first record
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if marker_re.match(line) or not line:
                break
        # Should now be at the start of a record, or end of the file
        while marker_re.match(line):
            length = len(line)
            # We expect the next line to be <accession>xxx</accession>
            # (possibly with leading spaces)
            # but allow it to be later on within the <entry>
            key = None
            while True:
                line = handle.readline()
                if key is None and start_acc_marker in line:
                    assert end_acc_marker in line, line
                    key = line[line.find(
                        start_acc_marker) + 11:].split(less_than, 1)[0]
                    length += len(line)
                elif end_entry_marker in line:
                    end_offset = handle.tell() - len(line) \
                        + line.find(end_entry_marker) + 8
                    break
                elif marker_re.match(line) or not line:
                    # Start of next record or end of file
                    raise ValueError("Didn't find end of record")
                else:
                    length += len(line)
            if not key:
                raise ValueError("Did not find <accession> line in bytes %i to %i"
                                 % (start_offset, end_offset))
            yield _bytes_to_string(key), start_offset, length
            # Find start of next record
            while not marker_re.match(line) and line:
                start_offset = handle.tell()
                line = handle.readline()
        assert not line, repr(line)

    def get_raw(self, offset):
        """Similar to the get method, but returns the record as a raw string."""
        handle = self._handle
        marker_re = self._marker_re
        end_entry_marker = _as_bytes("</entry>")
        handle.seek(offset)
        data = [handle.readline()]
        while True:
            line = handle.readline()
            i = line.find(end_entry_marker)
            if i != -1:
                data.append(line[:i + 8])
                break
            if marker_re.match(line) or not line:
                # End of file, or start of next record
                raise ValueError("Didn't find end of record")
            data.append(line)
        return _as_bytes("").join(data)

    def get(self, offset):
        # TODO - Can we handle this directly in the parser?
        # This is a hack - use get_raw for <entry>...</entry> and wrap it with
        # the apparently required XML header and footer.
        data = """<?xml version='1.0' encoding='UTF-8'?>
        <uniprot xmlns="http://uniprot.org/uniprot"
        xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
        xsi:schemaLocation="http://uniprot.org/uniprot
        http://www.uniprot.org/support/docs/uniprot.xsd">
        %s
        </uniprot>
        """ % _bytes_to_string(self.get_raw(offset))
        # TODO - For consistency, this function should not accept a string:
        return next(SeqIO.UniprotIO.UniprotIterator(data))


class IntelliGeneticsRandomAccess(SeqFileRandomAccess):
    """Random access to a IntelliGenetics file."""
    def __init__(self, filename, format, alphabet):
        SeqFileRandomAccess.__init__(self, filename, format, alphabet)
        self._marker_re = re.compile(_as_bytes("^;"))

    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        marker_re = self._marker_re
        semi_char = _as_bytes(";")
        while True:
            offset = handle.tell()
            line = handle.readline()
            length = len(line)
            if marker_re.match(line):
                # Now look for the first line which doesn't start ";"
                while True:
                    line = handle.readline()
                    if line[0:1] != semi_char and line.strip():
                        key = line.split()[0]
                        yield _bytes_to_string(key), offset, length
                        break
                    if not line:
                        raise ValueError("Premature end of file?")
                    length += len(line)
            elif not line:
                # End of file
                break

    def get_raw(self, offset):
        handle = self._handle
        handle.seek(offset)
        marker_re = self._marker_re
        lines = []
        line = handle.readline()
        semi_char = _as_bytes(";")
        while line.startswith(semi_char):
            lines.append(line)
            line = handle.readline()
        while line and not line.startswith(semi_char):
            lines.append(line)
            line = handle.readline()
        return _as_bytes("").join(lines)


class TabRandomAccess(SeqFileRandomAccess):
    """Random access to a simple tabbed file."""
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        tab_char = _as_bytes("\t")
        while True:
            start_offset = handle.tell()
            line = handle.readline()
            if not line:
                break  # End of file
            try:
                key = line.split(tab_char)[0]
            except ValueError as err:
                if not line.strip():
                    # Ignore blank lines
                    continue
                else:
                    raise err
            else:
                yield _bytes_to_string(key), start_offset, len(line)

    def get_raw(self, offset):
        """Like the get method, but returns the record as a raw string."""
        handle = self._handle
        handle.seek(offset)
        return handle.readline()


##########################
# Now the FASTQ indexers #
##########################

class FastqRandomAccess(SeqFileRandomAccess):
    """Random access to a FASTQ file (any supported variant).

    With FASTQ the records all start with a "@" line, but so can quality lines.
    Note this will cope with line-wrapped FASTQ files.
    """
    def __iter__(self):
        handle = self._handle
        handle.seek(0)
        id = None
        start_offset = handle.tell()
        line = handle.readline()
        if not line:
            # Empty file!
            return
        at_char = _as_bytes("@")
        plus_char = _as_bytes("+")
        if line[0:1] != at_char:
            raise ValueError("Problem with FASTQ @ line:\n%r" % line)
        while line:
            # assert line[0]=="@"
            # This record seems OK (so far)
            id = line[1:].rstrip().split(None, 1)[0]
            # Find the seq line(s)
            seq_len = 0
            length = len(line)
            while line:
                line = handle.readline()
                length += len(line)
                if line.startswith(plus_char):
                    break
                seq_len += len(line.strip())
            if not line:
                raise ValueError("Premature end of file in seq section")
            # assert line[0]=="+"
            # Find the qual line(s)
            qual_len = 0
            while line:
                if seq_len == qual_len:
                    if seq_len == 0:
                        # Special case, quality line should be just "\n"
                        line = handle.readline()
                        if line.strip():
                            raise ValueError("Expected blank quality line, not %r" % line)
                    # Should be end of record...
                    end_offset = handle.tell()
                    line = handle.readline()
                    if line and line[0:1] != at_char:
                        raise ValueError("Problem with line %r" % line)
                    break
                else:
                    line = handle.readline()
                    qual_len += len(line.strip())
                    length += len(line)
            if seq_len != qual_len:
                raise ValueError("Problem with quality section")
            yield _bytes_to_string(id), start_offset, length
            start_offset = end_offset
        # print("EOF")

    def get_raw(self, offset):
        """Similar to the get method, but returns the record as a raw string."""
        # TODO - Refactor this and the __init__ method to reduce code duplication?
        handle = self._handle
        handle.seek(offset)
        line = handle.readline()
        data = line
        at_char = _as_bytes("@")
        plus_char = _as_bytes("+")
        if line[0:1] != at_char:
            raise ValueError("Problem with FASTQ @ line:\n%r" % line)
        # Find the seq line(s)
        seq_len = 0
        while line:
            line = handle.readline()
            data += line
            if line.startswith(plus_char):
                break
            seq_len += len(line.strip())
        if not line:
            raise ValueError("Premature end of file in seq section")
        assert line[0:1] == plus_char
        # Find the qual line(s)
        qual_len = 0
        while line:
            if seq_len == qual_len:
                if seq_len == 0:
                    # Special case, quality line should be just "\n"
                    line = handle.readline()
                    if line.strip():
                        raise ValueError("Expected blank quality line, not %r" % line)
                    data += line
                # Should be end of record...
                line = handle.readline()
                if line and line[0:1] != at_char:
                    raise ValueError("Problem with line %r" % line)
                break
            else:
                line = handle.readline()
                data += line
                qual_len += len(line.strip())
        if seq_len != qual_len:
            raise ValueError("Problem with quality section")
        return data


###############################################################################

_FormatToRandomAccess = {"ace": SequentialSeqFileRandomAccess,
                         "embl": EmblRandomAccess,
                         "fasta": SequentialSeqFileRandomAccess,
                         "fastq": FastqRandomAccess,  # Class handles all three variants
                         "fastq-sanger": FastqRandomAccess,  # alias of the above
                         "fastq-solexa": FastqRandomAccess,
                         "fastq-illumina": FastqRandomAccess,
                         "genbank": GenBankRandomAccess,
                         "gb": GenBankRandomAccess,  # alias of the above
                         "ig": IntelliGeneticsRandomAccess,
                         "imgt": EmblRandomAccess,
                         "phd": SequentialSeqFileRandomAccess,
                         "pir": SequentialSeqFileRandomAccess,
                         "sff": SffRandomAccess,
                         "sff-trim": SffTrimedRandomAccess,
                         "swiss": SwissRandomAccess,
                         "tab": TabRandomAccess,
                         "qual": SequentialSeqFileRandomAccess,
                         "uniprot-xml": UniprotRandomAccess,
                         }
