# Copyright 2007-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "genbank" and "embl" file formats.

You are expected to use this module via the Bio.SeqIO functions.
Note that internally this module calls Bio.GenBank to do the actual
parsing of GenBank, EMBL and IMGT files.

See Also:
International Nucleotide Sequence Database Collaboration
http://www.insdc.org/

GenBank
http://www.ncbi.nlm.nih.gov/Genbank/

EMBL Nucleotide Sequence Database
http://www.ebi.ac.uk/embl/

DDBJ (DNA Data Bank of Japan)
http://www.ddbj.nig.ac.jp/

IMGT (use a variant of EMBL format with longer feature indents)
http://imgt.cines.fr/download/LIGM-DB/userman_doc.html
http://imgt.cines.fr/download/LIGM-DB/ftable_doc.html
http://www.ebi.ac.uk/imgt/hla/docs/manual.html

"""


import warnings
from datetime import datetime
from Bio import BiopythonWarning

from Bio.Seq import UnknownSeq
from Bio.GenBank.Scanner import GenBankScanner, EmblScanner, _ImgtScanner
from Bio import SeqIO
from Bio import SeqFeature
from .Interfaces import SequenceIterator, SequenceWriter


# NOTE
# ====
# The "brains" for parsing GenBank, EMBL and IMGT files (and any
# other flat file variants from the INSDC in future) is in
# Bio.GenBank.Scanner (plus the _FeatureConsumer in Bio.GenBank)
# However, all the writing code is in this file.


class GenBankIterator(SequenceIterator):
    """Parser for GenBank files."""

    def __init__(self, source):
        """Break up a Genbank file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.

        This gets called internally by Bio.SeqIO for the GenBank file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("GenBank/cor6_6.gb", "gb"):
        ...     print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        Equivalently,

        >>> with open("GenBank/cor6_6.gb") as handle:
        ...     for record in GenBankIterator(handle):
        ...         print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        """
        super().__init__(source, mode="t", fmt="GenBank")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = GenBankScanner(debug=0).parse_records(handle)
        return records


class EmblIterator(SequenceIterator):
    """Parser for EMBL files."""

    def __init__(self, source):
        """Break up an EMBL file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.

        This gets called internally by Bio.SeqIO for the EMBL file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("EMBL/epo_prt_selection.embl", "embl"):
        ...     print(record.id)
        ...
        A00022.1
        A00028.1
        A00031.1
        A00034.1
        A00060.1
        A00071.1
        A00072.1
        A00078.1
        CQ797900.1

        Equivalently,

        >>> with open("EMBL/epo_prt_selection.embl") as handle:
        ...     for record in EmblIterator(handle):
        ...         print(record.id)
        ...
        A00022.1
        A00028.1
        A00031.1
        A00034.1
        A00060.1
        A00071.1
        A00072.1
        A00078.1
        CQ797900.1

        """
        super().__init__(source, mode="t", fmt="EMBL")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = EmblScanner(debug=0).parse_records(handle)
        return records


class ImgtIterator(SequenceIterator):
    """Parser for IMGT files."""

    def __init__(self, source):
        """Break up an IMGT file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.
        """
        super().__init__(source, mode="t", fmt="IMGT")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = _ImgtScanner(debug=0).parse_records(handle)
        return records


class GenBankCdsFeatureIterator(SequenceIterator):
    """Parser for GenBank files, creating a SeqRecord for each CDS feature."""

    def __init__(self, source):
        """Break up a Genbank file into SeqRecord objects for each CDS feature.

        Argument source is a file-like object opened in text mode or a path to a file.

        Every section from the LOCUS line to the terminating // can contain
        many CDS features.  These are returned as with the stated amino acid
        translation sequence (if given).
        """
        super().__init__(source, mode="t", fmt="GenBank")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        return GenBankScanner(debug=0).parse_cds_features(handle)


class EmblCdsFeatureIterator(SequenceIterator):
    """Parser for EMBL files, creating a SeqRecord for each CDS feature."""

    def __init__(self, source):
        """Break up a EMBL file into SeqRecord objects for each CDS feature.

        Argument source is a file-like object opened in text mode or a path to a file.

        Every section from the LOCUS line to the terminating // can contain
        many CDS features.  These are returned as with the stated amino acid
        translation sequence (if given).
        """
        super().__init__(source, mode="t", fmt="EMBL")

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        return EmblScanner(debug=0).parse_cds_features(handle)


def _insdc_feature_position_string(pos, offset=0):
    """Build a GenBank/EMBL position string (PRIVATE).

    Use offset=1 to add one to convert a start position from python counting.
    """
    if isinstance(pos, SeqFeature.ExactPosition):
        return "%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.WithinPosition):
        return "(%i.%i)" % (
            pos.position + offset,
            pos.position + pos.extension + offset,
        )
    elif isinstance(pos, SeqFeature.BetweenPosition):
        return "(%i^%i)" % (
            pos.position + offset,
            pos.position + pos.extension + offset,
        )
    elif isinstance(pos, SeqFeature.BeforePosition):
        return "<%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.AfterPosition):
        return ">%i" % (pos.position + offset)
    elif isinstance(pos, SeqFeature.OneOfPosition):
        return "one-of(%s)" % ",".join(
            _insdc_feature_position_string(p, offset) for p in pos.position_choices
        )
    elif isinstance(pos, SeqFeature.AbstractPosition):
        raise NotImplementedError("Please report this as a bug in Biopython.")
    else:
        raise ValueError("Expected a SeqFeature position object.")


def _insdc_location_string_ignoring_strand_and_subfeatures(location, rec_length):
    if location.ref:
        ref = "%s:" % location.ref
    else:
        ref = ""
    assert not location.ref_db
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start.position == location.end.position
    ):
        # Special case, for 12:12 return 12^13
        # (a zero length slice, meaning the point between two letters)
        if location.end.position == rec_length:
            # Very special case, for a between position at the end of a
            # sequence (used on some circular genomes, Bug 3098) we have
            # N:N so return N^1
            return "%s%i^1" % (ref, rec_length)
        else:
            return "%s%i^%i" % (ref, location.end.position, location.end.position + 1)
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start.position + 1 == location.end.position
    ):
        # Special case, for 11:12 return 12 rather than 12..12
        # (a length one slice, meaning a single letter)
        return "%s%i" % (ref, location.end.position)
    elif isinstance(location.start, SeqFeature.UnknownPosition) or isinstance(
        location.end, SeqFeature.UnknownPosition
    ):
        # Special case for features from SwissProt/UniProt files
        if isinstance(location.start, SeqFeature.UnknownPosition) and isinstance(
            location.end, SeqFeature.UnknownPosition
        ):
            # warnings.warn("Feature with unknown location", BiopythonWarning)
            # return "?"
            raise ValueError("Feature with unknown location")
        elif isinstance(location.start, SeqFeature.UnknownPosition):
            # Treat the unknown start position as a BeforePosition
            return "%s<%i..%s" % (
                ref,
                location.nofuzzy_end,
                _insdc_feature_position_string(location.end),
            )
        else:
            # Treat the unknown end position as an AfterPosition
            return "%s%s..>%i" % (
                ref,
                _insdc_feature_position_string(location.start, +1),
                location.nofuzzy_start + 1,
            )
    else:
        # Typical case, e.g. 12..15 gets mapped to 11:15
        return (
            ref
            + _insdc_feature_position_string(location.start, +1)
            + ".."
            + _insdc_feature_position_string(location.end)
        )


def _insdc_location_string(location, rec_length):
    """Build a GenBank/EMBL location from a (Compound) FeatureLocation (PRIVATE).

    There is a choice of how to show joins on the reverse complement strand,
    GenBank used "complement(join(1,10),(20,100))" while EMBL used to use
    "join(complement(20,100),complement(1,10))" instead (but appears to have
    now adopted the GenBank convention). Notice that the order of the entries
    is reversed! This function therefore uses the first form. In this situation
    we expect the CompoundFeatureLocation and its parts to all be marked as
    strand == -1, and to be in the order 19:100 then 0:10.
    """
    try:
        parts = location.parts
        # CompoundFeatureLocation
        if location.strand == -1:
            # Special case, put complement outside the join/order/... and reverse order
            return "complement(%s(%s))" % (
                location.operator,
                ",".join(
                    _insdc_location_string_ignoring_strand_and_subfeatures(
                        p, rec_length
                    )
                    for p in parts[::-1]
                ),
            )
        else:
            return "%s(%s)" % (
                location.operator,
                ",".join(_insdc_location_string(p, rec_length) for p in parts),
            )
    except AttributeError:
        # Simple FeatureLocation
        loc = _insdc_location_string_ignoring_strand_and_subfeatures(
            location, rec_length
        )
        if location.strand == -1:
            return "complement(%s)" % loc
        else:
            return loc


class _InsdcWriter(SequenceWriter):
    """Base class for GenBank and EMBL writers (PRIVATE)."""

    MAX_WIDTH = 80
    QUALIFIER_INDENT = 21
    QUALIFIER_INDENT_STR = " " * QUALIFIER_INDENT
    QUALIFIER_INDENT_TMP = "     %s                "  # 21 if %s is empty
    FTQUAL_NO_QUOTE = (
        "anticodon",
        "citation",
        "codon_start",
        "compare",
        "direction",
        "estimated_length",
        "mod_base",
        "number",
        "rpt_type",
        "rpt_unit_range",
        "tag_peptide",
        "transl_except",
        "transl_table",
    )

    def _write_feature_qualifier(self, key, value=None, quote=None):
        if value is None:
            # Value-less entry like /pseudo
            self.handle.write("%s/%s\n" % (self.QUALIFIER_INDENT_STR, key))
            return

        if type(value) == str:
            value = value.replace(
                '"', '""'
            )  # NCBI says escape " as "" in qualifier values

        # Quick hack with no line wrapping, may be useful for testing:
        # self.handle.write('%s/%s="%s"\n' % (self.QUALIFIER_INDENT_STR, key, value))
        if quote is None:
            # Try to mimic unwritten rules about when quotes can be left out:
            if isinstance(value, int) or key in self.FTQUAL_NO_QUOTE:
                quote = False
            else:
                quote = True
        if quote:
            line = '%s/%s="%s"' % (self.QUALIFIER_INDENT_STR, key, value)
        else:
            line = "%s/%s=%s" % (self.QUALIFIER_INDENT_STR, key, value)
        if len(line) <= self.MAX_WIDTH:
            self.handle.write(line + "\n")
            return
        while line.lstrip():
            if len(line) <= self.MAX_WIDTH:
                self.handle.write(line + "\n")
                return
            # Insert line break...
            for index in range(
                min(len(line) - 1, self.MAX_WIDTH), self.QUALIFIER_INDENT + 1, -1
            ):
                if line[index] == " ":
                    break
            if line[index] != " ":
                # No nice place to break...
                index = self.MAX_WIDTH
            assert index <= self.MAX_WIDTH
            self.handle.write(line[:index] + "\n")
            line = self.QUALIFIER_INDENT_STR + line[index:].lstrip()

    def _wrap_location(self, location):
        """Split a feature location into lines (break at commas) (PRIVATE)."""
        # TODO - Rewrite this not to recurse!
        length = self.MAX_WIDTH - self.QUALIFIER_INDENT
        if len(location) <= length:
            return location
        index = location[:length].rfind(",")
        if index == -1:
            # No good place to split (!)
            warnings.warn("Couldn't split location:\n%s" % location, BiopythonWarning)
            return location
        return (
            location[: index + 1]
            + "\n"
            + self.QUALIFIER_INDENT_STR
            + self._wrap_location(location[index + 1 :])
        )

    def _write_feature(self, feature, record_length):
        """Write a single SeqFeature object to features table (PRIVATE)."""
        assert feature.type, feature
        location = _insdc_location_string(feature.location, record_length)
        f_type = feature.type.replace(" ", "_")
        line = (
            (self.QUALIFIER_INDENT_TMP % f_type)[: self.QUALIFIER_INDENT]
            + self._wrap_location(location)
            + "\n"
        )
        self.handle.write(line)
        # Now the qualifiers...
        # Note as of Biopython 1.69, this is an ordered-dict, don't sort it:
        for key, values in feature.qualifiers.items():
            if isinstance(values, (list, tuple)):
                for value in values:
                    self._write_feature_qualifier(key, value)
            else:
                # String, int, etc - or None for a /pseudo tpy entry
                self._write_feature_qualifier(key, values)

    @staticmethod
    def _get_annotation_str(record, key, default=".", just_first=False):
        """Get an annotation dictionary entry (as a string) (PRIVATE).

        Some entries are lists, in which case if just_first=True the first entry
        is returned.  If just_first=False (default) this verifies there is only
        one entry before returning it.
        """
        try:
            answer = record.annotations[key]
        except KeyError:
            return default
        if isinstance(answer, list):
            if not just_first:
                assert len(answer) == 1
            return str(answer[0])
        else:
            return str(answer)

    @staticmethod
    def _split_multi_line(text, max_len):
        """Return a list of strings (PRIVATE).

        Any single words which are too long get returned as a whole line
        (e.g. URLs) without an exception or warning.
        """
        # TODO - Do the line splitting while preserving white space?
        text = text.strip()
        if len(text) <= max_len:
            return [text]

        words = text.split()
        text = ""
        while words and len(text) + 1 + len(words[0]) <= max_len:
            text += " " + words.pop(0)
            text = text.strip()
        # assert len(text) <= max_len
        answer = [text]
        while words:
            text = words.pop(0)
            while words and len(text) + 1 + len(words[0]) <= max_len:
                text += " " + words.pop(0)
                text = text.strip()
            # assert len(text) <= max_len
            answer.append(text)
        assert not words
        return answer

    def _split_contig(self, record, max_len):
        """Return a list of strings, splits on commas (PRIVATE)."""
        # TODO - Merge this with _write_multi_line method?
        # It would need the addition of the comma splitting logic...
        # are there any other cases where that would be sensible?
        contig = record.annotations.get("contig", "")
        if isinstance(contig, (list, tuple)):
            contig = "".join(contig)
        contig = self.clean(contig)
        answer = []
        while contig:
            if len(contig) > max_len:
                # Split lines at the commas
                pos = contig[: max_len - 1].rfind(",")
                if pos == -1:
                    raise ValueError("Could not break up CONTIG")
                text, contig = contig[: pos + 1], contig[pos + 1 :]
            else:
                text, contig = contig, ""
            answer.append(text)
        return answer


class GenBankWriter(_InsdcWriter):
    """GenBank writer."""

    HEADER_WIDTH = 12
    QUALIFIER_INDENT = 21
    STRUCTURED_COMMENT_START = "-START##"
    STRUCTURED_COMMENT_END = "-END##"
    STRUCTURED_COMMENT_DELIM = " :: "
    LETTERS_PER_LINE = 60
    SEQUENCE_INDENT = 9

    def _write_single_line(self, tag, text):
        """Write single line in each GenBank record (PRIVATE).

        Used in the 'header' of each GenBank record.
        """
        assert len(tag) < self.HEADER_WIDTH
        if len(text) > self.MAX_WIDTH - self.HEADER_WIDTH:
            if tag:
                warnings.warn(
                    "Annotation %r too long for %r line" % (text, tag), BiopythonWarning
                )
            else:
                # Can't give such a precise warning
                warnings.warn("Annotation %r too long" % text, BiopythonWarning)
        self.handle.write(
            "%s%s\n" % (tag.ljust(self.HEADER_WIDTH), text.replace("\n", " "))
        )

    def _write_multi_line(self, tag, text):
        """Write multiple lines in each GenBank record (PRIVATE).

        Used in the 'header' of each GenBank record.
        """
        # TODO - Do the line splitting while preserving white space?
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_multi_line(text, max_len)
        self._write_single_line(tag, lines[0])
        for line in lines[1:]:
            self._write_single_line("", line)

    def _write_multi_entries(self, tag, text_list):
        # used for DBLINK and any similar later line types.
        # If the list of strings is empty, nothing is written.
        for i, text in enumerate(text_list):
            if i == 0:
                self._write_single_line(tag, text)
            else:
                self._write_single_line("", text)

    @staticmethod
    def _get_date(record):
        default = "01-JAN-1980"
        try:
            date = record.annotations["date"]
        except KeyError:
            return default
        # Cope with a list of one string:
        if isinstance(date, list) and len(date) == 1:
            date = date[0]
        if isinstance(date, datetime):
            date = date.strftime("%d-%b-%Y").upper()

        months = [
            "JAN",
            "FEB",
            "MAR",
            "APR",
            "MAY",
            "JUN",
            "JUL",
            "AUG",
            "SEP",
            "OCT",
            "NOV",
            "DEC",
        ]
        if not isinstance(date, str) or len(date) != 11:
            return default
        try:
            datetime(int(date[-4:]), months.index(date[3:6]) + 1, int(date[0:2]))
        except ValueError:
            date = default
        return date

    @staticmethod
    def _get_data_division(record):
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNK"
        if division in [
            "PRI",
            "ROD",
            "MAM",
            "VRT",
            "INV",
            "PLN",
            "BCT",
            "VRL",
            "PHG",
            "SYN",
            "UNA",
            "EST",
            "PAT",
            "STS",
            "GSS",
            "HTG",
            "HTC",
            "ENV",
            "CON",
            "TSA",
        ]:
            # Good, already GenBank style
            #    PRI - primate sequences
            #    ROD - rodent sequences
            #    MAM - other mammalian sequences
            #    VRT - other vertebrate sequences
            #    INV - invertebrate sequences
            #    PLN - plant, fungal, and algal sequences
            #    BCT - bacterial sequences [plus archaea]
            #    VRL - viral sequences
            #    PHG - bacteriophage sequences
            #    SYN - synthetic sequences
            #    UNA - unannotated sequences
            #    EST - EST sequences (expressed sequence tags)
            #    PAT - patent sequences
            #    STS - STS sequences (sequence tagged sites)
            #    GSS - GSS sequences (genome survey sequences)
            #    HTG - HTGS sequences (high throughput genomic sequences)
            #    HTC - HTC sequences (high throughput cDNA sequences)
            #    ENV - Environmental sampling sequences
            #    CON - Constructed sequences
            #    TSA - Transcriptome Shotgun Assembly
            #
            # (plus UNK for unknown)
            pass
        else:
            # See if this is in EMBL style:
            #    Division                 Code
            #    -----------------        ----
            #    Bacteriophage            PHG - common
            #    Environmental Sample     ENV - common
            #    Fungal                   FUN - map to PLN (plants + fungal)
            #    Human                    HUM - map to PRI (primates)
            #    Invertebrate             INV - common
            #    Other Mammal             MAM - common
            #    Other Vertebrate         VRT - common
            #    Mus musculus             MUS - map to ROD (rodent)
            #    Plant                    PLN - common
            #    Prokaryote               PRO - map to BCT (poor name)
            #    Other Rodent             ROD - common
            #    Synthetic                SYN - common
            #    Transgenic               TGN - ??? map to SYN ???
            #    Unclassified             UNC - map to UNK
            #    Viral                    VRL - common
            #
            # (plus XXX for submitting which we can map to UNK)
            embl_to_gbk = {
                "FUN": "PLN",
                "HUM": "PRI",
                "MUS": "ROD",
                "PRO": "BCT",
                "UNC": "UNK",
                "XXX": "UNK",
            }
            try:
                division = embl_to_gbk[division]
            except KeyError:
                division = "UNK"
        assert len(division) == 3
        return division

    def _get_topology(self, record):
        """Set the topology to 'circular', 'linear' if defined (PRIVATE)."""
        max_topology_len = len("circular")

        topology = self._get_annotation_str(record, "topology", default="")
        if topology and len(topology) <= max_topology_len:
            return topology.ljust(max_topology_len)
        else:
            return " " * max_topology_len

    def _write_the_first_line(self, record):
        """Write the LOCUS line (PRIVATE)."""
        locus = record.name
        if not locus or locus == "<unknown name>":
            locus = record.id
        if not locus or locus == "<unknown id>":
            locus = self._get_annotation_str(record, "accession", just_first=True)
        if len(locus) > 16:
            if len(locus) + 1 + len(str(len(record))) > 28:
                # Locus name and record length to long to squeeze in.
                # Per updated GenBank standard (Dec 15, 2018) 229.0
                # the Locus identifier can be any length, and a space
                # is added after the identifier to keep the identifier
                # and length fields separated
                warnings.warn(
                    "Increasing length of locus line to allow "
                    "long name. This will result in fields that "
                    "are not in usual positions.",
                    BiopythonWarning,
                )

        if len(locus.split()) > 1:
            raise ValueError("Invalid whitespace in %r for LOCUS line" % locus)
        if len(record) > 99999999999:
            # As of the GenBank release notes 229.0, the locus line can be
            # any length. However, long locus lines may not be compatible
            # with all software.
            warnings.warn(
                "The sequence length is very long. The LOCUS "
                "line will be increased in length to compensate. "
                "This may cause unexpected behavior.",
                BiopythonWarning,
            )

        # Get the molecule type
        mol_type = self._get_annotation_str(record, "molecule_type", None)
        if mol_type is None:
            raise ValueError("missing molecule_type in annotations")
        if mol_type and len(mol_type) > 7:
            # Deal with common cases from EMBL to GenBank
            mol_type = mol_type.replace("unassigned ", "").replace("genomic ", "")
            if len(mol_type) > 7:
                warnings.warn("Molecule type %r too long" % mol_type, BiopythonWarning)
                mol_type = "DNA"
        if mol_type in ["protein", "PROTEIN"]:
            mol_type = ""

        if mol_type == "":
            units = "aa"
        else:
            units = "bp"

        topology = self._get_topology(record)

        division = self._get_data_division(record)

        # Accommodate longer header, with long accessions and lengths
        if len(locus) > 16 and len(str(len(record))) > (11 - (len(locus) - 16)):
            name_length = locus + " " + str(len(record))

        # This is the older, standard 80 position header
        else:
            name_length = str(len(record)).rjust(28)
            name_length = locus + name_length[len(locus) :]
            assert len(name_length) == 28, name_length
            assert " " in name_length, name_length

        assert len(units) == 2
        assert len(division) == 3
        line = "LOCUS       %s %s    %s %s %s %s\n" % (
            name_length,
            units,
            mol_type.ljust(7),
            topology,
            division,
            self._get_date(record),
        )
        # Extra long header
        if len(line) > 80:
            splitline = line.split()
            if splitline[3] not in ["bp", "aa"]:
                raise ValueError(
                    "LOCUS line does not contain size units at "
                    "expected position:\n" + line
                )

            if not (
                splitline[4].strip() == ""
                or "DNA" in splitline[4].strip().upper()
                or "RNA" in splitline[4].strip().upper()
            ):
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "sequence type (DNA, RNA, ...):\n" + line
                )

            self.handle.write(line)

        # 80 position header
        else:
            assert len(line) == 79 + 1, repr(line)  # plus one for new line

            # We're bending the rules to allow an identifier over 16 characters
            # if we can steal spaces from the length field:
            # assert line[12:28].rstrip() == locus, \
            #     'LOCUS line does not contain the locus at the expected position:\n' + line
            # assert line[28:29] == " "
            # assert line[29:40].lstrip() == str(len(record)), \
            #     'LOCUS line does not contain the length at the expected position:\n' + line
            assert line[12:40].split() == [locus, str(len(record))], line

            # Tests copied from Bio.GenBank.Scanner
            if line[40:44] not in [" bp ", " aa "]:
                raise ValueError(
                    "LOCUS line does not contain size units at "
                    "expected position:\n" + line
                )
            if line[44:47] not in ["   ", "ss-", "ds-", "ms-"]:
                raise ValueError(
                    "LOCUS line does not have valid strand "
                    "type (Single stranded, ...):\n" + line
                )
            if not (
                line[47:54].strip() == ""
                or "DNA" in line[47:54].strip().upper()
                or "RNA" in line[47:54].strip().upper()
            ):
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "sequence type (DNA, RNA, ...):\n" + line
                )
            if line[54:55] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 55:\n" + line
                )
            if line[55:63].strip() not in ["", "linear", "circular"]:
                raise ValueError(
                    "LOCUS line does not contain valid "
                    "entry (linear, circular, ...):\n" + line
                )
            if line[63:64] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 64:\n" + line
                )
            if line[67:68] != " ":
                raise ValueError(
                    "LOCUS line does not contain space at position 68:\n" + line
                )
            if line[70:71] != "-":
                raise ValueError(
                    "LOCUS line does not contain - at position 71 in date:\n" + line
                )
            if line[74:75] != "-":
                raise ValueError(
                    "LOCUS line does not contain - at position 75 in date:\n" + line
                )

            self.handle.write(line)

    def _write_references(self, record):
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            data = str(number)
            # TODO - support more complex record reference locations?
            if ref.location and len(ref.location) == 1:
                molecule_type = record.annotations.get("molecule_type")
                if molecule_type and "protein" in molecule_type:
                    units = "residues"
                else:
                    units = "bases"
                data += "  (%s %i to %i)" % (
                    units,
                    ref.location[0].nofuzzy_start + 1,
                    ref.location[0].nofuzzy_end,
                )
            self._write_single_line("REFERENCE", data)
            if ref.authors:
                # We store the AUTHORS data as a single string
                self._write_multi_line("  AUTHORS", ref.authors)
            if ref.consrtm:
                # We store the consortium as a single string
                self._write_multi_line("  CONSRTM", ref.consrtm)
            if ref.title:
                # We store the title as a single string
                self._write_multi_line("  TITLE", ref.title)
            if ref.journal:
                # We store this as a single string - holds the journal name,
                # volume, year, and page numbers of the citation
                self._write_multi_line("  JOURNAL", ref.journal)
            if ref.medline_id:
                # This line type is obsolete and was removed from the GenBank
                # flatfile format in April 2005. Should we write it?
                # Note this has a two space indent:
                self._write_multi_line("  MEDLINE", ref.medline_id)
            if ref.pubmed_id:
                # Note this has a THREE space indent:
                self._write_multi_line("   PUBMED", ref.pubmed_id)
            if ref.comment:
                self._write_multi_line("  REMARK", ref.comment)

    def _write_comment(self, record):
        # This is a bit complicated due to the range of possible
        # ways people might have done their annotation...
        # Currently the parser uses a single string with newlines.
        # A list of lines is also reasonable.
        # A single (long) string is perhaps the most natural of all.
        # This means we may need to deal with line wrapping.
        lines = []
        if "structured_comment" in record.annotations:
            comment = record.annotations["structured_comment"]
            # Find max length of keys for equal padded printing
            padding = 0
            for key, data in comment.items():
                for subkey, subdata in data.items():
                    padding = len(subkey) if len(subkey) > padding else padding
            # Construct output
            for key, data in comment.items():
                lines.append(f"##{key}{self.STRUCTURED_COMMENT_START}")
                for subkey, subdata in data.items():
                    spaces = " " * (padding - len(subkey))
                    lines.append(
                        f"{subkey}{spaces}{self.STRUCTURED_COMMENT_DELIM}{subdata}"
                    )
                lines.append(f"##{key}{self.STRUCTURED_COMMENT_END}")
        if "comment" in record.annotations:
            comment = record.annotations["comment"]
            if isinstance(comment, str):
                lines += comment.split("\n")
            elif isinstance(comment, (list, tuple)):
                lines += list(comment)
            else:
                raise ValueError("Could not understand comment annotation")
        self._write_multi_line("COMMENT", lines[0])
        for line in lines[1:]:
            self._write_multi_line("", line)

    def _write_contig(self, record):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_contig(record, max_len)
        self._write_single_line("CONTIG", lines[0])
        for text in lines[1:]:
            self._write_single_line("", text)

    def _write_sequence(self, record):
        # Loosely based on code from Howard Salis
        # TODO - Force lower case?

        if isinstance(record.seq, UnknownSeq):
            # We have already recorded the length, and there is no need
            # to record a long sequence of NNNNNNN...NNN or whatever.
            if "contig" in record.annotations:
                self._write_contig(record)
            else:
                self.handle.write("ORIGIN\n")
            return

        # Catches sequence being None:
        data = self._get_seq_string(record).lower()
        seq_len = len(data)
        self.handle.write("ORIGIN\n")
        for line_number in range(0, seq_len, self.LETTERS_PER_LINE):
            self.handle.write(str(line_number + 1).rjust(self.SEQUENCE_INDENT))
            for words in range(
                line_number, min(line_number + self.LETTERS_PER_LINE, seq_len), 10
            ):
                self.handle.write(" %s" % data[words : words + 10])
            self.handle.write("\n")

    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_line(record)

        default = record.id
        if default.count(".") == 1 and default[default.index(".") + 1 :].isdigit():
            # Good, looks like accession.version and not something
            # else like identifier.start-end
            default = record.id.split(".", 1)[0]
        accession = self._get_annotation_str(
            record, "accession", default, just_first=True
        )
        acc_with_version = accession
        if record.id.startswith(accession + "."):
            try:
                acc_with_version = "%s.%i" % (
                    accession,
                    int(record.id.split(".", 1)[1]),
                )
            except ValueError:
                pass
        gi = self._get_annotation_str(record, "gi", just_first=True)

        descr = record.description
        if descr == "<unknown description>":
            descr = ""  # Trailing dot will be added later

        # The DEFINITION field must end with a period
        # see ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt [3.4.5]
        # and discussion https://github.com/biopython/biopython/pull/616
        # So let's add a period
        descr += "."
        self._write_multi_line("DEFINITION", descr)

        self._write_single_line("ACCESSION", accession)
        if gi != ".":
            self._write_single_line("VERSION", "%s  GI:%s" % (acc_with_version, gi))
        else:
            self._write_single_line("VERSION", "%s" % acc_with_version)

        # The NCBI initially expected two types of link,
        # e.g. "Project:28471" and "Trace Assembly Archive:123456"
        #
        # This changed and at some point the formatting switched to
        # include a space after the colon, e.g.
        #
        # LOCUS       NC_000011               1606 bp    DNA     linear   CON 06-JUN-2016
        # DEFINITION  Homo sapiens chromosome 11, GRCh38.p7 Primary Assembly.
        # ACCESSION   NC_000011 REGION: complement(5225466..5227071) GPC_000001303
        # VERSION     NC_000011.10  GI:568815587
        # DBLINK      BioProject: PRJNA168
        #             Assembly: GCF_000001405.33
        # ...
        #
        # Or,
        #
        # LOCUS       JU120277                1044 bp    mRNA    linear   TSA 27-NOV-2012
        # DEFINITION  TSA: Tupaia chinensis tbc000002.Tuchadli mRNA sequence.
        # ACCESSION   JU120277
        # VERSION     JU120277.1  GI:379775257
        # DBLINK      BioProject: PRJNA87013
        #             Sequence Read Archive: SRR433859
        # ...
        dbxrefs_with_space = []
        for x in record.dbxrefs:
            if ": " not in x:
                x = x.replace(":", ": ")
            dbxrefs_with_space.append(x)
        self._write_multi_entries("DBLINK", dbxrefs_with_space)
        del dbxrefs_with_space

        try:
            # List of strings
            # Keywords should be given separated with semi colons,
            keywords = "; ".join(record.annotations["keywords"])
            # with a trailing period:
            if not keywords.endswith("."):
                keywords += "."
        except KeyError:
            # If no keywords, there should be just a period:
            keywords = "."
        self._write_multi_line("KEYWORDS", keywords)

        if "segment" in record.annotations:
            # Deal with SEGMENT line found only in segmented records,
            # e.g. AH000819
            segment = record.annotations["segment"]
            if isinstance(segment, list):
                assert len(segment) == 1, segment
                segment = segment[0]
            self._write_single_line("SEGMENT", segment)

        self._write_multi_line("SOURCE", self._get_annotation_str(record, "source"))
        # The ORGANISM line MUST be a single line, as any continuation is the taxonomy
        org = self._get_annotation_str(record, "organism")
        if len(org) > self.MAX_WIDTH - self.HEADER_WIDTH:
            org = org[: self.MAX_WIDTH - self.HEADER_WIDTH - 4] + "..."
        self._write_single_line("  ORGANISM", org)
        try:
            # List of strings
            # Taxonomy should be given separated with semi colons,
            taxonomy = "; ".join(record.annotations["taxonomy"])
            # with a trailing period:
            if not taxonomy.endswith("."):
                taxonomy += "."
        except KeyError:
            taxonomy = "."
        self._write_multi_line("", taxonomy)

        if "references" in record.annotations:
            self._write_references(record)

        if (
            "comment" in record.annotations
            or "structured_comment" in record.annotations
        ):
            self._write_comment(record)

        handle.write("FEATURES             Location/Qualifiers\n")
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length)
        self._write_sequence(record)
        handle.write("//\n")


class EmblWriter(_InsdcWriter):
    """EMBL writer."""

    HEADER_WIDTH = 5
    QUALIFIER_INDENT = 21
    QUALIFIER_INDENT_STR = "FT" + " " * (QUALIFIER_INDENT - 2)
    QUALIFIER_INDENT_TMP = "FT   %s                "  # 21 if %s is empty
    # Note second spacer line of just FH is expected:
    FEATURE_HEADER = "FH   Key             Location/Qualifiers\nFH\n"

    LETTERS_PER_BLOCK = 10
    BLOCKS_PER_LINE = 6
    LETTERS_PER_LINE = LETTERS_PER_BLOCK * BLOCKS_PER_LINE
    POSITION_PADDING = 10

    def _write_contig(self, record):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_contig(record, max_len)
        for text in lines:
            self._write_single_line("CO", text)

    def _write_sequence(self, record):
        handle = self.handle  # save looking up this multiple times

        if isinstance(record.seq, UnknownSeq):
            # We have already recorded the length, and there is no need
            # to record a long sequence of NNNNNNN...NNN or whatever.
            if "contig" in record.annotations:
                self._write_contig(record)
            else:
                # TODO - Can the sequence just be left out as in GenBank files?
                handle.write("SQ   \n")
            return

        # Catches sequence being None
        data = self._get_seq_string(record).lower()
        seq_len = len(data)

        molecule_type = record.annotations.get("molecule_type")
        if molecule_type is not None and "DNA" in molecule_type:
            # TODO - What if we have RNA?
            a_count = data.count("A") + data.count("a")
            c_count = data.count("C") + data.count("c")
            g_count = data.count("G") + data.count("g")
            t_count = data.count("T") + data.count("t")
            other = seq_len - (a_count + c_count + g_count + t_count)
            handle.write(
                "SQ   Sequence %i BP; %i A; %i C; %i G; %i T; %i other;\n"
                % (seq_len, a_count, c_count, g_count, t_count, other)
            )
        else:
            handle.write("SQ   \n")

        for line_number in range(0, seq_len // self.LETTERS_PER_LINE):
            handle.write("    ")  # Just four, not five
            for block in range(self.BLOCKS_PER_LINE):
                index = (
                    self.LETTERS_PER_LINE * line_number + self.LETTERS_PER_BLOCK * block
                )
                handle.write(" %s" % data[index : index + self.LETTERS_PER_BLOCK])
            handle.write(
                str((line_number + 1) * self.LETTERS_PER_LINE).rjust(
                    self.POSITION_PADDING
                )
            )
            handle.write("\n")
        if seq_len % self.LETTERS_PER_LINE:
            # Final (partial) line
            line_number = seq_len // self.LETTERS_PER_LINE
            handle.write("    ")  # Just four, not five
            for block in range(self.BLOCKS_PER_LINE):
                index = (
                    self.LETTERS_PER_LINE * line_number + self.LETTERS_PER_BLOCK * block
                )
                handle.write(
                    (" %s" % data[index : index + self.LETTERS_PER_BLOCK]).ljust(11)
                )
            handle.write(str(seq_len).rjust(self.POSITION_PADDING))
            handle.write("\n")

    def _write_single_line(self, tag, text):
        assert len(tag) == 2
        line = tag + "   " + text
        if len(text) > self.MAX_WIDTH:
            warnings.warn("Line %r too long" % line, BiopythonWarning)
        self.handle.write(line + "\n")

    def _write_multi_line(self, tag, text):
        max_len = self.MAX_WIDTH - self.HEADER_WIDTH
        lines = self._split_multi_line(text, max_len)
        for line in lines:
            self._write_single_line(tag, line)

    def _write_the_first_lines(self, record):
        """Write the ID and AC lines (PRIVATE)."""
        if "." in record.id and record.id.rsplit(".", 1)[1].isdigit():
            version = "SV " + record.id.rsplit(".", 1)[1]
            accession = self._get_annotation_str(
                record, "accession", record.id.rsplit(".", 1)[0], just_first=True
            )
        else:
            version = ""
            accession = self._get_annotation_str(
                record, "accession", record.id, just_first=True
            )

        if ";" in accession:
            raise ValueError(
                "Cannot have semi-colon in EMBL accession, '%s'" % accession
            )
        if " " in accession:
            # This is out of practicality... might it be allowed?
            raise ValueError("Cannot have spaces in EMBL accession, '%s'" % accession)

        topology = self._get_annotation_str(record, "topology", default="")

        # Get the molecule type
        # TODO - record this explicitly in the parser?
        # Note often get RNA vs DNA discrepancy in real EMBL/NCBI files
        mol_type = record.annotations.get("molecule_type")
        if mol_type is None:
            raise ValueError("missing molecule_type in annotations")
        elif "DNA" in mol_type:
            units = "BP"
        elif "RNA" in mol_type:
            units = "BP"
        elif "PROTEIN" in mol_type.upper():
            mol_type = "PROTEIN"
            units = "AA"
        else:
            raise ValueError("failed to understand molecule_type '%s'" % mol_type)

        # Get the taxonomy division
        division = self._get_data_division(record)

        # TODO - Full ID line
        handle = self.handle
        # ID   <1>; SV <2>; <3>; <4>; <5>; <6>; <7> BP.
        # 1. Primary accession number
        # 2. Sequence version number
        # 3. Topology: 'circular' or 'linear'
        # 4. Molecule type
        # 5. Data class
        # 6. Taxonomic division
        # 7. Sequence length
        self._write_single_line(
            "ID",
            "%s; %s; %s; %s; ; %s; %i %s."
            % (accession, version, topology, mol_type, division, len(record), units),
        )
        handle.write("XX\n")
        self._write_single_line("AC", accession + ";")
        handle.write("XX\n")

    @staticmethod
    def _get_data_division(record):
        try:
            division = record.annotations["data_file_division"]
        except KeyError:
            division = "UNC"
        if division in [
            "PHG",
            "ENV",
            "FUN",
            "HUM",
            "INV",
            "MAM",
            "VRT",
            "MUS",
            "PLN",
            "PRO",
            "ROD",
            "SYN",
            "TGN",
            "UNC",
            "VRL",
            "XXX",
        ]:
            # Good, already EMBL style
            #    Division                 Code
            #    -----------------        ----
            #    Bacteriophage            PHG
            #    Environmental Sample     ENV
            #    Fungal                   FUN
            #    Human                    HUM
            #    Invertebrate             INV
            #    Other Mammal             MAM
            #    Other Vertebrate         VRT
            #    Mus musculus             MUS
            #    Plant                    PLN
            #    Prokaryote               PRO
            #    Other Rodent             ROD
            #    Synthetic                SYN
            #    Transgenic               TGN
            #    Unclassified             UNC (i.e. unknown)
            #    Viral                    VRL
            #
            # (plus XXX used for submitting data to EMBL)
            pass
        else:
            # See if this is in GenBank style & can be converted.
            # Generally a problem as the GenBank groups are wider
            # than those of EMBL. Note that GenBank use "BCT" for
            # both bacteria and archaea thus this maps to EMBL's
            # "PRO" nicely.
            gbk_to_embl = {"BCT": "PRO", "UNK": "UNC"}
            try:
                division = gbk_to_embl[division]
            except KeyError:
                division = "UNC"
        assert len(division) == 3
        return division

    def _write_keywords(self, record):
        # Put the keywords right after DE line.
        # Each 'keyword' can have multiple words and spaces, but we
        # must not split any 'keyword' between lines.
        # TODO - Combine short keywords onto one line
        for keyword in record.annotations["keywords"]:
            self._write_single_line("KW", keyword)
        self.handle.write("XX\n")

    def _write_references(self, record):
        # The order should be RN, RC, RP, RX, RG, RA, RT, RL
        number = 0
        for ref in record.annotations["references"]:
            if not isinstance(ref, SeqFeature.Reference):
                continue
            number += 1
            self._write_single_line("RN", "[%i]" % number)
            # TODO - support for RC line (needed in parser too)
            # TODO - support more complex record reference locations?
            if ref.location and len(ref.location) == 1:
                self._write_single_line(
                    "RP",
                    "%i-%i"
                    % (ref.location[0].nofuzzy_start + 1, ref.location[0].nofuzzy_end),
                )
            # TODO - record any DOI or AGRICOLA identifier in the reference object?
            if ref.pubmed_id:
                self._write_single_line("RX", "PUBMED; %s." % ref.pubmed_id)
            if ref.consrtm:
                self._write_single_line("RG", "%s" % ref.consrtm)
            if ref.authors:
                # We store the AUTHORS data as a single string
                self._write_multi_line("RA", ref.authors + ";")
            if ref.title:
                # We store the title as a single string
                self._write_multi_line("RT", '"%s";' % ref.title)
            if ref.journal:
                # We store this as a single string - holds the journal name,
                # volume, year, and page numbers of the citation
                self._write_multi_line("RL", ref.journal)
            self.handle.write("XX\n")

    def _write_comment(self, record):
        # This is a bit complicated due to the range of possible
        # ways people might have done their annotation...
        # Currently the parser uses a single string with newlines.
        # A list of lines is also reasonable.
        # A single (long) string is perhaps the most natural of all.
        # This means we may need to deal with line wrapping.
        comment = record.annotations["comment"]
        if isinstance(comment, str):
            lines = comment.split("\n")
        elif isinstance(comment, (list, tuple)):
            lines = comment
        else:
            raise ValueError("Could not understand comment annotation")
        # TODO - Merge this with the GenBank comment code?
        if not lines:
            return
        for line in lines:
            self._write_multi_line("CC", line)
        self.handle.write("XX\n")

    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        self._write_the_first_lines(record)

        # PR line (0 or 1 lines only), project identifier
        #
        # Assuming can't use 2 lines, we should prefer newer GenBank
        # DBLINK BioProject:... entries over the older GenBank DBLINK
        # Project:... lines.
        #
        # In either case, seems EMBL uses just "PR    Project:..."
        # regardless of the type of ID (old numeric only, or new
        # with alpha prefix), e.g. for CP002497 NCBI now uses:
        #
        # DBLINK      BioProject: PRJNA60715
        #             BioSample: SAMN03081426
        #
        # While EMBL uses:
        #
        # XX
        # PR   Project:PRJNA60715;
        # XX
        #
        # Sorting ensures (new) BioProject:... is before old Project:...
        for xref in sorted(record.dbxrefs):
            if xref.startswith("BioProject:"):
                self._write_single_line("PR", xref[3:] + ";")
                handle.write("XX\n")
                break
            if xref.startswith("Project:"):
                self._write_single_line("PR", xref + ";")
                handle.write("XX\n")
                break

        # TODO - DT lines (date)

        descr = record.description
        if descr == "<unknown description>":
            descr = "."
        self._write_multi_line("DE", descr)
        handle.write("XX\n")

        if "keywords" in record.annotations:
            self._write_keywords(record)

        # Should this be "source" or "organism"?
        self._write_multi_line("OS", self._get_annotation_str(record, "organism"))
        try:
            # List of strings
            taxonomy = "; ".join(record.annotations["taxonomy"]) + "."
        except KeyError:
            taxonomy = "."
        self._write_multi_line("OC", taxonomy)
        handle.write("XX\n")

        if "references" in record.annotations:
            self._write_references(record)

        if "comment" in record.annotations:
            self._write_comment(record)

        handle.write(self.FEATURE_HEADER)
        rec_length = len(record)
        for feature in record.features:
            self._write_feature(feature, rec_length)
        handle.write("XX\n")

        self._write_sequence(record)
        handle.write("//\n")


class ImgtWriter(EmblWriter):
    """IMGT writer (EMBL format variant)."""

    HEADER_WIDTH = 5
    QUALIFIER_INDENT = 25  # Not 21 as in EMBL
    QUALIFIER_INDENT_STR = "FT" + " " * (QUALIFIER_INDENT - 2)
    QUALIFIER_INDENT_TMP = "FT   %s                    "  # 25 if %s is empty
    FEATURE_HEADER = "FH   Key                 Location/Qualifiers\nFH\n"


def _genbank_convert_fasta(in_file, out_file):
    """Fast GenBank to FASTA (PRIVATE)."""
    # We don't need to parse the features...
    records = GenBankScanner().parse_records(in_file, do_features=False)
    return SeqIO.write(records, out_file, "fasta")


def _embl_convert_fasta(in_file, out_file):
    """Fast EMBL to FASTA (PRIVATE)."""
    # We don't need to parse the features...
    records = EmblScanner().parse_records(in_file, do_features=False)
    return SeqIO.write(records, out_file, "fasta")


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
