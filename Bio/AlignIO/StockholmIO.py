# Copyright 2006-2016 by Peter Cock.  All rights reserved.
# Revisions copyright 2015 by Ben Woodcroft.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.AlignIO support for "stockholm" format (used in the PFAM database).

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).

For example, consider a Stockholm alignment file containing the following::

    # STOCKHOLM 1.0
    #=GC SS_cons       .................<<<<<<<<...<<<<<<<........>>>>>>>..
    AP001509.1         UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGU
    #=GR AP001509.1 SS -----------------<<<<<<<<---..<<-<<-------->>->>..--
    AE007476.1         AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGU
    #=GR AE007476.1 SS -----------------<<<<<<<<-----<<.<<-------->>.>>----

    #=GC SS_cons       ......<<<<<<<.......>>>>>>>..>>>>>>>>...............
    AP001509.1         CUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU
    #=GR AP001509.1 SS -------<<<<<--------->>>>>--->>>>>>>>---------------
    AE007476.1         UUCUACAAGGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU
    #=GR AE007476.1 SS ------.<<<<<--------->>>>>.-->>>>>>>>---------------
    //

This is a single multiple sequence alignment, so you would probably load this
using the Bio.AlignIO.read() function:

    >>> from Bio import AlignIO
    >>> align = AlignIO.read("Stockholm/simple.sth", "stockholm")
    >>> print(align)
    Alignment with 2 rows and 104 columns
    UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-G...UGU AP001509.1
    AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-C...GAU AE007476.1
    >>> for record in align:
    ...     print("%s %i" % (record.id, len(record)))
    AP001509.1 104
    AE007476.1 104

In addition to the sequences themselves, this example alignment also includes
some GR lines for the secondary structure of the sequences.  These are
strings, with one character for each letter in the associated sequence:

    >>> for record in align:
    ...     print(record.id)
    ...     print(record.seq)
    ...     print(record.letter_annotations['secondary_structure'])
    AP001509.1
    UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGUCUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU
    -----------------<<<<<<<<---..<<-<<-------->>->>..---------<<<<<--------->>>>>--->>>>>>>>---------------
    AE007476.1
    AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGUUUCUACAAGGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU
    -----------------<<<<<<<<-----<<.<<-------->>.>>----------.<<<<<--------->>>>>.-->>>>>>>>---------------

Any general annotation for each row is recorded in the SeqRecord's annotations
dictionary.  Any per-column annotation for the entire alignment in in the
alignment's column annotations dictionary, such as the secondary structure
consensus in this example:

    >>> sorted(align.column_annotations.keys())
    ['secondary_structure']
    >>> align.column_annotations["secondary_structure"]
    '.................<<<<<<<<...<<<<<<<........>>>>>>>........<<<<<<<.......>>>>>>>..>>>>>>>>...............'

You can output this alignment in many different file formats
using Bio.AlignIO.write(), or the MultipleSeqAlignment object's format method:

    >>> print(format(align, "fasta"))
    >AP001509.1
    UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGUCUCUAC-A
    GGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU
    >AE007476.1
    AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGUUUCUACAA
    GGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU
    <BLANKLINE>

Most output formats won't be able to hold the annotation possible in a
Stockholm file:

    >>> print(format(align, "stockholm"))
    # STOCKHOLM 1.0
    #=GF SQ 2
    AP001509.1 UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGUCUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU
    #=GS AP001509.1 AC AP001509.1
    #=GS AP001509.1 DE AP001509.1
    #=GR AP001509.1 SS -----------------<<<<<<<<---..<<-<<-------->>->>..---------<<<<<--------->>>>>--->>>>>>>>---------------
    AE007476.1 AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGUUUCUACAAGGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU
    #=GS AE007476.1 AC AE007476.1
    #=GS AE007476.1 DE AE007476.1
    #=GR AE007476.1 SS -----------------<<<<<<<<-----<<.<<-------->>.>>----------.<<<<<--------->>>>>.-->>>>>>>>---------------
    #=GC SS_cons .................<<<<<<<<...<<<<<<<........>>>>>>>........<<<<<<<.......>>>>>>>..>>>>>>>>...............
    //
    <BLANKLINE>

Note that when writing Stockholm files, AlignIO does not break long sequences
up and interleave them (as in the input file shown above).  The standard
allows this simpler layout, and it is more likely to be understood by other
tools.

Finally, as an aside, it can sometimes be useful to use Bio.SeqIO.parse() to
iterate over the alignment rows as SeqRecord objects - rather than working
with Alignnment objects.

    >>> from Bio import SeqIO
    >>> for record in SeqIO.parse("Stockholm/simple.sth", "stockholm"):
    ...     print(record.id)
    ...     print(record.seq)
    ...     print(record.letter_annotations['secondary_structure'])
    AP001509.1
    UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGUCUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU
    -----------------<<<<<<<<---..<<-<<-------->>->>..---------<<<<<--------->>>>>--->>>>>>>>---------------
    AE007476.1
    AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGUUUCUACAAGGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU
    -----------------<<<<<<<<-----<<.<<-------->>.>>----------.<<<<<--------->>>>>.-->>>>>>>>---------------

Remember that if you slice a SeqRecord, the per-letter-annotations like the
secondary structure string here, are also sliced:

    >>> sub_record = record[10:20]
    >>> print(sub_record.seq)
    AUCGUUUUAC
    >>> print(sub_record.letter_annotations['secondary_structure'])
    -------<<<

Likewise with the alignment object, as long as you are not dropping any rows,
slicing specific columns of an alignment will slice any per-column-annotations:

    >>> align.column_annotations["secondary_structure"]
    '.................<<<<<<<<...<<<<<<<........>>>>>>>........<<<<<<<.......>>>>>>>..>>>>>>>>...............'
    >>> part_align = align[:,10:20]
    >>> part_align.column_annotations["secondary_structure"]
    '.......<<<'

You can also see this in the Stockholm output of this partial-alignment:

    >>> print(format(part_align, "stockholm"))
    # STOCKHOLM 1.0
    #=GF SQ 2
    AP001509.1 UCAACACUCU
    #=GS AP001509.1 AC AP001509.1
    #=GS AP001509.1 DE AP001509.1
    #=GR AP001509.1 SS -------<<<
    AE007476.1 AUCGUUUUAC
    #=GS AE007476.1 AC AE007476.1
    #=GS AE007476.1 DE AE007476.1
    #=GR AE007476.1 SS -------<<<
    #=GC SS_cons .......<<<
    //
    <BLANKLINE>

"""
from collections import OrderedDict

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import AlignmentIterator
from .Interfaces import SequentialAlignmentWriter


class StockholmWriter(SequentialAlignmentWriter):
    """Stockholm/PFAM alignment writer."""

    # These dictionaries should be kept in sync with those
    # defined in the StockholmIterator class.
    pfam_gr_mapping = {
        "secondary_structure": "SS",
        "surface_accessibility": "SA",
        "transmembrane": "TM",
        "posterior_probability": "PP",
        "ligand_binding": "LI",
        "active_site": "AS",
        "intron": "IN",
    }
    # These GC mappings are in addition to *_cons in GR mapping:
    pfam_gc_mapping = {"reference_annotation": "RF", "model_mask": "MM"}
    # Following dictionary deliberately does not cover AC, DE or DR
    pfam_gs_mapping = {"organism": "OS", "organism_classification": "OC", "look": "LO"}

    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file.

        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        count = len(alignment)

        self._length_of_sequences = alignment.get_alignment_length()
        self._ids_written = []

        if count == 0:
            raise ValueError("Must have at least one sequence")
        if self._length_of_sequences == 0:
            raise ValueError("Non-empty sequences are required")

        self.handle.write("# STOCKHOLM 1.0\n")
        self.handle.write("#=GF SQ %i\n" % count)
        for record in alignment:
            self._write_record(record)
        # This shouldn't be None... but just in case,
        if alignment.column_annotations:
            for k, v in sorted(alignment.column_annotations.items()):
                if k in self.pfam_gc_mapping:
                    self.handle.write("#=GC %s %s\n" % (self.pfam_gc_mapping[k], v))
                elif k in self.pfam_gr_mapping:
                    self.handle.write(
                        "#=GC %s %s\n" % (self.pfam_gr_mapping[k] + "_cons", v)
                    )
                else:
                    # It doesn't follow the PFAM standards, but should we record
                    # this data anyway?
                    pass
        self.handle.write("//\n")

    def _write_record(self, record):
        """Write a single SeqRecord to the file (PRIVATE)."""
        if self._length_of_sequences != len(record.seq):
            raise ValueError("Sequences must all be the same length")

        # For the case for stockholm to stockholm, try and use record.name
        seq_name = record.id
        if record.name is not None:
            if "accession" in record.annotations:
                if record.id == record.annotations["accession"]:
                    seq_name = record.name

        # In the Stockholm file format, spaces are not allowed in the id
        seq_name = seq_name.replace(" ", "_")

        if "start" in record.annotations and "end" in record.annotations:
            suffix = "/%s-%s" % (
                record.annotations["start"],
                record.annotations["end"],
            )
            if seq_name[-len(suffix) :] != suffix:
                seq_name = "%s/%s-%s" % (
                    seq_name,
                    record.annotations["start"],
                    record.annotations["end"],
                )

        if seq_name in self._ids_written:
            raise ValueError("Duplicate record identifier: %s" % seq_name)
        self._ids_written.append(seq_name)
        self.handle.write("%s %s\n" % (seq_name, record.seq))

        # The recommended placement for GS lines (per sequence annotation)
        # is above the alignment (as a header block) or just below the
        # corresponding sequence.
        #
        # The recommended placement for GR lines (per sequence per column
        # annotation such as secondary structure) is just below the
        # corresponding sequence.
        #
        # We put both just below the corresponding sequence as this allows
        # us to write the file using a single pass through the records.

        # AC = Accession
        if "accession" in record.annotations:
            self.handle.write(
                "#=GS %s AC %s\n"
                % (seq_name, self.clean(record.annotations["accession"]))
            )
        elif record.id:
            self.handle.write("#=GS %s AC %s\n" % (seq_name, self.clean(record.id)))

        # DE = description
        if record.description:
            self.handle.write(
                "#=GS %s DE %s\n" % (seq_name, self.clean(record.description))
            )

        # DE = database links
        for xref in record.dbxrefs:
            self.handle.write("#=GS %s DR %s\n" % (seq_name, self.clean(xref)))

        # GS = other per sequence annotation
        for key, value in record.annotations.items():
            if key in self.pfam_gs_mapping:
                data = self.clean(str(value))
                if data:
                    self.handle.write(
                        "#=GS %s %s %s\n"
                        % (seq_name, self.clean(self.pfam_gs_mapping[key]), data)
                    )
            else:
                # It doesn't follow the PFAM standards, but should we record
                # this data anyway?
                pass

        # GR = per row per column sequence annotation
        for key, value in record.letter_annotations.items():
            if key in self.pfam_gr_mapping and len(str(value)) == len(record.seq):
                data = self.clean(str(value))
                if data:
                    self.handle.write(
                        "#=GR %s %s %s\n"
                        % (seq_name, self.clean(self.pfam_gr_mapping[key]), data)
                    )
            else:
                # It doesn't follow the PFAM standards, but should we record
                # this data anyway?
                pass


class StockholmIterator(AlignmentIterator):
    """Loads a Stockholm file from PFAM into MultipleSeqAlignment objects.

    The file may contain multiple concatenated alignments, which are loaded
    and returned incrementally.

    This parser will detect if the Stockholm file follows the PFAM
    conventions for sequence specific meta-data (lines starting #=GS
    and #=GR) and populates the SeqRecord fields accordingly.

    Any annotation which does not follow the PFAM conventions is currently
    ignored.

    If an accession is provided for an entry in the meta data, IT WILL NOT
    be used as the record.id (it will be recorded in the record's
    annotations).  This is because some files have (sub) sequences from
    different parts of the same accession (differentiated by different
    start-end positions).

    Wrap-around alignments are not supported - each sequences must be on
    a single line.  However, interlaced sequences should work.

    For more information on the file format, please see:
    http://sonnhammer.sbc.su.se/Stockholm.html
    https://en.wikipedia.org/wiki/Stockholm_format
    http://bioperl.org/formats/alignment_formats/Stockholm_multiple_alignment_format.html

    For consistency with BioPerl and EMBOSS we call this the "stockholm"
    format.
    """

    # These dictionaries should be kept in sync with those
    # defined in the PfamStockholmWriter class.
    pfam_gr_mapping = {
        "SS": "secondary_structure",
        "SA": "surface_accessibility",
        "TM": "transmembrane",
        "PP": "posterior_probability",
        "LI": "ligand_binding",
        "AS": "active_site",
        "IN": "intron",
    }
    # These GC mappings are in addition to *_cons in GR mapping:
    pfam_gc_mapping = {"RF": "reference_annotation", "MM": "model_mask"}
    # Following dictionary deliberately does not cover AC, DE or DR
    pfam_gs_mapping = {"OS": "organism", "OC": "organism_classification", "LO": "look"}

    _header = None  # for caching lines between __next__ calls

    def __next__(self):
        """Parse the next alignment from the handle."""
        handle = self.handle

        if self._header is None:
            line = handle.readline()
        else:
            # Header we saved from when we were parsing
            # the previous alignment.
            line = self._header
            self._header = None

        if not line:
            # Empty file - just give up.
            raise StopIteration
        if line.strip() != "# STOCKHOLM 1.0":
            raise ValueError("Did not find STOCKHOLM header")

        # Note: If this file follows the PFAM conventions, there should be
        # a line containing the number of sequences, e.g. "#=GF SQ 67"
        # We do not check for this - perhaps we should, and verify that
        # if present it agrees with our parsing.

        seqs = {}
        ids = OrderedDict()  # Really only need an OrderedSet, but python lacks this
        gs = {}
        gr = {}
        gf = {}
        gc = {}
        passed_end_alignment = False
        while True:
            line = handle.readline()
            if not line:
                break  # end of file
            line = line.strip()  # remove trailing \n
            if line == "# STOCKHOLM 1.0":
                self._header = line
                break
            elif line == "//":
                # The "//" line indicates the end of the alignment.
                # There may still be more meta-data
                passed_end_alignment = True
            elif line == "":
                # blank line, ignore
                pass
            elif line[0] != "#":
                # Sequence
                # Format: "<seqname> <sequence>"
                assert not passed_end_alignment
                parts = [x.strip() for x in line.split(" ", 1)]
                if len(parts) != 2:
                    # This might be someone attempting to store a zero length sequence?
                    raise ValueError(
                        "Could not split line into identifier and sequence:\n" + line
                    )
                seq_id, seq = parts
                if seq_id not in ids:
                    ids[seq_id] = True
                seqs.setdefault(seq_id, "")
                seqs[seq_id] += seq.replace(".", "-")
            elif len(line) >= 5:
                # Comment line or meta-data
                if line[:5] == "#=GF ":
                    # Generic per-File annotation, free text
                    # Format: #=GF <feature> <free text>
                    feature, text = line[5:].strip().split(None, 1)
                    # Each feature key could be used more than once,
                    # so store the entries as a list of strings.
                    if feature not in gf:
                        gf[feature] = [text]
                    else:
                        gf[feature].append(text)
                elif line[:5] == "#=GC ":
                    # Generic per-Column annotation, exactly 1 char per column
                    # Format: "#=GC <feature> <exactly 1 char per column>"
                    feature, text = line[5:].strip().split(None, 2)
                    if feature not in gc:
                        gc[feature] = ""
                    gc[feature] += text.strip()  # append to any previous entry
                    # Might be interleaved blocks, so can't check length yet
                elif line[:5] == "#=GS ":
                    # Generic per-Sequence annotation, free text
                    # Format: "#=GS <seqname> <feature> <free text>"
                    try:
                        seq_id, feature, text = line[5:].strip().split(None, 2)
                    except ValueError:
                        # Free text can sometimes be empty, which a one line split throws an error for.
                        # See https://github.com/biopython/biopython/issues/2982 for more details
                        seq_id, feature = line[5:].strip().split(None, 1)
                        text = ""
                    # if seq_id not in ids:
                    #    ids.append(seq_id)
                    if seq_id not in gs:
                        gs[seq_id] = {}
                    if feature not in gs[seq_id]:
                        gs[seq_id][feature] = [text]
                    else:
                        gs[seq_id][feature].append(text)
                elif line[:5] == "#=GR ":
                    # Generic per-Sequence AND per-Column markup
                    # Format: "#=GR <seqname> <feature> <exactly 1 char per column>"
                    seq_id, feature, text = line[5:].strip().split(None, 2)
                    # if seq_id not in ids:
                    #    ids.append(seq_id)
                    if seq_id not in gr:
                        gr[seq_id] = {}
                    if feature not in gr[seq_id]:
                        gr[seq_id][feature] = ""
                    gr[seq_id][feature] += text.strip()  # append to any previous entry
                    # Might be interleaved blocks, so can't check length yet
            # Next line...

        assert len(seqs) <= len(ids)
        # assert len(gs)   <= len(ids)
        # assert len(gr)   <= len(ids)

        self.ids = ids.keys()
        self.sequences = seqs
        self.seq_annotation = gs
        self.seq_col_annotation = gr

        if ids and seqs:

            if (
                self.records_per_alignment is not None
                and self.records_per_alignment != len(ids)
            ):
                raise ValueError(
                    "Found %i records in this alignment, told to expect %i"
                    % (len(ids), self.records_per_alignment)
                )

            alignment_length = len(list(seqs.values())[0])
            records = []  # Alignment obj will put them all in a list anyway
            for seq_id in ids:
                seq = seqs[seq_id]
                if alignment_length != len(seq):
                    raise ValueError(
                        "Sequences have different lengths, or repeated identifier"
                    )
                name, start, end = self._identifier_split(seq_id)
                record = SeqRecord(
                    Seq(seq),
                    id=seq_id,
                    name=name,
                    description=seq_id,
                    annotations={"accession": name},
                )
                # Accession will be overridden by _populate_meta_data if an explicit
                # accession is provided:
                record.annotations["accession"] = name

                if start is not None:
                    record.annotations["start"] = start
                if end is not None:
                    record.annotations["end"] = end

                self._populate_meta_data(seq_id, record)
                records.append(record)
            for k, v in gc.items():
                if len(v) != alignment_length:
                    raise ValueError(
                        "%s length %i, expected %i" % (k, len(v), alignment_length)
                    )
            alignment = MultipleSeqAlignment(records)

            for k, v in sorted(gc.items()):
                if k in self.pfam_gc_mapping:
                    alignment.column_annotations[self.pfam_gc_mapping[k]] = v
                elif k.endswith("_cons") and k[:-5] in self.pfam_gr_mapping:
                    alignment.column_annotations[self.pfam_gr_mapping[k[:-5]]] = v
                else:
                    # Ignore it?
                    alignment.column_annotations["GC:" + k] = v

            # TODO - Introduce an annotated alignment class?
            # For now, store the annotation a new private property:
            alignment._annotations = gr

            return alignment
        else:
            raise StopIteration

    def _identifier_split(self, identifier):
        """Return (name, start, end) string tuple from an identifier (PRIVATE)."""
        if "/" in identifier:
            name, start_end = identifier.rsplit("/", 1)
            if start_end.count("-") == 1:
                try:
                    start, end = start_end.split("-")
                    return name, int(start), int(end)
                except ValueError:
                    # Non-integers after final '/' - fall through
                    pass
        return identifier, None, None

    def _get_meta_data(self, identifier, meta_dict):
        """Take an itentifier and returns dict of all meta-data matching it (PRIVATE).

        For example, given "Q9PN73_CAMJE/149-220" will return all matches to
        this or "Q9PN73_CAMJE" which the identifier without its /start-end
        suffix.

        In the example below, the suffix is required to match the AC, but must
        be removed to match the OS and OC meta-data::

            # STOCKHOLM 1.0
            #=GS Q9PN73_CAMJE/149-220  AC Q9PN73
            ...
            Q9PN73_CAMJE/149-220               NKA...
            ...
            #=GS Q9PN73_CAMJE OS Campylobacter jejuni
            #=GS Q9PN73_CAMJE OC Bacteria

        This function will return an empty dictionary if no data is found.
        """
        name, start, end = self._identifier_split(identifier)
        if name == identifier:
            identifier_keys = [identifier]
        else:
            identifier_keys = [identifier, name]
        answer = {}
        for identifier_key in identifier_keys:
            try:
                for feature_key in meta_dict[identifier_key]:
                    answer[feature_key] = meta_dict[identifier_key][feature_key]
            except KeyError:
                pass
        return answer

    def _populate_meta_data(self, identifier, record):
        """Add meta-date to a SecRecord's annotations dictionary (PRIVATE).

        This function applies the PFAM conventions.
        """
        seq_data = self._get_meta_data(identifier, self.seq_annotation)
        for feature in seq_data:
            # Note this dictionary contains lists!
            if feature == "AC":  # ACcession number
                assert len(seq_data[feature]) == 1
                record.annotations["accession"] = seq_data[feature][0]
            elif feature == "DE":  # DEscription
                record.description = "\n".join(seq_data[feature])
            elif feature == "DR":  # Database Reference
                # Should we try and parse the strings?
                record.dbxrefs = seq_data[feature]
            elif feature in self.pfam_gs_mapping:
                record.annotations[self.pfam_gs_mapping[feature]] = ", ".join(
                    seq_data[feature]
                )
            else:
                # Ignore it?
                record.annotations["GS:" + feature] = ", ".join(seq_data[feature])

        # Now record the per-letter-annotations
        seq_col_data = self._get_meta_data(identifier, self.seq_col_annotation)
        for feature in seq_col_data:
            # Note this dictionary contains strings!
            if feature in self.pfam_gr_mapping:
                record.letter_annotations[self.pfam_gr_mapping[feature]] = seq_col_data[
                    feature
                ]
            else:
                # Ignore it?
                record.letter_annotations["GR:" + feature] = seq_col_data[feature]


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
