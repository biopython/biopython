# Copyright 2006-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for alignment files in the Stockholm file format.

You are expected to use this module via the Bio.Align functions.

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
using the Bio.Align.read() function:

    >>> from Bio import Align
    >>> alignment = Align.read("Stockholm/simple.sth", "stockholm")
    # >>> print(align)
    # Alignment with 2 rows and 104 columns
    # UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-G...UGU AP001509.1
    # AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-C...GAU AE007476.1
    >>> for record in alignment.sequences:
    ...     print("%s %i" % (record.id, len(record)))
    AP001509.1 104
    AE007476.1 104

In addition to the sequences themselves, this example alignment also includes
some GR lines for the secondary structure of the sequences.  These are
strings, with one character for each letter in the associated sequence:

    >>> for record in alignment:
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
from collections import defaultdict

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Stockholm file format."""

    # These dictionaries should be kept in sync with those
    # defined in the AlignmentIterator class.
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

    def write_header(self, alignments):
        """Use this to write the file header."""
        stream = self.stream

    def write_alignment(self, alignment):
        """Use this to write the alignment to an open file.

        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        stream = self.stream
        count = len(alignment)

        n, self._length_of_sequences = alignment.shape()
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
                    self.handle.write(f"#=GC {self.pfam_gc_mapping[k]} {v}\n")
                elif k in self.pfam_gr_mapping:
                    self.handle.write(f"#=GC {self.pfam_gr_mapping[k]}_cons {v}\n")
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
            suffix = f"/{record.annotations['start']}-{record.annotations['end']}"
            if seq_name[-len(suffix) :] != suffix:
                seq_name = "%s/%s-%s" % (
                    seq_name,
                    record.annotations["start"],
                    record.annotations["end"],
                )

        if seq_name in self._ids_written:
            raise ValueError(f"Duplicate record identifier: {seq_name}")
        self._ids_written.append(seq_name)
        self.handle.write(f"{seq_name} {record.seq}\n")

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
                f"#=GS {seq_name} AC {self.clean(record.annotations['accession'])}\n"
            )
        elif record.id:
            self.handle.write(f"#=GS {seq_name} AC {self.clean(record.id)}\n")

        # DE = description
        if record.description:
            self.handle.write(f"#=GS {seq_name} DE {self.clean(record.description)}\n")

        # DE = database links
        for xref in record.dbxrefs:
            self.handle.write(f"#=GS {seq_name} DR {self.clean(xref)}\n")

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



class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for PFAM alignment files in the Stockholm format.

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

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="Clustal")
        stream = self.stream
        try:
            line = next(stream)
        except StopIteration:
            raise ValueError("Empty file.") from None
        if line.strip() != "# STOCKHOLM 1.0":
            raise ValueError("Did not find STOCKHOLM header")

    def _create_alignment(self, aligned_sequences, records, strands, gc, gf, gr):
        gf = dict(gf)
        for key, value in gf.items():
            if len(value) > 1:
                raise Exception(key)
        for k, v in gc.items():
            if len(v) != length:
                raise ValueError(f"{k} length is {len(v)}, expected {length}")
        coordinates = Alignment.infer_coordinates(aligned_sequences)
        for i, (record, strand) in enumerate(zip(records, strands)):
            if strand == "-":
                n = len(record.seq)
                coordinates[i, :] = n - coordinates[i, :]
        alignment = Alignment(records, coordinates)

        for k, v in sorted(gc.items()):
            if k in self.pfam_gc_mapping:
                alignment.column_annotations[self.pfam_gc_mapping[k]] = v
            elif k.endswith("_cons") and k[:-5] in self.pfam_gr_mapping:
                alignment.column_annotations[self.pfam_gr_mapping[k[:-5]]] = v
            else:
                # Ignore it?
                alignment.column_annotations["GC:" + k] = v

        alignment.annotations = gr

        return alignment

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        # Note: If this file follows the PFAM conventions, there should be
        # a line containing the number of sequences, e.g. "#=GF SQ 67"
        # We do not check for this - perhaps we should, and verify that
        # if present it agrees with our parsing.

        records = []
        aligned_sequences = []
        strands = []
        length = None
        gs = {}
        gr = {}
        gf = defaultdict(list)
        gc = {}
        passed_end_alignment = False
        for line in stream:
            line = line.strip()  # remove trailing \n
            if not line:
                continue
            elif line == "# STOCKHOLM 1.0":
                yield self._create_alignment(aligned_sequences, records, strands, gc, gf, gr)
                # Starting a new alignment
                passed_end_alignment = False
            elif line == "//":
                # The "//" line indicates the end of the alignment.
                # There may still be more meta-data
                passed_end_alignment = True
            elif not line.startswith("#"):
                # Sequence
                # Format: "<seqname> <sequence>"
                assert not passed_end_alignment
                try:
                    seqname, aligned_sequence = line.split(None, 1)
                except ValueError:
                    # This might be someone attempting to store a zero length sequence?
                    raise ValueError(
                        "Could not split line into sequence name and aligned sequence:\n" + line
                    ) from None
                if length is None:
                    length = len(aligned_sequence)
                elif length != len(aligned_sequence):
                    raise ValueError(f"Sequence {seqname} length is inconsistent")
                aligned_sequence = aligned_sequence.replace(".", "-")
                sequence = aligned_sequence.replace("-", "")
                aligned_sequences.append(aligned_sequence)
                name, start, end = self._identifier_split(seqname)
                if start is None:
                    seq = Seq(sequence)
                else:
                    if start < end:
                        strands.append("+")
                    else:
                        start, end = end, start
                        strands.append("-")
                        sequence = reverse_complement(sequence)
                    start -= 1  # 0-based index
                    if start + len(sequence) != end:
                        raise ValueError(f"Start and end of sequence {name} are not consistent with sequence length")
                    seq = Seq({start: sequence}, end)
                record = SeqRecord(seq, id=name)
                records.append(record)
            elif line.startswith("#=GF "):
                # Generic per-File annotation, free text
                # Format: #=GF <feature> <free text>
                feature, text = line[5:].strip().split(None, 1)
                # Each feature key could be used more than once,
                # so store the entries as a list of strings.
                gf[feature].append(text)
            elif line.startswith("#=GC "):
                # Generic per-Column annotation, exactly 1 char per column
                # Format: "#=GC <feature> <exactly 1 char per column>"
                feature, text = line[5:].strip().split(None, 2)
                if feature not in gc:
                    gc[feature] = ""
                gc[feature] += text.strip()  # append to any previous entry
                # Might be interleaved blocks, so can't check length yet
            elif line.startswith("#=GS "):
                # Generic per-Sequence annotation, free text
                # Format: "#=GS <seqname> <feature> <free text>"
                try:
                    seqname, feature, text = line[5:].strip().split(None, 2)
                except ValueError:
                    # Free text can sometimes be empty, which a one line split throws an error for.
                    # See https://github.com/biopython/biopython/issues/2982 for more details
                    seqname, feature = line[5:].strip().split(None, 1)
                    text = ""
                gs.setdefault(seqname, {})
                gs[seqname].setdefault(feature, [])
                gs[seqname][feature].append(text)
            elif line[:5] == "#=GR ":
                # Generic per-Sequence AND per-Column markup
                # Format: "#=GR <seqname> <feature> <exactly 1 char per column>"
                seqname, feature, text = line[5:].strip().split(None, 2)
                gr.setdefault(seqname, {})
                gr[seqname].setdefault(feature, "")
                gr[seqname][feature] += text.strip()  # append to any previous entry
                # Might be interleaved blocks, so can't check length yet
        yield _create_alignment(aligned_sequences, records, strands, gc, gf, gr)

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
