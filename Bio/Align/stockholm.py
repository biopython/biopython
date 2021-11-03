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
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for PFAM alignment files in the Stockholm format.

    The file may contain multiple concatenated alignments, which are loaded
    and returned incrementally.

    This parser will detect if the Stockholm file follows the PFAM
    conventions for sequence specific meta-data (lines starting #=GS
    and #=GR) and populates the SeqRecord fields accordingly.

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
    """

    gf_mapping = {
        "AC": "accession",
        "ID": "identification",
        "DE": "definition",
        "AU": "author",
        "SE": "source of seed",
        "SS": "source of structure",
        "GA": "gathering method",
        "BM": "build method",
        "SM": "search method",
        "TC": "trusted cutoff",
        "TP": "type",
        "NC": "noise cutoff",
        "PI": "previous identifier",
        "CC": "comment",
        "CL": "clan",
        "WK": "wikipedia",
        "CB": "calibration method",
        "**": "**",  # Found in Rfam
    }

    gr_mapping = {
        "SS": "secondary_structure",
        "PP": "posterior_probability",
        "CSA": "Catalytic Site Atlas",  # used in CATH
        # These features are included in the Stockholm file format
        # documentation, but currently not used in the PFAM, RFAM, and CATH
        # databases:
        "SA": "surface_accessibility",
        "TM": "transmembrane",
        "LI": "ligand_binding",
        "AS": "active_site",
        "pAS": "active site - Pfam predicted",
        "sAS": "active site - from SwissProt",
        "IN": "intron",
    }

    gc_mapping = {"RF": "reference_coordinate_annotation",
                  "seq_cons": "consensus_sequence",
                  "scorecons": "consensus_score",  # used in CATH
                  "scorecons_70": "consensus_score_70",  # used in CATH
                  "scorecons_80": "consensus_score_80",  # used in CATH
                  "scorecons_90": "consensus_score_90",  # used in CATH
                  # This feature is included in the Stockholm file format
                  # documentation, but currently not used in the PFAM, RFAM,
                  # and CATH databases:
                  "MM": "model_mask",
                 }
    # Add *_cons from GR mapping:
    for key, value in gr_mapping.items():
        gc_mapping[key + "_cons"] = "consensus_" + value

    # These GC keywords are used in Rfam:
    for keyword in ("RNA_elements",
                    "RNA_structural_element",
                    "RNA_structural_elements",
                    "RNA_ligand_AdoCbl",
                    "RNA_ligand_AqCbl",
                    "RNA_ligand_FMN",
                    "RNA_ligand_Guanidinium",
                    "RNA_ligand_SAM",
                    "RNA_ligand_THF_1",
                    "RNA_ligand_THF_2",
                    "RNA_ligand_TPP",
                    "RNA_ligand_preQ1",
                    "RNA_motif_k_turn",
                    "Repeat_unit",
                    "2L3J_B_SS",
                    "CORE",
                    "PK",
                    "PK_SS",
                    "cons",
                   ):
        gc_mapping[keyword] = keyword
    gs_mapping = {"AC": "accession",
                  # "DE": description,  # handled separately
                  # "DR": "database_references",  # handled separately
                  "OS": "organism",
                  # These two features are included in the Stockholm file
                  # format documentation, but currently not used in the PFAM,
                  # RFAM, and CATH databases:
                  "OC": "organism_classification",
                  "LO": "look",
                 }

    def __init__(self, source):
        """Create an AlignmentIterator object.

        Arguments:
         - source   - input data or file name

        """
        super().__init__(source, mode="t", fmt="Stockholm")

    def _add_per_file_annotations(self, alignment, gf, rows):
        for key, value in gf.items():
            if key == "WK":
                lines = iter(value)
                references = []
                for line in lines:
                    reference = ""
                    while line.endswith("/"):
                        reference += line[:-1]
                        line = next(lines)
                    reference += line
                    references.append(reference)
                value = references
            elif key in ("SM", "CC", "**"):
                value = " ".join(value)
            elif key == "SQ":
                assert len(value) == 1
                if int(value.pop()) != rows:
                    raise ValueError("Inconsistent number of sequences in alignment")
                continue
            elif key == "AU":
                pass
            else:
                assert len(value) == 1, (key, value)
                value = value.pop()
            alignment.annotations[self.gf_mapping[key]] = value

    def _add_per_column_annotations(self, alignment, gc, columns, skipped_columns):
        if gc:
            alignment.column_annotations = {}
            for key, value in gc.items():
                if skipped_columns:
                    value = "".join(letter for index, letter in enumerate(value) if index not in skipped_columns)
                if len(value) != columns:
                    raise ValueError(f"{key} length is {len(value)}, expected {columns}")
                alignment.column_annotations[self.gc_mapping[key]] = value

    def _add_per_sequence_annotations(self, alignment, gs):
        for seqname, annotations in gs.items():
            for record in alignment.sequences:
                if record.id == seqname:
                    break
            else:
                raise ValueError(f"Failed to find seqname {seqname}")
            for key, value in annotations.items():
                if key == "DE":
                    record.description = value
                elif key == "DR":
                    record.dbxrefs = value
                else:
                    key = self.gs_mapping[key]
                    record.annotations[key] = value

    def _add_per_sequence_and_per_column_annotations(self, alignment, gr):
        for seqname, letter_annotations in gr.items():
            for record in alignment.sequences:
                if record.id == seqname:
                    break
            else:
                raise ValueError(f"Failed to find seqname {seqname}")
            for keyword, letter_annotation in letter_annotations.items():
                feature = self.gr_mapping[keyword]
                if keyword == "CSA":
                    letter_annotation = letter_annotation.replace("-", "")
                else:
                    letter_annotation = letter_annotation.replace(".", "")
                record.letter_annotations[feature] = letter_annotation

    def parse(self, stream):
        """Parse the next alignment from the stream."""
        if stream is None:
            raise StopIteration

        for line in stream:
            line = line.strip()
            if not line:
                continue
            elif line == "# STOCKHOLM 1.0":
                # Starting a new alignment
                records = []
                aligned_sequences = []
                references = []
                reference_comments = []
                database_references = []
                nested_domains = []
                gf = defaultdict(list)
                gc = {}
                gs = defaultdict(lambda: {"DR": []})
                gr = defaultdict(dict)
                length = None
            elif line == "//":
                # The "//" line indicates the end of the alignment.
                # There may still be more meta-data
                skipped_columns = []
                coordinates = Alignment.infer_coordinates(aligned_sequences, skipped_columns)
                skipped_columns = set(skipped_columns)
                alignment = Alignment(records, coordinates)
                alignment.annotations = {}
                if references:
                    alignment.annotations["references"] = []
                    for reference in references:
                        reference = dict(reference)
                        reference["title"] = " ".join(reference["title"])
                        reference["author"] = " ".join(reference["author"])
                        reference["location"] = " ".join(reference["location"])
                        alignment.annotations["references"].append(reference)
                if database_references:
                    alignment.annotations["database_references"] = database_references
                if nested_domains:
                    alignment.annotations["nested_domains"] = nested_domains
                rows, columns = alignment.shape
                self._add_per_file_annotations(alignment, gf, rows)
                self._add_per_column_annotations(alignment, gc, columns, skipped_columns)
                self._add_per_sequence_annotations(alignment, gs)
                self._add_per_sequence_and_per_column_annotations(alignment, gr)
                yield alignment
            elif not line.startswith("#"):
                # Sequence
                # Format: "<seqname> <sequence>"
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
                    raise ValueError(f"Aligned sequence {seqname} consists of {len(aligned_sequence)} letters, expected {length} letters)")
                aligned_sequence = aligned_sequence.replace(".", "-")
                sequence = aligned_sequence.replace("-", "")
                aligned_sequences.append(aligned_sequence)
                seq = Seq(sequence)
                record = SeqRecord(seq, id=seqname)
                records.append(record)
            elif line.startswith("#=GF "):
                # Generic per-File annotation, free text
                # Format: #=GF <feature> <free text>
                feature, text = line[5:].strip().split(None, 1)
                if feature == "RN":
                    assert text.startswith("[")
                    assert text.endswith("]")
                    number = int(text[1:-1])
                    reference = defaultdict(list)
                    reference["number"] = number
                    if reference_comments:
                        reference["comment"] = " ".join(reference_comments)
                        reference_comments = []
                    references.append(reference)
                elif feature == "RM":
                    assert not reference["medline"]
                    reference["medline"] = text
                elif feature == "RT":
                    reference["title"].append(text)
                elif feature == "RA":
                    reference["author"].append(text)
                elif feature == "RL":
                    reference["location"].append(text)
                elif feature == "RC":
                    reference_comments.append(text)
                elif feature == "DR":
                    database_reference = {"reference": text}
                    database_references.append(database_reference)
                elif feature == "DC":
                    assert "comment" not in database_reference
                    database_reference["comment"] = text
                elif feature == "NE":
                    nested_domain = {"accession": text}
                    nested_domains.append(nested_domain)
                elif feature == "NL":
                    assert "location" not in nested_domain
                    nested_domain["location"]  = text
                else:
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
                if feature == "DR":
                    gs[seqname][feature].append(text)
                else:
                    assert feature not in gs[seqname]
                    gs[seqname][feature] = text
            elif line[:5] == "#=GR ":
                # Generic per-Sequence AND per-Column markup
                # Format: "#=GR <seqname> <feature> <exactly 1 char per column>"
                terms = line[5:].split(None, 2)
                assert terms[0] == seqname
                feature = terms[1]
                gr[seqname][feature] = terms[2].strip()


class AlignmentWriter(interfaces.AlignmentWriter):
    """Alignment file writer for the Stockholm file format."""

    gr_mapping = {value: key for key, value in AlignmentIterator.gr_mapping.items()}
    gc_mapping = {value: key for key, value in AlignmentIterator.gc_mapping.items()}
    gs_mapping = {value: key for key, value in AlignmentIterator.gs_mapping.items()}

    def write_alignment(self, alignment):
        """Use this to write the alignment to an open file.

        Note that sequences and their annotation are recorded
        together (rather than having a block of annotation followed
        by a block of aligned sequences).
        """
        stream = self.stream
        count = len(alignment)

        n, self._length_of_sequences = alignment.shape
        self._ids_written = []

        if count == 0:
            raise ValueError("Must have at least one sequence")
        if self._length_of_sequences == 0:
            raise ValueError("Non-empty sequences are required")

        stream.write("# STOCKHOLM 1.0\n")
        stream.write("#=GF SQ %i\n" % count)
        for record in alignment:
            self._write_record(record)
        # This shouldn't be None... but just in case,
        if alignment.column_annotations:
            for k, v in sorted(alignment.column_annotations.items()):
                if k in self.gc_mapping:
                    stream.write(f"#=GC {self.gc_mapping[k]} {v}\n")
                elif k in self.gr_mapping:
                    stream.write(f"#=GC {self.gr_mapping[k]}_cons {v}\n")
                else:
                    # It doesn't follow the PFAM standards, but should we record
                    # this data anyway?
                    pass
        stream.write("//\n")

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
        stream.write(f"{seq_name} {record.seq}\n")

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
            stream.write(
                f"#=GS {seq_name} AC {self.clean(record.annotations['accession'])}\n"
            )
        elif record.id:
            stream.write(f"#=GS {seq_name} AC {self.clean(record.id)}\n")

        # DE = description
        if record.description:
            stream.write(f"#=GS {seq_name} DE {self.clean(record.description)}\n")

        # DE = database links
        for xref in record.dbxrefs:
            stream.write(f"#=GS {seq_name} DR {self.clean(xref)}\n")

        # GS = other per sequence annotation
        for key, value in record.annotations.items():
            if key in self.pfam_gs_mapping:
                data = self.clean(str(value))
                if data:
                    stream.write(
                        "#=GS %s %s %s\n"
                        % (seq_name, self.clean(self.pfam_gs_mapping[key]), data)
                    )
            else:
                # It doesn't follow the PFAM standards, but should we record
                # this data anyway?
                pass

        # GR = per row per column sequence annotation
        for key, value in record.letter_annotations.items():
            if key in self.gr_mapping and len(str(value)) == len(record.seq):
                data = self.clean(str(value))
                if data:
                    stream.write(
                        "#=GR %s %s %s\n"
                        % (seq_name, self.clean(self.gr_mapping[key]), data)
                    )
            else:
                # It doesn't follow the PFAM standards, but should we record
                # this data anyway?
                pass


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
