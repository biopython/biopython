# Copyright 2006-2016 by Peter Cock.  All rights reserved.
# Copyright 2021 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.Align support for alignment files in the Stockholm file format.

You are expected to use this module via the Bio.Align functions.

For example, consider this alignment from PFAM for the HAT helix motif::

    # STOCKHOLM 1.0
    #=GF ID   HAT
    #=GF AC   PF02184.18
    #=GF DE   HAT (Half-A-TPR) repeat
    #=GF AU   SMART;
    #=GF SE   Alignment kindly provided by SMART
    #=GF GA   21.00 21.00;
    #=GF TC   21.00 21.00;
    #=GF NC   20.90 20.90;
    #=GF BM   hmmbuild HMM.ann SEED.ann
    #=GF SM   hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq
    #=GF TP   Repeat
    #=GF CL   CL0020
    #=GF RN   [1]
    #=GF RM   9478129
    #=GF RT   The HAT helix, a repetitive motif implicated in RNA processing.
    #=GF RA   Preker PJ, Keller W;
    #=GF RL   Trends Biochem Sci 1998;23:15-16.
    #=GF DR   INTERPRO; IPR003107;
    #=GF DR   SMART; HAT;
    #=GF DR   SO; 0001068; polypeptide_repeat;
    #=GF CC   The HAT (Half A TPR) repeat is found in several RNA processing
    #=GF CC   proteins [1].
    #=GF SQ   3
    #=GS CRN_DROME/191-222     AC P17886.2
    #=GS CLF1_SCHPO/185-216    AC P87312.1
    #=GS CLF1_SCHPO/185-216    DR PDB; 3JB9 R; 185-216;
    #=GS O16376_CAEEL/201-233  AC O16376.2
    CRN_DROME/191-222                KEIDRAREIYERFVYVH.PDVKNWIKFARFEES
    CLF1_SCHPO/185-216               HENERARGIYERFVVVH.PEVTNWLRWARFEEE
    #=GR CLF1_SCHPO/185-216    SS    --HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH
    O16376_CAEEL/201-233             KEIDRARSVYQRFLHVHGINVQNWIKYAKFEER
    #=GC SS_cons                     --HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH
    #=GC seq_cons                    KEIDRARuIYERFVaVH.P-VpNWIKaARFEEc
    //

Parsing this file using Bio.Align stores the alignment, its annotations, as
well as the sequences and their annotations::

    >>> from Bio.Align import stockholm
    >>> alignments = stockholm.AlignmentIterator("Stockholm/example.sth")
    >>> alignment = next(alignments)
    >>> alignment.shape
    (3, 33)
    >>> alignment[0]
    'KEIDRAREIYERFVYVH-PDVKNWIKFARFEES'

Alignment meta-data are stored in alignment.annotations::

    >>> alignment.annotations["accession"]
    'PF02184.18'
    >>> alignment.annotations["references"][0]["title"]
    'The HAT helix, a repetitive motif implicated in RNA processing.'

Annotations of alignment columns are stored in alignment.column_annotations::

    >>> alignment.column_annotations["consensus secondary structure"]
    '--HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH'

Sequences and their annotations are stored in alignment.sequences::

   >>> alignment.sequences[0].id
   'CRN_DROME/191-222'
   >>> alignment.sequences[0].seq
   Seq('KEIDRAREIYERFVYVHPDVKNWIKFARFEES')
   >>> alignment.sequences[1].letter_annotations["secondary structure"]
   '--HHHHHHHHHHHHHHS--HHHHHHHHHHHHH'

Slicing specific columns of an alignment will slice any per-column-annotations:

    >>> alignment.column_annotations["consensus secondary structure"]
    '--HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH'
    >>> part_alignment = alignment[:,10:20]
    >>> part_alignment.column_annotations["consensus secondary structure"]
    'HHHHHHS.--'
"""
import textwrap
from collections import defaultdict

from Bio.Align import Alignment
from Bio.Align import interfaces
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class AlignmentIterator(interfaces.AlignmentIterator):
    """Alignment iterator for alignment files in the Stockholm format.

    The file may contain multiple concatenated alignments, which are loaded
    and returned incrementally.

    Alignment meta-data (lines starting with #=GF) are stored in the dictionary
    alignment.annotations. Column annotations (lines starting with #=GC) are
    stored in the dictionary alignment.column_annotations. Sequence names are
    stored in record.id. Sequence record meta-data (lines starting with #=GS)
    are stored in the dictionary record.annotations. Sequence letter
    annotations (lines starting with #=GR) are stored in the dictionary
    record.letter_annotations.

    Wrap-around alignments are not supported - each sequence must be on
    a single line.

    For more information on the file format, please see:
    http://sonnhammer.sbc.su.se/Stockholm.html
    https://en.wikipedia.org/wiki/Stockholm_format
    """

    fmt = "Stockholm"

    gf_mapping = {
        "ID": "identifier",
        "AC": "accession",
        "DE": "definition",
        "AU": "author",
        "SE": "source of seed",
        "SS": "source of structure",
        "GA": "gathering method",
        "TC": "trusted cutoff",
        "NC": "noise cutoff",
        "BM": "build method",
        "SM": "search method",
        "TP": "type",
        "PI": "previous identifier",
        "CC": "comment",
        "CL": "clan",
        "WK": "wikipedia",
        "CB": "calibration method",
        "**": "**",  # Found in Rfam
    }

    gr_mapping = {
        "SS": "secondary structure",
        "PP": "posterior probability",
        "CSA": "Catalytic Site Atlas",  # used in CATH
        # These features are included in the Stockholm file format
        # documentation, but currently not used in the PFAM, RFAM, and CATH
        # databases:
        "SA": "surface accessibility",
        "TM": "transmembrane",
        "LI": "ligand binding",
        "AS": "active site",
        "pAS": "active site - Pfam predicted",
        "sAS": "active site - from SwissProt",
        "IN": "intron",
    }

    gc_mapping = {
        "RF": "reference coordinate annotation",
        "seq_cons": "consensus sequence",
        "scorecons": "consensus score",  # used in CATH
        "scorecons_70": "consensus score 70",  # used in CATH
        "scorecons_80": "consensus score 80",  # used in CATH
        "scorecons_90": "consensus score 90",  # used in CATH
        # This feature is included in the Stockholm file format
        # documentation, but currently not used in the PFAM, RFAM,
        # and CATH databases:
        "MM": "model mask",
    }
    # Add *_cons from GR mapping:
    for key, value in gr_mapping.items():
        gc_mapping[key + "_cons"] = "consensus " + value

    # These GC keywords are used in Rfam:
    for keyword in (
        "RNA_elements",
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
        gc_mapping[keyword] = keyword.replace("_", " ")
    gs_mapping = {
        "AC": "accession",
        # "DE": description,  # handled separately
        # "DR": "database_references",  # handled separately
        "OS": "organism",
        # These two features are included in the Stockholm file
        # format documentation, but currently not used in the PFAM,
        # RFAM, and CATH databases:
        "OC": "organism classification",
        "LO": "look",
    }

    @staticmethod
    def _store_per_file_annotations(alignment, gf, rows):
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
            try:
                alignment.annotations[AlignmentIterator.gf_mapping[key]] = value
            except KeyError:
                pass

    @staticmethod
    def _store_per_column_annotations(alignment, gc, columns, skipped_columns):
        if gc:
            alignment.column_annotations = {}
            for key, value in gc.items():
                if skipped_columns:
                    value = "".join(
                        letter
                        for index, letter in enumerate(value)
                        if index not in skipped_columns
                    )
                if len(value) != columns:
                    raise ValueError(
                        f"{key} length is {len(value)}, expected {columns}"
                    )
                alignment.column_annotations[
                    AlignmentIterator.gc_mapping.get(key, key)
                ] = value

    @staticmethod
    def _store_per_sequence_annotations(alignment, gs):
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
                    record.annotations[
                        AlignmentIterator.gs_mapping.get(key, key)
                    ] = value

    @staticmethod
    def _store_per_sequence_and_per_column_annotations(alignment, gr):
        for seqname, letter_annotations in gr.items():
            for record in alignment.sequences:
                if record.id == seqname:
                    break
            else:
                raise ValueError(f"Failed to find seqname {seqname}")
            for key, letter_annotation in letter_annotations.items():
                feature = AlignmentIterator.gr_mapping.get(key, key)
                if key == "CSA":
                    letter_annotation = letter_annotation.replace("-", "")
                else:
                    letter_annotation = letter_annotation.replace(".", "")
                record.letter_annotations[feature] = letter_annotation

    def _read_next_alignment(self, stream):
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
                # Reached the end of the alignment.
                skipped_columns = []
                coordinates = Alignment.infer_coordinates(
                    aligned_sequences, skipped_columns
                )
                skipped_columns = set(skipped_columns)
                alignment = Alignment(records, coordinates)
                for index in sorted(skipped_columns, reverse=True):
                    del operations[index]  # noqa: F821
                alignment.operations = operations  # noqa: F821
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
                    alignment.annotations["database references"] = database_references
                if nested_domains:
                    alignment.annotations["nested domains"] = nested_domains
                rows, columns = alignment.shape
                AlignmentIterator._store_per_file_annotations(alignment, gf, rows)
                AlignmentIterator._store_per_column_annotations(
                    alignment, gc, columns, skipped_columns
                )
                AlignmentIterator._store_per_sequence_annotations(alignment, gs)
                AlignmentIterator._store_per_sequence_and_per_column_annotations(
                    alignment, gr
                )
                return alignment
            elif not line.startswith("#"):
                # Sequence
                # Format: "<seqname> <sequence>"
                try:
                    seqname, aligned_sequence = line.split(None, 1)
                except ValueError:
                    # This might be someone attempting to store a zero length sequence?
                    raise ValueError(
                        "Could not split line into sequence name and aligned sequence:\n"
                        + line
                    ) from None
                if length is None:
                    length = len(aligned_sequence)
                    operations = bytearray(b"M" * length)
                elif length != len(aligned_sequence):
                    raise ValueError(
                        f"Aligned sequence {seqname} consists of {len(aligned_sequence)} letters, expected {length} letters)"
                    )
                for i, letter in enumerate(aligned_sequence):
                    if letter == "-":
                        assert operations[i] != ord("I")
                        operations[i] = ord("D")  # deletion
                    elif letter == ".":
                        assert operations[i] != ord("D")
                        operations[i] = ord("I")  # insertion
                aligned_sequence = aligned_sequence.replace(".", "-")
                sequence = aligned_sequence.replace("-", "")
                aligned_sequences.append(aligned_sequence)
                seq = Seq(sequence)
                record = SeqRecord(seq, id=seqname, description="")
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
                    nested_domain["location"] = text
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

    gf_mapping = {value: key for key, value in AlignmentIterator.gf_mapping.items()}
    gs_mapping = {value: key for key, value in AlignmentIterator.gs_mapping.items()}
    gr_mapping = {value: key for key, value in AlignmentIterator.gr_mapping.items()}
    gc_mapping = {value: key for key, value in AlignmentIterator.gc_mapping.items()}

    fmt = "Stockholm"

    def format_alignment(self, alignment):
        """Return a string with a single alignment in the Stockholm format."""
        rows, columns = alignment.shape

        if rows == 0:
            raise ValueError("Must have at least one sequence")
        if columns == 0:
            raise ValueError("Non-empty sequences are required")

        try:
            alignment_annotations = alignment.annotations
        except AttributeError:
            alignment_annotations = {}
        lines = []
        lines.append("# STOCKHOLM 1.0\n")
        # #=GF Above the alignment; alignment.annotations
        for key, feature in self.gf_mapping.items():
            if key == "comment":
                # write this last
                continue
            value = alignment_annotations.get(key)
            if value is not None:
                feature = self.gf_mapping[key]
                if key in ("author", "wikipedia"):
                    for item in value:
                        lines.append(f"#=GF {feature}   {item}\n")
                else:
                    lines.append(f"#=GF {feature}   {value}\n")
        nested_domains = alignment_annotations.get("nested domains")
        if nested_domains is not None:
            for nested_domain in nested_domains:
                accession = nested_domain.get("accession")
                if accession is not None:
                    lines.append(f"#=GF NE   {accession}\n")
                location = nested_domain.get("location")
                if location is not None:
                    lines.append(f"#=GF NL   {location}\n")
        references = alignment_annotations.get("references")
        if references is not None:
            for reference in references:
                comment = reference.get("comment")
                lines.append(AlignmentWriter._format_long_text("#=GF RC   ", comment))
                lines.append(f"#=GF RN   [{reference['number']}]\n")
                lines.append(f"#=GF RM   {reference['medline']}\n")
                title = reference["title"]
                lines.append(AlignmentWriter._format_long_text("#=GF RT   ", title))
                lines.append(f"#=GF RA   {reference['author']}\n")
                lines.append(f"#=GF RL   {reference['location']}\n")
        database_references = alignment_annotations.get("database references")
        if database_references is not None:
            for database_reference in database_references:
                lines.append(f"#=GF DR   {database_reference['reference']}\n")
                comment = database_reference.get("comment")
                if comment is not None:
                    lines.append(f"#=GF DC   {comment}\n")
        key = "comment"
        value = alignment_annotations.get(key)
        if value is not None:
            prefix = "#=GF %s   " % self.gf_mapping[key]
            lines.append(AlignmentWriter._format_long_text(prefix, value))
        for key in alignment_annotations:
            if key in self.gf_mapping:
                continue
            if key == "nested domains":
                continue
            if key == "references":
                continue
            if key == "database references":
                continue
            raise ValueError(
                "Unknown annotation %s found in alignment.annotations" % key
            )
        lines.append("#=GF SQ   %i\n" % rows)
        # #=GS Above the alignment or just below the corresponding sequence;
        #    record.annotations
        # #=GR Just below the corresponding sequence;
        #    record.letter_annotations
        width = max(len(record.id) for record in alignment.sequences)
        start = max(width, 20) + 12
        for record in alignment.sequences:
            name = record.id.ljust(width)
            for key, value in record.annotations.items():
                feature = self.gs_mapping.get(key, key)
                lines.append(f"#=GS {name}  {feature} {value}\n")
            if record.description:
                lines.append(f"#=GS {name}  DE {record.description}\n")
            for value in record.dbxrefs:
                lines.append(f"#=GS {name}  DR {value}\n")

        try:
            operations = alignment.operations
        except AttributeError:
            operations = bytes(b"M" * columns)
        else:
            assert len(operations) == columns
        for aligned_sequence, record in zip(alignment, alignment.sequences):
            aligned_sequence = "".join(
                "." if letter == "-" and operation == ord("I") else letter
                for operation, letter in zip(operations, aligned_sequence)
            )
            lines.extend(
                AlignmentWriter._format_record(width, start, aligned_sequence, record)
            )
        # #=GC Below the alignment;
        #    alignment.column_annotations
        if alignment.column_annotations:
            for key, value in alignment.column_annotations.items():
                feature = self.gc_mapping.get(key, key)
                line = f"#=GC {feature} ".ljust(start) + value + "\n"
                lines.append(line)
        lines.append("//\n")
        return "".join(lines)

    @staticmethod
    def _format_long_text(prefix, text):
        """Format the text as wrapped lines (PRIVATE)."""
        if text is None:
            return ""
        return (
            textwrap.fill(
                text,
                width=79,
                break_long_words=False,
                initial_indent=prefix,
                subsequent_indent=prefix,
            )
            + "\n"
        )

    @staticmethod
    def _format_record(width, start, aligned_sequence, record):
        """Format lines for a single SeqRecord (PRIVATE)."""
        name = record.id.ljust(start)
        line = name + aligned_sequence + "\n"
        yield line
        indices = [
            index for index, letter in enumerate(aligned_sequence) if letter in ".-"
        ]
        indices.reverse()
        name = record.id.ljust(width)
        for key, value in record.letter_annotations.items():
            feature = AlignmentWriter.gr_mapping.get(key, key)
            j = 0
            values = bytearray(b"." * len(aligned_sequence))
            for i, letter in enumerate(aligned_sequence):
                if letter not in ".-":
                    values[i] = ord(value[j])
                    j += 1
            value = values.decode()
            line = f"#=GR {name}  {feature} ".ljust(start) + value + "\n"
            yield line


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
