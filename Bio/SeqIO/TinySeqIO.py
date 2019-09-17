# Copyright 2019 by Peter Cock. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "tinyseq" file format, NCBI Tiny Sequence XML.

This is a very simple XML format, expanding on the FASTA format with
explicit fields for accession, organism, sequence type, in addition to
an essentially free text description.

See also: https://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd
"""

from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import SequentialSequenceWriter


class TinySeqWriter(SequentialSequenceWriter):
    """Write SeqRecord objects as NCBI TinySeq XML."""

    def write_header(self):
        """Write the XML header including opening <TSeqSet> tag."""
        SequentialSequenceWriter.write_header(self)
        self.handle.write(
            '<?xml version="1.0"?>\n'
            '<!DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd">\n'
            "<TSeqSet>\n"
        )

    def write_record(self, record):
        """Write one record from <TSeq> to </TSeq>."""
        # TODO - Use Python standard library to handle encoding values nicely...

        # Sequence type - mandatory
        alpha = Alphabet._get_base_alphabet(record.seq.alphabet)
        if isinstance(alpha, Alphabet.NucleotideAlphabet):
            seq_type = "nucleotide"
        elif isinstance(alpha, Alphabet.ProteinAlphabet):
            seq_type = "protein"
        else:
            raise ValueError("Need a Nucleotide or Protein alphabet")

        # Identifiers - several options here
        gi_tag = ""  # <TSeq_gi>11321596</TSeq_gi>\n
        accver_tag = ""  # <TSeq_accver>NM_002253.1</TSeq_accver>\n
        sid_tag = ""  # <TSeq_sid>ref|NM_002253.1|</TSeq_sid>\n
        local_tag = (
            "  <TSeq_local>%s</TSeq_local>\n" % record.id
            if record.id and record.id != "<unknown id>"
            else ""
        )

        # Taxonomy - following is following GenBank/EMBL parser output
        org = record.annotations.get("organism", None)
        taxid = None
        if record.features and record.features[0].type == "source":
            source = record.features[0]
            try:
                # This is usally the scientific name, while the GenBank/EMBL
                # header often adds an informal name in brackets
                org = source.qualifiers["organism"][0]
            except KeyError:
                pass
            for ref in source.qualifiers.get("db_xref", []):
                if ref.startswith("taxon:"):
                    try:
                        taxid = int(ref[6:])
                    except ValueError:
                        pass
        org_tag = "  <TSeq_orgname>%s</TSeq_orgname>\n" % org if org else ""
        taxid_tag = "  <TSeq_taxid>%i</TSeq_taxid>\n" % taxid if taxid else ""

        # Defline - mandatory
        if "\n" in record.description:
            raise ValueError("Not sure if should allow new lines in TinySeq defline...")
        def_tag = "  <TSeq_defline>%s</TSeq_defline>\n" % (
            record.description
            if record.description and record.description != "<unknown description>"
            else ""
        )

        self.handle.write(
            "<TSeq>\n"
            '  <TSeq_seqtype value="%s"/>\n%s%s%s%s%s%s'
            "  <TSeq_length>%i</TSeq_length>\n"
            "  <TSeq_sequence>%s</TSeq_sequence>\n"
            "</TSeq>\n"
            % (
                seq_type,
                gi_tag,
                sid_tag,
                local_tag,
                taxid_tag,
                org_tag,
                def_tag,
                len(record),
                record.seq,
            )
        )

    def write_footer(self):
        """Write the closing </TSeqSet> tag."""
        SequentialSequenceWriter.write_footer(self)
        self.handle.write("</TSeqSet>\n")


def _parse_value(text):
    """Extract key:value from '<TSeq_key>value</TSeq_key>' (PRIVATE)."""
    assert text.startswith("<TSeq_")
    start = text.split(">", 1)[0] + ">"
    end = "</" + start[1:]
    assert text.endswith(end), "%r.endswith(%r)" % (text, end)
    return start[6:-1], text[len(start) : -len(end)]


def TinySeqIterator(handle):
    """Iterate over TinySeq XML as SeqRecord objects.

    Niave line-centric parser, will not cope with all valid XML!
    """
    line = handle.readline()
    if not line:
        raise ValueError("Empty file.")
    if line == '<?xml version="1.0"?>\n':
        # Good
        line = handle.readline()
        if not line:
            raise ValueError("Premature end of file")
    elif line.strip() == "<TSeqSet>":
        # OK
        pass
    else:
        raise ValueError(
            'TinySeq XML files should start <?xml version="1.0"?>, or <TSeqSet>, but not: %r'
            % line
        )

    if line.startswith("<!DOCTYPE"):
        if line not in (
            '<!DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "https://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd">\n',
            '<!DOCTYPE TSeqSet PUBLIC "-//NCBI//NCBI TSeq/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_TSeq.dtd">\n',
        ):
            raise ValueError("Unexpected DOCTYPE for TinySeq: %r" % line)
        line = handle.readline()
    if line.strip() != "<TSeqSet>":
        raise ValueError("Missing TinySeq opening <TSeqSet> tag: %r" % line)

    while True:
        line = handle.readline()
        while line and not line.strip():
            line = handle.readline()
        if not line:
            raise ValueError("Premature end of file mid-record")
        line = line.strip()
        if line == "</TSeqSet>":
            break
        if line != "<TSeq>":
            raise ValueError(
                "TinySeq records should start with <TSeq> line, not: %r" % line
            )
        local = sid = accver = defline = taxid = org = seq = alpha = length = None
        while True:
            line = handle.readline().strip()
            if line == '<TSeq_seqtype value="nucleotide"/>':
                alpha = Alphabet.generic_nucleotide
            elif line == '<TSeq_seqtype value="protein"/>':
                alpha = Alphabet.generic_protein
            elif line.startswith("<TSeq_"):
                key, value = _parse_value(line)
                if key == "local":
                    local = value
                elif key == "sid":
                    sid = value
                elif key == "accver":
                    accver = value
                elif key == "orgname":
                    org = value
                elif key == "taxid":
                    taxid = int(value)
                elif key == "defline":
                    defline = value
                elif key == "sequence":
                    seq = value.strip()
                elif key == "length":
                    length = int(value)
                elif key == "gi":
                    # NCBI are deprecating this
                    pass
                else:
                    raise ValueError("Unexpected attribute %s from: %r" % (key, line))
            elif line == "</TSeq>":
                break
            else:
                raise ValueError("Unexpected line: %r" % line)
        if alpha is None or length is None or defline is None or seq is None:
            raise ValueError("Missing mandatory values in TinySeq record")
        if length != len(seq):
            raise ValueError(
                "Mismatch, declared length %i versus actual %i" % (length, len(seq))
            )
        record = SeqRecord(
            Seq(seq, alpha), id=local or sid or accver, description=defline
        )
        record.name = record.id  # better choice?
        # How to record all the identification attributes?
        if org:
            record.annotations["organism"] = org
        if taxid:
            # Follow key used for BioSQL, SwissProt and SeqXml
            record.annotations["ncbi_taxid"] = taxid
        yield record

    line = handle.readline()
    while line:
        if line.strip():
            raise ValueError("Unexpected text after </TSeqSet> tag: %r" % line)
        line = handle.readline()
