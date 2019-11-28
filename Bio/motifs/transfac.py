# Copyright 2003 by Bartek Wilczynski.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Parsing TRANSFAC files."""


from Bio import motifs


class Motif(motifs.Motif, dict):
    """Store the information for one TRANSFAC motif.

    This class inherits from the Bio.motifs.Motif base class, as well
    as from a Python dictionary. All motif information found by the parser
    is stored as attributes of the base class when possible; see the
    Bio.motifs.Motif base class for a description of these attributes. All
    other information associated with the motif is stored as (key, value)
    pairs in the dictionary, where the key is the two-letter fields as found
    in the TRANSFAC file. References are an exception: These are stored in
    the .references attribute.

    These fields are commonly found in TRANSFAC files::

        AC:    Accession number
        AS:    Accession numbers, secondary
        BA:    Statistical basis
        BF:    Binding factors
        BS:    Factor binding sites underlying the matrix
               [sequence; SITE accession number; start position for matrix
               sequence; length of sequence used; number of gaps inserted;
               strand orientation.]
        CC:    Comments
        CO:    Copyright notice
        DE:    Short factor description
        DR:    External databases
               [database name: database accession number]
        DT:    Date created/updated
        HC:    Subfamilies
        HP:    Superfamilies
        ID:    Identifier
        NA:    Name of the binding factor
        OC:    Taxonomic classification
        OS:    Species/Taxon
        OV:    Older version
        PV:    Preferred version
        TY:    Type
        XX:    Empty line; these are not stored in the Record.

    References are stored in an .references attribute, which is a list of
    dictionaries with the following keys::

        RN:    Reference number
        RA:    Reference authors
        RL:    Reference data
        RT:    Reference title
        RX:    PubMed ID

    For more information, see the TRANSFAC documentation.
    """

    multiple_value_keys = {"BF", "OV", "HP", "BS", "HC", "DT", "DR"}
    # These keys can occur multiple times for one motif

    reference_keys = {"RX", "RA", "RT", "RL"}
    # These keys occur for references


class Record(list):
    """Store the information in a TRANSFAC matrix table.

    The record inherits from a list containing the individual motifs.

    Attributes:
     - version - The version number, corresponding to the 'VV' field
       in the TRANSFAC file;

    """

    def __init__(self):
        """Initialize."""
        self.version = None

    def __str__(self):
        """Turn the TRANSFAC matrix into a string."""
        return write(self)


def read(handle, strict=True):
    """Parse a transfac format handle into a Record object."""
    annotations = {}
    references = []
    counts = None
    record = Record()
    for line in handle:
        line = line.strip()
        if not line:
            continue
        key_value = line.split(None, 1)
        key = key_value[0].strip()
        if strict:
            if len(key) != 2:
                raise ValueError("The key value of a TRANSFAC motif line "
                                 "should have 2 characters: "
                                 '"{0:s}"'.format(line))
        if len(key_value) == 2:
            value = key_value[1].strip()
            if strict:
                if not line.partition("  ")[1]:
                    raise ValueError("A TRANSFAC motif line should have 2 "
                                     "spaces between key and value columns: "
                                     '"{0:s}"'.format(line))
        if key == "VV":
            record.version = value
        elif key in ("P0", "PO"):  # Old TRANSFAC files use PO instead of P0
            counts = {}
            if value.split()[:4] != ["A", "C", "G", "T"]:
                raise ValueError('A TRANSFAC matrix "{0:s}" line should be '
                                 'followed by "A C G T": '
                                 '"{0:s}"'.format(key, line))
            length = 0
            for c in "ACGT":
                counts[c] = []
            for line in handle:
                line = line.strip()
                key_value = line.split(None, 1)
                key = key_value[0].strip()
                if len(key_value) == 2:
                    value = key_value[1].strip()
                    if strict:
                        if not line.partition("  ")[1]:
                            raise ValueError("A TRANSFAC motif line should "
                                             "have 2 spaces between key and "
                                             "value columns: "
                                             '"{0:s}"'.format(line))
                try:
                    i = int(key)
                except ValueError:
                    break
                if length == 0 and i == 0:
                    if strict:
                        raise ValueError("A TRANSFAC matrix should start with "
                                         '"01" as first row of the matrix, '
                                         'but this matrix uses "00": '
                                         '"{0:s}"'.format(line))
                else:
                    length += 1
                if i != length:
                    raise ValueError("The TRANSFAC matrix row number does not "
                                     "match the position in the matrix: "
                                     '"{0:s}"'.format(line))
                if strict:
                    if len(key) == 1:
                        raise ValueError("A TRANSFAC matrix line should have a "
                                         "2 digit key at the start of the line "
                                         '("{0:02d}"), but this matrix uses '
                                         '"{0:d}": "{1:s}".'.format(i, line))
                    if len(key_value) != 2:
                        raise ValueError("A TRANSFAC matrix line should have "
                                         "a key and a value: "
                                         '"{0:s}"'.format(line))
                values = value.split()[:4]
                if len(values) != 4:
                    raise ValueError("A TRANSFAC matrix line should have a "
                                     "value for each nucleotide "
                                     '(A, C, G and T): "{0:s}"'.format(line))
                for c, v in zip("ACGT", values):
                    counts[c].append(float(v))
        if line == "XX":
            pass
        elif key == "RN":
            index, separator, accession = value.partition(";")
            if index[0] != "[":
                raise ValueError('The index "{0:s}" in a TRANSFAC RN line '
                                 "should start with a "
                                 '"[": "{0:s}"'.format(index, line))
            if index[-1] != "]":
                raise ValueError('The index "{0:s}" in a TRANSFAC RN line '
                                 "should end with a "
                                 '"]": "{0:s}"'.format(index, line))
            index = int(index[1:-1])
            if len(references) != index - 1:
                raise ValueError('The index "{0:d}" of the TRANSFAC RN line '
                                 "does not match the current number of seen "
                                 'references "{1:d}": "{2:s}"'.format(index, len(references) + 1, line))
            reference = {key: value}
            references.append(reference)
        elif key == "//":
            if counts is not None:
                motif = Motif(alphabet="ACGT", counts=counts)
                motif.update(annotations)
                motif.references = references
                record.append(motif)
            annotations = {}
            references = []
        elif key in Motif.reference_keys:
            reference[key] = value
        elif key in Motif.multiple_value_keys:
            if key not in annotations:
                annotations[key] = []
            annotations[key].append(value)
        else:
            annotations[key] = value
    return record


def write(motifs):
    """Write the representation of a motif in TRANSFAC format."""
    blocks = []
    try:
        version = motifs.version
    except AttributeError:
        pass
    else:
        if version is not None:
            block = """\
VV  %s
XX
//
""" % version
            blocks.append(block)
    multiple_value_keys = Motif.multiple_value_keys
    sections = (("AC", "AS",),  # Accession
                ("ID",),        # ID
                ("DT", "CO"),   # Date, copyright
                ("NA",),        # Name
                ("DE",),        # Short factor description
                ("TY",),        # Type
                ("OS", "OC"),   # Organism
                ("HP", "HC"),   # Superfamilies, subfamilies
                ("BF",),        # Binding factors
                ("P0",),        # Frequency matrix
                ("BA",),        # Statistical basis
                ("BS",),        # Factor binding sites
                ("CC",),        # Comments
                ("DR",),        # External databases
                ("OV", "PV",),  # Versions
                )
    for motif in motifs:
        lines = []
        for section in sections:
            blank = False
            for key in section:
                if key == "P0":
                    # Frequency matrix
                    length = motif.length
                    if length == 0:
                        continue
                    sequence = motif.degenerate_consensus
                    letters = sorted(motif.alphabet)
                    line = "      ".join(["P0"] + letters)

                    lines.append(line)
                    for i in range(length):
                        line = " ".join(
                            ["%02.d"] +
                            ["%6.20g" for l in letters]) + \
                            "      %s"
                        line = line % tuple(
                            [i + 1] +
                            [motif.counts[l][i] for l in letters] +
                            [sequence[i]]
                        )
                        lines.append(line)
                    blank = True
                else:
                    try:
                        value = motif.get(key)
                    except AttributeError:
                        value = None
                    if value is not None:
                        if key in multiple_value_keys:
                            for v in value:
                                line = "%s  %s" % (key, v)
                                lines.append(line)
                        else:
                            line = "%s  %s" % (key, value)
                            lines.append(line)
                        blank = True
                if key == "PV":
                    # References
                    try:
                        references = motif.references
                    except AttributeError:
                        pass
                    else:
                        keys = ("RN", "RX", "RA", "RT", "RL")
                        for reference in references:
                            for key in keys:
                                value = reference.get(key)
                                if value is None:
                                    continue
                                line = "%s  %s" % (key, value)
                                lines.append(line)
                                blank = True
            if blank:
                line = "XX"
                lines.append(line)
        # Finished this motif; glue the lines together
        line = "//"
        lines.append(line)
        block = "\n".join(lines) + "\n"
        blocks.append(block)
    # Finished all motifs; glue the blocks together
    text = "".join(blocks)
    return text
