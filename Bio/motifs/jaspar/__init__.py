# Copyright 2013 by Anthony Mathelier and David Arenillas. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""JASPAR2014 module."""

from Bio.Seq import Seq
import re
import math

from Bio import motifs
from Bio import Align


class Motif(motifs.Motif):
    """A subclass of Bio.motifs.Motif used to represent a JASPAR profile.

    Additional metadata information are stored if available. The metadata
    availability depends on the source of the JASPAR motif (a 'pfm' format
    file, a 'jaspar' format file or a JASPAR database).
    """

    def __init__(
        self,
        matrix_id,
        name,
        alphabet="ACGT",
        alignment=None,
        counts=None,
        collection=None,
        tf_class=None,
        tf_family=None,
        species=None,
        tax_group=None,
        acc=None,
        data_type=None,
        medline=None,
        pazar_id=None,
        comment=None,
    ):
        """Construct a JASPAR Motif instance."""
        motifs.Motif.__init__(self, alphabet, alignment, counts)
        self.name = name
        self.matrix_id = matrix_id
        self.collection = collection
        self.tf_class = tf_class
        self.tf_family = tf_family
        # May have multiple so species is a list.
        # The species are actually specified as
        # taxonomy IDs.
        self.species = species
        self.tax_group = tax_group
        self.acc = acc  # May have multiple so acc is a list.
        self.data_type = data_type
        self.medline = medline
        self.pazar_id = pazar_id
        self.comment = comment

    @property
    def base_id(self):
        """Return the JASPAR base matrix ID."""
        (base_id, __) = split_jaspar_id(self.matrix_id)
        return base_id

    @property
    def version(self):
        """Return the JASPAR matrix version."""
        (__, version) = split_jaspar_id(self.matrix_id)
        return version

    def __str__(self):
        """Return a string represention of the JASPAR profile.

        We choose to provide only the filled metadata information.
        """
        tf_name_str = f"TF name\t{self.name}\n"
        matrix_id_str = f"Matrix ID\t{self.matrix_id}\n"
        the_string = "".join([tf_name_str, matrix_id_str])
        if self.collection:
            collection_str = f"Collection\t{self.collection}\n"
            the_string = "".join([the_string, collection_str])
        if self.tf_class:
            tf_class_str = f"TF class\t{self.tf_class}\n"
            the_string = "".join([the_string, tf_class_str])
        if self.tf_family:
            tf_family_str = f"TF family\t{self.tf_family}\n"
            the_string = "".join([the_string, tf_family_str])
        if self.species:
            species_str = f"Species\t{','.join(self.species)}\n"
            the_string = "".join([the_string, species_str])
        if self.tax_group:
            tax_group_str = f"Taxonomic group\t{self.tax_group}\n"
            the_string = "".join([the_string, tax_group_str])
        if self.acc:
            acc_str = f"Accession\t{self.acc}\n"
            the_string = "".join([the_string, acc_str])
        if self.data_type:
            data_type_str = f"Data type used\t{self.data_type}\n"
            the_string = "".join([the_string, data_type_str])
        if self.medline:
            medline_str = f"Medline\t{self.medline}\n"
            the_string = "".join([the_string, medline_str])
        if self.pazar_id:
            pazar_id_str = f"PAZAR ID\t{self.pazar_id}\n"
            the_string = "".join([the_string, pazar_id_str])
        if self.comment:
            comment_str = f"Comments\t{self.comment}\n"
            the_string = "".join([the_string, comment_str])
        matrix_str = f"Matrix:\n{self.counts}\n\n"
        the_string = "".join([the_string, matrix_str])
        return the_string

    def __hash__(self):
        """Return the hash key corresponding to the JASPAR profile.

        :note: We assume the unicity of matrix IDs

        """
        return self.matrix_id.__hash__()

    def __eq__(self, other):
        """Return True if matrix IDs are the same."""
        return self.matrix_id == other.matrix_id


class Record(list):
    """Represent a list of jaspar motifs.

    Attributes:
     - version: The JASPAR version used

    """

    def __init__(self):
        """Initialize the class."""
        self.version = None

    def __str__(self):
        """Return a string of all motifs in the Record."""
        return "\n".join(str(the_motif) for the_motif in self)

    def to_dict(self):
        """Return the list of matrices as a dictionary of matrices."""
        dic = {}
        for motif in self:
            dic[motif.matrix_id] = motif
        return dic


def read(handle, format):
    """Read motif(s) from a file in one of several different JASPAR formats.

    Return the record of PFM(s).
    Call the appropriate routine based on the format passed.
    """
    format = format.lower()
    if format == "pfm":
        record = _read_pfm(handle)
        return record
    elif format == "sites":
        record = _read_sites(handle)
        return record
    elif format == "jaspar":
        record = _read_jaspar(handle)
        return record
    else:
        raise ValueError("Unknown JASPAR format %s" % format)


def write(motifs, format):
    """Return the representation of motifs in "pfm" or "jaspar" format."""
    letters = "ACGT"
    lines = []
    if format == "pfm":
        motif = motifs[0]
        counts = motif.counts
        for letter in letters:
            terms = [f"{value:6.2f}" for value in counts[letter]]
            line = f"{' '.join(terms)}\n"
            lines.append(line)
    elif format == "jaspar":
        for m in motifs:
            counts = m.counts
            try:
                matrix_id = m.matrix_id
            except AttributeError:
                matrix_id = None
            line = f">{matrix_id} {m.name}\n"
            lines.append(line)
            for letter in letters:
                terms = [f"{value:6.2f}" for value in counts[letter]]
                line = f"{letter} [{' '.join(terms)}]\n"
                lines.append(line)
    else:
        raise ValueError("Unknown JASPAR format %s" % format)

    # Finished; glue the lines together
    text = "".join(lines)

    return text


def _read_pfm(handle):
    """Read the motif from a JASPAR .pfm file (PRIVATE)."""
    alphabet = "ACGT"
    counts = {}

    for letter, line in zip(alphabet, handle):
        words = line.split()
        # if there is a letter in the beginning, ignore it
        if words[0] == letter:
            words = words[1:]
        counts[letter] = [float(x) for x in words]

    motif = Motif(matrix_id=None, name=None, alphabet=alphabet, counts=counts)
    motif.mask = "*" * motif.length
    record = Record()
    record.append(motif)

    return record


def _read_sites(handle):
    """Read the motif from JASPAR .sites file (PRIVATE)."""
    alphabet = "ACGT"
    instances = []

    for line in handle:
        if not line.startswith(">"):
            break
        # line contains the header ">...."
        # now read the actual sequence
        line = next(handle)
        instance = ""
        for c in line.strip():
            if c.isupper():
                instance += c
        instance = Seq(instance)
        instances.append(instance)

    alignment = Align.Alignment(instances)
    motif = Motif(matrix_id=None, name=None, alphabet=alphabet, alignment=alignment)
    motif.mask = "*" * motif.length
    record = Record()
    record.append(motif)

    return record


def _read_jaspar(handle):
    """Read motifs from a JASPAR formatted file (PRIVATE).

    Format is one or more records of the form, e.g.::

      - JASPAR 2010 matrix_only format::

                >MA0001.1 AGL3
                A  [ 0  3 79 40 66 48 65 11 65  0 ]
                C  [94 75  4  3  1  2  5  2  3  3 ]
                G  [ 1  0  3  4  1  0  5  3 28 88 ]
                T  [ 2 19 11 50 29 47 22 81  1  6 ]

      - JASPAR 2010-2014 PFMs format::

                >MA0001.1 AGL3
                0	3	79	40	66	48	65	11	65	0
                94	75	4	3	1	2	5	2	3	3
                1	0	3	4	1	0	5	3	28	88
                2	19	11	50	29	47	22	81	1	6

    """
    alphabet = "ACGT"
    counts = {}

    record = Record()

    head_pat = re.compile(r"^>\s*(\S+)(\s+(\S+))?")
    row_pat_long = re.compile(r"\s*([ACGT])\s*\[\s*(.*)\s*\]")
    row_pat_short = re.compile(r"\s*(.+)\s*")

    identifier = None
    name = None
    row_count = 0
    nucleotides = ["A", "C", "G", "T"]
    for line in handle:
        line = line.strip()

        head_match = head_pat.match(line)
        row_match_long = row_pat_long.match(line)
        row_match_short = row_pat_short.match(line)

        if head_match:
            identifier = head_match.group(1)
            if head_match.group(3):
                name = head_match.group(3)
            else:
                name = identifier
        elif row_match_long:
            (letter, counts_str) = row_match_long.group(1, 2)
            words = counts_str.split()
            counts[letter] = [float(x) for x in words]
            row_count += 1
            if row_count == 4:
                record.append(Motif(identifier, name, alphabet=alphabet, counts=counts))
                identifier = None
                name = None
                counts = {}
                row_count = 0
        elif row_match_short:
            words = row_match_short.group(1).split()
            counts[nucleotides[row_count]] = [float(x) for x in words]
            row_count += 1
            if row_count == 4:
                record.append(Motif(identifier, name, alphabet=alphabet, counts=counts))
                identifier = None
                name = None
                counts = {}
                row_count = 0

    return record


def calculate_pseudocounts(motif):
    """Calculate pseudocounts.

    Computes the root square of the total number of sequences multiplied by
    the background nucleotide.
    """
    alphabet = motif.alphabet
    background = motif.background

    # It is possible to have unequal column sums so use the average
    # number of instances.
    total = 0
    for i in range(motif.length):
        total += sum(motif.counts[letter][i] for letter in alphabet)

    avg_nb_instances = total / motif.length
    sq_nb_instances = math.sqrt(avg_nb_instances)

    if background:
        background = dict(background)
    else:
        background = dict.fromkeys(sorted(alphabet), 1.0)

    total = sum(background.values())
    pseudocounts = {}

    for letter in alphabet:
        background[letter] /= total
        pseudocounts[letter] = sq_nb_instances * background[letter]

    return pseudocounts


def split_jaspar_id(id):
    """Split a JASPAR matrix ID into its component.

    Components are base ID and version number, e.g. 'MA0047.2' is returned as
    ('MA0047', 2).
    """
    id_split = id.split(".")

    base_id = None
    version = None
    if len(id_split) == 2:
        base_id = id_split[0]
        version = id_split[1]
    else:
        base_id = id

    return (base_id, version)
