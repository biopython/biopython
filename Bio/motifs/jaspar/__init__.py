# Copyright 2013 by Anthony Mathelier and David Arenillas. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_dna as dna
import re
import math

from Bio._py3k import range

from Bio import motifs


class Motif(motifs.Motif):
    """
    A subclass of Bio.motifs.Motif used to represent a JASPAR profile with
    additional metadata information if available. The metadata availability
    depends on the source of the JASPAR motif (a 'pfm' format file, a 'jaspar'
    format file or a JASPAR database).

    """
    def __init__(self, matrix_id, name, alphabet=dna, instances=None,
                 counts=None, collection=None, tf_class=None, tf_family=None,
                 species=None, tax_group=None, acc=None, data_type=None,
                 medline=None, pazar_id=None, comment=None):
        """
        Construct a JASPAR Motif instance.

        """
        motifs.Motif.__init__(self, alphabet, instances, counts)
        self.name = name
        self.matrix_id = matrix_id
        self.collection = collection
        self.tf_class = tf_class
        self.tf_family = tf_family
        self.species = species      # May have multiple so species is a list.
                                    # The species are actually specified as
                                    # taxonomy IDs.
        self.tax_group = tax_group
        self.acc = acc              # May have multiple so acc is a list.
        self.data_type = data_type
        self.medline = medline
        self.pazar_id = pazar_id
        self.comment = comment

    @property
    def base_id(self):
        """
        Return the JASPAR base matrix ID
        """
        (base_id, version) = split_jaspar_id(self.matrix_id)
        return base_id

    @property
    def version(self):
        """
        Return the JASPAR matrix version
        """
        (base_id, version) = split_jaspar_id(self.matrix_id)
        return version

    def __str__(self):
        """
        Return a string represention of the JASPAR profile. We choose to
        provide only the filled metadata information.

        """
        tf_name_str = "TF name\t{0}\n".format(self.name)
        matrix_id_str = "Matrix ID\t{0}\n".format(self.matrix_id)
        the_string = "".join([tf_name_str, matrix_id_str])
        if self.collection:
            collection_str = "Collection\t{0}\n".format(self.collection)
            the_string = "".join([the_string, collection_str])
        if self.tf_class:
            tf_class_str = "TF class\t{0}\n".format(self.tf_class)
            the_string = "".join([the_string, tf_class_str])
        if self.tf_family:
            tf_family_str = "TF family\t{0}\n".format(self.tf_family)
            the_string = "".join([the_string, tf_family_str])
        if self.species:
            species_str = "Species\t{0}\n".format(",".join(self.species))
            the_string = "".join([the_string, species_str])
        if self.tax_group:
            tax_group_str = "Taxonomic group\t{0}\n".format(self.tax_group)
            the_string = "".join([the_string, tax_group_str])
        if self.acc:
            acc_str = "Accession\t{0}\n".format(self.acc)
            the_string = "".join([the_string, acc_str])
        if self.data_type:
            data_type_str = "Data type used\t{0}\n".format(self.data_type)
            the_string = "".join([the_string, data_type_str])
        if self.medline:
            medline_str = "Medline\t{0}\n".format(self.medline)
            the_string = "".join([the_string, medline_str])
        if self.pazar_id:
            pazar_id_str = "PAZAR ID\t{0}\n".format(self.pazar_id)
            the_string = "".join([the_string, pazar_id_str])
        if self.comment:
            comment_str = "Comments\t{0}\n".format(self.comment)
            the_string = "".join([the_string, comment_str])
        matrix_str = "Matrix:\n{0}\n\n".format(self.counts)
        the_string = "".join([the_string, matrix_str])
        return the_string

    def __hash__(self):
        """
        Return the hash key corresponding to the JASPAR profile

        :note: We assume the unicity of matrix IDs

        """
        return self.matrix_id.__hash__()

    def __eq__(self, other):
        return self.matrix_id == other.matrix_id


class Record(list):
    """
    Represents a list of jaspar motifs

    Attribute:
        o version: The JASPAR version used

    """

    def __init__(self):
        self.version = None

    def __str__(self):
        return "\n".join(str(the_motif) for the_motif in self)

    def to_dict(self):
        """
        Return the list of matrices as a dictionnary of matrices

        """

        dic = {}
        for motif in self:
            dic[motif.matrix_id] = motif
        return dic


def read(handle, format):
    """
    Read motif(s) from a file in one of several different JASPAR formats.
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
    """Returns the representation of the motifs in "pfm" or "jaspar" format
    """
    letters = "ACGT"
    lines = []
    if format == 'pfm':
        motif = motifs[0]
        counts = motif.counts
        for letter in letters:
            terms = ["{0:6.2f}".format(value) for value in counts[letter]]
            line = "{0}\n".format(" ".join(terms))
            lines.append(line)
    elif format == 'jaspar':
        for m in motifs:
            counts = m.counts
            line = ">{0} {1}\n".format(m.matrix_id, m.name)
            lines.append(line)
            for letter in letters:
                terms = ["{0:6.2f}".format(value) for value in counts[letter]]
                line = "{0} [{1}]\n".format(letter, " ".join(terms))
                lines.append(line)
    else:
        raise ValueError("Unknown JASPAR format %s" % format)

    # Finished; glue the lines together
    text = "".join(lines)

    return text


def _read_pfm(handle):
    """
    Reads the motif from a JASPAR .pfm file
    """
    alphabet = dna
    counts = {}

    letters = "ACGT"
    for letter, line in zip(letters, handle):
        words = line.split()
        #if there is a letter in the beginning, ignore it
        if words[0] == letter:
            words = words[1:]
        counts[letter] = [float(x) for x in words]

    motif = Motif(matrix_id=None, name=None, alphabet=alphabet, counts=counts)
    motif.mask = "*" * motif.length
    record = Record()
    record.append(motif)

    return record


def _read_sites(handle):
    """
    Reads the motif from JASPAR .sites file
    """

    alphabet = dna
    instances = []

    for line in handle:
        if not line.startswith(">"):
            break
        # line contains the header ">...."
        # now read the actual sequence
        line = next(handle)
        instance = ""
        for c in line.strip():
            if c == c.upper():
                instance += c
        instance = Seq(instance, alphabet)
        instances.append(instance)

    instances = motifs.Instances(instances, alphabet)
    motif = Motif(
        matrix_id=None, name=None, alphabet=alphabet, instances=instances
    )
    motif.mask = "*" * motif.length
    record = Record()
    record.append(motif)

    return record


def _read_jaspar(handle):
    """
    Read motifs from a JASPAR formatted file

    Format is one or more records of the form, e.g.:
    >MA0001.1 AGL3
    A  [ 0  3 79 40 66 48 65 11 65  0 ]
    C  [94 75  4  3  1  2  5  2  3  3 ]
    G  [ 1  0  3  4  1  0  5  3 28 88 ]
    T  [ 2 19 11 50 29 47 22 81  1  6 ]

    """

    alphabet = dna
    counts = {}

    record = Record()

    head_pat = re.compile(r"^>\s*(\S+)(\s+(\S+))?")
    row_pat = re.compile(r"\s*([ACGT])\s*\[\s*(.*)\s*\]")

    identifier = None
    name = None
    row_count = 0
    for line in handle:
        line.rstrip('\r\n')

        head_match = head_pat.match(line)
        row_match = row_pat.match(line)

        if head_match:
            identifier = head_match.group(1)
            if head_match.group(2):
                name = head_match.group(2)
            else:
                name = identifier
        elif row_match:
            (letter, counts_str) = row_match.group(1, 2)

            words = counts_str.split()

            counts[letter] = [float(x) for x in words]

            row_count += 1

            if row_count == 4:
                record.append(Motif(identifier, name, alphabet=alphabet,
                                    counts=counts))

                identifier = None
                name = None
                counts = {}
                row_count = 0

    return record

def calculate_pseudocounts(motif):
    alphabet = motif.alphabet
    background = motif.background

    # It is possible to have unequal column sums so use the average
    # number of instances.
    total = 0
    for i in range(motif.length):
        total += sum(float(motif.counts[letter][i]) for letter in alphabet.letters)

    avg_nb_instances = total / motif.length
    sq_nb_instances = math.sqrt(avg_nb_instances)

    if background:
        background = dict(background)
    else:
        background = dict.fromkeys(sorted(alphabet.letters), 1.0)

    total = sum(background.values())
    pseudocounts = {}

    for letter in alphabet.letters:
        background[letter] /= total
        pseudocounts[letter] = sq_nb_instances * background[letter]

    return pseudocounts

def split_jaspar_id(id):
    """
    Utility function to split a JASPAR matrix ID into its component base ID
    and version number, e.g. 'MA0047.2' is returned as ('MA0047', 2).
    """

    id_split = id.split('.')

    base_id = None
    version = None
    if len(id_split) == 2:
        base_id = id_split[0]
        version = id_split[1]
    else:
        base_id = id

    return (base_id, version)
