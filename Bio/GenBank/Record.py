"""Hold GenBank data in a straightforward format.

classes:
o Record - All of the information in a GenBank record.
o Reference - hold reference data for a record.
o Feature - Hold the information in a Feature Table.
o Qualifier - Qualifiers on a Feature.
"""
class Record:
    """Hold GenBank information in a format similar to the original record.

    The Record class is meant to make data easy to get to when you are
    just interested in looking at GenBank data.

    Attributes:
    o locus - The name specified after the LOCUS keyword in the GenBank
    record. This may be the accession number, or a clone id or something else.
    o size - The size of the record.
    o residue_type - The type of residues making up the sequence in this
    record. Normally something like RNA, DNA or PROTEIN, but may be as
    esoteric as 'ss-RNA circular'.
    o data_file_division - The division this record is stored under in
    GenBank (ie. PLN -> plants; PRI -> humans, primates; BCT -> bacteria...)
    o date - The date of submission of the record, in a form like '28-JUL-1998'
    o accession - list of all accession numbers for the sequence.
    o nid - Nucleotide identifier number.
    o version - The accession number + version (ie. AB01234.2)
    o gi - The NCBI gi identifier for the record.
    o keywords - A list of keywords related to the record.
    o segment - If the record is one of a series, this is info about which
    segment this record is (something like '1 of 6').
    o source - The source of material where the sequence came from.
    o organism - The genus and species of the organism (ie. 'Homo sapiens')
    o taxonomy - A listing of the taxonomic classification of the organism,
    starting general and getting more specific.
    o references - A list of Reference objects.
    o comment - Text with any kind of comment about the record.
    o features - A listing of Features making up the feature table.
    o base_counts - A string with the counts of bases for the sequence.
    o origin - A string specifying info about the origin of the sequence.
    o sequence - A string with the sequence itself.
    """
    def __init__(self):
        self.locus = ''
        self.size = ''
        self.residue_type = ''
        self.data_file_division = ''
        self.date = ''
        self.definition = ''
        self.accession = []
        self.nid = ''
        self.version = ''
        self.gi = ''
        self.keywords = []
        self.segment = ''
        self.source = ''
        self.organism = ''
        self.taxonomy = []
        self.references = []
        self.comment = ''
        self.features = []
        self.base_counts = ''
        self.origin = ''
        self.sequence = ''

class Reference:
    """Hold information from a GenBank reference.

    Attributes:
    o number - The number of the reference in the listing of references.
    o bases - The bases in the sequence the reference refers to.
    o authors - String with all of the authors.
    o title - The title of the reference.
    o journal - Information about the journal where the reference appeared.
    o medline_id - The medline id for the reference.
    o pubmed_id - The pubmed_id for the reference.
    o remark - Free-form remarks about the reference.
    """
    def __init__(self):
        self.number = ''
        self.bases = ''
        self.title = ''
        self.journal = ''
        self.medline_id = ''
        self.pubmed_id = ''
        self.remark = ''

class Feature:
    """Hold information about a Feature in the Feature Table of GenBank record.

    Attributes:
    o key - The key name of the featue (ie. source)
    o location - The string specifying the location of the feature.
    o qualfiers - A listing Qualifier objects in the feature.
    """
    def __init__(self):
        self.key = ''
        self.location = ''
        self.qualifiers = []

class Qualifier:
    """Hold information about a qualifier in a GenBank feature.

    Attributes:
    o key - The key name of the qualifier (ie. /organim=)
    o value - The value of the qualifier ("Dictyostelium discoideum").
    """
    def __init__(self):
        self.key = ''
        self.value = ''
