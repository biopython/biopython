from xml.sax import handler

from Martel import Dispatch

from Bio import SeqRecord, StdHandler, Std, DBXRef, Seq
from Bio import Alphabet, SeqFeature
from Bio.Alphabet import IUPAC

alphabet_table = {
    "iupac-protein": IUPAC.protein,
    "iupac-dna": IUPAC.unambiguous_dna,
    "iupac-rna": IUPAC.unambiguous_rna,
    "iupac-extended-protein": IUPAC.extended_protein,
    "iupac-ambiguous-dna": IUPAC.ambiguous_dna,
    "iupac-ambiguous-rna": IUPAC.ambiguous_rna,
    "protein": Alphabet.generic_protein,
    "dna": Alphabet.generic_dna,
    "rna": Alphabet.generic_rna,
    "unknown": Alphabet.single_letter_alphabet,
    }


# Convert from the internal Feature data structure used by the parser
# into the standard Biopytho form
def convert_std_feature(feature):
    f = SeqFeature.SeqFeature()
    loc = feature.location
    return feature

class BuildSeqRecord(Dispatch.Dispatcher):
    def __init__(self):
        Dispatch.Dispatcher.__init__(self)
        self.acquire(StdHandler.Handle_dbid(self.add_dbid),
                     prefix = Std.NS)
        self.acquire(StdHandler.Handle_dbxref(self.add_dbxref_dbids),
                     prefix = Std.NS)
        self.acquire(StdHandler.Handle_description(self.add_description),
                     prefix = Std.NS)
        self.acquire(StdHandler.Handle_sequence(self.add_sequence),
                     prefix = Std.NS)
        self.acquire(StdHandler.Handle_features(self.add_features),
                     prefix = Std.NS)
        self.acquire(StdHandler.Handle_dbxref(self.add_dbxref),
                     prefix = Std.NS)

    def start_record(self, tag, attrs):
        self.dbname = None
        self.id_text = None
        self.name_text = '<unknown name>'
        self.description = None
        self.alphabet = None
        self.seq = None
        self.features = None
        self.dbxrefs = []

    def add_dbid(self, text, attrs):
        if attrs.get("type") == "primary":
            self.dbname = attrs.get("dbname", "unknown")
            self.id_text = text
        # use the first accession/secondary id as the name
        # this should be equivalent to what Biopython does
        elif attrs.get("type") in ["accession", "secondary"]:
            self.name_text = text

    def add_dbxref_dbids(self, dbname_style, dbname, idtype, dbid, negate):
        """Handle setting name and id attributes from the dbxref ids.

        Likely we'll either have a dbid or dbxref dbids to use. We default
        to using the dbid if it exists.
        """
        # first deal with the primary id: SeqFeature.id
        # set the id if we haven't yet set an id (take the first id we get)
        # and if we have a primary id
        if (self.id_text is None and idtype == "primary"):
            self.id_text = dbid

        # now deal with secondary ids: SeqFeature.name
        if idtype == "secondary":
            self.name_text = dbid


    def add_description(self, text):
        self.description = text

    def add_sequence(self, (alphabet, seq, gapchar, stopchar)):
        alphabet = alphabet_table.get(alphabet,
                                      Alphabet.single_letter_alphabet)
        self.seq = Seq.Seq(seq, alphabet)

    def add_dbxref(self, dbname_style, dbname, idtype, dbid, negate):
        """Store all id cross references.
        """
        self.dbxrefs.append(DBXRef.from_parser(dbname_style, dbname, idtype,
                                               dbid, negate))

    def add_features(self, features):
        # Brad -- I can't understand this assertion -- there is no
        # self.features on the first call to add features and then you'll
        # expect to have some on future calls
        # assert self.features is None
        #print [feature.location for feature in features]
        self.features = map(convert_std_feature, features)

    def end_record(self, tag):
        self.document = SeqRecord.SeqRecord(
            seq = self.seq,
            id = self.id_text,
            name = self.name_text,
            description = self.description,
            dbxrefs = self.dbxrefs,
            features = self.features,
            )
        
make_builder = BuildSeqRecord
