# Define the IUPAC Alphabets you know and love

from Bio import Alphabet
from Bio.Data import IUPACData

##################### Protein

# From the IUPAC definition at:
#   http://www.chem.qmw.ac.uk/iupac/AminoAcid/A2021.html#AA21
class IUPACProtein(Alphabet.ProteinAlphabet):
    letters = IUPACData.protein_letters

protein = IUPACProtein()

# This could be considered the base class for the standard IUPAC
# protein, except that some encodings will use "X" to mean "unknown
# character", which causes a collision.  If you use X for
# selenocysteines, then you'll need a new alphabet.

class ExtendedIUPACProtein(Alphabet.ProteinAlphabet):
    letters = IUPACData.extended_protein_letters
    # B = "Asx";  aspartic acid or asparagine
    # X = "Sec";  selenocysteine  Note: IUPAC is moving to use 'U' for this
    # Z = "Glx";  glutamic acid or glutamine (or substances such as
    #         4-carboxyglutamic acid and 5-oxoproline that yield glutamic
    #         acid on acid hydrolysis of peptides)

extended_protein = ExtendedIUPACProtein()

##################### DNA

# The next two are the IUPAC definitions, from:
#   http://www.chem.qmw.ac.uk/iubmb/misc/naseq.html
class IUPACAmbiguousDNA(Alphabet.DNAAlphabet):
    letters = IUPACData.ambiguous_dna_letters

ambiguous_dna = IUPACAmbiguousDNA()

class IUPACUnambiguousDNA(IUPACAmbiguousDNA):
    letters = IUPACData.unambiguous_dna_letters

unambiguous_dna = IUPACUnambiguousDNA()


# Also from the URL, but not part of the standard
class ExtendedIUPACDNA(Alphabet.DNAAlphabet):
    letters = IUPACData.extended_dna_letters
    #   B == 5-bromouridine
    #   D == 5,6-dihydrouridine
    #   S == thiouridine
    #   W == wyosine

extended_dna = ExtendedIUPACDNA()

##################### RNA

class IUPACAmbiguousRNA(Alphabet.RNAAlphabet):
    letters = IUPACData.ambiguous_rna_letters

ambiguous_rna = IUPACAmbiguousRNA()

class IUPACUnambiguousRNA(IUPACAmbiguousRNA):
    letters = IUPACData.unambiguous_rna_letters

unambiguous_rna = IUPACUnambiguousRNA()

# are there extended forms?
#class ExtendedIUPACRNA(Alphabet.RNAAlphabet):
#    letters = extended_rna_letters
#    #   B == 5-bromouridine
#    #   D == 5,6-dihydrouridine
#    #   S == thiouridine
#    #   W == wyosine


# We need to load the property resolution information, but we need to
# wait until after the systems have been loaded. (There's a nasty loop
# where, eg, translation objects need an alphabet, which need to be
# assocated with translators.)

from Bio.PropertyManager import default_manager

def _bootstrap(manager, klass, property):
    assert manager is default_manager
    del default_manager.class_resolver[IUPACProtein]
    del default_manager.class_resolver[ExtendedIUPACProtein]
    del default_manager.class_resolver[IUPACAmbiguousDNA]
    del default_manager.class_resolver[IUPACUnambiguousDNA]
    del default_manager.class_resolver[ExtendedIUPACDNA]
    del default_manager.class_resolver[IUPACAmbiguousRNA]
    del default_manager.class_resolver[IUPACUnambiguousRNA]

    from Bio.Encodings import IUPACEncoding

    return manager.resolve_class(klass, property)

default_manager.class_resolver[IUPACProtein] = _bootstrap
default_manager.class_resolver[ExtendedIUPACProtein] = _bootstrap
default_manager.class_resolver[IUPACAmbiguousDNA] = _bootstrap
default_manager.class_resolver[IUPACUnambiguousDNA] = _bootstrap
default_manager.class_resolver[ExtendedIUPACDNA] = _bootstrap
default_manager.class_resolver[IUPACAmbiguousRNA] = _bootstrap
default_manager.class_resolver[IUPACUnambiguousRNA] = _bootstrap
