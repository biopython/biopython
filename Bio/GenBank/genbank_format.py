"""Martel based parser to read GenBank formatted files.

This is a huge regular regular expression for GenBank, built using
the 'regular expressiona on steroids' capabilities of Martel.

Notes:
Just so I remember -- the new end of line syntax is:
  New regexp syntax - \R
     \R    means "\n|\r\n?"
     [\R]  means "[\n\r]"

This helps us have endlines be consistent across platforms.

Documentation for GenBank format that I found:

o GenBank/EMBL feature tables are described at:
http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html

o There are also descriptions of different GenBank lines at:
http://www.ibc.wustl.edu/standards/gbrel.txt
"""
# standard library
import string
     
# Martel
import Martel
from Martel import RecordReader

# --- first set up some helper constants and functions

# - useful constants for dealing with the blank space in GenBank documents
# this is useful since blank space can be significant in GenBank flat files.
INDENT = 12
FEATURE_KEY_INDENT = 5
FEATURE_QUALIFIER_INDENT = 21

blank_space = Martel.Rep1(Martel.Str(" "))
small_indent_space = Martel.Str(" " * 2)
big_indent_space = Martel.Str(" " * FEATURE_KEY_INDENT)
qualifier_space = Martel.Str(" " * FEATURE_QUALIFIER_INDENT)

# - useful functions
def define_block(identifier, block_tag, block_data):
    """Define a Martel grouping which can parse a block of text.

    Many of the GenBank lines we'll want to process are grouped into
    a block like:

    IDENTIFIER   Blah blah blah

    Where blah blah blah can wrap for multiple lines. This function makes
    it easy to consistently define a definition for these blocks.

    Arguments:
    o identifier - The identifier that begins the block (like DEFINITION).
    o block_tag - A callback tag for the entire block.
    o block_data - A callback tag for the data in the block (ie. the
    stuff you are interested in).
    """
    diff = INDENT - len(identifier)
    assert diff > 0, diff

    return Martel.Group(block_tag,
                        Martel.Str(identifier + " " * diff) +
                        Martel.ToEol(block_data) +
                        Martel.Rep(Martel.Str(" " * INDENT) + 
                                   Martel.ToEol(block_data)))


# first line
# LOCUS       AC007323    86436 bp    DNA             PLN       19-JAN-2000
locus = Martel.Group("locus",
                     Martel.Re("[\w]+"))
size = Martel.Group("size",
                    Martel.Rep1(Martel.Integer()))

# deal with the different kinds of residues we can have
valid_residue_prefixes = ["ss-", "ds-", "ms-"]
valid_residue_types = ["DNA", "RNA", "mRNA", "tRNA", "rRNA", "uRNA",
                       "snRNA", "PROTEIN"]

residue_prefixes = map(Martel.Str, valid_residue_prefixes)
residue_types = map(Martel.Str, valid_residue_types)

residue_type = Martel.Group("residue_type",
                            Martel.Opt(Martel.Alt(*residue_prefixes)) +
                            Martel.Opt(Martel.Alt(*residue_types)) +
                            Martel.Opt(blank_space +
                                       Martel.Str("circular")))
date = Martel.Group("date",
                    Martel.Re("[-\w]+"))

# the PLN, etc stuff indicates data file divisions
valid_divisions = ["PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "RNA",
                   "VRL", "PHG", "SYN", "UNA", "EST", "PAT", "STS", "GSS",
                   "HTG", "HTC"]
divisions = map(Martel.Str, valid_divisions)
data_file_division = Martel.Group("data_file_division",
                                  Martel.Alt(*divisions))

locus_line = Martel.Group("locus_line",
                          Martel.Str("LOCUS") +
                          blank_space +
                          locus +
                          blank_space +
                          size +
                          blank_space +
                          Martel.Re("bp|aa") +
                          blank_space +
                          Martel.Opt(residue_type +
                                     blank_space) +
                          data_file_division +
                          blank_space +
                          date +
                          Martel.AnyEol()) 

# definition line
# DEFINITION  Genomic sequence for Arabidopsis thaliana BAC T25K16 from
#             chromosome I, complete sequence.
definition_block = define_block("DEFINITION", "definition_block", "definition")

# accession line
# ACCESSION   AC007323
accession = Martel.Group("accession",
                         Martel.Re("[\w]+"))

accession_block = Martel.Group("accession_block",
                               Martel.Str("ACCESSION") +
                               Martel.Rep1(blank_space +
                                           Martel.Rep1(accession +
                                                Martel.Opt(Martel.Str(" "))) +
                                           Martel.AnyEol()))

# accession_block = define_block("ACCESSION", "accession_block", "accession")

# NID         g44010
nid = Martel.Group("nid",
                   Martel.Re("[\w\d]+"))
nid_line = Martel.Group("nid_line",
                        Martel.Str("NID") +
                        blank_space +
                        nid +
                        Martel.AnyEol())

# version and GI line
# VERSION     AC007323.5  GI:6587720
version = Martel.Group("version",
                       Martel.Re("[\w\d\.]+"))

gi = Martel.Group("gi",
                  Martel.Re("[\d]+"))

version_line = Martel.Group("version_line",
                            Martel.Str("VERSION") +
                            blank_space +
                            version +
                            blank_space +
                            Martel.Str("GI:") +
                            gi +
                            Martel.AnyEol())

# keywords line
# KEYWORDS    antifreeze protein homology; cold-regulated gene; cor6.6 gene;
#             KIN1 homology.
keywords_block = define_block("KEYWORDS", "keywords_block", "keywords")

# SEGMENT     1 of 6
segment = Martel.Group("segment",
                       Martel.Integer("segment_num") + \
                       Martel.Str(" of ") + \
                       Martel.Integer("segment_total"))
segment_line = Martel.Group("segment_line",
                            Martel.Str("SEGMENT     ") + segment + \
                            Martel.AnyEol())

# SOURCE      thale cress.
source_block = define_block("SOURCE", "source_block", "source")

# ORGANISM  Arabidopsis thaliana
#           Eukaryota; Viridiplantae; Embryophyta; Tracheophyta; Spermatophyta;
#           Magnoliophyta; eudicotyledons; core eudicots; Rosidae; eurosids II;
#            Brassicales; Brassicaceae; Arabidopsis.
organism = Martel.Group("organism",
                        Martel.ToEol())

taxonomy = Martel.Group("taxonomy",
                        Martel.Rep1(blank_space +
                                    Martel.ToEol()))

organism_block = Martel.Group("organism_block",
                              Martel.Str("  ORGANISM") +
                              blank_space +
                              organism +
                              taxonomy)

# REFERENCE   1  (bases 1 to 86436)
#   AUTHORS   Thomashow,M.F.
#   TITLE     Direct Submission
#   JOURNAL   Submitted (01-FEB-1991) M.F. Thomashow, Dept. Crop and Soil
#             Sciences, Dept. Microbiology, Michigan State University, East
#             Lansing, Michigan 48824, USA
reference_num = Martel.Group("reference_num",
                             Martel.Re("[\d]+"))

# can have normal references, like that shown above, or references like:
# REFERENCE   1  (sites)
# with no base information or even:
# REFERENCE   2  (bases 1 to 105654; 110423 to 111122)
reference_bases = Martel.Group("reference_bases",
                               Martel.Str("(") +
                               Martel.Re("[;\w\d \R]+") +
                               Martel.Str(")"))
reference_line = Martel.Group("reference_line",
                              Martel.Str("REFERENCE") +
                              blank_space +
                              reference_num +
                              Martel.Opt(blank_space +
                                         reference_bases) +
                              Martel.AnyEol())

authors_block = define_block("  AUTHORS", "authors_block", "authors")
title_block = define_block("  TITLE", "title_block", "title")
journal_block = define_block("  JOURNAL", "journal_block", "journal")

#  MEDLINE   92119220
medline_line = Martel.Group("medline_line",
                            Martel.Str("  MEDLINE   ") +
                            Martel.Integer("medline_id") +
                            Martel.AnyEol())

# PUBMED   10617197
pubmed_line = Martel.Group("pubmed_line",
                           Martel.Str("   PUBMED   ") +
                           Martel.Integer("pubmed_id") +
                           Martel.AnyEol())

# REMARK    This sequence is of BAC F10O3 from Arabidopsis thaliana chromosome
remark_block = define_block("  REMARK", "remark_block", "remark")

# an entire reference for the sequence
reference = Martel.Group("reference",
                         reference_line +
                         authors_block +
                         Martel.Opt(title_block) +
                         journal_block +
                         Martel.Opt(medline_line) +
                         Martel.Opt(pubmed_line) +
                         Martel.Opt(remark_block))

# COMMENT     On Dec 16, 1999 this sequence version replaced gi:5729683.
comment_block = define_block("COMMENT", "comment_block", "comment")

# start on the feature table. Eeek -- This is the part I was afraid of
# most!

# the header, so that we know we are heading into some features
# FEATURES             Location/Qualifiers
features_line = Martel.Group("features_line",
                             Martel.Str("FEATURES") +
                             blank_space +
                             Martel.Str("Location/Qualifiers") +
                             Martel.AnyEol())

# --- now we need to read in the features one at a time
# -- first, set up the feature keys and locations
# a listing of valid feature keys
# XXX Because of a bug in Martel, if one name is a substring of
# XXX another, put the longer name first.  Eg, "primer_bind" before "bind"
feature_key_names = (
    "allele",           # Obsolete; see variation feature key
    "attenuator",       # Sequence related to transcription termination
    "C_region",         # Span of the C immunological feature
    "CAAT_signal",      # 'CAAT box' in eukaryotic promoters
    "CDS",              # Sequence coding for amino acids in protein (includes
                        #   stop codon)
    "conflict",         # Independent sequence determinations differ
    "D-loop",           # Displacement loop
    "D_segment",        # Span of the D immunological feature
    "enhancer",         # Cis-acting enhancer of promoter function
    "exon",             # Region that codes for part of spliced mRNA
    "GC_signal",        # 'GC box' in eukaryotic promoters
    "gene",             # Region that defines a functional gene, possibly
                        #   including upstream (promotor, enhancer, etc)
                        #   and downstream control elements, and for which
                        #   a name has been assigned.
    "iDNA",             # Intervening DNA eliminated by recombination
    "intron",           # Transcribed region excised by mRNA splicing
    "J_segment",         # Span of the J immunological feature
    "LTR",              # Long terminal repeat
    "mat_peptide",      # Mature peptide coding region (does not include
                        #   stop codon)
    "misc_binding",     # Miscellaneous binding site
    "misc_difference",  # Miscellaneous difference feature
    "misc_feature",     # Region of biological significance that cannot
                        #   be described by any other feature
    "misc_recomb",      # Miscellaneous recombination feature
    "misc_RNA",         # Miscellaneous transcript feature not defined by
                        #   other RNA keys
    "misc_signal",      # Miscellaneous signal
    "misc_structure",   # Miscellaneous DNA or RNA structure
    "modified_base",    # The indicated base is a modified nucleotide
    "mRNA",             # Messenger RNA
    "mutation",         # Obsolete: see variation feature key
    "N_region",         # Span of the N immunological feature
    "old_sequence",     # Presented sequence revises a previous version
    "polyA_signal",     # Signal for cleavage & polyadenylation
    "polyA_site",       # Site at which polyadenine is added to mRNA
    "precursor_RNA",    # Any RNA species that is not yet the mature
                        #   RNA product
    "prim_transcript",  # Primary (unprocessed) transcript
    "primer_bind",      # Non-covalent primer binding site
    "primer",           # Primer binding region used with PCR  XXX not in 
                        #   http://www.ncbi.nlm.nih.gov/collab/FT/index.html
    "promoter",         # A region involved in transcription initiation
    "protein_bind",     # Non-covalent protein binding site on DNA or RNA
    "RBS",              # Ribosome binding site
    "rep_origin",       # Replication origin for duplex DNA
    "repeat_region",    # Sequence containing repeated subsequences
    "repeat_unit",      # One repeated unit of a repeat_region
    "rRNA",             # Ribosomal RNA
    "S_region",         # Span of the S immunological feature
    "satellite",        # Satellite repeated sequence
    "scRNA",            # Small cytoplasmic RNA
    "sig_peptide",      # Signal peptide coding region
    "snRNA",            # Small nuclear RNA
    "source",           # Biological source of the sequence data
                        #   represented by a GenBank record. Mandatory
                        #   feature, one or more per record.  For organisms
                        #   that have been incorporated within the NCBI
                        #   taxonomy database, an associated
                        #   /db_xref="taxon:NNNN" qualifier will be present
                        #   (where NNNNN is the numeric identifier assigned
                        #   to the organism within the NCBI taxonomy
                        #   database).
    "stem_loop",        # Hair-pin loop structure in DNA or RNA
    "STS",              # Sequence Tagged Site; operationally unique
                        #   sequence that identifies the combination of
                        #   primer spans used in a PCR assay
    "TATA_signal",      # 'TATA box' in eukaryotic promoters
    "terminator",       # Sequence causing transcription termination
    "transit_peptide",  # Transit peptide coding region
    "transposon",       # Transposable element (TN)
    "tRNA",             # Transfer RNA
    "unsure",           # Authors are unsure about the sequence in this region
    "V_region",         # Span of the V immunological feature
    "V_segment",        # Variable segment of immunoglobulin light and heavy
                        #   chains, and T-cell receptor alpha, beta, and
                        #   gamma chains
    "variation",        # A related population contains stable mutation
    "-10_signal",       # 'Pribnow box' in prokaryotic promoters
    "-35_signal",       # '-35 box' in prokaryotic promoters
    "3'clip",           # 3'-most region of a precursor transcript removed'
                        #   in processing
    "3'UTR",            # 3' untranslated region (trailer)'
    "5'clip",           # 5'-most region of a precursor transcript removed'
                        # in processing
    "5'UTR",            # 5' untranslated region (leader)'
    "-"                 # (hyphen)      Placeholder
)
valid_feature_keys = map(Martel.Str, feature_key_names)

feature_key = Martel.Group("feature_key",
                           Martel.Alt(*valid_feature_keys))

# handle lots of different kinds of locations
# complement(10..20)
# join(10..20,30..40)
# 10..20
# we can have an optional reference to another accession number, ie:
# J00194:(100..202)
# can also have a version ie. A10000.1
location_ref = Martel.Group("location_ref",
                            Martel.Re("[_\d\w\.]+") +
                            Martel.Str(":"))
"""
location_part  = Martel.Group("location_part",
                              Martel.Rep1(Martel.Re("[\<\>\(\)\^\.\,\d]") |
                                          Martel.Str("complement") |
                                          Martel.Str("join") |
                                          Martel.Str("order") |
                                          Martel.Str("replace") |
                                          (Martel.Str('"') +
                                           Martel.Opt(Martel.Re("\w")) +
                                           Martel.Str('"')) |
                                          location_ref))

location = Martel.Group("location",
                        Martel.Rep1(blank_space +
                                    Martel.Rep1(location_part +
                                                Martel.AnyEol())))
"""

location = Martel.Group("location",
                      Martel.ToEol("feature_location") + \
                      Martel.Rep(qualifier_space + \
                                 Martel.Re("(?!/)") + \
                                 Martel.ToEol("feature_location")))



feature_key_line = Martel.Group("feature_key_line",
                                big_indent_space +
                                feature_key +
                                location)

# -- now set up all of the info we can have for qualifiers
# a listing of valid qualifier keys
# For now use a simple list and get all of the matching text.
# In the future could allow more specific matches for each key.
# XXX Because of a bug in Martel, if one name is a substring of
# XXX another, put the longer name first.  Eg, "clone_lib" before "clone"
feature_qualifier_names = (
    "allele",         # Name of the allele for the a given gene
    "anticodon",      # Location of the anticodon of tRNA and the amino
                      #   acid for which it codes
    "bound_moiety",   # Moiety bound
    "cell_line",      # Cell line from which the sequence was obtained
    "cell_type",      # Cell type from which the sequence was obtained
    "chromosome",     # Chromosome (e.g. Chromosome number) from which
                      #   the sequence was obtained
    "chloroplast",    # Organelle type from which the sequence was obtained
    "chromoplast",    # Organelle type from which the sequence was obtained
    "citation",       # Reference to a citation providing the claim of or
                      #   evidence for a feature
    "clone_lib",      # Clone library from which the sequence was obtained
    "clone",          # Clone from which the sequence was obtained
    "codon_start",    # Indicates the first base of the first complete codon
                      #   in a CDS (as 1 or 2 or 3)
    "codon",          # Specifies a codon that is different from any found
                      #   in the reference genetic code
    "cons_splice",    # Identifies intron splice sites that do not conform to
                      #   the 5'-GT... AG-3' splice site consensus
    "country",        # Country of origin for DNA sample, intended for
                      #   epidemiological or population studies.
    "cultivar",       # Variety of plant from which sequence was obtained
    "cyanelle",       # Organelle type from which the sequence was obtained
    "db_xref",        # A database cross-reference; pointer to related
                      #   information in another database. A description of
                      #   all cross-references can be found at:
                      #   http://www.ncbi.nlm.nih.gov/collab/db_xref.html
    "dev_stage",      # If the sequence was obtained from an organism in
                      #   a specific developmental stage, it is specified
                      #   with this qualifier
    "direction",      # Direction of DNA replication
    "EC_number",      # Enzyme Commission number for the enzyme product
                      #   of the sequence
    "evidence",       # Value indicating the nature of supporting evidence
    "exception",      # Indicates that the amino acid or RNA sequence
                      #   will not translate or agree with the DNA sequence
                      #   according to standard biological rules
    "focus",          # Defines the source feature of primary biological
                      #   interest for records that have multiple source
                      #   features originating from different organisms
    "frequency",      # Frequency of the occurrence of a feature
    "function",       # Function attributed to a sequence
    "gene",           # Symbol of the gene corresponding to a sequence
                      #   region (usable with all features)
    "germline",       # If the sequence shown is DNA and a member of the
                      #   immunoglobulin family, this qualifier is used to
                      #   denote that the sequence is from unrearranged DNA
    "haplotype",      # Haplotype of organism from which the sequence was
                      #   obtained
    "insertion_seq",  # Insertion sequence element from which the sequence
                      #   was obtained
    "isolate",        # Individual isolate from which the sequence was
                      #   obtained
    "kinetoplast",    # Organelle type from which the sequence was obtained
    "label",          # A label used to permanently identify a feature
    "lab_host",       # Laboratory host used to propagate the organism
                      #   from which the sequence was obtained
    "macronuclear",   # If the sequence shown is DNA and from an organism
                      #   which undergoes chromosomal differentiation
                      #   between macronuclear and micronuclear stages,
                      #   this qualifier is used to denote that the
                      #   sequence is from macronuclear DNA.
    "map",            # Map position of the feature in free-format text
    "mitochondrion",  # Organelle type from which the sequence was obtained
    "mod_base",       # Abbreviation for a modified nucleotide base
    "note",           # Any comment or additional information
    "number",         # A number indicating the order of genetic elements
                      #   (e.g., exons or introns) in the 5 to 3 direction
    "organelle",      # Type of membrane-bound intracellular structure from
                      #  which the sequence was obtained
    "organism",       # Name of the organism that is the source of the
                      #   sequence data in the record. 
    "partial",        # Differentiates between complete regions and
                      #   partial ones
    "PCR_conditions", # Description of reaction conditions and components
                      #   for PCR
    "phenotype",      # Phenotype conferred by the feature
    "plasmid",        # Name of plasmid from which sequence was obtained
    "pop_variant",    # Population variant from which the sequence was obtained
    "product",        # Name of a product encoded by a coding region (CDS)
                      #   feature
    "protein_id",     # Protein Identifier, issued by International
                      #   collaborators.  This qualifier consists of a stable
                      #   ID portion (3+5 format with 3 position letters and
                      #   5 numbers) plus a version number after the decimal
                      #   point
    "proviral",       # If the sequence shown is viral and integrated into
                      #   another organism's genome, this qualifier is used
                      #   to denote that.
    "pseudo",         # Indicates that this feature is a non-functional
                      #   version of the element named by the feature key
    "rearranged",     # If the sequence shown is DNA and a member of the
                      #   immunoglobulin family, this qualifier is used to
                      #   denote that the sequence is from rearranged DNA
    "replace",        # Indicates that the sequence identified a feature's
                      #   intervals is replaced by the  sequence shown in
                      #   "text"
    "rpt_family",     # Type of repeated sequence; Alu or Kpn, for example
    "rpt_type",       # Organization of repeated sequence
    "rpt_unit",       # Identity of repeat unit that constitutes a
                      #   repeat_region
    "sequenced_mol",  # Molecule from which the sequence was obtained
    "serotype",       # Variety of a species (usually bacteria or virus)
                      #   characterized by its antigenic properties
    "sex",            # Sex of the organism from which the sequence
                      #   was obtained
    "specific_host",  # Natural host from which the sequence was obtained
    "specimen_voucher", # An identifier of the individual or collection
                      #   of the source organism and the place where it
                      #   is currently stored, usually an institution.
    "standard_name",  # Accepted standard name for this feature
    "strain",         # Strain from which the sequence was obtained
    "sub_clone",      # Sub-clone from which the sequence was obtained
    "sub_species",    # Sub-species name of organism  from which the
                      #   sequence was obtained
    "sub_strain",     # Sub_strain from which the sequence was obtained
    "tissue_lib",     # Tissue library from which the sequence was obtained
    "tissue_type",    # Tissue type from which the sequence was obtained
    "translation",    # Amino acid translation of a coding region
    "transl_except",  # Translational exception: single codon, the
                      #   translation of which does not conform to the
                      #   reference genetic code
    "transl_table",   # Definition of genetic code table used if other
                      #   than universal genetic code table
    "transposon",     # Transposable element from which the sequence
                      #   was obtained
    "type",           # Name of a strain if different from that in the
                      #   SOURCE field  (XXX not in
                      #   http://www.ncbi.nlm.nih.gov/collab/FT/index.html )
    "usedin",         # Indicates that feature is used in a compound
                      #   feature in another entry
    "variety",        # Variety from which sequence was obtained
    "virion",         # Viral genomic sequence as it is encapsidated
                      #   (distinguished from its proviral form integrated
                      #   in a host cell's chromosome) 
    
)

feature_qualifiers = map(Martel.Str, feature_qualifier_names)

qualifier_key = Martel.Group("qualifier_key",
                             Martel.Opt(blank_space) +
                             Martel.Str("/") +
                             Martel.Alt(*feature_qualifiers) +
                             Martel.Opt(Martel.Str("=")))

# this fails on really annoying records that have / in the first
# line not signalling a keyword.
# qualifier_value = Martel.Group("qualifier_value",
#                               Martel.ToEol() +
#                               Martel.Rep(qualifier_space +
#                                          Martel.AnyBut("/") +
#                                          Martel.ToEol()))

qualifier_value = Martel.Group("qualifier_value",
                               Martel.ToEol() +
                               Martel.Rep(qualifier_space +
                                          (Martel.AnyBut("/") |
                                           (Martel.Str("/") +
                                            Martel.Opt(blank_space) +
                                            Martel.Re("[\w\-\,]+") +
                                            blank_space)) + 
                                          Martel.ToEol()))


qualifier = Martel.Group("qualifier",
                         qualifier_key +
                         qualifier_value)
feature = Martel.Group("feature",
                       feature_key_line +
                       Martel.Rep(qualifier))


# BASE COUNT    28300 a  15069 c  15360 g  27707 t
base_count = Martel.Group("base_count",
                          Martel.Re("[\w\d ]+"))
base_count_line = Martel.Group("base_count_line",
                               Martel.Str("BASE COUNT") +
                               blank_space +
                               base_count +
                               Martel.AnyEol())

# ORIGIN      
#       1 ggacaaggcc aaggatgctg ctgctgcagc tggagcttcc gcgcaacaag taaacagata
origin_line = Martel.Group("origin_line",
                           Martel.Str("ORIGIN") +
                           (Martel.ToEol("origin_name") |
                            Martel.AnyEol()))

base_number = Martel.Group("base_number",
                           Martel.Re("[\d]+"))
sequence = Martel.Group("sequence",
                        Martel.ToEol())
sequence_line = Martel.Group("sequence_line",
                             blank_space +
                             Martel.Opt(base_number) +
                             sequence)

sequence_entry = Martel.Group("sequence_entry",
                              origin_line +
                              Martel.Rep1(sequence_line))

# all done!
# //
record_end = Martel.Group("record_end",
                          Martel.Str("//") +
                          Martel.Rep1(Martel.AnyEol()))

record = Martel.Group("genbank_record",
                      locus_line + \
                      definition_block + \
                      accession_block + \
                      Martel.Opt(nid_line) + \
                      Martel.Opt(version_line) + \
                      keywords_block + \
                      Martel.Opt(segment_line) + \
                      source_block + \
                      organism_block + \
                      Martel.Rep(reference) + \
                      Martel.Opt(comment_block) + \
                      features_line + \
                      Martel.Rep1(feature) + \
                      base_count_line + \
                      sequence_entry + \
                      record_end)

record_format = Martel.ParseRecords("genbank_file", record,
                                    RecordReader.StartsWith, ("LOCUS",) )


# if you download a big mess of GenBank files, it'll have a header
# in that case you should be using 'format' instead of the standard
# 'record'
header = Martel.Re("""\
(?P<filename>[^ ]+) +Genetic Sequence Data Bank
 *(?P<release_day>\d+) (?P<release_month>\w+) (?P<release_year>\d+)

 *(?P<data_bank_name>[^\R]+)

 *(?P<data_bank_name>[^\R]+)

 *(?P<num_loci>\d+) loci, *(?P<num_bases>\d+) bases, from *(?P<num_reports>\d+) reported sequences


""")

format = Martel.HeaderFooter("genbank",
                             header, RecordReader.CountLines, (10,),
                             record, RecordReader.EndsWith, ("//",),
                             None, None, None,
                             )

multirecord = Martel.ParseRecords("genbank", record,
                                  RecordReader.EndsWith, ("//",))


                          



