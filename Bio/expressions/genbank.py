"""Martel based parser to read GenBank formatted files.

This is a huge regular regular expression for GenBank, built using
the 'regular expressions on steroids' capabilities of Martel.

Documentation for GenBank format that I found:

o GenBank/EMBL feature tables are described at:
http://www.ebi.ac.uk/embl/Documentation/FT_definitions/feature_table.html

o There are also descriptions of different GenBank lines at:
http://www.ibc.wustl.edu/standards/gbrel.txt
"""
# Martel
import Martel
from Martel import RecordReader

# identify certain items as important for format converters
from Bio import Std

# --- first set up some helper constants and functions

# - useful constants for dealing with the blank space in GenBank documents
# this is useful since blank space can be significant in GenBank flat files.
INDENT = 12
FEATURE_KEY_INDENT = 5
FEATURE_QUALIFIER_INDENT = 21

blank_space = Martel.Spaces()
small_indent_space = Martel.Str(" " * 2)
big_indent_space = Martel.Str(" " * FEATURE_KEY_INDENT)
qualifier_space = Martel.Str(" " * FEATURE_QUALIFIER_INDENT) | \
                  Martel.Str("\t" + " " * (FEATURE_QUALIFIER_INDENT - 8))

# - useful functions
def define_block(identifier, block_tag, block_data, std_block_tag = None,
                 std_tag = None):
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
    o std_block_tag - A Bio.Std Martel tag used to register the entire
    block as having being a "standard" type of information.
    o std_tag - A Bio.Std Martel tag used to register just the information
    in the block as being "standard"
    """
    diff = INDENT - len(identifier)
    assert diff > 0, diff
   
    # if no std_tag info is defined, just make std_tag a no-op function
    if std_tag is None:
        def do_nothing(martel_info):
            return martel_info
        std_tag = do_nothing

    identifier_and_text = Martel.Str(identifier) + \
                          Martel.Rep(Martel.Str(" ")) + \
                          std_tag(Martel.UntilEol(block_data)) + \
                          Martel.AnyEol()
    indented_text = Martel.Str(" "*INDENT) + \
                    std_tag(Martel.UntilEol(block_data)) + \
                    Martel.AnyEol()
    block_info = Martel.Group(
        block_tag,
        identifier_and_text +
        Martel.Rep(Martel.Alt(Martel.AnyEol(), indented_text))
        )
    # tag the info as some standard Martel element if specified
    if std_block_tag is not None:
        block_info = std_block_tag(block_info)

    return block_info


# first line
# LOCUS       AC007323    86436 bp    DNA             PLN       19-JAN-2000
locus = Martel.Group("locus",
                     Martel.Re(r"[\w\-]+"))
size = Martel.Group("size",
                    Martel.Rep1(Martel.Integer()))

# deal with the different kinds of residues we can have
valid_residue_prefixes = ["ss-", "ds-", "ms-"]
valid_residue_types = ["DNA", "RNA", "mRNA", "tRNA", "rRNA", "uRNA",
                       "scRNA", "snRNA", "snoRNA", "PROTEIN"]

residue_prefixes = map(Martel.Str, valid_residue_prefixes)
residue_types = map(Martel.Str, valid_residue_types)

residue_type = Martel.Group("residue_type",
                            Martel.Opt(Martel.Alt(*residue_prefixes)) +
                            Martel.Opt(Martel.Alt(*residue_types)) +
                            Martel.Opt(Martel.Opt(blank_space) +
                                       Martel.Alt(Martel.Str("circular"),
                                                  Martel.Str("linear"))))

date = Martel.Group("date",
                    Martel.Re("[-\w]+"))

# the PLN, etc stuff indicates data file divisions
valid_divisions = ["PRI", "ROD", "MAM", "VRT", "INV", "PLN", "BCT", "RNA",
                   "VRL", "PHG", "SYN", "UNA", "EST", "PAT", "STS", "GSS",
                   "HTG", "HTC", "CON", "ENV"]
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
definition_block = define_block("DEFINITION", "definition_block",
                                "definition", Std.description_block,
                                Std.description)

# accession line
# ACCESSION   AC007323
# or
# ACCESSION   NC_004353 REGION: 1..1281640
accession = Martel.Group("accession",
                         Martel.Re("[\w]+"))

region = Martel.Group("region",
                      Martel.Re("[\d]+..[\d]+"))

accession_block = Martel.Group("accession_block",
                               Martel.Str("ACCESSION") +
                               Martel.Rep1(blank_space +
                                   Martel.Rep1(accession +
                                       Martel.Opt(
                                           Martel.Opt(Martel.Str(" ")) +
                                           Martel.Str("REGION:") +
                                           Martel.Opt(Martel.Str(" ")) +
                                           region) +
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

# PID         g6754304
pid = Martel.Group("pid",
                   Martel.Re("[\w\d]+"))
pid_line = Martel.Group("pid_line",
                        Martel.Str("PID") +
                        blank_space +
                        pid +
                        Martel.AnyEol())

# version and GI line
# VERSION     AC007323.5  GI:6587720
version = Martel.Group("version",
                       Std.dbid(Martel.Re("[\w\d\.]+"),
                                {"type" : "primary", "dbname" : "genbank"}))

gi = Martel.Group("gi",
                  Std.dbid(Martel.Re("[\d]+"),
                           {"type" : "secondary", "dbname" : "genbank"}))

version_line = Martel.Group("version_line",
                            Martel.Str("VERSION") +
                            blank_space +
                            version +
                            Martel.Opt(blank_space +
                                       Martel.Str("GI:") +
                                       gi) +
                            Martel.AnyEol())

# DBSOURCE    REFSEQ: accession NM_010510.1
db_source_block = define_block("DBSOURCE", "db_source_block", "db_source")

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
consrtm_block = define_block("  CONSRTM", "consrtm_block", "consrtm")
title_block   = define_block("  TITLE",   "title_block",   "title")
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
                         Martel.Opt(authors_block) +
                         Martel.Opt(consrtm_block) +
                         Martel.Opt(title_block) +
                         journal_block +
                         Martel.Opt(medline_line) +
                         Martel.Opt(pubmed_line) +
                         Martel.Opt(remark_block))

# COMMENT     On Dec 16, 1999 this sequence version replaced gi:5729683.
comment_block = define_block("COMMENT", "comment_block", "comment")
# PRIMARY
primary_line = Martel.Group("primary_line",
                            Martel.Str("PRIMARY") +
                            blank_space +
                            Martel.Str("TPA_SPAN") +
                            blank_space +
                            Martel.Str("PRIMARY_IDENTIFIER") +
                            blank_space +
                            Martel.Str("PRIMARY_SPAN") +
                            blank_space +
                            Martel.Str("COMP") +
                            Martel.ToEol())

primary_ref_line =Martel.Group("primary_ref_line",
                               blank_space +
                               Martel.Re(r"\d+\-\d+") +
                               blank_space +
                               Martel.Re("[\S]+") +
                               blank_space +
                               Martel.Re("\d+\-\d+")+
                               Martel.Opt(blank_space +  Martel.Str("c"))+
                               Martel.ToEol())
                              
primary =  Martel.Group("primary",primary_line +
                        Martel.Rep1(primary_ref_line))
                                       
# start on the feature table. Eeek -- This is the part I was afraid of
# most!

# the header, so that we know we are heading into some features
# FEATURES             Location/Qualifiers
features_line = Martel.Group("features_line",
                             Martel.Str("FEATURES") +
                             blank_space +
                             Martel.Str("Location/Qualifiers") +
                             Martel.AnyEol())

# feature key names are basically words, but sometimes have additional
# characters ("-10_signal", "3'UTR"...)
feature_key = Martel.Group("feature_key",
                            Martel.Re("[\w'-]+"))
"""
location = Martel.Group("location",
                      Martel.ToEol("feature_location") + \
                      Martel.Rep(qualifier_space + \
                                 Martel.Re("(?!/)") + \
                                 Martel.ToEol("feature_location")))
"""

location = Martel.Group("location",
                        Std.feature_location(Martel.UntilEol()) +
                        Martel.AnyEol() +
                        Martel.Rep(qualifier_space +
                                   Martel.AssertNot(Martel.Str("/")) +
                                   Std.feature_location(Martel.UntilEol()) +
                                   Martel.AnyEol())
                        )

feature_key_line = Martel.Group("feature_key_line",
                                big_indent_space +
                                Std.feature_name(feature_key) +
                                location)

# qualifiers escape quotes using double quotes
quote = Martel.Str('"')
quoted_chars = Std.feature_qualifier_description(Martel.Re(r'([^"\R]|"")*'))

quoted_string = (quote + quoted_chars +
                 Martel.Rep(Martel.AnyEol() + qualifier_space + quoted_chars) +
                 quote + Martel.AnyEol())

unquoted_string = Martel.AssertNot(quote) + \
                  Std.feature_qualifier_description(Martel.UntilEol()) + \
                  Martel.AnyEol()

qualifier = Std.feature_qualifier(
    qualifier_space +
    Martel.Str("/") +
    Std.feature_qualifier_name(Martel.Word("feature_qualifier_name")) +
    (Martel.AnyEol() | #  '/pseudo'
     (Martel.Str("=") +
     Martel.Group("feature_qualifier_description",
      (unquoted_string |  #  '/evidence=experimental'
       quoted_string))))   #  '/translation="AAAAAAAA....
                          #   AAAAAAAAAAAAAAAAAAAA'
    )

feature = Std.feature(feature_key_line +
                      Martel.Rep(qualifier))

feature_block = Std.feature_block(Martel.Rep1(feature),
                            {"location-style" : "genbank"})

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
sequence = Std.sequence(Martel.Group("sequence",
                        Martel.Re("[\w]+")))
sequence_plus_spaces = Martel.Group("sequence_plus_spaces",
                                    Martel.Rep1(Martel.Str(" ") +
                                    Martel.Opt(sequence)) +
                                    Martel.Opt(Martel.Str(" ")))
sequence_line = Martel.Group("sequence_line",
                             blank_space +
                             Martel.Opt(base_number) +
                             sequence_plus_spaces +
                             Martel.AnyEol())

sequence_entry = Std.sequence_block(Martel.Group("sequence_entry",
                                    origin_line +
                                    Martel.Rep1(sequence_line)))

# CONTIG
# this is the contig information for RefSeq records

contig_location = Martel.Group("contig_location",
                    Martel.ToEol("feature_location") + \
                    Martel.Rep(Martel.Str(" " * INDENT) + \
                               Martel.Re("(?!/)") + \
                               Martel.ToEol("feature_location")))

contig_block = Martel.Group("contig_block",
                            Martel.Str("CONTIG") +
                            blank_space +
                            contig_location)

# all done!
# //
record_end = Martel.Group("record_end",
                          Martel.Str("//") +
                          Martel.Rep1(Martel.AnyEol()))

record = Std.record(Martel.Group("genbank_record",
                      locus_line + \
                      definition_block + \
                      accession_block + \
                      Martel.Opt(nid_line) + \
                      Martel.Opt(pid_line) + \
                      Martel.Opt(version_line) + \
                      Martel.Opt(db_source_block) + \
                      keywords_block + \
                      Martel.Opt(segment_line) + \
                      source_block + \
                      organism_block + \
                      Martel.Rep(reference) + \
                      Martel.Opt(primary) +\
                      Martel.Opt(comment_block) + \
                      features_line + \
                      feature_block + \
                      Martel.Alt(Martel.Opt(base_count_line) +
                                 sequence_entry,
                                 contig_block) + \
                      record_end))

# if you download a big mess of GenBank files, it'll have a header
# in that case you should be using 'ncbi_format' instead of the standard
# 'format'
header = Martel.Re("""\
(?P<filename>[^ ]+) +Genetic Sequence Data Bank
 *(?P<release_day>\d+) (?P<release_month>\w+) (?P<release_year>\d+)

 *(?P<data_bank_name>[^\R]+)

 *(?P<data_bank_name>[^\R]+)

 *(?P<num_loci>\d+) loci, *(?P<num_bases>\d+) bases, from *(?P<num_reports>\d+) reported sequences


""")

ncbi_format = Martel.HeaderFooter("genbank", {"format" : "ncbi_genbank"},
                             header, RecordReader.CountLines, (10,),
                             record, RecordReader.EndsWith, ("//",),
                             None, None, None,
                             )

format = Martel.ParseRecords("genbank", {"format" : "genbank"},
                             record, RecordReader.StartsWith, ("LOCUS ",))
