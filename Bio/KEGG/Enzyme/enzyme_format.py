"""Martel based parser to read KEGG Ligand/Enzyme files.

This is a huge regular expression for KEGG Ligand/Enzyme,
built using the 'regular expressions on steroids' capabilities of
Martel.

A description of the format can be found in the 'ligand.doc' file
from the Ligand distribution, available from:

 http://www.genome.ad.jp/kegg

(Note that as of KEGG release 19.0, the 'enzyme' file that comes
as part of the Ligand distribution does not conform 100% to this
description. In particular, some entries contain data items that
are not ordered as they should be according the description. This
Biopython format description has been weakened to accomodate parsing
of these irregular entries.)
 
"""

# Martel
from Martel import *
from Martel import RecordReader

# --- First set up some helper constants and functions
INDENT = 12

blank_space = Rep1(Str1(" "))
point = Str1(".")

def _define_block_seq(block_tag, block_group, item):
    """Define a Martel grouping that can parse a block of groups.

    Many of the KEGG entries are formatted as:

    TAG         blah blah blah blah blah blah blah blah
                blah blah
                zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

    Where 'blah blah...' and 'zzzz...' are identically formatted
    items that can wrap over multiple lines. This function defines
    a consistent grouping for these blocks.

    Arguments:
    o block_tag - A string of length <= INDENT giving the block tag
                  (i.e. 'TAG' in the example above).
    o block_group - A callbakc tag for the entire block.
    o item - A Martel group/expression for a single instance of the
             data item in this block. (Must NOT include the final
             newline character).
    """
    return Group(block_group,
                 Str1(block_tag + " " * (INDENT - len(block_tag))) +
                 item +
                 Rep(AnyEol() +
                     Str1(" " * INDENT) + item) +
                 AnyEol())


# --- Create regular expressions for each data item

# ENTRY
entry_tag = "ENTRY"
entry = Group("entry",
              Str1("EC ") +
              Rep(Rep1(Integer()) + point) +
              Rep1(Integer()))

entry_line = Group("entry_line",
                   Str1(entry_tag + " " * (INDENT - len(entry_tag))) +
                   entry +
                   AnyEol())


# NAME
name_tag = "NAME"
partial_name = Re("[\S]+")
composite_name = Rep1(Opt(Str1(" ")) + partial_name + Opt(Str1(" ")))
wrapped_name = composite_name + Rep(AnyEol() +
                                    Str1(" " * INDENT + "$") +
                                    composite_name)

name = Group("name", wrapped_name)
name_block = _define_block_seq("NAME", "name_block", name)


# CLASS (renamed 'classname' to avoid conflict with Python)
class_name = Group("classname", composite_name)
class_block = _define_block_seq("CLASS", "class_block", class_name)


# SYSNAME
sysname = Group("sysname", wrapped_name)
sysname_block = _define_block_seq("SYSNAME", "sysname_block", sysname)


# REACTION
#
# Note: It would be nice to do something fancier with this block,
#       like doing a complete parse of the reaction equations.
#       Unfortunately this would'nt work in general as the KEGG
#       format specification allows this field to be an arbitrary
#       string.
#
reaction_block = _define_block_seq("REACTION",
                                   "reaction_block",
                                   Group("reaction", Re(".*")))


# SUBSTRATE
substrate = Group("substrate", wrapped_name)
substrate_block = _define_block_seq("SUBSTRATE",
                                    "substrate_block",
                                    substrate)


# PRODUCT
product = Group("product", wrapped_name)
product_block = _define_block_seq("PRODUCT",
                                  "product_block",
                                  product)


# INHIBITOR
inhibitor = Group("inhibitor", wrapped_name)
inhibitor_block = _define_block_seq("INHIBITOR",
                                    "inhibitor_block",
                                    inhibitor)


# COFACTOR
cofactor = Group("cofactor", wrapped_name) 
cofactor_block = _define_block_seq("COFACTOR",
                                   "cofactor_block",
                                   cofactor)


# EFFECTOR
effector = Group("effector", wrapped_name)
effector_block = _define_block_seq("EFFECTOR",
                                   "effector_block",
                                   effector)


# COMMENT
comment_block = _define_block_seq("COMMENT",
                                  "comment_block",
                                  Group("comment", Re(".*")))


# PATHWAY
pathway_db = Group("pathway_db", Re("[\w]+"))
pathway_id = Group("pathway_id", Re("[\w]+"))
pathway_desc = Group("pathway_desc", Re(".*"))
pathway_line = pathway_db + Str1(":") + \
               blank_space + \
               pathway_id + \
               blank_space + \
               pathway_desc + \
               Rep(AnyEol() +
                   Str1(" " * (INDENT + 1)) +
                   blank_space +
                   pathway_desc)

pathway_block = _define_block_seq("PATHWAY", "pathway_block", pathway_line)


# GENES
organism = Group("organism", Re("[\w]{3}"))
gene_id = Group("gene_id", Re("[\S]+"))
genes = organism + \
        Str1(":") + \
        Rep1(blank_space + gene_id) + \
        Rep(AnyEol() +
            Str1(" " * (INDENT + 1)) +
            Rep1(blank_space + gene_id))

genes_block = _define_block_seq("GENES", "genes_block", genes)


# DISEASE
disease_db = Group("disease_db", Re("[\w]+"))
disease_id = Group("disease_id", Re("[\w]+"))
partial_desc = Group("disease_desc", Re(".*"))
disease_desc = partial_desc + \
               Rep(AnyEol() +
                   Str1(" " * (INDENT + 1)) +
                   blank_space +
                   partial_desc)
disease_line = disease_db + Str1(":") + \
               blank_space + \
               disease_id + \
               blank_space + \
               disease_desc

disease_block = _define_block_seq("DISEASE", "disease_block", disease_line)


# MOTIF
partial_motif = Re("\S.*")
wrapped_motif = partial_motif + \
                Rep(AnyEol() +
                    Str1(" " * (INDENT + 1)) +
                    blank_space +
                    partial_motif)

motif_db = Group("motif_db", Re("[\w]+"))
motif_id = Group("motif_id", Re("[\w]+"))
motif = Group("motif", wrapped_motif)
motif_line = motif_db + Str1(":") + \
             blank_space + \
             motif_id + \
             blank_space + \
             motif

motif_block = _define_block_seq("MOTIF", "motif_block", motif_line)


# STRUCTURES
structure_db = Group("structure_db", Re("[\w]+"))
structure_id = Group("structure_id", Re("[\w]+"))
structures_line = structure_db + \
                  Str1(":") + \
                  blank_space + \
                  Rep1(structure_id + blank_space) + \
                  Rep(AnyEol() +
                      Str1(" " * (INDENT + 1)) +
                      blank_space + 
                      Rep1(structure_id + blank_space))

structures_block = _define_block_seq("STRUCTURES",
                                     "structures_block",
                                     structures_line)


# DBLINKS
dblinks_db = Group("dblinks_db", Rep1(AnyBut(":")))
dblinks_id = Group("dblinks_id", Re("[\S]+"))
dblinks_line = dblinks_db + \
               Str1(":") + \
               blank_space + \
               Rep1(dblinks_id + Opt(blank_space)) + \
               Rep(AnyEol() +
                   Str1(" " * (INDENT + 1)) +
                   blank_space +
                   Rep1(dblinks_id + blank_space))

dblinks_block = _define_block_seq("DBLINKS", "dblinks_line", dblinks_line)


# --- Assemble the full record format

record_end = Group("record_end",
                   Str1("///") +
                   Rep1(AnyEol()))

# The irregularities in the following format are present to avoid
# conflict with a few inconsistently formatted records in the KEGG
# database.
record = Group("KEGG_Enzyme_record",
               entry_line +
               Rep1(name_block) +
               Rep(Alt(class_block,
                       sysname_block,
                       reaction_block,
                       substrate_block,
                       product_block,
                       inhibitor_block,
                       cofactor_block,
                       effector_block,
                       comment_block,
                       pathway_block,
                       genes_block,
                       disease_block,
                       motif_block,
                       structures_block,
                       dblinks_block)) +
               record_end)

record_format = ParseRecords("KEGG_Enzyme_file", {},
                             record, RecordReader.EndsWith, ("///",))


