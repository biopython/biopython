"""Martel based parser to read KEGG Ligand/Compound files.

This is a huge regular regular expression for KEGG Ligand/Compound,
built using the 'regular expressions on steroids' capabilities of
Martel.

A description of the format can be found in the 'ligand.doc' file
from the Ligand distribution, available from:

 http://www.genome.ad.jp/kegg
 
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
entry = Group("entry", Re("C.*"))

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


# FORMULA
formula_tag = "FORMULA"
formula = Group("formula", Re(".+"))

formula_line = Group("formula_line",
                     Str1(formula_tag + " " * (INDENT - len(formula_tag))) +
                     formula +
                     AnyEol())


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


# ENZYME
enzyme_id = Group("enzyme_id", Rep1(Integer() + Opt(Str1("."))))
enzyme_role = Group("enzyme_role", Str("R", "C", "I", "E"))
enzyme_entry = enzyme_id + \
               blank_space + \
               Str1("(") + enzyme_role + Str1(")")
enzyme_line = Rep1(enzyme_entry + Opt(blank_space))

enzyme_block = _define_block_seq("ENZYME", "enzyme_block", enzyme_line)


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
dblinks_line = Opt(dblinks_db +
                   Str1(":") +
                   blank_space) + \
               dblinks_id + \
               Rep(Opt(AnyEol() +
                       Str1(" " * (INDENT + 1))) +
                   Rep1(blank_space + dblinks_id))

dblinks_block = _define_block_seq("DBLINKS", "dblinks_line", dblinks_line)


# --- Assemble the full record format

record_end = Group("record_end",
                   Str1("///") +
                   Rep1(AnyEol()))

record = Group("KEGG_Compound_record",
               entry_line +
               name_block +
               Opt(formula_line) +
               Opt(pathway_block) +
               Opt(enzyme_block) +
               Opt(structures_block) +
               Opt(dblinks_block) +
               Opt(formula_line) +
               record_end)

record_format = ParseRecords("KEGG_Compound_file", {},
                             record, RecordReader.EndsWith, ("///",))


