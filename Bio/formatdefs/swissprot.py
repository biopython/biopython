# Define the various SWISS-PROT formats

from Bio import register_format, link_format

import sequence

register_format(
    name = "swissprot",
#    filter = "Bio.expressions.swissprot.filter",
)

register_format(
    name = "swissprot/38",
    abbrev = "sprot38",
    expression = "Bio.expressions.swissprot.sprot38.format",
)

register_format(
    name = "swissprot/39",
    abbrev = "sprot39",
    expression = "Bio.expressions.swissprot.sprot38.format",
)


link_format("swissprot", "swissprot/38")
#link_format("swissprot", "swissprot/39", before = "swissprot/38")

#link_format("swissprot", "swissprot/40")  # default inserts at end

link_format("sequence", "swissprot")
