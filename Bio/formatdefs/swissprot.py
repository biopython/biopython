# Define the various SWISS-PROT formats

from Bio.config.FormatRegistry import FormatObject, FormatGroup


sprot38 = FormatObject(
    name = "swissprot/38",
    abbrev = "sprot38",
    expression = "Bio.expressions.swissprot.sprot38.format",
)

sprot40 = FormatObject(
    name = "swissprot/40",
    abbrev = "sprot40",
    expression = "Bio.expressions.swissprot.sprot40.format",
)

ipi = FormatObject(
        name = "ipi",
        abbrev = "ipi",
        expression = "Bio.expressions.swissprot.ipi.format"
)


swissprot = FormatGroup(
    name = "swissprot",
#    filter = "Bio.expressions.swissprot.filter",
)
swissprot.add(sprot38)
swissprot.add_before(sprot40, sprot38)
swissprot.add_before(ipi, sprot40)

from Bio.formatdefs import sequence
sequence.sequence.add(swissprot)
