from Bio.config.FormatRegistry import FormatObject, FormatGroup

embl65 = FormatObject(
    name = "embl/65",
    abbrev = "embl65",
    expression = "Bio.expressions.embl.embl65.format",
    )

embl = FormatGroup(
    name = "embl",
    )
embl.add(embl65)

from Bio.formatdefs import sequence
sequence.sequence.add(embl)
