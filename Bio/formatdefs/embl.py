from Bio import register_format, link_format

import sequence

register_format(
    name = "embl",
    )

register_format(
    name = "embl/65",
    abbrev = "embl65",
    expression = "Bio.expressions.embl.embl65.format",
    )

link_format("embl", "embl/65")
link_format("sequence", "embl")
