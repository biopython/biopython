from Bio import register_format, link_format

import sequence

register_format(
    name = "genbank",
    expression = "Bio.GenBank.genbank_format.record_format"
    )

link_format("sequence", "genbank")
