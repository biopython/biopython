from Bio import register_format, link_format

register_format(
    name = "genbank",
    )

# This is a sequence of genbank records
register_format(
    name = "genbank-records",
    abbrev = "genbank_records",
    expression = "Bio.expressions.genbank.multirecord",
    )

# This is the format as released by NCBI.  It includes
# the header information.
register_format(
    name = "genbank-release",
    abbrev = "genbank_release",
    expression = "Bio.expressions.genbank.format",
    )


link_format("sequence", "genbank")
link_format("genbank", "genbank-release")
link_format("genbank", "genbank-release")

