from Bio.config.FormatRegistry import FormatObject, FormatGroup


# This is a sequence of genbank records
genbank_records = FormatObject(
    name = "genbank-records",
    abbrev = "genbank_records",
    expression = "Bio.expressions.genbank.multirecord",
    )

# This is the format as released by NCBI.  It includes
# the header information.
genbank_release = FormatObject(
    name = "genbank-release",
    abbrev = "genbank_release",
    expression = "Bio.expressions.genbank.format",
    )


genbank = FormatGroup(
    name = "genbank",
    )
genbank.add(genbank_records)
genbank.add(genbank_release)

from Bio.formatdefs import sequence
sequence.sequence.add(genbank)
