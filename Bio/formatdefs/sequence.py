from Bio.config.FormatRegistry import FormatObject, FormatGroup

fasta = FormatObject(
    name="fasta",
    expression="Bio.expressions.fasta.format",
    )

empty = FormatGroup(
    name="empty"
    )

sequence = FormatGroup(
    name="sequence",
    )
sequence.add(fasta)
