from Bio import register_format, link_format

register_format(
    name = "sequence",
)

register_format(
    name = "fasta",
    expression = "Bio.expressions.fasta.format",
)

link_format("sequence", "fasta")

register_format(
    name = "empty",
)
