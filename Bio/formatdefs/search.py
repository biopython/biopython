from Bio import register_format, link_format

register_format(
    name = "search",
    multirecord = 0,
)

register_format(
    name = "blast",
    filter = "Bio.expressions.blast.ncbiblast.blast_filter",
    multirecord = 0,
)
link_format("search", "blast")

register_format(
    name = "blastp",
    multirecord = 0,
)
link_format("blast", "blastp")

register_format(
    name = "blastn",
    multirecord = 0,
)
link_format("blast", "blastn")

register_format(
    name = "blastx",
    multirecord = 0,
)
link_format("blast", "blastx")

register_format(
    name = "tblastn",
    multirecord = 0,
)
link_format("blast", "tblastn")

register_format(
    name = "tblastx",
    multirecord = 0,
)
link_format("blast", "tblastx")

### NCBI BLASTs
register_format(
    name = "ncbi-blastp",
    abbrev = "ncbi_blastp",
    expression = "Bio.expressions.blast.ncbiblast.blastp",
    filter = "Bio.expressions.blast.ncbiblast.blastp_appheader",
    multirecord = 0,
)
link_format("blastp", "ncbi-blastp")

register_format(
    name = "ncbi-blastn",
    abbrev = "ncbi_blastn",
    expression = "Bio.expressions.blast.ncbiblast.blastn",
    filter = "Bio.expressions.blast.ncbiblast.blastn_appheader",
    multirecord = 0,
)
link_format("blastn", "ncbi-blastn")

register_format(
    name = "ncbi-blastx",
    abbrev = "ncbi_blastx",
    expression = "Bio.expressions.blast.ncbiblast.blastx",
    filter = "Bio.expressions.blast.ncbiblast.blastx_appheader",
    multirecord = 0,
)
link_format("blastx", "ncbi-blastx")

register_format(
    name = "ncbi-tblastn",
    abbrev = "ncbi_tblastn",
    expression = "Bio.expressions.blast.ncbiblast.tblastn",
    filter = "Bio.expressions.blast.ncbiblast.tblastn_appheader",
    multirecord = 0,
)
link_format("tblastn", "ncbi-tblastn")

register_format(
    name = "ncbi-tblastx",
    abbrev = "ncbi_tblastx",
    expression = "Bio.expressions.blast.ncbiblast.tblastx",
    filter = "Bio.expressions.blast.ncbiblast.tblastx_appheader",
    multirecord = 0,
)
link_format("tblastx", "ncbi-tblastx")


### WU-BLASTs
register_format(
    name = "wu-blastp",
    abbrev = "wu_blastp",
    expression = "Bio.expressions.blast.wublast.blastp",
    filter = "Bio.expressions.blast.wublast.blastp_appheader",
    multirecord = 0,
)
link_format("blastp", "wu-blastp")

register_format(
    name = "wu-blastn",
    abbrev = "wu_blastn",
    expression = "Bio.expressions.blast.wublast.blastn",
    filter = "Bio.expressions.blast.wublast.blastn_appheader",
    multirecord = 0,
)
link_format("blastn", "wu-blastn")

register_format(
    name = "wu-blastx",
    abbrev = "wu_blastx",
    expression = "Bio.expressions.blast.wublast.blastx",
    filter = "Bio.expressions.blast.wublast.blastx_appheader",
    multirecord = 0,
)
link_format("blastx", "wu-blastx")

