from Bio.config.FormatRegistry import FormatObject, FormatGroup

search = FormatGroup(
    name = "search",
    multirecord = 0,
)

blast = FormatGroup(
    name = "blast",
    filter = "Bio.expressions.blast.ncbiblast.blast_filter",
    multirecord = 0,
)
search.add(blast)

blastp = FormatGroup(
    name = "blastp",
    multirecord = 0,
)
blast.add(blastp)

blastn = FormatGroup(
    name = "blastn",
    multirecord = 0,
)
blast.add(blastn)

blastx = FormatGroup(
    name = "blastx",
    multirecord = 0,
)
blast.add(blastx)

tblastn = FormatGroup(
    name = "tblastn",
    multirecord = 0,
)
blast.add(tblastn)

tblastx = FormatGroup(
    name = "tblastx",
    multirecord = 0,
)
blast.add(tblastx)

### NCBI BLASTs
ncbi_blastp = FormatObject(
    name = "ncbi-blastp",
    abbrev = "ncbi_blastp",
    expression = "Bio.expressions.blast.ncbiblast.blastp",
    filter = "Bio.expressions.blast.ncbiblast.blastp_appheader",
    multirecord = 0,
)
blastp.add(ncbi_blastp)

ncbi_blastn = FormatObject(
    name = "ncbi-blastn",
    abbrev = "ncbi_blastn",
    expression = "Bio.expressions.blast.ncbiblast.blastn",
    filter = "Bio.expressions.blast.ncbiblast.blastn_appheader",
    multirecord = 0,
)
blastn.add(ncbi_blastn)

ncbi_blastx = FormatObject(
    name = "ncbi-blastx",
    abbrev = "ncbi_blastx",
    expression = "Bio.expressions.blast.ncbiblast.blastx",
    filter = "Bio.expressions.blast.ncbiblast.blastx_appheader",
    multirecord = 0,
)
blastx.add(ncbi_blastx)

ncbi_tblastn = FormatObject(
    name = "ncbi-tblastn",
    abbrev = "ncbi_tblastn",
    expression = "Bio.expressions.blast.ncbiblast.tblastn",
    filter = "Bio.expressions.blast.ncbiblast.tblastn_appheader",
    multirecord = 0,
)
tblastn.add(ncbi_tblastn)

ncbi_tblastx = FormatObject(
    name = "ncbi-tblastx",
    abbrev = "ncbi_tblastx",
    expression = "Bio.expressions.blast.ncbiblast.tblastx",
    filter = "Bio.expressions.blast.ncbiblast.tblastx_appheader",
    multirecord = 0,
)
tblastx.add(ncbi_tblastx)


### WU-BLASTs
wu_blastp = FormatObject(
    name = "wu-blastp",
    abbrev = "wu_blastp",
    expression = "Bio.expressions.blast.wublast.blastp",
    filter = "Bio.expressions.blast.wublast.blastp_appheader",
    multirecord = 0,
)
blastp.add(wu_blastp)

wu_blastn = FormatObject(
    name = "wu-blastn",
    abbrev = "wu_blastn",
    expression = "Bio.expressions.blast.wublast.blastn",
    filter = "Bio.expressions.blast.wublast.blastn_appheader",
    multirecord = 0,
)
blastn.add(wu_blastn)

wu_blastx = FormatObject(
    name = "wu-blastx",
    abbrev = "wu_blastx",
    expression = "Bio.expressions.blast.wublast.blastx",
    filter = "Bio.expressions.blast.wublast.blastx_appheader",
    multirecord = 0,
)
blastx.add(wu_blastx)
