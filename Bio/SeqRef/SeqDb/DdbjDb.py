from .NcbiDb import NcbiNucleotideDb


class DdbjDb(NcbiNucleotideDb):
    name = "DDBJ"
    base_url = "http://getentry.ddbj.nig.ac.jp"
    entry_url = "http://getentry.ddbj.nig.ac.jp/getentry/na/{}"
    fetch_url = "http://getentry.ddbj.nig.ac.jp/getentry/na/{id}/?format={format}"
    fetch_file_format_map = dict(fasta="fasta", ddbj="flatfile")
    fetch_file_format_default = "ddbj"
