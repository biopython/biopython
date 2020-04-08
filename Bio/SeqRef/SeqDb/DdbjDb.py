from .NcbiDb import NcbiNucleotideDb


class DdbjDb(NcbiNucleotideDb):
    name = "DDBJ"
    base_url = "http://getentry.ddbj.nig.ac.jp"
    entry_url = "http://getentry.ddbj.nig.ac.jp/getentry/na/{}"
