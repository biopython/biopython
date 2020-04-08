from ._commons import SeqDb


class EbiDB(SeqDb):
    name = "EBI-ENA"
    base_url = "https://www.ebi.ac.uk/ena"
    entry_url = "https://www.ebi.ac.uk/ena/data/view/{}"
