from ._commons import _SeqDb


class EbiEnaDB(_SeqDb):
    name = "EBI-ENA"
    base_url = "https://www.ebi.ac.uk/ena"
    entry_url = "old: https://www.ebi.ac.uk/ena/data/view/{0} new: https://www.ebi.ac.uk/ena/browser/view/{0}"
