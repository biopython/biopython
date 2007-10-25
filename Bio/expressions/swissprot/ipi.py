"""Expression for IPI format.

IPI is nearly swissprot, but contains some differents which makes the
Swissprot parsers choke.
"""

import warnings
warnings.warn("Bio.expressions was deprecated, as it does not work with recent versions of mxTextTools. If you want to continue to use this module, please get in contact with the Biopython developers at biopython-dev@biopython.org to avoid permanent removal of this module from Biopython", DeprecationWarning)



from Bio import Std
import Martel
from Martel import Time
import sprot40

# The ID line contains a versioned period number
ID_exp = Martel.Group("ID",
                  Martel.Str("ID   ") + \
                  Std.dbid(Martel.Group("entry_name", Martel.Re("[\w.]+")),
                      {"type": "primary", "dbname": "sp"}) + \
                  Martel.Spaces() + \
                  Martel.Word("data_class_table") + \
                  Martel.Str(";") + Martel.Spaces() + \
                  Martel.Word("molecule_type") + \
                  Martel.Str(";") + Martel.Spaces() + \
                  Martel.Digits("sequence_length") + \
                  Martel.Str(" AA.") + \
                  Martel.AnyEol()
                  )

# The DT formatted lines look different, and there is not
# a third DT line for annotations
# DT   04-MAR-2003 (IPI Human rel. 2.17, Created)
# DT   04-MAR-2003 (IPI Human rel. 2.17, Last sequence update)

DT_created_exp = (Martel.Str("DT   ") +
                  Time.make_expression("%(DD)-%(Jan)-%(YYYY)") + \
                  Martel.Str(" (IPI Human rel. ") + \
                  Martel.Float("release") + \
                  Martel.Str(", Created)") + Martel.AnyEol())

DT_seq_update_exp = (Martel.Str("DT   ") +
                  Time.make_expression("%(DD)-%(Jan)-%(YYYY)") + \
                  Martel.Str(" (IPI Human rel. ") + \
                  Martel.Float("release") + \
                  Martel.Str(", Last sequence update)") + Martel.AnyEol())

DT_ann_update_exp = (Martel.Str("DT   ") +
                  Time.make_expression("%(DD)-%(Jan)-%(YYYY)") + \
                  Martel.Str(" (IPI Human rel. ") + \
                  Martel.Float("release") + \
                  Martel.Str(", Last annotation update)") + Martel.AnyEol())


replacements = [
    ("ID", ID_exp),
    ("DT_created", DT_created_exp),
    ("DT_seq_update", DT_seq_update_exp),
    ("DT_ann_update", Martel.Opt(DT_ann_update_exp))
    ]

record = Martel.replace_groups(sprot40.record, replacements)


format_expression = Martel.replace_groups(
    sprot40.format_expression, replacements)

format = Martel.replace_groups(sprot40.format, replacements)
