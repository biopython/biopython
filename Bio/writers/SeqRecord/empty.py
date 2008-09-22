"""Part of an old unused and undocumented sequence writing framework (DEPRECATED)."""
from Bio import Writer

class WriteEmpty(Writer.Writer):
    pass

make_writer = WriteEmpty
