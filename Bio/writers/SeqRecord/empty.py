from Bio import Writer

class WriteEmpty(Writer.Writer):
    pass

make_writer = WriteEmpty
