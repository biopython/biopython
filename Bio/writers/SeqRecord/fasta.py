from Bio import Writer

class WriteFasta(Writer.Writer):
    def __init__(self, outfile, seqwidth = 72):
        Writer.Writer.__init__(self, outfile)
        assert seqwidth > 0, seqwidth
        self.seqwidth = seqwidth
        
    def write(self, record):
        self.outfile.write(">%s %s\n" % (record.id, record.description))
        seq = record.seq
        seqwidth = self.seqwidth
        for i in range(0, len(seq), seqwidth):
            self.outfile.write(seq[i:i+seqwidth])
            self.outfile.write("\n")

make_writer = WriteFasta
