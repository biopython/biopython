import warnings
warnings.warn("Bio.Writer was deprecated, as any Biopython modules that used "\
               +"it have been deprecated. If you do use this module, please "\
               +"get in touch via the Biopython mailing list or bugzilla.")

class Writer:
    def __init__(self, outfile):
        self.outfile = outfile
    def writeHeader(self):
        pass
    def write(self, record):
        pass
    def writeFooter(self):
        pass
