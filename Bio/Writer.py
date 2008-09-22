"""Part of an old unused and undocumented sequence writing framework (DEPRECATED)."""

import warnings
warnings.warn("Bio.Writer and Bio.writer.* are deprecated. If you do use"\
              +" these modules, please get in touch via the mailing list or"\
              +" bugzilla to avoid their permanent removal from Biopython.", \
              DeprecationWarning)

class Writer:
    def __init__(self, outfile):
        self.outfile = outfile
    def writeHeader(self):
        pass
    def write(self, record):
        pass
    def writeFooter(self):
        pass
