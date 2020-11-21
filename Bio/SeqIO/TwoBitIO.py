from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import _twoBitIO
from .Interfaces import SequenceIterator


class TwoBitIterator(SequenceIterator):
    """Parser for UCSC twoBit (.2bit) files."""

    def __init__(self, source):
        super().__init__(source, mode="b", fmt="twoBit")
        isByteSwapped, names, sequences = _twoBitIO.TwoBitIterator(self.stream)
        self.isByteSwapped = isByteSwapped
        self.names = names
        self.sequences = [Seq(sequence) for sequence in sequences]
        # wait to close the file until the TwoBitIterator goes out of scope:
        self.should_close_stream = False

    def parse(self, stream):
        for name, sequence in zip(self.names, self.sequences):
            sequence = Seq(sequence)
            record = SeqRecord(sequence, id=name)
            yield record

    def __getitem__(self, key):
        try:
            index = self.names.index(key)
        except ValueError:
            raise KeyError(key) from None
        sequence = self.sequences[index]
        return SeqRecord(sequence, id=key)

    def keys(self):
        return self.names

    def __len__(self):
        return len(self.names)
