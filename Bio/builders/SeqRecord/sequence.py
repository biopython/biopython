from xml.sax import handler

from Bio import SeqRecord, StdHandler, Dispatch, Std

class BuildSeqRecord(Dispatch.Dispatcher):
    def __init__(self):
        Dispatch.Dispatcher.__init__(self)
        self.acquire(StdHandler.Handle_dbid(self.add_dbid),
                     prefix = Std.NS)
        self.acquire(StdHandler.Handle_description(self.add_description),
                     prefix = Std.NS)
        self.acquire(StdHandler.Handle_sequence(self.add_sequence),
                     prefix = Std.NS)

    def start_record(self, tag, attrs):
        self.dbname = None
        self.id_text = None
        self.description = None
        self.alphabet = None
        self.seq = None
        

    def add_dbid(self, text, attrs):
        if attrs.get("type") == "primary":
            self.dbname = attrs.get("dbname", "unknown")
            self.id_text = text

    def add_description(self, text):
        self.description = text

    def add_sequence(self, alphabet, seq):
        self.alphabet = alphabet
        self.seq = seq

    def end_record(self, tag):
        self.document = SeqRecord.SeqRecord(
            seq = self.seq,    # Need a real Seq record!
            id = self.id_text,
            description = self.description,
            )
        
make_builder = BuildSeqRecord
