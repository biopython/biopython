"""Code to interact with Primer-related programs from EMBOSS (DEPRECATED).

Bio.Emboss.Primer has been deprecated, please use Bio.Emboss.Primer3 or
Bio.Emboss.PrimerSearch instead.

To parse primersearch output into a PrimerSearch.OutputRecord, use
    from Bio.Emboss import PrimerSearch
    handle = open('myprimersearchoutputfile.txt')
    record = PrimerSearch.read(handle)

To parse primer3 output into a Primer3.Record, use
    from Bio.Emboss import Primer3
    handle = open('myprimer3outputfile.txt')
    record = Primer3.read(handle)
"""

import warnings
warnings.warn("""\
Bio.Emboss.Primer has been deprecated.
Please use Bio.Emboss.Primer3 or Bio.Emboss.PrimerSearch instead.

To parse primersearch output into a PrimerSearch.OutputRecord, use
    from Bio.Emboss import PrimerSearch
    handle = open('myprimersearchoutputfile.txt')
    record = PrimerSearch.read(handle)

To parse primer3 output into a Primer3.Record, use
    from Bio.Emboss import Primer3
    handle = open('myprimer3outputfile.txt')
    record = Primer3.read(handle)
""", DeprecationWarning)

# standard library
import string
from xml.sax import handler

# Martel
import Martel

# biopython stuff
from Bio.ParserSupport import AbstractConsumer
from Bio.ParserSupport import EventGenerator

import primersearch_format
import primer3_format

# --- primersearch

class PrimerSearchInputRecord:
    """Represent the input file into the primersearch program.

    This makes it easy to add primer information and write it out to the
    simple primer file format.
    """ 
    def __init__(self):
        self.primer_info = []

    def __str__(self):
        output = ""
        for name, primer1, primer2 in self.primer_info:
            output += "%s %s %s\n" % (name, primer1, primer2)
        return output

    def add_primer_set(self, primer_name, first_primer_seq, 
                       second_primer_seq):
        """Add primer information to the record.
        """
        self.primer_info.append((primer_name, first_primer_seq,
                                 second_primer_seq))
           
class PrimerSearchParser:
    """Parse primersearch output into a PrimerSearchOutputRecord.
    """
    def __init__(self, debug_level = 0):
        self._scanner = _PrimerSearchScanner(debug_level)

    def parse(self, handle):
        self._consumer = _PrimerSearchRecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class PrimerSearchOutputRecord:
    """Represent the information from a primersearch job.

    amplifiers is a dictionary where the keys are the primer names and
    the values are a list of PrimerSearchAmplifier objects.
    """
    def __init__(self):
        self.amplifiers = {}

class PrimerSearchAmplifier:
    """Represent a single amplification from a primer.
    """
    def __init__(self):
        self.hit_info = ""
        self.length = 0

class _PrimerSearchRecordConsumer(AbstractConsumer):
    """Get output from primersearch into a PrimerSearchOutputRecord
    """
    def __init__(self):
        self.data = PrimerSearchOutputRecord()
        self._cur_primer = None
        self._cur_amplifier = None

    def _add_last_amplifier(self):
        # add on the last amplifier
        if self._cur_primer is not None and self._cur_amplifier is not None:
            self.data.amplifiers[self._cur_primer].append(self._cur_amplifier)

    def primer_name(self, name):
        self._add_last_amplifier() 
        self.data.amplifiers[name] = []
        self._cur_primer = name

    def amplifier(self, amplifier_num):
        self._add_last_amplifier()
        self._cur_amplifier = PrimerSearchAmplifier()

    def amplifier_sequence(self, sequence_info):
        self._cur_amplifier.hit_info = sequence_info

    def amplifier_length(self, amplifier_info):
        self._cur_amplifier.length = int(amplifier_info)

    def end_record(self):
        self._add_last_amplifier()

class _PrimerSearchScanner:
    """Scan output from the primersearch program.
    """
    def __init__(self, debug = 0):
        self.interest_tags = ["primer_name", "amplifier", 
                              "amplifier_sequence", "amplifier_length",
                              "end_record"]

        expression = Martel.select_names(primersearch_format.record,
                                            self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        self._parser.setContentHandler(EventGenerator(consumer,
                                                      self.interest_tags,
                                                      _strip_and_combine))
        self._parser.setErrorHandler(handler.ErrorHandler())
        self._parser.parseFile(handle)

        consumer.end_record()

# --- primer3

class Primer3Record:
    """Represent information from a primer3 run finding primers.

    Members:

    primers   A list of primers that are generated (usually 5)
    """
    def __init__(self):
        self.comments = ""
        self.primers = []

class Primer3Primers:
    """A primer set designed by Primer3.

    Members:

    size
    forward_seq
    forward_start
    forward_length
    forward_tm
    forward_gc
    reverse_seq
    reverse_start
    reverse_length
    reverse_tm
    reverse_gc
    """
    def __init__(self):
        self.size = 0
        self.forward_seq = ""
        self.forward_start = 0
        self.forward_length = 0
        self.forward_tm = 0.0
        self.forward_gc = 0.0
        self.reverse_seq = ""
        self.reverse_start = 0
        self.reverse_length = 0
        self.reverse_tm = 0.0
        self.reverse_gc = 0.0

class Primer3Parser:
    """Parse primer3 output into a Primer3Record.
    """
    def __init__(self, debug_level = 0):
        self._scanner = _Primer3Scanner(debug_level)

    def parse(self, handle):
        self._consumer = _Primer3RecordConsumer()
        self._scanner.feed(handle, self._consumer)
        return self._consumer.data

class _Primer3RecordConsumer(AbstractConsumer):
    """Get output from prime3 into a Primer3Record
    """
    def __init__(self):
        self.data = Primer3Record()
        self._cur_primer = None

    def _add_last_primer(self):
        # add on the last amplifier
        if self._cur_primer is not None:
            self.data.primers.append(self._cur_primer)

    def comments(self, comment):
        self.data.comments = comment

    def start_primer(self, junk):
        self._add_last_primer()
        self._cur_primer = Primer3Primers()

    def single_primer_line(self, junk):
        self.start_primer(junk)

    def product_size(self, size):
        self._cur_primer.size = int(size)

    def forward_start(self, start):
        self._cur_primer.forward_start = int(start)

    def forward_length(self, length):
        self._cur_primer.forward_length = int(length)

    def forward_tm(self, tm):
        self._cur_primer.forward_tm = float(tm)

    def forward_gc(self, gc):
        self._cur_primer.forward_gc = float(gc)

    def forward_seq(self, seq):
        self._cur_primer.forward_seq = seq

    def reverse_start(self, start):
        self._cur_primer.reverse_start = int(start)

    def reverse_length(self, length):
        self._cur_primer.reverse_length = int(length)

    def reverse_tm(self, tm):
        self._cur_primer.reverse_tm = float(tm)

    def reverse_gc(self, gc):
        self._cur_primer.reverse_gc = float(gc)

    def reverse_seq(self, seq):
        self._cur_primer.reverse_seq = seq

    def internal_start(self, start):
        self._cur_primer.internal_start = int(start)

    def internal_length(self, length):
        self._cur_primer.internal_length = int(length)

    def internal_tm(self, tm):
        self._cur_primer.internal_tm = float(tm)

    def internal_gc(self, gc):
        self._cur_primer.internal_gc = float(gc)

    def internal_seq(self, seq):
        self._cur_primer.internal_seq = seq

    def end_record(self):
        self._add_last_primer()

class _Primer3Scanner:
    """Scan output from the primer3 program.
    """
    def __init__(self, debug = 0):
        self.interest_tags = ["comments", "single_primer_line",
                              "start_primer", "product_size",
                              "forward_start", "forward_length",
                              "forward_tm", "forward_gc", "forward_seq",
                              "reverse_start", "reverse_length",
                              "reverse_tm", "reverse_gc", "reverse_seq",
                              "internal_start", "internal_length",
                              "internal_tm", "internal_gc", "internal_seq",
                              "end_record"]

        expression = Martel.select_names(primer3_format.record,
                                         self.interest_tags)
        self._parser = expression.make_parser(debug_level = debug)

    def feed(self, handle, consumer):
        self._parser.setContentHandler(EventGenerator(consumer,
                                                      self.interest_tags,
                                                      _strip_and_combine))
        self._parser.setErrorHandler(handler.ErrorHandler())
        self._parser.parseFile(handle)

        consumer.end_record()

def _strip_and_combine(line_list):
    """Combine multiple lines of content separated by spaces.
    """
    # first strip out extra whitespace
    stripped_line_list = map(string.strip, line_list)
    # now combine everything with spaces
    return ' '.join(stripped_line_list)
