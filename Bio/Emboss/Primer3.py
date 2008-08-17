"""Code to interact with the primer3 program.
"""

# --- primer3

class Record:
    """Represent information from a primer3 run finding primers.

    Members:

    primers   A list of primers that are generated (usually 5)
    """
    def __init__(self):
        self.comments = ""
        self.primers = []

class Primers:
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

def read(handle):
    """Parse primer3 output into a Primer3Record.
    """
    handle = iter(handle)
    record = Record()
    primer = None

    # Skip empty lines at the top of the file
    for line in handle:
        if line.strip():
            break

    # Read the comment lines
    while True:
        if line.strip() and not line[0]=='#':
            break
        record.comments += line
        try:
            line = handle.next()
        except StopIteration:
            return record

    # Read the primers
    while True:
        if not line.strip():
            pass
        elif line[5:19]=="PRODUCT SIZE: ":
            primer = Primers()
            primer.size = int(line[19:])
            record.primers.append(primer)
        elif line[5:19]=="FORWARD PRIMER":
            words = line.split()
            if not primer or primer.size==0:
                primer = Primers()
                record.primers.append(primer)
            primer.forward_start = int(words[2])
            primer.forward_length = int(words[3])
            primer.forward_tm = float(words[4])
            primer.forward_gc = float(words[5])
            primer.forward_seq = words[6]
        elif line[5:19]=="REVERSE PRIMER":
            words = line.split()
            if not primer or primer.size==0:
                primer = Primers()
                record.primers.append(primer)
            primer.reverse_start = int(words[2])
            primer.reverse_length = int(words[3])
            primer.reverse_tm = float(words[4])
            primer.reverse_gc = float(words[5])
            primer.reverse_seq = words[6]
        elif line[5:19]=="INTERNAL OLIGO":
            words = line.split()
            if not primer or primer.size==0:
                primer = Primers()
                record.primers.append(primer)
            primer.internal_start = int(words[2])
            primer.internal_length = int(words[3])
            primer.internal_tm = float(words[4])
            primer.internal_gc = float(words[5])
            primer.internal_seq = words[6]
        try:
            line = handle.next()
        except StopIteration:
            break

    return record
