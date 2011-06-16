"""Code to interact with the primersearch program from EMBOSS.
"""


class InputRecord(object):
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
           
class OutputRecord(object):
    """Represent the information from a primersearch job.

    amplifiers is a dictionary where the keys are the primer names and
    the values are a list of PrimerSearchAmplifier objects.
    """
    def __init__(self):
        self.amplifiers = {}

class Amplifier(object):
    """Represent a single amplification from a primer.
    """
    def __init__(self):
        self.hit_info = ""
        self.length = 0

def read(handle):
    """Get output from primersearch into a PrimerSearchOutputRecord
    """
    record = OutputRecord()

    for line in handle:
        if not line.strip():
            continue
        elif line.startswith("Primer name"):
            name = line.split()[-1]
            record.amplifiers[name] = []
        elif line.startswith("Amplimer"):
            amplifier = Amplifier()
            record.amplifiers[name].append(amplifier)
        elif line.startswith("\tSequence: "):
            amplifier.hit_info = line.replace("\tSequence: ", "")
        elif line.startswith("\tAmplimer length: "):
            length = line.split()[-2]
            amplifier.length = int(length)
        else:
            amplifier.hit_info += line

    for name in record.amplifiers:
        for amplifier in record.amplifiers[name]:
            amplifier.hit_info = amplifier.hit_info.rstrip()

    return record
