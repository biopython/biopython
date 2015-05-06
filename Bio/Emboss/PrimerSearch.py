# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Code to interact with the primersearch program from EMBOSS.
"""

__docformat__ = "restructuredtext en"


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
        self.seq_id = ""
        self.seq_description = ""
        self.binding_sites = {}


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
            seq_id = line.replace("\tSequence: ", "")
            amplifier.seq_id = seq_id.strip()
            amplifier.hit_info = seq_id
        elif line.startswith("\tAmplimer length: "):
            length = line.split()[-2]
            amplifier.length = int(length)
        elif line.endswith("mismatches\n"):
            words = line.strip().split()
            misses = int(words[-2])
            location = words[-4]
            # PrimerSearch notation for distance from end is [<index>]
            if location.startswith("["):
                # PrimerSearch is 1-indexed
                location = -int(location[1:-1]) + 1
            else:
                location = int(location) - 1
            seq = words[0]
            amplifier.binding_sites[seq] = (location, misses)
            amplifier.hit_info += line


        else:
            amplifier.seq_description = line.strip()
            amplifier.hit_info += line

    for name in record.amplifiers:
        for amplifier in record.amplifiers[name]:
            amplifier.hit_info = amplifier.hit_info.rstrip()

    return record
