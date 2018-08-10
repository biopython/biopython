# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Code to interact with the primersearch program from EMBOSS."""


class InputRecord(object):
    """Represent the input file into the primersearch program.

    This makes it easy to add primer information and write it out to the
    simple primer file format.
    """

    def __init__(self):
        """Initialize the class."""
        self.primer_info = []

    def __str__(self):
        output = ""
        for name, primer1, primer2 in self.primer_info:
            output += "%s %s %s\n" % (name, primer1, primer2)
        return output

    def add_primer_set(self, primer_name, first_primer_seq,
                       second_primer_seq):
        """Add primer information to the record."""
        self.primer_info.append((primer_name, first_primer_seq,
                                 second_primer_seq))


class OutputRecord(object):
    """Represent the information from a primersearch job.

    amplifiers is a dictionary where the keys are the primer names and
    the values are a list of PrimerSearchAmplifier objects.
    """

    def __init__(self):
        """Initialize the class."""
        self.amplifiers = {}


class Amplifier(object):
    """Represent a single amplification from a primer."""

    def __init__(self, primer_name):
        """Initialize the class."""
        self.primer_name = primer_name
        self.hit_info = ""
        self.length = 0


def parse(handle):
    current_aplifier = None
    current_primer_name = None
    for line in handle:
        print('AAAAAAA')
        if not line.strip():
            continue

        elif line.startswith("Primer name"):
            current_primer_name = line.split()[-1]
        elif line.startswith("Amplimer"):
            if current_aplifier is not None:
                current_aplifier.hit_info = current_aplifier.hit_info.rstrip()
                yield current_aplifier
            current_aplifier = Amplifier(current_primer_name)
        elif line.startswith("\tSequence: "):
            current_aplifier.hit_info = line.replace("\tSequence: ", "")
        elif line.startswith("\tAmplimer length: "):
            length = line.split()[-2]
            current_aplifier.length = int(length)
        else:
            current_aplifier.hit_info += line
    if current_aplifier is not None:
        current_aplifier.hit_info = current_aplifier.hit_info.rstrip()
        yield current_aplifier


def read(handle):
    """Get output from primersearch into a PrimerSearchOutputRecord."""
    record = OutputRecord()

    for amplifier in parse(handle):
        if amplifier.primer_name in record.amplifiers:
            record.amplifiers[amplifier.primer_name].append(amplifier)
        else:
            record.amplifiers[amplifier.primer_name] = [amplifier]

    return record
