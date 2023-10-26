# Copyright 2008 Michiel de Hoon. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Code to interact with the primersearch program from EMBOSS."""


class InputRecord:
    """Represent the input file into the primersearch program.

    This makes it easy to add primer information and write it out to the
    simple primer file format.
    """

    def __init__(self):
        """Initialize the class."""
        self.primer_info = []

    def __str__(self):
        """Summarize the primersearch input record as a string."""
        output = ""
        for name, primer1, primer2 in self.primer_info:
            output += f"{name} {primer1} {primer2}\n"
        return output

    def add_primer_set(self, primer_name, first_primer_seq, second_primer_seq):
        """Add primer information to the record."""
        self.primer_info.append((primer_name, first_primer_seq, second_primer_seq))


class OutputRecord:
    """Represent the information from a primersearch job.

    amplifiers is a dictionary where the keys are the primer names and
    the values are a list of PrimerSearchAmplifier objects.
    """

    def __init__(self):
        """Initialize the class."""
        self.amplifiers = {}


class Amplifier:
    """Represent a single amplification from a primer."""

    def __init__(self):
        """Initialize the class."""
        self.hit_info = ""
        self.length = 0


def read(handle):
    """Get output from primersearch into a PrimerSearchOutputRecord."""
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
