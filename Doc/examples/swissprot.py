# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Example of connecting with exPASy and parsing SwissProt records."""

# biopython

from Bio import ExPASy
from Bio import SwissProt

# 'O23729', 'O23730', 'O23731', Chalcone synthases from Orchid

ids = ["O23729", "O23730", "O23731"]

for id in ids:
    handle = ExPASy.get_sprot_raw(id)
    record = SwissProt.read(handle)
    print(f"description: {record.description}")
    for ref in record.references:
        print(f"authors: {ref.authors}")
        print(f"title: {ref.title}")

    print(f"classification: {record.organism_classification}")
    print("")
