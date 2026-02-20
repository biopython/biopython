# Copyright 2025 by Oz Nova.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""
Benchmark Bio.PDB.MMCIFParser and Bio.PDB.FastMMCIFParser using available cif files
"""

import os
from statistics import mean, stdev
import time

from Bio.PDB import MMCIFParser, FastMMCIFParser


ITERATIONS = 10


def _main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_dir = os.path.join(script_dir, "..", "..", "Tests", "PDB")
    files = sorted(
        os.path.join(pdb_dir, f)
        for f in os.listdir(pdb_dir)
        if f.endswith(".cif") and "_" not in f
    )

    print(f"Benchmarking with {len(files)} files.")
    print()
    print(f'{"File":<12} {"MMCIF (s)":<12} {"FastMMCIF (s)":<12}')
    print(f'{"---":<12} {"---":<12} {"---":<12}')

    parsers = (MMCIFParser(QUIET=True), FastMMCIFParser(QUIET=True))
    times = ([], [])
    for test_file in files:
        pdb_id = os.path.splitext(os.path.basename(test_file))[0]
        for out, parser in zip(times, parsers):
            start_time = time.time()
            for _ in range(ITERATIONS):
                _ = parser.get_structure(pdb_id, test_file)
            end_time = time.time()
            elapsed = end_time - start_time
            out.append(elapsed / ITERATIONS)

        print(f"{pdb_id:<12} {times[0][-1]:<12.4f} {times[1][-1]:<12.4f}")

    print()
    print(f'{"TOTAL":<12} {sum(times[0]):<12.4f} {sum(times[1]):<12.4f}')
    print(f'{"MEAN":<12} {mean(times[0]):<12.4f} {mean(times[1]):<12.4f}')
    print(f'{"STDEV":<12} {stdev(times[0]):<12.4f} {stdev(times[1]):<12.4f}')
    print()


if __name__ == "__main__":
    _main()
