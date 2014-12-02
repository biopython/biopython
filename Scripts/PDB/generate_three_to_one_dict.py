# File download/unzip written 2012 by Lenna X. Peterson (arklenna@gmail.com)
# Dictionary extraction written 2011 by Hongbo Zhu
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Download PDB Chemical Component Dictionary and generate dict.

Download and parse PDB Chemical Component Dictionary,
then write out dict for to_one_letter_code.
"""

from __future__ import print_function

import gzip
import inspect
import os
import warnings

from Bio._py3k import urlopen

url = "ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"

# extract name of gzip file
gzname = os.path.basename(url)
# extract name of cif file (split by sep, remove last, rejoin)
cifname = os.extsep.join(gzname.split(os.extsep)[:-1])

url_handle = urlopen(url)

with open(gzname, 'wb') as gzh:
    print("Downloading file... (approx. 29 MB)")
    while True:
        data = url_handle.read(1024)
        if len(data) == 0:
            break
        gzh.write(data)

# size as of 13 April 2012
if os.path.getsize(gzname) < 29944258:
    warnings.warn("ERROR: Downloaded file is too small",
                  RuntimeWarning)

fh = gzip.open(gzname, 'rb')

# write extracted file to disk (not necessary)
# with open(cifname, 'wb') as cifh:
#     print("Extracting file...")
#     cifh.write(fh.read())

# The following code written by Hongbo Zhu
# generate three_to_one_dict
# two records in PDB Chemical Component Dictionary are parsed to
# generate the dictionary:
# _chem_comp.one_letter_code
# _chem_comp.three_letter_code

three_to_one_buf = []      # all three-letter codes
three_to_one_buf_noq = []  # only those with non-'?' one-letter codes

current_line = 'to_one_letter_code = {'
current_line_noq = 'to_one_letter_code = {'

found_one = False    # found one-letter code
found_three = False  # found three-letter code

counter = 0
counter_noq = 0

line = fh.readline()

while line:
    if line.startswith('_chem_comp.one_letter_code'):
        one = line.strip().split()[-1]
        found_one = True
    if line.startswith('_chem_comp.three_letter_code'):
        three = '%-3s' % (line.strip().split()[-1],)  # make it three-letter
        found_three = True

    if found_one and found_three:
        if counter % 5 == 0:
            three_to_one_buf.append('%s\n' % (current_line,))
            current_line = '    '

        current_line = '%s\'%s\':\'%s\',' % (current_line, three, one)
        counter += 1

        if one != '?':
            if counter_noq % 5 == 0:
                three_to_one_buf_noq.append('%s\n' % (current_line_noq,))
                current_line_noq = '    '

            current_line_noq = '%s\'%s\':\'%s\',' % (current_line_noq, three, one)
            counter_noq += 1

        found_one = False
        found_three = False

    line = fh.readline()

if len(current_line) < 5:
    three_to_one_buf[-1] = three_to_one_buf[:-1]  # remove the last comma
    three_to_one_buf.append('}')
else:
    three_to_one_buf.append('%s }' % (current_line[:-1]))

if len(current_line_noq) < 5:
    three_to_one_buf_noq[-1] = three_to_one_buf_noq[:-1]
    three_to_one_buf_noq.append('}')
else:
    three_to_one_buf_noq.append('%s }' % (current_line_noq[:-1]))

# Find path of current script
_scriptPath = os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0])
# Path to SCOP module
_rafPath = os.path.normpath(os.path.join(_scriptPath, "..", "..", "Bio", "SCOP"))
_threeAllPath = os.path.join(_rafPath, 'three_to_one_all.py')
_threePath = os.path.join(_rafPath, 'three_to_one_dict.py')

# with open(_threeAllPath, 'w') as fh:
#     fh.writelines(three_to_one_buf)
with open(_threePath, 'w') as fh:
    fh.writelines(three_to_one_buf_noq)
