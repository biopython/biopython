# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Turn an mmCIF file into a dictionary."""

from __future__ import print_function

import shlex

from Bio.File import as_handle


class MMCIF2Dict(dict):
    """Parse a mmCIF file and return a dictionary."""

    def __init__(self, filename):
        """Parse a mmCIF file and return a dictionary.

        Arguments:
         - file - name of the PDB file OR an open filehandle

        """
        with as_handle(filename) as handle:
            loop_flag = False
            key = None
            tokens = self._tokenize(handle)
            token = next(tokens)
            self[token[0:5]] = token[5:]
            i = 0
            n = 0
            for token in tokens:
                if token == "loop_":
                    loop_flag = True
                    keys = []
                    i = 0
                    n = 0
                    continue
                elif loop_flag:
                    if token.startswith("_"):
                        if i > 0:
                            loop_flag = False
                        else:
                            self[token] = []
                            keys.append(token)
                            n += 1
                            continue
                    else:
                        self[keys[i % n]].append(token)
                        i += 1
                        continue
                if key is None:
                    key = token
                else:
                    self[key] = token
                    key = None

    # Private methods

    def _tokenize(self, handle):
        for line in handle:
            if line.startswith("#"):
                continue
            elif line.startswith(";"):
                token = line[1:].strip()
                for line in handle:
                    line = line.strip()
                    if line == ';':
                        break
                    token += line
                yield token
            else:
                tokens = shlex.split(line)
                for token in tokens:
                    yield token
