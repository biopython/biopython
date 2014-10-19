# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Turn an mmCIF file into a dictionary."""

from __future__ import print_function

from Bio._py3k import input as _input

import shlex


class MMCIF2Dict(dict):

    def __init__(self, filename):
        with open(filename) as handle:
            loop_flag = False
            key = None
            tokens = self._tokenize(handle)
            token = next(tokens)
            self[token[0:5]]=token[5:]
            for token in tokens:
                if token=="loop_":
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
                        self[keys[i%n]].append(token)
                        i+=1
                        continue
                if key is None:
                    key = token
                else:
                    self[key] = token
                    key = None

    def _tokenize(self, handle):
        for line in handle:
            if line.startswith("#"):
                continue
            elif line.startswith(";"):
                token = line[1:].strip()
                for line in handle:
                    line = line.strip()
                    if line==';':
                        break
                    token += line
                yield token
            else:
                tokens = shlex.split(line)
                for token in tokens:
                    yield token


if __name__=="__main__":

    import sys

    if len(sys.argv)!=2:
        print("Usage: python MMCIF2Dict filename.")

    filename=sys.argv[1]

    mmcif_dict = MMCIF2Dict(filename)

    entry = ""
    print("Now type a key ('q' to end, 'k' for a list of all keys):")
    while(entry != "q"):
        entry = _input("MMCIF dictionary key ==> ")
        if entry == "q":
            sys.exit()
        if entry == "k":
            for key in mmcif_dict:
                print(key)
            continue
        try:
            value=mmcif_dict[entry]
            if isinstance(value, list):
                for item in value:
                    print(item)
            else:
                print(value)
        except KeyError:
            print("No such key found.")
