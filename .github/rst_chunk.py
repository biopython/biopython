#!/usr/bin/env python
import os
import sys

if len(sys.argv) == 2:
    prefix = sys.argv[1]
else:
    prefix = "chapter"
marker = ".. _%s" % prefix

handle = sys.stdout
for line in sys.stdin:
    if line.rstrip().startswith(marker) and line.rstrip()[-1] == ":":
        if handle != sys.stdout:
            handle.close()
        filename = "%s_%s.rst" % (prefix, line.rstrip()[len(marker) + 1:-1])
        sys.stderr.write("Starting %s\n" % filename)
        handle = open(filename, "w")
    handle.write(line)
if handle != sys.stdout:
    handle.close()
