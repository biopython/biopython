#!/usr/bin/env python
import re
import sys

# re_link = re.compile(r"`\[(.+)\] <#(.+)>`__")
re_link = re.compile(r"`\[([A-Za-z0-9_:\-.]+)\] <#([A-Za-z0-9_:\-.]+)>`__")

assert re_link.findall(r"Chapter \ `[chapter:quick_start] <#chapter:quick_start>`__ before\n")
assert re_link.findall(r"(see Section `[sec:appendix-handles] <#sec:appendix-handles>`__):")
assert re_link.findall(r"Section `[sec:Bio.SeqIO-and-StringIO] <#sec:Bio.SeqIO-and-StringIO>`__):")
assert re_link.findall(r"Figure `[fig:three_track_cl2] <#fig:three_track_cl2>`__.")
assert list(re_link.finditer(r"in Chapter \ `[chapter:seq_annot] <#chapter:seq_annot>`__. This aims to"))
assert len(list(re_link.finditer(r"functions (`[eq:OP] <#eq:OP>`__) and (`[eq:NOP] <#eq:NOP>`__)."))) == 2

re_section = re.compile(r"`[0-9.]+ <#([A-Za-z0-9_:\-.]+)>`__")

assert re_section.findall(r"the label’s color (used in Section `1.1.9 <#sec:gd_nice_example>`__).")

def fix_line(line):
    # e.g. Chapter \ `[chapter:quick_start] <#chapter:quick_start>`__
    line = line.replace("\xa0\\ ", " ")
    for match in re_link.finditer(line):
        old = match.group()
        ref = match.group(1)
        assert ref == match.group(2), old
        new = r":ref:`%s`" % ref
        line = line.replace(old, new)
        print("%s -> %s" % (old, new))
    for match in re_section.finditer(line):
        old = match.group()
        ref = match.group(1)
        assert old.endswith("<#%s>`__" % ref), old
        new = r":ref:`%s`" % ref
        line = line.replace(old, new)
        print("%s -> %s" % (old, new))
    return line

def fix_file(filename):
    with open(filename) as handle:
        lines = list(handle)
    with open(filename, "w") as handle:
        for line in lines:
            handle.write(fix_line(line))

for f in sys.argv[1:]:
    sys.stderr.write("Fixing %s\n" % f)
    fix_file(f)
