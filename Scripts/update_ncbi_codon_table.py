# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Helper script to update Codon tables from the NCBI.

These tables are based on parsing the NCBI file:
ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt

This script is used to update Bio/Data/CodonTable.py

Note that the NCBI sometimes revise the older tables,
so don't just add new tables - replace all of them
and check for any differences in the old tables.
"""

import re

INDENT = len("register_ncbi_table(")


def line_wrap(text, indent=0, max_len=78):
    if len(text) <= max_len:
        return text
    line = text[:max_len]
    assert " " in line, line
    line, rest = line.rsplit(" ", 1)
    rest = " " * indent + rest + text[max_len:]
    assert len(line) < max_len
    if indent + len(rest) <= max_len:
        return line + "\n" + rest
    else:
        return line + "\n" + line_wrap(rest, max_len)


print("##########################################################################")
print("# Start of auto-generated output from Scripts/update_ncbi_codon_table.py #")
print("##########################################################################")
print("")
print("")

for line in open("gc.prt").readlines():
    if line[:2] == " {":
        names = []
        id = None
        aa = None
        start = None
        bases = []
    elif line[:6] == "  name":
        names.append(re.search('"([^"]*)"', line).group(1))
    elif line[:8] == "    name":
        names.append(re.search('"(.*)$', line).group(1))
    elif line == ' Mitochondrial; Mycoplasma; Spiroplasma" ,\n':
        names[-1] = names[-1] + " Mitochondrial; Mycoplasma; Spiroplasma"
    elif line[:4] == "  id":
        id = int(re.search('(\d+)', line).group(1))
    elif line[:10] == "  ncbieaa ":
        aa = line[12:12 + 64]
    elif line[:10] == "  sncbieaa":
        start = line[12:12 + 64]
    elif line[:9] == "  -- Base":
        bases.append(line[12:12 + 64])
    elif line[:2] == " }":
        assert names != [] and id is not None and aa is not None
        assert start is not None and bases != []
        if len(names) == 1:
            names.append(None)
        print("register_ncbi_table(name=%r," % names[0])
        print(" " * INDENT + "alt_name=%r, id=%d," % (names[1], id))
        print(" " * INDENT + "table={")
        s = "    "
        for i in range(64):
            if aa[i] != "*":
                t = " '%s%s%s': '%s'," % (bases[0][i], bases[1][i],
                                          bases[2][i], aa[i])
                if len(s) + len(t) > 75:
                    print(s)
                    s = "    " + t
                else:
                    s += t
        print("%s }," % s)
        codons = [bases[0][i] + bases[1][i] + bases[2][i]
                  for i in range(64) if aa[i] == "*"]
        print(line_wrap(" " * INDENT + "stop_codons=%r," % codons, indent=INDENT + 13))
        codons = [bases[0][i] + bases[1][i] + bases[2][i]
                  for i in range(64) if start[i] == "M"]
        print(line_wrap(" " * INDENT + "start_codons=%r)" % codons, indent=INDENT + 14))
        print("")
    elif line[:2] == "--" or line in ("\n", "}\n", 'Genetic-code-table ::= {\n'):
        pass
    else:
        raise Exception("Unparsed: " + repr(line))

print("")
print("########################################################################")
print("# End of auto-generated output from Scripts/update_ncbi_codon_table.py #")
print("########################################################################")
