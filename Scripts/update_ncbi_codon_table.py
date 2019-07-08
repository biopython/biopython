# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Helper script to update Codon tables from the NCBI.

These tables are based on parsing the NCBI file:
ftp://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt

More detailed information about the tables are here:
https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

This script is used to update Bio/Data/CodonTable.py

Note that the NCBI sometimes revise the older tables,
so don't just add new tables - replace all of them
and check for any differences in the old tables.
"""

import re

INDENT = len("register_ncbi_table(")
INDENT2 = INDENT + len("table={")


def line_wrap(text, indent=0, max_len=78, string=False):
    """Return a wrapped line if length is larger max_len.

    The new parameter 'string' allows to wrap quoted text which is delimited
    by single quotes. It adds " '" to the end of the line and "'" to the start
    of the next line.
    """
    split_len = max_len if not string else max_len - 2
    if len(text) <= max_len:
        return text
    line = text[:split_len]
    assert " " in line, line
    line, rest = line.rsplit(" ", 1)
    # New:
    if string:
        line += " '"
        rest = "'" + rest
    rest = " " * indent + rest + text[split_len:]
    assert len(line) < max_len
    if indent + len(rest) <= max_len:
        return line + "\n" + rest
    else:
        return line + "\n" + line_wrap(rest, indent, max_len, string)


print("""
##########################################################################
# Start of auto-generated output from Scripts/update_ncbi_codon_table.py #
##########################################################################

"""[1:])

version = ''
for line in open("gc.prt").readlines():
    if not version and line.startswith("--  Version"):
        version = line.split('Version', 1)[1].strip()
        print("# Data from NCBI genetic code table version %s\n" % version)
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
        id = int(re.search(r'(\d+)', line).group(1))
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
        # Use %r instead of %s to include the quotes of the string!
        print(line_wrap("register_ncbi_table(name=%r," % names[0],
              indent=INDENT, string=True))
        print(" " * INDENT + "alt_name=%r, id=%d," % (names[1], id))
        s = " " * INDENT + "table={"
        for i in range(64):
            if aa[i] != "*":
                t = "'%s%s%s': '%s', " % (bases[0][i], bases[1][i],
                                          bases[2][i], aa[i])
                if len(s) + len(t) > 75:
                    print(s.rstrip())
                    s = " " * INDENT2 + t
                else:
                    s += t
        print("%s}," % s[:-2])  # remove ', ' from last entry before '}'
        codons = [bases[0][i] + bases[1][i] + bases[2][i]
                  for i in range(64) if start[i] == "*"]
        print(line_wrap(" " * INDENT + "stop_codons=%r,"
                        % codons, indent=INDENT + 13))
        codons = [bases[0][i] + bases[1][i] + bases[2][i]
                  for i in range(64) if start[i] == "M"]
        print(line_wrap(" " * INDENT + "start_codons=%r)"
                        % codons, indent=INDENT + 14))
        print("")
    elif line[:2] == "--" or line in ("\n", "}\n",
                                      'Genetic-code-table ::= {\n'):
        pass
    else:
        raise Exception("Unparsed: " + repr(line))

print("""
########################################################################
# End of auto-generated output from Scripts/update_ncbi_codon_table.py #
########################################################################""")
