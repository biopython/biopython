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
        aa = line[12:12+64]
    elif line[:10] == "  sncbieaa":
        start = line[12:12+64]
    elif line[:9] == "  -- Base":
        bases.append(line[12:12+64])
    elif line[:2] == " }":
        assert names != [] and id is not None and aa is not None
        assert start is not None and bases != []
        if len(names) == 1:
            names.append(None)
        print("register_ncbi_table(name=%s," % repr(names[0]))
        print("                    alt_name=%s, id=%d," % \
              (repr(names[1]), id))
        print("                    table={")
        s = "    "
        for i in range(64):
            if aa[i] != "*":
                t = " '%s%s%s': '%s'," % (bases[0][i], bases[1][i],
                                          bases[2][i], aa[i])
                if len(s) + len(t) > 75:
                    print(s)
                    s = "    " + t
                else:
                    s = s + t
        print("%s }," % s)
        s = "                    stop_codons=["
        for i in range(64):
            if aa[i] == "*":
                t = "'%s%s%s'," % (bases[0][i], bases[1][i], bases[2][i])
                if len(s) + len(t) > 75:
                    s_with_spaces = s.replace("','", "', '")
                    print(s_with_spaces)
                    s = "                                    " + t
                else:
                    s = s + t
        s_with_spaces = s.replace("','", "', '")
        print("%s ]," % s_with_spaces)
        s = "                    start_codons=["
        for i in range(64):
            if start[i] == "M":
                t = "'%s%s%s'," % (bases[0][i], bases[1][i], bases[2][i])
                if len(s) + len(t) > 75:
                    s_with_spaces = s.replace("','", "', '")
                    print(s_with_spaces)
                    s = "                                    " + t
                else:
                    s = s + t
        s_with_spaces = s.replace("','", "', '")
        print("%s ]" % s_with_spaces)
        print("                    )")
        # Two blank lines between each function call
        print("")
        print("")
    elif line[:2] == "--" or line == "\n" or line == "}\n" or \
         line == 'Genetic-code-table ::= {\n':
        pass
    else:
        raise Exception("Unparsed: " + repr(line))


print("########################################################################")
print("# End of auto-generated output from Scripts/update_ncbi_codon_table.py #")
print("########################################################################")
