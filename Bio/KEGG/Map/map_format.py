# Copyright 2001 Tarjei Mikkelsen. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Martel based parser to read KEGG Map files from the Ligand database.

This is a regular expression for KEGG Map files, built using the 'regular
expressions on steroids' capabilities of Martel.

A description of the format can be found in the 'ligand.doc' file
from the Ligand distribution, available from:

 http://www.genome.ad.jp/kegg
 
"""

# Martel
from Martel import *
from Martel import RecordReader


reaction_id = Group("reaction_id", Re("R\d{5}"))
ec_id = Group("ec_id", Rep(Rep1(Integer()) + Str1(".")) + Alt(Rep1(Integer()),
                                                          Str1("-")))
stoch = Group("stoch", Rep1(Integer()))
reactant = Group("reactant", AnyBut(" +") + Rep1(AnyBut(" <=>\n")) +
                 Rep(Str1(" ") + AnyBut(" +") + Rep1(AnyBut(" <=>\n"))))
plus = Group("plus", Rep1(Str1(" ")) + Str1("+") + Rep1(Str1(" ")))
arrows = Group("arrows", Str1("<=>"))
end = Group("end", AnyEol())

reaction = Group("reaction",
                 reaction_id + Str1(":") + Opt(ec_id + Str1(":")) + Rep(Str1(" ")) +
                 Rep1(Opt(stoch + Str1(" ")) + Rep(Str1(" ")) + reactant + Opt(plus)) + 
                 Rep(Str1(" ")) + arrows + Rep(Str1(" ")) +
                 Rep1(Opt(stoch + Str1(" ")) + Rep(Str1(" ")) + reactant + Opt(plus)) +
                 Rep(Str1(" ")) + end)

record_format = ParseRecords("KEGG_Map_File",
                             reaction, RecordReader.CountLines, (1,))



