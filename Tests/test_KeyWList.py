# Copyright 1999 by Jeffrey Chang.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import os
from Bio.SwissProt import KeyWList

filename = os.path.join("SwissProt", "keywlist.txt")
handle = open(filename)
records = KeyWList.parse(handle)

record = records.next()
assert record["ID"]=="2Fe-2S."
assert record["AC"]=="KW-0001"
assert record["DE"]=="Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron atoms complexed to 2 inorganic sulfides and 4 sulfur atoms of cysteines from the protein."
assert record["SY"]=="Fe2S2; [2Fe-2S] cluster; [Fe2S2] cluster; Fe2/S2 (inorganic) cluster; Di-mu-sulfido-diiron; 2 iron, 2 sulfur cluster binding."
assert len(record["GO"])==1
assert record["GO"]==["GO:0051537; 2 iron, 2 sulfur cluster binding"]
assert len(record["HI"])==2
assert record["HI"][0]=="Ligand: Iron; Iron-sulfur; 2Fe-2S."
assert record["HI"][1]=="Ligand: Metal-binding; 2Fe-2S."
assert record["CA"]=="Ligand."

record = records.next()
assert record["IC"]=="Molecular function."
assert record["AC"]=="KW-9992"
assert record["DE"]=="Keywords assigned to proteins due to their particular molecular function."

record = records.next()
assert record["ID"]=="Zymogen."
assert record["AC"]=="KW-0865"
assert record["DE"]=="The enzymatically inactive precursor of mostly proteolytic enzymes."
assert record["SY"]=="Proenzyme."
assert len(record["HI"])==1
assert record["HI"][0]=="PTM: Zymogen."
assert record["CA"]=="PTM."
