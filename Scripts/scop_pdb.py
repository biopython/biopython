#!/usr/bin/env python

# Copyright 2001 by Gavin E. Crooks.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


import getopt
import sys
import types
import urllib

from Bio.SCOP import *

def usage() :
    print \
"""Extract a SCOP domain's ATOM and HETATOM records from the relevant PDB file.

For example:
  scop_pdb.py astral-rapid-access-1.55.raf dir.cla.scop.txt_1.55 d3hbib_

A result file, d3hbib_.ent, will be generated in the working directory.

The required RAF file can be found at [http://astral.stanford.edu/raf.html],
and the SCOP CLA file at [http://scop.berkeley.edu/parse/index.html].

Note: Errors will occur if the PDB file has been altered since the creation
of the SCOP CLA and ASTRAL RAF files.
 
Usage: scop_pdb [-h] [-i file] [-o file] [-p pdb_url_prefix]
                 raf_url cla_url [sid] [sid] [sid] ...

 -h        -- Print this help message.

 -i file   -- Input file name. Each line should start with an sid (Scop domain
              identifier). Blank lines, and lines starting with '#' are
              ignored. If file is '-' then data is read from stdin. If not
              given then sids are taken from the command line. 

 -o file   -- Output file name. If '-' then data is written to stdout. If not
              given then data is written to files named sid+'.ent'.

 -p pdb_url-- A URL for PDB files. The token '%s' will be replaced with the
              4 character PDB ID. If the pdb_url is not given then the latest
              PDB file is retrieved directly from rcsb.org.


  raf_url  -- The URL or filename of an ASTRAL Rapid Access File sequence map.
              See [http://astral.stanford.edu/raf.html]

  cla_url  -- The URL or filename of a SCOP parsable CLA file.
              See [http://scop.berkeley.edu/parse/index.html]

  sid      -- A SCOP domain identifier. e.g. d3hbib_ 
"""



default_pdb_url = "http://www.rcsb.org/pdb/cgi/export.cgi/somefile.pdb?" \
                      "format=PDB&pdbId=%s&compression=None"
#default_pdb_url = "file://usr/local/db/pdb/data/010331/snapshot/all/pdb%s.ent"

def open_pdb(pdbid, pdb_url=None) :
    if pdb_url ==None: pdb_url = default_pdb_url
    url = pdb_url % pdbid
    fn, header = urllib.urlretrieve(url)
    return open(fn)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hp:o:i:",
             ["help", "usage","pdb=","output=","input="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    input= None
    in_handle = None
    output = None
    pdb_url = None
    cla_url = None
    raf_url = None
    
    for o, a in opts:
        if o in ("-h", "--help","--usage"):
            usage()
            sys.exit()
        elif o in ("-o", "--output"):
            output = a
        elif o in ("-i", "--input"):
            input = a
        elif o in ("-p", "--pdb"):
            pdb_url = a

    if len(args) <2 :
        print >>sys.stderr, \
             "Not enough arguments. Try --help for more details."
        sys.exit(2)

    raf_url = args[0]
    cla_url = args[1]
    
    (raf_filename, headers) = urllib.urlretrieve(raf_url)
    seqMapIndex = Raf.SeqMapIndex(raf_filename)

    (cla_filename, headers) = urllib.urlretrieve(cla_url)
    claIndex = Cla.Index(cla_filename)

    if input == None :
        sids = args[2:]
    elif input == '-' :
        sids = sys.stdin.xreadlines()
    else :
        in_handle = open(input)
        sids = in_handle.xreadlines()

    try:
        for sid in sids :
            if not sid or sid[0:1]=='#': continue
            id = sid[0:7]
            pdbid=id[1:5]
            s = pdbid[0:1]
            if s=='0' or s=='s' :
                print >>sys.stderr,"No coordinates for domain "+id
                continue

            if output == None :
                filename = id+".ent"
                out_handle = open(filename, "w+")
            elif output == '-' :
                out_handle = sys.stdout
            else :
                out_handle = open(output, "w+")
                
            try:
                try:
                    claRec = claIndex[id]
                    residues = claRec.residues
                    seqMap = seqMapIndex.getSeqMap(residues)
                    pdbid = residues.pdbid

                    f = open_pdb(pdbid, pdb_url) 
                    try:
                        seqMap.getAtoms(f, out_handle)
                    finally :
                        f.close()
                except (IOError, KeyError, RuntimeError), e:
                    print >>sys.stderr, "I cannot do SCOP domain ",id,":",e
            finally :
                out_handle.close()
    finally :
        if in_handle != None :
            in_handle.close()
                
                
if __name__ == "__main__":
    main()


















