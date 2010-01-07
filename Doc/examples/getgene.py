#!/usr/bin/env python
# Created: Sun Oct 15 16:16:20 2000
# Last changed: Time-stamp: <01/02/15 09:01:27 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: getgene.py

""" Example code to index a non-reduntant protein database of
    SwissProt + TrEMBL for fast lookup and retrieval.
    
    To build the database and index it:
        cd /opt/bio/data/
        wget -N -nd -r -l1 -A'.dat.Z' ftp://expasy.cbr.nrc.ca/databases/sp_tr_nrdb/
        zcat *.dat.Z > nr.dat
        ./getgene.py --index nr.dat
        setenv PYPHY '/opt/bio/data'
        
    To retrieve entries from the command line:
    ./getgene.py EFTU_ECOLI

    To use from a python script:
    from getgene import DB_Index

    db_index = DB_Index()

    # retrieve a complete entry:
    db_index.Get('EFTU_ECOLI')

    # get organism, lineage and gene
    db_index.Get_OS_OC_GN('EFTU_ECOLI')
"""    
    


import string, re
import os, sys
import gdbm

class DB_Index:
    def __init__(self, open = 1):
        if open: self.Open()
            
    def Create(self, infile, outfile):
        db = gdbm.open(outfile, 'n')
        fid = open(infile)

        db['datafile'] = os.path.abspath(infile)
        
        while 1:
            line = fid.readline()
            if not line or not len(line): break

            if line[:3] == 'ID ':
                id = string.split(line)[1]
                start = fid.tell() - len(line)

            elif line[:3] == 'AC ':
                acc = string.split(line)[1]
                if acc[-1] ==';': acc = acc[:-1]

            elif line[:2] =='//':
                stop = fid.tell()
                try:
                    value = '%d %d' % (start, stop)
                    db[id] = value
                    db[acc] = value
                    id, acc, start, stop = None, None, None, None
                except:
                    print 'AARRGGGG', start, stop, type(start), type(stop)
                    print id, acc
                    
        db.close()
        fid.close()

    def Open(self, indexfile = None):
        if not indexfile:
            indexfile = os.path.join(os.environ['PYPHY'],'nr.dat.indexed')

        self.db = gdbm.open(indexfile)
        self.datafile = self.db['datafile']
        self.fid = open(self.datafile)

    def Close(self):
        self.db.close()
        
    def Get(self, id):
        try:
            values = self.db[id]
        except:
            return None
        start, stop= map(int,string.split(values))
        self.fid.seek(start)
        txt = self.fid.read(stop - start)
        return txt
        
    def Get_Organism(self, id):
        entry = self.Get(id)
        if not entry: return None
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
                return OS
            if line[0:2] =="//": break
        return OS

    def FixOS(self, os):
        os = string.split(os,',')[0]
        os = string.split(os,'(')[0]
        return string.strip(os)
    
    def Get_Taxonomy(self, id):
        entry = self.Get(id)
        if not entry: return None
        OC = ""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OC   ':
                OC = OC + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:2] =="//": break
        return OC
    
    def Get_Kingdom(self, id):
        res = self.Get_Taxonomy(id)
        #print id, res
        if not res: return "U"
        kd = string.strip(string.split(res,";")[0])
        if kd == "Eubacteria" or kd == "Prokaryota" or kd == "Bacteria": return "B"
        elif kd == "Eukaryota" or kd =="Eukaryotae": return "E"
        elif kd == "Archaebacteria" or kd == "Archaea": return "A"
        elif kd == "Viridae" or kd == "Viruses": return "V"
        else:
            print kd, "UNKNOWN"
            return "U"
        
    def Get_Gene(self, id):
        entry = self.Get(id)
        if not entry: return None
        GN = ''
        for line in string.split(entry, '\n'):
            if line[0:5] == 'GN   ':
                GN = string.strip(line[5:])
                if GN[-1] ==".": GN = GN[0:-1]
                return GN
            if line[0:2] =="//": break
        return GN


    def Get_OS_OC_GN(self, id):
        entry = self.Get(id)
        if not entry: return None, None, None
        OS, OC, GN = "","",""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            if line[0:5] == 'OC   ':
                OC = OC + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:5] == 'GN   ':
                GN = string.strip(line[5:])
                if GN[-1] ==".": GN = GN[0:-1]
            if line[0:2] =="//": break
        return OS, OC, GN
    
    def Get_OS_OC_OG(self, id):
        entry = self.Get(id)
        if not entry: return None, None, None
        OS, OC, OG = "","",""
        for line in string.split(entry, '\n'):
            if line[0:5] == 'OS   ':
                OS = string.strip(line[5:])
                if OS[-1] ==".": OS = OS[0:-1]
            if line[0:5] == 'OC   ':
                OC = OC + string.strip(line[5:])
                if OC[-1] ==".": OC = OC[0:-1]
            if line[0:5] == 'OG   ':
                OG = string.strip(line[5:])
                if OG[-1] ==".": OG = OG[0:-1]
            if line[0:2] =="//": break
        return OS, OC, OG

    def Get_SQ(self, id, fasta = 1):
        entry = self.Get(id)
        if not entry: return ""
        SQ = ""
        record = 0
        for line in string.split(entry, '\n'):
            if record: SQ = SQ + string.strip(line[5:])
            if line[0:5] == 'SQ   ': record = 1
            if line[0:2] =="//": break

        SQ = re.sub('[ \n]','',SQ)
        if fasta: SQ = '>%s\n%s' % (id, re.sub('(.{60})','\\1\n',SQ))
        return SQ

    def Get_XX(self, id, xx):
        entry = self.Get(id)
        if not entry: return ""
        XX = ""
        for line in string.split(entry, '\n'):
            if line[0:5] == '%s   ' % xx:
                XX = XX + string.strip(line[5:])
                if XX[-1] ==".": XX = XX[0:-1]
            if line[0:2] =="//": break
        return XX
        
    def Get_Keywords(self, id):
        entry = self.Get(id)
        if not entry: return []
        keywords = []
        for line in string.split(entry, '\n'):
            if line[0:5] == 'KW   ':
                for i in string.split(string.strip(line[5:]),';'):
                    kw = string.strip(i)
                    if len(kw) < 2: continue
                    if kw[-1] == '.': kw = kw[:-1]
                    keywords.append(kw)
            if line[0:2] =="//": break
        return keywords



def help(exit = 0):
    name = os.path.basename(sys.argv[0])
    print 'Usage: %s <db> <gene ID>' % name
    print '  or   %s --index <db.dat>' % name
    if exit: sys.exit(0)

if __name__ == '__main__':
    pyphy_home = os.environ.get('PYPHY', None)

    if len(sys.argv) == 1: help(exit = 1)
    db_index = DB_Index(open = 0)
    func = db_index.Get
    for arg in sys.argv[1:]:
        if arg == '--index':
            sys.argv.remove(arg)
            infile = sys.argv[1]
            outfile = os.path.basename(infile) + '.indexed'
            db_index.Create(infile, outfile)
            sys.exit(0)

        elif arg[:4] == '-Get':
            sys.argv.remove(arg)
            func = getattr(db_index, arg[1:])

        elif arg == '-h' or arg == '--help': help(exit = 1)

    db = 'nr.dat'
    if len(sys.argv) == 2:
        # shortcut, mostly we use nr.dat so dont bother to name it
        ids = sys.argv[1:]
    else:
        try:
             db = sys.argv[1]
             ids = sys.argv[2:]
        except:
            help(exit = 1)
        
    dbfile = os.path.join(pyphy_home, db + '.indexed')
    db_index.Open(dbfile)
    for id in ids:
        #print db_index.Get(id)
        print func(id)
        


