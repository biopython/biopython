#!/usr/bin/python
#
# PDBList.py
#
# A tool for tracking changes in the PDB Protein Structure Database.
#
# (c) 2003 Kristian Rother
# This work was supported by the German Ministry of Education
# and Research (BMBF). Project http://www.bcbio.de
# 
# Contact the author
#    homepage : http://www.rubor.de/bioinf
#    email    : kristian.rother@charite.de
#
#
# This Code is released under the conditions of the Biopython license.
# It may be distributed freely with respect to the original author.
# Any maintainer of the BioPython code may change this notice
# when appropriate.
#
# Last modified on Tue, Oct 21st 2003, Berlin
#

__doc__="Access the PDB over the internet (for example to download structures)."

import urllib,string,re,os,sys

class PDBList:
    """
    This class provides quick access to the structure lists on the
    PDB server or its mirrors. The structure lists contain
    four-letter PDB codes, indicating that structures are
    new, have been modified or are obsolete. The lists are released
    on a weekly basis.

    It also provides a function to retrieve PDB files from the server.
    To use it properly, prepare a directory /pdb or the like,
    where PDB files are stored.

    If You want to use this module from inside a proxy, add
    the proxy variable to Your environment, e.g. in Unix
    export HTTP_PROXY='http://realproxy.charite.de:888'    
    (This can also be added to ~/.bashrc)
    """
    
    PDB_REF="""
    The Protein Data Bank: a computer-based archival file for macromolecular structures.
    F.C.Bernstein, T.F.Koetzle, G.J.B.Williams, E.F.Meyer Jr, M.D.Brice, J.R.Rodgers, O.Kennard, T.Shimanouchi, M.Tasumi
    J. Mol. Biol. 112 pp. 535-542 (1977)
    http://www.pdb.org/.
    """

    def __init__(self,server='ftp://ftp.pdb.mdc-berlin.de',pdb='/pdb'):
        """Initialize the class with the default server or a custom one."""
        # remote pdb server
        self.pdb_server = server

        # local pdb file tree
        self.local_pdb = pdb

        
    def get_recent_filenames(self):
        """Returns names of the newest three weekly files (added,mod,obsolete).
        
        Reads the directories with changed entries from the PDB server and
        returns a tuple of three URL's to the files of new, modified and
        obsolete entries from the most recent list. The directory with the
        largest numerical name is used.
        Returns None if something goes wrong.
        
        Contents of the data/status dir (20031013 would be used);
drwxrwxr-x   2 1002     sysadmin     512 Oct  6 18:28 20031006
drwxrwxr-x   2 1002     sysadmin     512 Oct 14 02:14 20031013
-rw-r--r--   1 1002     sysadmin    1327 Mar 12  2001 README
        """     
        url = urllib.urlopen(self.pdb_server+'/pub/pdb/data/status/')
        file = url.readlines()
        maxdate = 0
        for l in file:
            try:
                # check if this is a valid date
                date = int(l[54:62])
                if date > maxdate: maxdate = date
            except:
                pass
        if maxdate>0:
            return self.get_filenames_by_date(str(maxdate))
        return None

    def get_filenames_by_date(self,date):
        """Returns URL's of three weekly files (added,mod,obsolete).
        
        Returns names of added, modified and obsolete pdb status files
        for a given date, e.g. '20031013'.
        """
        path = self.pdb_server+'/pub/pdb/data/status/%s/'%(date)
        return (path+'added.pdb',path+'modified.pdb',path+'obsolete.pdb')

    def get_list(self,url):
        """Retrieves a list of pdb codes from the given URL.
        
        Returns a list of pdb codes in the pdb status file with the
        given URL. The URLs are created by get_recent_filenames() or
        get_filenames_by_date(date).
        
        Typical contents of the list files parsed by this method;
-rw-r--r--   1 rcsb     rcsb      330156 Oct 14  2003 pdb1cyq.ent
-rw-r--r--   1 rcsb     rcsb      333639 Oct 14  2003 pdb1cz0.ent
        """
        url = urllib.urlopen(url)
        file = url.readlines()
        list = []
        for l in file:
            try:
                if l[61:65] == '.ent': list.append(l[57:61]) 
            except:
                pass
        return list

    def get_all_obsolete(self):
        """Returns a list of all obsolete entries ever in the PDB.

        Returns a list of all obsolete pdb codes that have ever been
        in the PDB.
        
        Gets and parses the file from the PDB server in the format
        (the first pdb_code column is the one used).
 LIST OF OBSOLETE COORDINATE ENTRIES AND SUCCESSORS
OBSLTE     30-SEP-03 1Q1D      1QZR
OBSLTE     26-SEP-03 1DYV      1UN2    
        """
        url = urllib.urlopen(self.pdb_server+'/pub/pdb/data/status/obsolete.dat')
        file = url.readlines()
        obsolete = []
        for l in file:
            if l[:6] == 'OBSLTE':
                pdb_code = l[21:25]
                obsolete.append(string.lower(pdb_code))

        return obsolete

    def changed_this_week(self):
        """Returns 3 lists of new/modified/obsolete PDB entries for weekly updates.
        
        Returns all three lists (new, modified, obsolete) pdb codes
        for this week.
        Uses get_recent_status() and get_list() for that.
        """
        urls = self.get_recent_filenames()
        tw = []
        tw.append(self.get_list(urls[0]))
        tw.append(self.get_list(urls[1]))
        tw.append(self.get_list(urls[2]))
        return tw

    def retrieve_pdb_file(self,pdb_code, compression='.Z', uncompress="gunzip", write=1):
        """Retrieves a PDB structure file from the PDB server and
        stores it in a local file tree (if 'write' is set true).
        The PDB structure is returned as a single string.
        The compression should be '.Z' or '.gz'. 'uncompress' is
        the command called to uncompress the files.
        """
        # get the structure
        code=string.lower(pdb_code)
        url = self.pdb_server+'/pub/pdb/data/structures/divided/pdb/%s/pdb%s.ent%s'%(code[1:3],code,compression)
        lines = urllib.urlopen(url).read()

        # save the structure
        path = self.local_pdb + os.sep + code[1:3]
        filename = path + os.sep+"pdb%s.ent%s"%(code,compression)
        if not os.access(path,os.F_OK):
            os.mkdir(path)
        open(filename,'w').write(lines)

        # uncompress the file
        os.system("%s %s" % (uncompress, filename))

        return lines
            

    def update_pdb(self):
        """
        I guess this is the 'most wanted' function from this module.
        It gets the weekly lists of new and modified pdb entries and
        automatically downloads the according PDB files.
        You can call this module as a weekly cronjob.
        """
        changes  = self.changed_this_week()
        new      = changes[0]
        modified = changes[1]
        to_download = new + modified

        for pdb_code in to_download:
            self.retrieve_pdb_file(pdb_code)
            try:
                print 'retrieving %s'%(pdb_code)            
            except:
                print 'error %s'%(pdb_code)
                # you can insert here some more log notes that
                # something has gone wrong.            

        #
        # delete the obsolete files
        #    this part could easily misbehave, so i commented it out.
        #
        # obsolete = changes[2]
        # for pdb_code in obsolete:
        #     file = self.local_pdb + os.sep + pdb_code[1:3] + os.sep + 'pdb%s.ent'%(pdb_code)
        #     os.remove(file)


if __name__ == '__main__':
    doc = """PDBList.py
    (c) Kristian Rother 2003, Contributed to BioPython

    Standalone usage

    PDBList.py update - write weekly PDB updates to local /pdb tree.
    PDBList.py        - simple usage examples.
    
    """
    print doc

    if len(sys.argv)>1:
        # update PDB 
        if sys.argv[1] == 'update':
            pl = PDBList()
            pl.update_pdb()
            sys.exit(0)

    #
    # usage example
    #
    
    # 1. create object
    pl = PDBList()

    # 2. get all obsolete structure codes
    print "\nAll obsolete structures from the PDB server:"
    obsolete = pl.get_all_obsolete()
    print string.join(obsolete,'    ')

    # 3. get the weekly updated lists
    changes  = pl.changed_this_week()
    print "\nThis weeks new structures:"    
    new      = changes[0]
    print string.join(new,'    ')

    print "\nThis weeks modified structures:"    
    modified = changes[1]
    print string.join(modified,'    ')

    print "\nThis weeks obsolete structures:"    
    obsolete = changes[2]
    print string.join(obsolete,'    ')
    
    
