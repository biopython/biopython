#!/usr/bin/env python
# Created: Thu Jul 13 14:07:25 2000
# Last changed: Time-stamp: <00/07/13 17:53:37 thomas>
# Thomas.Sicheritz@molbio.uu.se, http://evolution.bmc.uu.se/~thomas
# File: xbb_blast.py

import string, re, regsub
import posixpath, posix
import os, sys, glob
from threading import *
import commands
from Tkinter import *
import Pmw
sys.path.insert(0, '.')


class xbb_blast_thread(Thread):
    def __init__(self, command, callback):
        Thread.__init__(self, name="Blaster")
        self.result = ''
        self.callback = callback
        self.command = command


    def run(self):
        res = commands.getoutput(self.command)
        self.callback(res)
        

class blast_popup(Toplevel):
    def __init__(self, getseq = None, callback = None):
        Toplevel.__init__(self)
        self.callback = callback
        self.getseq = getseq
        self.blastdb = os.environ.get('BLASTDB', './')
        self.get_blast_databases()
        self.blast_aa_progs = ('blastp', 'blastx')
        self.blast_nt_progs = ('blastn', 'tblastx')
        self.blast_dbs = self.aa_dbs + self.nt_dbs
        self.db = StringVar()
        self.blast_program = StringVar()
        self.blast_program.set('blastp')
        self.method_menu = Pmw.OptionMenu(self,
                                          labelpos = 'w',
                                          label_text = 'Program:',
                                          menubutton_textvariable = self.blast_program,
                                          items = self.blast_aa_progs + self.blast_nt_progs,
                                          menubutton_width = 10,
                                          command = self.check_combination,
	)

	self.method_menu.pack(anchor = 'w', padx = 10, pady = 10, side = LEFT)
        self.db_menu = Pmw.OptionMenu(self,
                                      labelpos = 'w',
                                      label_text = 'Database:',
                                      menubutton_textvariable = self.db,
                                      menubutton_width = 10,
                                      command = self.check_combination,
                                      )
	self.db_menu.pack(anchor = 'w', padx = 10, pady = 10, side = LEFT)
        self.db_menu.setitems(self.blast_dbs,0)


        self.ok = Button(self, text = 'Doit', command = self.run)
        self.ok.pack(padx = 10, pady = 10, side = LEFT)

        self.check_combination(None)

        
    def run(self, event = None):
        seq = self.getseq()
        com = 'echo %s | blastall -d %s -p %s' % (seq, self.db.get(), self.blast_program.get())
        print com
        BI = xbb_blast_thread(com, self.callback)
        BI.start()
        
    def get_old_blast_databases(self):
        self.old_aa_dbs = glob.glob(self.blastdb + '/*.atb')
        self.old_nt_dbs = glob.glob(self.blastdb + '/*.ntb')

    def get_blast_databases(self):
        self.aa_dbs = map(lambda x: string.split(os.path.basename(x),'.')[0],glob.glob(self.blastdb + '/*.pin'))
        self.nt_dbs = map(lambda x: string.split(os.path.basename(x),'.')[0],glob.glob(self.blastdb + '/*.nin'))

    def check_combination(self, event):
        color = 'red'
        prog = self.blast_program.get()
        db = self.db.get()
        p_is_aa = prog in self.blast_aa_progs
        p_is_nt = prog in self.blast_nt_progs
        db_is_nt = db in self.nt_dbs
        db_is_aa = db in self.aa_dbs

        if p_is_nt and p_is_aa: color = 'green'
        elif p_is_aa and db_is_aa: color = 'green'
        elif p_is_nt and db_is_nt: color = 'green'
        else: color ='red'
        self.ok.configure(bg = color)
        if color == 'red':
            self.ok.configure(state = 'disabled')
        else:
            self.ok.configure(state = 'normal')
            
        
    def check_progs(self, event):
        if self.blast_dbs == self.aa_dbs:
            if not self.blast_program in self.blast_aa_progs:
                self.method_menu.setitems(items = self.blast_aa_progs)
        else:
            if not self.blast_program in self.blast_nt_progs:
                self.method_menu.setitems(items = self.blast_nt_progs)
        

if __name__ == '__main__':
    os.environ['BLASTDB'] = '/opt/bio/data'
    test = blast_popup()
    
