#!/usr/bin/env python
# Created: Thu Jul 13 12:09:58 2000
# Last changed: Time-stamp: <00/07/13 12:52:49 thomas>
# Thomas.Sicheritz@molbio.uu.se, http://evolution.bmc.uu.se/~thomas
# File: xbb_utils.py

import string, re, regsub
import posixpath, posix
import os, sys  # os.system, sys.argv
sys.path.insert(0, '.')
from Tkinter import *
import Pmw

class NotePad(Toplevel):
    def __init__(self, master= None):
        Toplevel.__init__(self, master)
        self.menubar = Menu(self)
        self.filemenu = Menu(self.menubar)
        self.filemenu.add_command(label = "Save", command = self.save)
        self.filemenu.add_separator()
        self.filemenu.add_command(label = "Dismiss", command = self.destroy)
    
	self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.configure(menu = self.menubar)
        self.tid = Pmw.ScrolledText(self)
        self.tid.pack(fill = BOTH, expand = 1)


    def text_id(self): return self.tid
    def insert(self,start, txt):
        self.tid.insert(start, txt)
        
    def save(self):
        fd = SaveFileDialog(self)
        file = fd.go(key="test")
        if file:
            fid = open(file, 'w')
            fid.write(self.tid.get(0.0,END))
            fid.close()
