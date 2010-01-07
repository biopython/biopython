#!/usr/bin/env python
# Created: Thu Jul 13 12:09:58 2000
# Last changed: Time-stamp: <00/12/03 12:10:58 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xbb_utils.py

import sys
sys.path.insert(0, '.')
from Tkinter import *
from FileDialog import SaveFileDialog

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
        self.yscroll = Scrollbar(self,orient=VERTICAL)
        self.tid = Text(self, yscrollcommand = self.yscroll.set)
        self.yscroll.configure(command = self.tid.yview)
        self.tid.pack(side = LEFT, fill = BOTH, expand = 1)
        self.yscroll.pack(side = RIGHT, fill = Y)
        

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
