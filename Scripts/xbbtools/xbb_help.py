#!/usr/bin/env python
# Created: Tue Sep  4 09:05:16 2001
# Last changed: Time-stamp: <01/09/04 09:42:37 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xbb_help.py

from Tkinter import *
from ScrolledText import ScrolledText

class xbbtools_help(Toplevel):
    def __init__(self, *args):
        Toplevel.__init__(self)
        self.tid = ScrolledText(self)
        self.tid.pack(fill = BOTH, expand = 1)
        self.Styles()
        self.Show()

    def Styles(self):
        for c in ['red', 'blue', 'magenta', 'yellow', 'green', 'red4', 'green4', 'blue4']:
            self.tid.tag_configure(c, foreground = c)

        self.tid.tag_config('underline', underline =1)
        self.tid.tag_config('italic', font = ('Courier', 6, 'italic'))
        self.tid.tag_config('bold', font = ('Courier', 8, 'bold'))
        self.tid.tag_config('title', font = ('Courier', 12, 'bold'))
        self.tid.tag_config('small', font = ('Courier', 6, ''))
        self.tid.tag_config('highlight', background = 'gray')
        

    def Show(self):
        t = self.tid
        t.insert(END, "XBBtools Help\n", 'title')
        t.insert(END, """
Copyright 2001 by Thomas Sicheritz-Ponten.  All rights reserved.
This code is part of the Biopython distribution and governed by its
license.  Please see the LICENSE file that should have been included
as part of this package.\n
""", 'italic')
        t.insert(END, 'thomas@biopython.org\n\n', 'blue')
        t.insert(END, '* Goto Field\n', 'bold')
        t.insert(END, '\tinserting one position moves cursor to position\n')
        t.insert(END, "\tinserting two positions, sperated by ':' ")
        t.insert(END, 'highlights', 'highlight')
        t.insert(END, ' selected range\n')
        t.insert(END, '\n')
        t.insert(END, '* Search\n', 'bold')
        t.insert(END, '\tambiguous dna values are\n')
        t.insert(END, """
                A: A
                C: C
                G: G
                T: T
                M: AC
                R: AG
                W: AT
                S: CG
                Y: CT
                K: GT
                V: ACG
                H: ACT
                D: AGT
                B: CGT
                X: GATC
                N: GATC

                """, 'small')

        
