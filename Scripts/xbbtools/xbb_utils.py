#!/usr/bin/env python
# Copyright 2000 by Thomas Sicheritz-Ponten.
# Copyright 2016 by Markus Piotrowski.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Created: Thu Jul 13 12:09:58 2000
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas

"""Utility code for graphical Xbbtools tool."""

import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog


class NotePad(tk.Toplevel):
    """Top level window for results (translations, BLAST searches...)."""

    def __init__(self, master=None):
        """Set up notepad window."""
        tk.Toplevel.__init__(self, master)
        self.menubar = tk.Menu(self)
        self.filemenu = tk.Menu(self.menubar)
        self.filemenu.add_command(label="Save", command=self.save)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Dismiss", command=self.destroy)

        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.configure(menu=self.menubar)
        self.yscroll = ttk.Scrollbar(self, orient="vertical")
        self.tid = tk.Text(self, width=88, yscrollcommand=self.yscroll.set)
        self.yscroll.configure(command=self.tid.yview)
        self.tid.pack(side="left", fill="both", expand=1)
        self.yscroll.pack(side="right", fill="y")

    def text_id(self):
        """Get reference to notepad window."""
        return self.tid

    def insert(self, start, txt):
        """Add text to notepad window."""
        self.tid.insert(start, txt)

    def save(self):
        """Save text from notepad to file."""
        filename = filedialog.asksaveasfilename()
        if filename:
            with open(filename, "w") as fid:
                fid.write(self.tid.get(0.0, "end"))
