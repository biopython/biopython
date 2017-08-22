#!/usr/bin/env python
# Copyright 2000 by Thomas Sicheritz-Ponten.
# Copyrigth 2016 by Markus Piotrowski.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Created: Sun Dec  3 13:38:52 2000
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas

"""Search code for graphical Xbbtools tool."""

import re


try:  # Python 2
    import Tkinter as tk
    import ttk
    import tkColorChooser as colorchooser
except ImportError:  # Python 3
    import tkinter as tk
    import tkinter.ttk as ttk
    from tkinter import colorchooser

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import reverse_complement

import xbb_widget


class DNAsearch(object):
    def __init__(self):
        self.init_alphabet()
        self.sequence = ''

    def init_alphabet(self):
        self.alphabet = ambiguous_dna_values
        other = ''.join(self.alphabet)
        self.alphabet['N'] = self.alphabet['N'] + other
        for key in self.alphabet:
            if key == 'N':
                continue
            if key in self.alphabet[key]:
                continue
            self.alphabet[key] = self.alphabet[key] + key

    def SetSeq(self, seq):
        self.sequence = seq

    def SetPattern(self, pattern):
        self.pattern = pattern
        self.rx_pattern = self.IUPAC2regex(pattern)
        self.rx = re.compile(self.rx_pattern)

    def IUPAC2regex(self, s):
        rx = ''
        for i in s:
            r = self.alphabet.get(i, i)
            if len(r) > 1:
                rx = '%s[%s]' % (rx, r)
            else:
                rx += r
        return rx

    def _Search(self, start=0):
        pos = self.rx.search(self.sequence, start)
        return pos

    def Search(self, start=0):
        pos = self.rx.search(self.sequence, start)
        if pos:
            return pos.start()
        else:
            return -1

    def SearchAll(self):
        pos = -1
        positions = []
        while True:
            m = self._Search(pos + 1)
            if not m:
                break
            pos = m.start()
            if pos == -1:
                break

            positions.append(pos)
        return positions


class XDNAsearch(tk.Toplevel, DNAsearch):
    def __init__(self, seq='', master=None, highlight=0):
        DNAsearch.__init__(self)
        self.master = master
        self.highlight = highlight
        self.colors = []
        self.init_graphics()
        self.sequence = seq
        self.cur_pos = 0

    def init_graphics(self):
        tk.Toplevel.__init__(self, self.master)
        self.frame = ttk.Frame(self)
        self.frame.pack(fill=tk.BOTH, expand=1)

        self.search_entry = ttk.Entry(self.frame)
        self.search_entry.pack(fill=tk.BOTH, expand=1)

        f2 = ttk.Frame(self.frame)
        f2.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        f = f2
        self.forward = ttk.Button(f, text='Search +', command=self.do_search)
        self.forward.pack(side=tk.LEFT)
        self.forward = ttk.Button(
            f, text='Search -',
            command=lambda x=self.do_search: x(other_strand=1))
        self.forward.pack(side=tk.LEFT)
        self.cancel = ttk.Button(f, text='Cancel', command=self.exit)
        self.cancel.pack(side=tk.LEFT)
        self.current_color = 'cyan'
        self.colorb = ttk.Button(f, text='Color', command=self.change_color)
        self.colorb.pack(side=tk.LEFT)
        self.config_color(self.current_color)

    def config_color(self, color=None):
        if not self.highlight:
            return
        if not color:
            color = colorchooser.askcolor()[1]
            if not color:
                color = 'cyan'
        self.current_color = color
        self.current_tag = 'searched_%s' % self.current_color
        self.master.tag_config(self.current_tag, background=self.current_color)
        self.master.tag_config(self.current_tag + 'R',
                               background=self.current_color, underline=1)
        self.colors.append(color)

    def change_color(self):
        self.config_color()

    def get_pattern(self):
        pattern = self.search_entry.get()
        return pattern

    def do_search(self, other_strand=0):
        pattern = self.get_pattern()
        if other_strand:
            pattern = reverse_complement(pattern)
        self.SetPattern(pattern)
        pos = self.Search(self.cur_pos)
        self.cur_pos = pos + 1
        w = self.master
        if pos != -1:
            if self.highlight:
                start, stop = pos, pos + len(self.pattern)
                if other_strand:
                    w.tag_add(self.current_tag + 'R', '1.%d' % start,
                              '1.%s' % stop)
                else:
                    w.tag_add(self.current_tag, '1.%d' % start, '1.%s' % stop)
                w.see('1.%d' % start)

    def exit(self):
        for c in self.colors:
            self.master.tag_remove('searched_%s' % c, 1.0, tk.END)
            self.master.tag_remove('searched_%sR' % c, 1.0, tk.END)
        self.destroy()
        del(self)


if __name__ == '__main__':
    win = tk.Tk()
    xbbtools = xbb_widget.xbb_widget()

    seq = 'ATGGTGTGTGTGTACGATCGCCCCCCCCAGTCGATCGATGCATCGTA'
    xbbtools.insert_sequence(('Test_seq', seq))
    xbbtools.search()

    win.mainloop()
