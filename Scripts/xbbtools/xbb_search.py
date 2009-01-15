#!/usr/bin/env python
# Created: Sun Dec  3 13:38:52 2000
# Last changed: Time-stamp: <01/09/04 09:51:21 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xbb_search.py

import re
import os, sys, commands
sys.path.insert(0, '.')
from Tkinter import *
from tkColorChooser import askcolor
from Bio.Data.IUPACData import ambiguous_dna_values
import re

from Bio.Seq import reverse_complement

class DNAsearch:
    def __init__(self):
        self.init_alphabet()
        self.sequence = ''
        
    def init_alphabet(self):
        self.alphabet = ambiguous_dna_values
        other = ''.join(self.alphabet.keys())
        self.alphabet['N'] = self.alphabet['N'] + other
        for key in self.alphabet.keys():
            if key == 'N': continue
            if key in self.alphabet[key]: continue
            self.alphabet[key] = self.alphabet[key] + key

    def SetSeq(self, seq): self.sequence = seq
    
    def SetPattern(self, pattern):
        self.pattern = pattern
        self.rx_pattern = self.IUPAC2regex(pattern)
        self.rx = re.compile(self.rx_pattern)
    
    def IUPAC2regex(self, s):
        rx = ''
        for i in s:
            r = self.alphabet.get(i,i)
            if len(r) > 1:
                rx = '%s[%s]' % (rx, r)
            else:
                rx += r
        return rx
    
    def _Search(self, start = 0):
        pos = self.rx.search(self.sequence, start)
        return pos
    
    def Search(self, start = 0):
        pos = self.rx.search(self.sequence, start)
        if pos:
            return pos.start()
        else:
            return -1
        

    def SearchAll(self):
        pos = -1
        positions = []
        while 1:
            m   = self._Search(pos+1)
            if not m: break
            pos = m.start()
            if pos == -1:
                break
            
            positions.append(pos)
        return positions
    
        
class XDNAsearch(Toplevel, DNAsearch):
    def __init__(self, seq= '', master= None, highlight = 0):
        DNAsearch.__init__(self)
        self.master = master
        self.highlight = highlight
        self.colors = []
        self.init_graphics()
        self.sequence = seq
        self.cur_pos = 0
        
    def init_graphics(self):
        Toplevel.__init__(self, self.master)
        self.frame = Frame(self)
        self.frame.pack(fill = BOTH, expand = 1)

        self.search_entry = Entry(self.frame)
        self.search_entry.pack(fill = BOTH, expand = 1)

        f2 = Frame(self.frame)
        f2.pack(side = TOP, fill = BOTH, expand = 1)
        
        f = f2
        self.forward = Button(f, text = 'Search +', command = self.do_search)
        self.forward.pack(side = LEFT)
        self.forward = Button(f, text = 'Search -',
                              command = lambda x=self.do_search: x(other_strand=1))
        self.forward.pack(side = LEFT)
        self.cancel = Button(f, text = 'Cancel', command = self.exit)
        self.cancel.pack(side = LEFT)
        self.current_color = 'cyan'
        self.colorb = Button(f, text = 'Color', command = self.change_color, foreground = self.current_color)
        self.colorb.pack(side = LEFT)
        self.config_color(self.current_color)

        

        
    def config_color(self, color = None):
        if not self.highlight: return
        if not color:
            try:
                color = askcolor()[1]
            except:
                color = 'cyan'
        self.current_color = color
        self.current_tag = 'searched_%s' % self.current_color
        self.master.tag_config(self.current_tag, background=self.current_color)
        self.master.tag_config(self.current_tag+'R', background=self.current_color, underline = 1)
        self.colors.append(color)
        
    def change_color(self):
        self.config_color()
        self.colorb.configure(foreground = self.current_color)
        self.colorb.update()
            
    def get_pattern(self):
        pattern = self.search_entry.get()
        return pattern
        
    def do_search(self, other_strand = 0):
        pattern = self.get_pattern()
        if other_strand: pattern = reverse_complement(pattern)
        self.SetPattern(pattern)
        pos = self.Search(self.cur_pos)
        self.cur_pos = pos +1
        w = self.master
        if pos != -1:
            if self.highlight:
                start, stop = pos, pos + len(self.pattern)
                if other_strand:
                    w.tag_add(self.current_tag+'R', '1.%d' % start, '1.%s' % stop)
                else:
                    w.tag_add(self.current_tag, '1.%d' % start, '1.%s' % stop)
                w.see('1.%d' % start)


    
    def exit(self):
        for c in self.colors:
            self.master.tag_remove('searched_%s' % c, 1.0, END)
            self.master.tag_remove('searched_%sR' % c, 1.0, END)
        self.destroy()
        del(self)
    
    def showcolor(self):
        pass

    

if __name__ == '__main__':
    seq = 'ATGGTGTGTGTGTACGATCGCCCCCCCCAGTCGATCGATGCATCGTA'
    win = Tk()
    xtest = XDNAsearch(seq = seq, master = win)

    win.mainloop()
    
    
        
