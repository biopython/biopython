#!/usr/bin/env python
# Created: Wed Jun 21 10:28:14 2000
# Last changed: Time-stamp: <01/09/04 09:42:06 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xbb_widget.py

# Copyright 2000 by Thomas Sicheritz-Ponten.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

import re
import sys
import time

from Tkinter import *
from tkFileDialog import askopenfilename, asksaveasfilename

sys.path.insert(0, '.')
from xbb_utils import *
from xbb_translations import xbb_translations
from xbb_blast import BlastIt
from xbb_search import XDNAsearch
from xbb_help import xbbtools_help
from Bio.Data import CodonTable
from Bio.SeqUtils import quick_FASTA_reader


class xbb_widget:
    def __init__(self, parent = None):
        self.is_a_master = (parent == None)
        self.parent = parent
            
        self.init_variables()
        self.init_colors()
        # master frame
        self.main_frame = Frame(parent)
        if not parent:
            self.init_optionsdb()
            self.parent = self.main_frame.master
            
        self.main_frame.pack(fill = BOTH, expand = 1)
        
        # sequence info (GC%, positins etc.)
        self.info_frame = Frame(self.main_frame)
        self.info_frame.pack(fill = BOTH, expand = 1)
        
        self.create_menu(self.info_frame)
        self.create_seqinfo(self.info_frame)

        # sequence field and fast buttons
        self.seq_frame = Frame(self.main_frame)
        self.seq_frame.pack(fill = BOTH, expand = 1)
        
        self.create_buttons(self.seq_frame)
        self.create_seqfield(self.seq_frame)

        self.create_bindings()
        self.blastit = 'xbb_blast.py'
        
    def init_variables(self):
        self.seqwidth = 60
        self.translation_tables = {}
        for i, table in CodonTable.unambiguous_dna_by_id.iteritems():
            self.translation_tables[table.names[0]] = i
        self.translator = xbb_translations()

    def init_colors(self):
        self.colorsbg = {'frame':'brown',
                         'canvas':'brown',
                         'button':'darkgreen',
                         'radiobutton':'darkgrey',
                         'checkbutton':'darkgrey',
                         'label':'dimgrey',
                         'text':'bisque1',
                         'entry':'bisque1',
                         'menu':'darkgreen',
                         'menubutton':'darkgreen',
                         'seqinfo':'dimgrey'}
        self.colorsfg = {'button':'lightblue',
                         'radiobutton':'lightblue',
                         'checkbutton':'lightblue',
                         'label':'green3',
                         'text':'black',
                         'entry':'black',
                         'menu':'black',
                         'menubutton':'lightblue'}
        self.colorsNT = {'A':'green',
                         'C':'lightblue',
                         'G':'orange',
                         'T':'tomato'
                         }

        self.colorsPMWbg = {
            'ComboBox':'darkgreen',
            'ComboBox.Label':'darkgreen',
            }
        
        self.colorsPMWfg = {
            'ComboBox.Label':'lightblue',
            }


    def init_optionsdb(self):
        # does anybody know a better way of defining colors ?
        # how would one implement Tk's -class ?
        tk = self.main_frame.master
        for k,v in self.colorsbg.items():
            name = '*' + k[0].upper() + k[1:] + '.background'
            tk.option_add(name, v)

        for k,v in self.colorsfg.items():
            name = '*' + k[0].upper() + k[1:] + '.foreground'
            tk.option_add(name, v)
            
        for k,v in self.colorsPMWbg.items():
            name = '*' + k[0].upper() + k[1:] + '.background'
            tk.option_add(name, v)

        for k,v in self.colorsPMWfg.items():
            name = '*' + k[0].upper() + k[1:] + '.foreground'
            tk.option_add(name, v)
            
    def create_menu(self, parent):
        self.menubar = Menu(self.main_frame)
        
        # File menu
        self.file_menu = Menu(self.menubar)
        menu = self.file_menu
        menu.add_command(label='Exit', command = self.exit)
        self.menubar.add_cascade(label="File", menu=self.file_menu)

        # Edit menu
        self.edit_menu = Menu(self.menubar)
        menu = self.edit_menu
        menu.add_command(label='Complement', command = self.complement)
        menu.add_command(label='Antiparallel', command = self.antiparallel)
        menu.add_command(label='Reverse', command = self.reverse)
        menu.add_command(label='Fix sequence', command = self.fix_sequence)
        menu.add_command(label='Search', command = self.search)
        self.menubar.add_cascade(label="Edit", menu=self.edit_menu)

        # Translation menu
        self.translation_menu = Menu(self.menubar)
        menu = self.translation_menu
        menu.add_command(label='+1 Frame', command = self.translate)
        menu.add_command(label='6 Frames', command = self.gcframe)
        menu.add_command(label='Extract to FASTA', command = self.extract)

        self.current_codon_table = StringVar()
        self.current_codon_table.set('Standard')
        self.current_codon_table_id = 1
        
        keys = self.translation_tables.keys()
        keys.remove('Standard')
        keys.sort()
        keys = ['Standard'] + keys

        self.gencode_menu = Menu(self.translation_menu)
        menu = self.gencode_menu
        for table in keys:
            menu.add_radiobutton(label=table, command = self.set_codon_table, variable = self.current_codon_table)
        self.translation_menu.add_cascade(label="Genetic Codes", menu=self.gencode_menu)


        self.menubar.add_cascade(label="Translations", menu=self.translation_menu)

        # Tools menu
        self.tools_menu = Menu(self.menubar)
        menu = self.tools_menu
        menu.add_command(label='Blast', command = self.blast)
        menu.add_command(label='Stats', command = self.statistics)
        self.menubar.add_cascade(label="Tools", menu=self.tools_menu)
        
        # Help menu
        self.help_menu = Menu(self.menubar, name = 'help')
        menu = self.help_menu
        menu.add_command(label='Help', command = xbbtools_help)
        self.menubar.add_cascade(label="Help", menu=self.help_menu)

        self.parent.config(menu = self.menubar)

    def set_codon_table(self):
        self.current_codon_table_id = self.translation_tables[self.current_codon_table.get()]
        
    def exit(self, *args):
        # depending on if this widget is the first created or a child widget
        if self.is_a_master:
            sys.exit(0)
        else:
            self.main_frame.destroy()
            
    def create_seqinfo(self, parent):
        # all the sequence information in the top labels
        self.seq_info1 = Frame(parent, relief = RIDGE,
                               borderwidth = 5, height = 30)
        self.seq_info1.pack(fill = BOTH, expand = 1, side = TOP)

        self.position_ids = {}
        d = self.position_ids
        d['id'] = Label(self.seq_info1, width = 10)
        d['from_id'] = Label(self.seq_info1, width = 10)
        d['to_id'] = Label(self.seq_info1, width = 10)
        d['length_id'] = Label(self.seq_info1, width = 10)
        d['label'] = Label(self.seq_info1, width = 10)
        for i in ['id', 'from_id', 'to_id', 'length_id', 'label']:
            d[i].pack(side = LEFT, fill = BOTH, expand = 1)


        
        self.seq_info2 = Frame(parent, relief = RIDGE,
                               borderwidth = 5, height = 30)
        self.seq_info2.pack(fill = BOTH, expand = 1, side = TOP)
        self.statistics_ids = {}
        d = self.statistics_ids
        d['length_id'] = Label(self.seq_info2, width = 10)
        d['length_id'].pack(side = LEFT, fill = BOTH, expand = 1)
        for nt in ['A','C','G','T']:
            d[nt] = Label(self.seq_info2, width = 10, fg = self.colorsNT[nt])
            d[nt].pack(side = LEFT, fill = BOTH, expand = 1)
            

        
    def create_buttons(self, parent):
        self.button_frame = Frame(parent)
        self.button_frame.pack(fill = Y, side = LEFT)
        self.buttons = {}
        for text, func in [('Open', self.open),
                           ('Export', self.export),
                           ('GC Frame', self.gcframe),
                           ('Blast', self.blast),
                           ('Exit', self.exit)]:
            b_id = Button(self.button_frame, text = text,
                          command = func, width = 7)
            b_id.pack(side = TOP, pady = 5, padx = 10)
            self.buttons[text] = b_id

        f = Frame(self.button_frame)
        l = Label(f, text = 'Goto:', bg = self.colorsbg['frame'], fg = self.colorsfg['button'])
        l.pack(side = LEFT)
        l.bind('<Button-1>', self.goto)
        
        self.goto_entry = Entry(f, width = 5)
        self.goto_entry.pack(side = RIGHT, pady = 5, padx = 4)
        self.goto_entry.bind('<Return>', self.goto)
        f.pack(side = BOTTOM)
        
    def create_seqfield(self, parent):
        self.sequence_id = Text(parent, wrap = 'char',
                                width = self.seqwidth)
        self.sequence_id.pack(fill = BOTH, expand = 1, side = RIGHT)

    def create_bindings(self):
        self.sequence_id.bind('<Motion>', self.position)
        self.sequence_id.bind('<Leave>', lambda x,s = self:
                              s.position_ids['id'].configure(text = ''))
        self.sequence_id.bind('<1>', self.zero)
        self.sequence_id.bind('<B1-Motion>', self.count_selection)
        self.sequence_id.bind('<Double-Button-1>', self.select_all)
        
    def zero(self, event):
        p = self.position_ids
        for i in ['from_id', 'to_id', 'length_id']:
            self.position_ids[i].configure(text = '')

    def get_length(self):
        self.sequence_length = len(self.sequence_id.get(1.0,END))
        return self.sequence_length
    
    def select_all(self, event):
        self.select(1, self.get_length())
        self.count_selection(None)
        
    def select(self, a, b):
        w = self.sequence_id
        w.selection_own()
        w.tag_add('sel', '1.%d' % (a - 1), '1.%d' % b)
        self.count_selection(None)

    def get_selection_or_sequence(self):
        w = self.sequence_id
        seq = self.get_selection()
        if not len(seq):
            seq = self.sequence_id.get(1.0,END)

        seq = re.sub('[^A-Z]','',seq)    
        return seq
    
    def get_selection(self):
        w = self.sequence_id
        #print w.selection_own()
        #w.selection_own()
        try:
            return w.selection_get()
            #return string.upper(w.get(sel.first, sel.last))
        except:
            return ''
        
    def get_self_selection(self):
        w = self.sequence_id
        #w.selection_own()
        try:
            return w.selection_get()
            #return string.upper(w.get(sel.first, sel.last))
            #return string.upper(w.selection_own_get())
        except:
            return ''
        
    def count_selection(self, event):
        w = self.sequence_id
        w.selection_own()
        try:
            a = int(w.index('sel.first').split('.')[1]) +1
            b = int(w.index('sel.last').split('.')[1])
            length = b - a + 1

            self.position_ids['from_id'].configure(text = 'Start:%d'% a)
            self.position_ids['to_id'].configure(text = 'Stop:%d'% b)
            self.position_ids['length_id'].configure(text = '%d nt' % length)

            self.statistics_ids['length_id'].configure(text = 'Length=%d' % length)
            seq = self.get_self_selection()
            for nt in ['A','C','G','T']:
                n = seq.count(nt)
                self.statistics_ids[nt].configure(text = '%s=%d' % (nt,n))
                
            
        except:
            pass
        
    def position(self, event):
        x = event.x
        y = event.y
        pos = self.sequence_id.index('@%d,%d' % (x,y)).split('.')
        pos = int(pos[1]) + 1
        self.position_ids['id'].configure(text = str(pos))
        
    def open(self, file = None):
        if not file:
            file = askopenfilename()
        if not file: return
        genes = quick_FASTA_reader(file)
        self.insert_sequence(genes[0])

    def insert_sequence(self, (name, sequence)):
        self.sequence_id.delete(0.0, END)
        self.sequence_id.insert(END, sequence.upper())
        self.fix_sequence()
        self.update_label(name)

    def fix_sequence(self):
        seq = self.sequence_id.get(1.0,END)
        seq = seq.upper()
        seq = re.sub('[^A-Z]','',seq)
        self.sequence_id.delete(0.0,END)
        self.sequence_id.insert(END, seq)
        
    def update_label(self, header):
        name = header.split(' ')[0]
        name = name.split(',')[0]
        self.position_ids['label'].configure(text = name)
        
    def export(self):
        seq = self.get_self_selection()
        print seq, len(seq)
        
    def gcframe(self):
        seq = self.get_selection_or_sequence()
        if not seq: return
        np = NotePad()
        tid = np.text_id()
        tid.insert(END, self.translator.gcframe(seq, self.current_codon_table_id))

    def translate(self, frame = 1):
        seq = self.get_selection_or_sequence()
        if not seq: return
        np = NotePad()
        tid = np.text_id()
        tid.insert(END, self.translator.frame_nice(seq, frame, self.current_codon_table_id))

    def extract(self, frame = 1):
        seq = self.get_selection_or_sequence()
        if not seq: return
        aa_seq = self.translator.frame(seq, frame, self.current_codon_table_id)
        print '>%s<' % aa_seq
        aa_seq = re.sub('(.{50})','\\1\n',str(aa_seq))
        np = NotePad()
        tid = np.text_id()
        tid.insert(END,'>frame%d\n%s' % (frame,aa_seq))

    def statistics(self):
        seq = self.get_selection_or_sequence()
        if not seq: return
        seq = seq.upper()
        aa = {'A':0,'C':0,'G':0,'T':0,'N':0}
        for nt in seq:
            if nt not in aa: nt = 'N'
            aa[nt] = aa[nt] + 1

        GC = (100.0*(aa['G'] + aa['C']))/len(seq)
        
        np = NotePad()
        tid = np.text_id()

        tid.insert(END,"""%s

Length = %d
A=%d C=%d G=%d T=%d other=%d
GC=%f

""" % (time.strftime('%y %b %d, %X\n', time.localtime(time.time())),
               len(seq), aa['A'], aa['C'], aa['G'], aa['T'], aa['N'], GC)
                   )
        
    def blast(self):
        seq = self.get_selection_or_sequence()
        self.blaster = BlastIt(seq, self.parent)

    def reverse(self):
        w = self.sequence_id
        w.selection_own()
        try:
            start, stop = w.tag_ranges(SEL)
        except:
            start, stop = 1.0, self.sequence_id.index(END)

        seq = w.get(start, stop)
        seq = map(None,re.sub('[^A-Z]','',seq))
        seq.reverse()
        seq = ''.join(seq)

        w.delete(start, stop)
        w.insert(start, seq)
        w.tag_remove(SEL, 1.0, start)
        w.tag_add(SEL, start, stop)
        w.tag_remove(SEL, stop, END)

    def complement(self):
        w = self.sequence_id
        w.selection_own()
        try:
            start, stop = w.tag_ranges(SEL)
        except:
            start, stop = 1.0, self.sequence_id.index(END)

        seq = w.get(start, stop)
        seq = re.sub('[^A-Z]','',seq)    

        #print 'seq >%s<' % seq
        complementary = self.translator.complement(seq)
        w.delete(start, stop)
        w.insert(start, complementary)
        w.tag_remove(SEL, 1.0, start)
        w.tag_add(SEL, start, stop)
        w.tag_remove(SEL, stop, END)
             
    def antiparallel(self):
        w = self.sequence_id
        w.selection_own()
        try:
            start, stop = w.tag_ranges(SEL)
        except:
            start, stop = 1.0, self.sequence_id.index(END)

        seq = w.get(start, stop)
        seq = re.sub('[^A-Z]','',seq)    

        antip = self.translator.antiparallel(seq)
        w.delete(start, stop)
        w.insert(start, antip)
        w.tag_remove(SEL, 1.0, start)
        w.tag_add(SEL, start, stop)
        w.tag_remove(SEL, stop, END)

    def search(self):
        seq = self.get_selection_or_sequence()
        searcher = XDNAsearch(seq, master = self.sequence_id, highlight = 1)

    def goto(self, *args):
        pos = self.goto_entry.get()
        try:
            pos = int(pos) -1
        except:
            try:
                start, stop = pos.split(':')
                start = int(start)-1
                stop = int(stop)
                self.mark(start, stop)
                return
            except:
                import traceback
                traceback.print_exc()

                self.goto_entry.delete(0,END)
                return
            
        self.sequence_id.focus()
        self.sequence_id.mark_set('insert','1.%d' % pos)

    def mark(self, start, stop):
        self.sequence_id.focus()
        self.sequence_id.mark_set('insert','1.%d' % start)
        self.sequence_id.tag_add(SEL, '1.%d' % start, '1.%d' % stop)
        
if __name__ == '__main__':
    xbbtools = xbb_widget()
    xbbtools.main_frame.option_add('*frame.background', 'dimgrey')
    xbbtools.open('test.fas')
