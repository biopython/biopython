#!/usr/bin/env python
# Created: Thu Jul 13 14:07:25 2000
# Last changed: Time-stamp: <00/12/03 13:27:24 thomas>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
# File: xbb_blast.py

from __future__ import print_function

import glob
import os
import sys
from threading import *

try:
    from Tkinter import * # Python 2
except ImportError:
    from tkinter import * # Python 3

import Pmw
sys.path.insert(0, '.')

from xbb_utils import NotePad
import xbb_blastbg


class BlastIt:
    def __init__(self, seq, parent=None):
        self.seq = seq
        self.parent = parent
        self.toplevel = Toplevel(parent)
        Pmw.initialise(parent)
        self.GetBlasts()
        self.Choices()

    def GetBlasts(self):
        pin, nin = [], []
        try:
            pin.extend(glob.glob(os.environ['BLASTDB'] + '/*.pin'))
        except:
            pass
        pin.extend(glob.glob('*.pin'))

        try:
            nin.extend(glob.glob(os.environ['BLASTDB'] + '/*.nin'))
        except:
            pass
        nin.extend(glob.glob('*.nin'))

        self.pin = [os.path.splitext(x)[0] for x in pin]
        self.nin = [os.path.splitext(x)[0] for x in nin]

    def Choices(self):
        self.GetBlasts()
        self.cf = Frame(self.toplevel)
        self.cf.pack(side=TOP, expand=1, fill=X)
        self.dbs = Pmw.ComboBox(self.cf,
                                label_text='Blast Databases:',
                                labelpos='nw',
                                scrolledlist_items=self.nin + self.pin,
                                selectioncommand=self.Validate
                                )
        self.blasts = Pmw.ComboBox(self.cf,
                                   label_text='Blast Programs:',
                                   labelpos='nw',
                                   scrolledlist_items=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'],
                                   selectioncommand=self.Validate
                                   )
        self.dbs.pack(side=LEFT, expand=1, fill=X)
        self.blasts.pack(side=LEFT, expand=1, fill=X)

        self.alternative_f = Frame(self.cf)
        self.alternative = Entry(self.alternative_f)
        self.alternative_f.pack(side=TOP, fill=X, expand=1)
        self.alternative.pack(side=LEFT, fill=X, expand=1)
        self.ok = Button(self.alternative_f, text='Run',
                         command=self._Run)
        self.ok.pack(side=RIGHT)

        self.dbs.selectitem(0)
        self.blasts.selectitem(0)
        self.Validate()

    def Validate(self, *args):
        db = self.dbs.get()
        prog = self.blasts.get()
        color = 'red'
        if (prog in ['blastn', 'tblastx', 'tblastn']) == (db in self.nin):
            color = 'green'
        elif (prog in ['blastp', 'blastx']) == (db in self.pin):
            color = 'green'

        self.dbs.component('entry').configure(bg=color)
        self.blasts.component('entry').configure(bg=color)

    def _Run(self):
        alternative_command = self.alternative.get()
        if len(alternative_command.strip()):
            self.command = alternative_command.strip()
        else:
            db = self.dbs.get()
            prog = self.blasts.get()
            self.command = 'echo %s | nice blastall -d %s -p %s' % (self.seq, db, prog)

        self.Run()

    def Update(self):
        self.notepad.update()
        self.notepad.after(1, self.Update)

    def oldRun(self):
        self.notepad = NotePad()
        self.notepad.menubar.configure(bg='red')
        self.notepad.bind('<Destroy>', self.Exit)

        self.Update()

        print(self.command)
        self.pipe = posix.popen(self.command)
        while True:
            try:
                char = self.pipe.read(1)
                self.notepad.insert(END, char)
                self.notepad.update()
            except:
                break
            if not char:
                break

        try:
            self.notepad.menubar.configure(bg='green')
        except:
            pass

    def Run(self):
        self.notepad = NotePad()
        tid = self.notepad.tid
        self.notepad.menubar.configure(bg='red')

        self.toplevel.destroy()
        blastbg = xbb_blastbg.BlastDisplayer(self.command, tid)
        blastbg.RunCommand()

        # indicate the finished run by changing color
        try:
            self.notepad.menubar.configure(bg='green4')
        except:
            pass

    def Exit(self, *args):

        try:
            self.pipe.close()
            del(pipe)
        except:
            pass
        self.notepad.destroy()

        sys.exit(0)


if __name__ == '__main__':
    seq = sys.argv[1]
    win = Tk()
    test = BlastIt(seq)
    win.mainloop()
