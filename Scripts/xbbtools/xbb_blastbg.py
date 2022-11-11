#!/usr/bin/env python
# Copyright 2000 by Thomas Sicheritz-Ponten.
# Copyright 2016 by Markus Piotrowski.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Created: Sat Dec  2 16:02:17 2000
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas

"""BLAST code for graphical Xbbtools tool."""


import os
import tempfile
import threading

from tkinter import Tk, messagebox

from Bio.Blast.Applications import (
    NcbiblastnCommandline,
    NcbiblastpCommandline,
    NcbiblastxCommandline,
    NcbitblastnCommandline,
    NcbitblastxCommandline,
)


class BlastDisplayer:
    """A class for running and displaying a BLAST search."""

    def __init__(self, command_data, text_id=None):
        """Take command data and notpad id."""
        self.command_data = command_data
        self.tid = text_id

    def RunCommand(self):
        """Run the BLAST search."""
        self.fh_in, self.infile = tempfile.mkstemp()
        self.fh_out, self.outfile = tempfile.mkstemp()

        with open(self.infile, "w+") as f:
            f.write(">Name\n")
            f.write(self.command_data[0])

        blast_program = self.command_data[1]
        database = self.command_data[2]

        # Check if user supplied additional options and extract them
        if self.command_data[3]:
            option = self.command_data[3]
            options = {}
            for x in range(0, len(option.split()) - 1, 2):
                options[option.split()[x]] = option.split()[x + 1]
        else:
            options = {}

        args, kwargs = (
            blast_program,
            {"query": self.infile, "db": database, "out": self.outfile},
        )

        if blast_program.endswith("blastn"):
            blast_cmd = NcbiblastnCommandline(args, **kwargs)
        elif blast_program.endswith("blastp"):
            blast_cmd = NcbiblastpCommandline(args, **kwargs)
        elif blast_program.endswith("blastx"):
            blast_cmd = NcbiblastxCommandline(args, **kwargs)
        elif blast_program.endswith("tblastn"):
            blast_cmd = NcbitblastnCommandline(args, **kwargs)
        elif blast_program.endswith("tblastx"):
            blast_cmd = NcbitblastxCommandline(args, **kwargs)
        else:
            return

        if options:
            try:
                for key in options:
                    blast_cmd.set_parameter(key, options[key])
            except ValueError as e:
                messagebox.showerror("xbb tools", "Commandline error:\n\n" + str(e))
                self.tid.destroy()
                return

        self.worker = BlastWorker(blast_cmd)
        self.worker.start()

        self.UpdateResults()

    def UpdateResults(self):
        """Write BLAST result data into notepad."""
        # open the oufile and displays new appended text
        self.tid.insert("end", "BLAST is running...")
        while True:
            self.tid.update()
            if self.worker.finished:
                self.tid.delete("1.0", "end")
                break
        with open(self.outfile) as fid:
            try:
                txt = fid.read()
                self.tid.insert("end", txt)
                self.tid.update()
            except Exception:
                # text widget is destroyed, we assume the search
                # has been cancelled
                pass
        self.Exit()

    def Exit(self):
        """Clean up on exit."""
        if os.path.exists(self.outfile):
            os.close(self.fh_out)
            os.remove(self.outfile)
        if os.path.exists(self.infile):
            os.close(self.fh_in)
            os.remove(self.infile)

        del self.worker


class BlastWorker(threading.Thread):
    """Allows multiple BLAST searches by threading."""

    def __init__(self, blast_command):
        """Initialize the worker."""
        self.com = blast_command
        threading.Thread.__init__(self)
        self.finished = 0

    def run(self):
        """Start worker."""
        try:
            self.com()
        except Exception as e:
            messagebox.showwarning("BLAST error", "BLAST error:\n\n" + str(e))
        self.finished = 1


if __name__ == "__main__":
    from xbb_blast import BlastIt

    win = Tk()
    win.title("Dummy windows for BLAST test")
    test = BlastIt("ATGACAAAGCTAATTATTCACTTGGTTTCAGACTCTTCTGTGCAAACTGC")
    win.mainloop()
