#!/usr/bin/env python
# Copyright 2000 by Thomas Sicheritz-Ponten.
# Copyright 2016 by Markus Piotrowski.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# Created: Thu Jul 13 14:07:25 2000
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas

"""BLAST code for graphical Xbbtools tool."""


import glob
import os
import sys
import subprocess

import tkinter as tk
import tkinter.ttk as ttk
from tkinter import filedialog
from tkinter import messagebox


from xbb_utils import NotePad
import xbb_blastbg


class BlastIt:
    """Local BLAST integration for xbbtools."""

    nin, pin = [], []
    blast_ok = False
    blast_path = ""

    def __init__(self, seq, parent=None):
        """Set up new top-level window for BLAST search."""
        self.seq = seq
        self.parent = parent
        self.toplevel = tk.Toplevel(parent)
        self.toplevel.title("BLAST parameters")
        if not self.get_blast_databases() or not self.get_blast_binaries():
            return
        self.Choices()
        self.dbs.bind("<<ComboboxSelected>>", self.Validate)
        self.blasts.bind("<<ComboboxSelected>>", self.Validate)

    def get_blast_databases(self):
        """Try to locate the BLAST databases and put into lists."""
        if not (BlastIt.nin and BlastIt.pin):
            pin, nin = [], []

            try:
                pin.extend(glob.glob(os.environ["BLASTDB"] + "/*.pin"))
            except KeyError:
                pass
            pin.extend(glob.glob("C:*.pin"))

            try:
                nin.extend(glob.glob(os.environ["BLASTDB"] + "/*.nin"))
            except KeyError:
                pass

            # If no system variable BLASTDB exists, give user the chance to
            # locate his database folder:
            if not (nin and pin):
                database_folder = filedialog.askdirectory(
                    title="Please locate your BLAST database(s) folder:"
                )
                nin.extend(glob.glob(database_folder + "/*.nin"))
                pin.extend(glob.glob(database_folder + "/*.pin"))
                if not (nin and pin):
                    messagebox.showerror(
                        "xbb tools", "This folder does not contain any BLAST databases!"
                    )
                    self.toplevel.destroy()
                    return False

            self.pin = [os.path.splitext(x)[0] for x in pin]
            self.nin = [os.path.splitext(x)[0] for x in nin]

            BlastIt.pin = self.pin
            BlastIt.nin = self.nin

        return True

    def get_blast_binaries(self):
        """Test if BLAST binaries are in PATH or let user locate them."""
        if not BlastIt.blast_ok:
            # Test if blast binaries are in path
            if subprocess.call(
                ["blastn", "-version"]
            ):  # Return of non-zero means error
                self.blast_path = filedialog.askdirectory(
                    title="Please locate your BLAST program folder:"
                )
                if subprocess.call(
                    [os.path.join(self.blast_path, "blastn"), "-version"]
                ):
                    messagebox.showerror(
                        "xbb tools",
                        "Wrong folder or missing BLAST"
                        " binaries!\n  To run BLAST you must install the "
                        " standalone BLAST binaries.",
                    )
                    self.toplevel.destroy()
                    return False
                else:
                    BlastIt.blast_ok = True

            else:  # BLAST binaries are in PATH
                BlastIt.blast_ok = True
                self.blast_path = ""
        BlastIt.blast_path = self.blast_path
        self.toplevel.lift()
        return True

    def database_readable(self, db_paths):
        """Return the name of the blast database without path and extension."""
        db_names = [entry.split(os.sep)[-1].split(".")[0] for entry in db_paths]
        return db_names

    def convert_dbname_to_dbpath(self, db_name):
        """Return the full path for a given blast database name."""
        database_path = ""
        for database in self.nin:
            if database.endswith(db_name):
                database_path = database
                break
        for database in self.pin:
            if database.endswith(db_name):
                database_path = database
                break
        return database_path

    def Choices(self):
        """Set up window to select BLAST program and database."""
        self.blast_string = tk.StringVar()
        self.blast_string.set("blastn")
        self.cf = ttk.Frame(self.toplevel)
        self.cf.pack(side="top", expand=1, fill="x")
        self.dbs_frame = ttk.LabelFrame(self.cf, text="Databases")
        self.dbs_frame.pack(side="left", padx=5, pady=5, expand=1, fill="x")
        nin_values = self.database_readable(self.nin)
        pin_values = self.database_readable(self.pin)
        self.dbs = ttk.Combobox(
            self.dbs_frame, exportselection=0, values=nin_values + pin_values
        )
        self.dbs.current(0)

        self.blast_frame = ttk.LabelFrame(self.cf, text="BLAST programs")
        self.blast_frame.pack(side="left", padx=5, pady=5, expand=1, fill="x")
        self.blasts = ttk.Combobox(
            self.blast_frame,
            exportselection=0,
            textvariable=self.blast_string,
            values=["blastn", "blastp", "blastx", "tblastn", "tblastx"],
        )

        self.dbs.pack(side="left", padx=5, pady=5, expand=1, fill="x")
        self.blasts.pack(side="left", padx=5, pady=5, expand=1, fill="x")

        self.option_f = ttk.LabelFrame(self.cf, text="Command line options")
        self.option_f.pack(side="left", padx=5, pady=5, expand=1, fill="x")
        self.option = ttk.Entry(self.option_f)
        self.option.pack(side="left", padx=5, pady=5, fill="x", expand=1)
        self.ok = ttk.Button(self.cf, text="Run", command=self._Run, state="disabled")
        self.ok.pack(side="right")

        self.Validate()

    def Validate(self, *args):
        """Check everything and enable/disable 'Run' button."""
        db = self.convert_dbname_to_dbpath(self.dbs.get())
        prog = self.blasts.get()
        if (prog in ["blastn", "tblastx", "tblastn"]) == (db in self.nin):
            self.ok.config(state="normal")
        elif (prog in ["blastp", "blastx"]) == (db in self.pin):
            self.ok.config(state="normal")
        else:
            self.ok.config(state="disabled")

    def _Run(self):
        """Initialise options for Blast commandline (PRIVATE)."""
        command_options = self.option.get()
        options = ""
        if len(command_options.strip()):
            options = command_options.strip()

        db = self.convert_dbname_to_dbpath(self.dbs.get())
        prog = self.blast_path + self.blasts.get()
        self.command_data = [self.seq, prog, db, options]

        self.Run()

    def Run(self):
        """Open new notepad and initialize running BLAST."""
        self.notepad = NotePad()
        tid = self.notepad.tid

        self.toplevel.destroy()
        blastbg = xbb_blastbg.BlastDisplayer(self.command_data, tid)
        blastbg.RunCommand()


if __name__ == "__main__":
    try:
        seq = sys.argv[1]
    except IndexError:  # Started script without providing a sequence
        seq = "ATGACAAAGCTAATTATTCACTTGGTTTCAGACTCTTCTGTGCAAACTGC"
    win = tk.Tk()
    win.title("Dummy windows for BLAST test")
    test = BlastIt(seq)
    win.mainloop()
