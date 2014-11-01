# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Parser for FSSP files, used in a database of protein fold classifications.

This is a module to handle FSSP files. For now it parses only the header,
summary and alignment sections.

See: Holm and Sander (1996) The FSSP database: fold classification based on
structure-structure alignment of proteins.

functions: read_fssp(file_handle): reads an fssp file into the records. Returns a
tuple of two instances.
mult_align: returns a Biopython alignment object
"""
from __future__ import print_function

import re
from . import fssp_rec
from Bio.Align import Generic
from Bio import Alphabet
fff_rec = fssp_rec.fff_rec
header_records = {
    'database': re.compile('^DATABASE'),
    'pdbid': re.compile('^PDBID'),
    'header': re.compile('^HEADER'),
    'compnd': re.compile('^COMPND'),
    'author': re.compile('^AUTHOR'),
    'source': re.compile('^SOURCE'),
    'seqlength': re.compile('^SEQLENGTH'),
    'nalign': re.compile('^NALIGN')
}

summary_title = re.compile('## +SUMMARY')
summary_rec = re.compile(' *[0-9]+: +[1-9][0-9a-z]{3,3}')
alignments_title= re.compile('## +ALIGNMENTS')
alignments_rec = re.compile(' *[0-9]+ +-{0,1}[0-9]+')
equiv_title = re.compile('## +EQUIVALENCES')


class FSSPHeader(object):
    def __init__(self):
        self.database = None
        self.pdbid = ''
        self.header = ''
        self.compnd = ''
        self.source = ''
        self.author = []
        self.seqlength = 0
        self.nalign = 0

    def fill_header(self, inline):
        for i in header_records:
            if header_records[i].match(inline):
                if i == 'database' or i == 'seqlength' or i == 'nalign':
                    setattr(self, i, int(inline.split()[1]))
                elif i == 'compnd' or i == 'author':
                    setattr(self, i, inline.split()[1:])
                elif i == 'source' or i == 'header':
                    attr = inline[inline.find(' ')+1:].strip()
                    setattr(self, i, attr)
                else:
                    setattr(self, i, inline.split()[1])


class PosAlign(object):
    def __init__(self, inStr):
        inStr = inStr.strip()
        if len(inStr) != 1 and len(inStr) != 2:
            raise ValueError('PosAlign: length not 2 chars' + inStr)
        if inStr == '..':
            self.aa = '-'
            self.gap = 1
        else:
            self.gap = 0
            self.aa = inStr[0]
            if self.aa == self.aa.lower():
                self.aa = 'C'
            if len(inStr) == 2:
                self.ss = inStr[1].upper()
            else:
                self.ss = '0'

    def __repr__(self):
        if self.gap:
            outstring = '..'
        else:
            outstring = self.aa+self.ss.lower()
        return outstring

    __str__ = __repr__


class FSSPSumRec(object):
    """ Contains info from an FSSP summary record"""
    def __init__(self, in_str):
        self.raw = in_str
        in_rec = in_str.strip().split()
        # print(in_rec)
        self.nr = int(in_rec[0][:-1])
        self.pdb1 = in_rec[1][:4]
        if len(in_rec[1]) == 4:
            self.chain1 = '0'
        elif len(in_rec[1]) == 5:
            self.chain1=in_rec[1][4]
        else:
            raise ValueError('Bad PDB ID 1')
        self.pdb2 = in_rec[2][:4]
        if len(in_rec[2]) == 4:
            self.chain2='0'
        elif len(in_rec[2]) == 5:
            self.chain2=in_rec[2][4]
        else:
            raise ValueError('Bad PDB ID 2')
        self.zscore = float(in_rec[3])
        self.rmsd = float(in_rec[4])
        self.lali = float(in_rec[5])
        self.lseq2 = float(in_rec[6])
        self.pID = float(in_rec[7])
        self.revers = int(in_rec[8])
        self.permut = int(in_rec[9])
        self.nfrag = int(in_rec[10])
        self.topo = in_rec[11]
        self.doc = ''
        for i in in_rec[12:]:
            self.doc = self.doc + i + ' '
        self.doc = self.doc.rstrip() + '\n'

    def __repr__(self):
        return self.raw
    __str__ = __repr__


class FSSPAlignRec(object):
    def __init__(self, in_fff_rec):
        # print(in_fff_rec)
        self.abs_res_num = int(in_fff_rec[fssp_rec.align.abs_res_num])
        self.pdb_res_num = in_fff_rec[fssp_rec.align.pdb_res_num].strip()
        self.chain_id = in_fff_rec[fssp_rec.align.chain_id]
        if self.chain_id == ' ':
            self.chain_id = '0'
        self.res_name = in_fff_rec[fssp_rec.align.res_name]
        if self.res_name == self.res_name.lower():
            self.res_name = 'C'
        self.ss1 = in_fff_rec[fssp_rec.align.ss1]
        self.turn3 = in_fff_rec[fssp_rec.align.turn3]
        self.turn4 = in_fff_rec[fssp_rec.align.turn4]
        self.turn5 = in_fff_rec[fssp_rec.align.turn5]
        self.pos_align_dict = {}
        self.PosAlignList = []

    def add_align_list(self, align_list):
        for i in align_list:
            self.PosAlignList.append(PosAlign(i))

    def pos_align_list2dict(self):
        j = 1
        for i in self.PosAlignList:
            self.pos_align_dict[j] = i
            j = j + 1


class FSSPAlignDict(dict):
    def __init__(self):
        # The following two dictionaries are pointers to records in self
        # The first dictionary is a "pdb_residue_number: self_key"
        # The second dictionary is a "absolute_residue_number: self_key"
        self.pdb_res_dict = {}
        self.abs_res_dict = {}
        self.data = {}

    def build_resnum_list(self):
        for i in self:
            self.abs_res_dict[self[i].abs_res_num] = i
            self.pdb_res_dict[self[i].pdb_res_num] = i

    # Given an absolute residue number & chain, returns the relevant fssp
    # record
    def abs(self, num):
        return self[self.abs_res_dict[num]]

    # Given an PDB residue number & chain, returns the relevant fssp
    # record
    def pdb(self, num):
        return self[self.pdb_res_dict[num]]

    # Returns a sequence string
    def sequence(self, num):
        s = ''
        for i in sorted(self.abs_res_dict):
            s += self.abs(i).pos_align_dict[num].aa
        return s

    def fasta_mult_align(self):
        mult_align_dict = {}
        for j in self.abs(1).pos_align_dict:
            mult_align_dict[j] = ''
        for fssp_rec in self.values():
            for j in fssp_rec.pos_align_dict:
                mult_align_dict[j] += fssp_rec.pos_align_dict[j].aa
        out_str = ''
        for i in sorted(mult_align_dict):
            out_str += '> %d\n' % i
            k = 0
            for j in mult_align_dict[i]:
                k += 1
                if k % 72 == 0:
                    out_str += '\n'
                out_str += j
            out_str += '\n'
        return out_str


class FSSPSumDict(dict):
    pass


#
# Process a fssp file into its constituents. Return a 2-tuple containing
# a list of FSSPSumRecs and a dictionary of alignment records.
#
def read_fssp(fssp_handle):
    header = FSSPHeader()
    sum_dict = FSSPSumDict()
    align_dict = FSSPAlignDict()
    curline = fssp_handle.readline()
    while not summary_title.match(curline):
        # Still in title
        header.fill_header(curline)
        curline = fssp_handle.readline()

    if not summary_title.match(curline):
        raise ValueError('Bad FSSP file: no summary record found')
    curline = fssp_handle.readline()  # Read the title line, discard
    curline = fssp_handle.readline()  # Read the next line
    # Process the summary records into a list
    while summary_rec.match(curline):
        cur_sum_rec = FSSPSumRec(curline)
        sum_dict[cur_sum_rec.nr] = cur_sum_rec
        curline = fssp_handle.readline()

    # Outer loop: process everything up to the EQUIVALENCES title record
    while not equiv_title.match(curline):
        while (not alignments_title.match(curline) and
               not equiv_title.match(curline)):
            curline = fssp_handle.readline()
        if not alignments_title.match(curline):
            if equiv_title.match(curline):
                # print("Reached equiv_title")
                break
            else:
                raise ValueError('Bad FSSP file: no alignments title record found')

        if equiv_title.match(curline):
            break
        # If we got to this point, this means that we have matched an
        # alignments title. Parse the alignment records in a loop.
        curline = fssp_handle.readline()  # Read the title line, discard
        curline = fssp_handle.readline()  # Read the next line
        while alignments_rec.match(curline):
            align_rec = FSSPAlignRec(fff_rec(curline))
            key = align_rec.chain_id + align_rec.res_name + str(align_rec.pdb_res_num)
            align_list = curline[fssp_rec.align.start_aa_list:].strip().split()
            if key not in align_dict:
                align_dict[key] = align_rec
            align_dict[key].add_align_list(align_list)
            curline = fssp_handle.readline()
            if not curline:
                print('EOFEOFEOF')
                raise EOFError
    for i in align_dict.values():
        i.pos_align_list2dict()
        del i.PosAlignList
    align_dict.build_resnum_list()
    return (header, sum_dict, align_dict)
