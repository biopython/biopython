# Copyright 2011 by Wibowo Arindrarto (w.arindrarto@gmail.com)
# This code is part of the Biopython distribution and governed by its
# license. Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO parser for the ABIF format.

ABIF is the format used by Applied Biosystem's sequencing machines to store
sequencing results. 

For more details on the format specification, visit:
http://www.appliedbiosystem.com/support/software_community/ABIF_File_Format.pdf

"""

__docformat__ = "epytext en"

import datetime
import struct

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# dictionary for deciding which tags goes into SeqRecord annotation
# each key is tag_name + tag_number
# if a tag entry needs to be added, just add its key and its key
# for the annotations dictionary here
_EXTRACT = {
            'TUBE1': 'sample well',
            'MODL1': 'machine model',
            'DySN1': 'dye',
            'GTyp1': 'polymer',
            'RUND1': 'run start date',
            'RUND2': 'run finish date',
            'RUND3': 'data collection start date',
            'RUND4': 'data collection finish date',
            'RUNT1': 'run start time',
            'RUNT2': 'run finish time',
            'RUNT3': 'data collection start time',
            'RUNT4': 'data collection finish time',
           }

# dictionary for data unpacking format
_BYTEFMT = {
            1: 'b',     # byte
            2: 's',     # char
            3: 'H',     # word
            4: 'h',     # short
            5: 'i',     # long
            6: '2i',    # rational, legacy unsupported
            7: 'f',     # float
            8: 'd',     # double
            10: 'h2B',  # date
            11: '4B',   # time
            12: '2i2b', # thumb
            13: 'B',    # bool
            14: '2h',   # point, legacy unsupported
            15: '4h',   # rect, legacy unsupported
            16: '2i',   # vPoint, legacy unsupported
            17: '4i',   # vRect, legacy unsupported
            18: 's',    # pString
            19: 's',    # cString
            20: '2i',   # tag, legacy unsupported
           }

class _Dir(object):
    """Class representing directory content. (PRIVATE)"""
    def __init__(self, tag_entry, handle):
        """Instantiates the _Dir object.

        tag_entry - tag name, tag number, element type code, element size,
                    number of elements, data size, data offset,
                    directory handle, and tag start position
        handle - the abif file object from which the tags would be unpacked
        """
        self.elem_code = tag_entry[2]
        self.elem_num = tag_entry[4]
        self.data_size = tag_entry[5]
        self.data_offset = tag_entry[6]
        self.tag_offset = tag_entry[8]

        # if data size <= 4 bytes, data is stored inside tag
        # so offset needs to be changed
        if self.data_size <= 4:
            self.data_offset = self.tag_offset + 20

        self.tag_data = self._unpack_tag(handle)

    def _unpack_tag(self, handle):
        """"Returns tag data. (PRIVATE)
        
        handle - the abif file object from which the tags would be unpacked
        """ 
        if self.elem_code in _BYTEFMT:
            
            # because '>1s' unpack differently from '>s'
            num = '' if self.elem_num == 1 else str(self.elem_num)
            fmt = '>' + num + _BYTEFMT[self.elem_code]
            fmt_size = struct.calcsize(fmt)
            start = self.data_offset

            handle.seek(start)
            data = struct.unpack(fmt, handle.read(fmt_size))

            # no need to use tuple if len(data) == 1
            if self.elem_code not in [10, 11] and len(data) == 1:
                data = data[0]

            # account for different data types
            if self.elem_code == 10:
                return str(datetime.date(*data))
            elif self.elem_code == 11:
                return str(datetime.time(*data))
            elif self.elem_code == 13:
                return bool(data)
            elif self.elem_code == 18:
                return data[1:]
            elif self.elem_code == 19:
                return data[:-1]
            else:
                return data
        else:
            return None

def AbifIterator(handle, alphabet=IUPAC.unambiguous_dna, trim=False):
    """Iterator for the Abif file format.
    """
    file_id = handle.name.replace('.ab1','')
    annot = {}

    if not handle.read(4) == 'ABIF':
        raise IOError('%s is not a valid ABIF file.' % file_id)

    handle.seek(0)
    header = struct.unpack('>4sH4sI2H3I', handle.read(30))

    for entry in _abif_parse_header(header, handle):
        key = entry.tag_name + str(entry.tag_num)
        # PBAS2 is base-called sequence
        if key == 'PBAS2': 
            seq = entry.tag_data
            ambigs = 'KYWMRS'
            if len(set(seq).intersection(ambigs)) > 0:
                alphabet = IUPAC.ambiguous_dna
        # PCON2 is quality values of base-called sequence
        elif key == 'PCON2':
            qual = [ord(val) for val in entry.tag_data]
        # SMPL1 is sample id entered before sequencing run
        elif key == 'SMPL1':
            sample_id = entry.tag_data
        else:
            if key in _EXTRACT:
                annot[_EXTRACT[key]] = entry.tag_data
            else:
                raise KeyError('The %s tag can not be found.' % key)
    
    record = SeqRecord(Seq(seq, alphabet),
                       id=file_id, name=sample_id,
                       description='',
                       annotations=annot,
                       letter_annotations={'phred_quality': qual}
                      )

    if not trim:
        yield record
    else:
        yield _abif_trim(record)

def _AbifTrimIterator(handle):
    """Iterator for the Abif file format that yields trimmed SeqRecord objects.
    """
    return AbifIterator(handle, alphabet=IUPAC.unambiguous_dna, trim=True)

def _abif_parse_header(header, handle):
    """Generator that returns directory contents.
    """
    # header structure:
    # file type, file version, tag name, tag number,
    # element type code, element size, number of elements
    # data size, data offset, handle (not file handle)
    head_elem_size = header[5]
    head_elem_num = header[6]
    head_offset = header[8]
    index = 0

    while index < head_elem_num:
        start = head_offset + index * head_elem_size
        # add directory offset to tuple
        # to handle directories with data size <= 4 bytes
        handle.seek(start)
        dir_content = struct.unpack('>4sI2H4I', handle.read(28)) + (start,)
        index += 1
        yield _Dir(dir_content, handle)


def _abif_trim(seq_record):
    """Trims the sequence using Richard Mott's modified trimming algorithm.

    seq_record - SeqRecord object to be trimmed.

    Trimmed bases are determined from their segment score, which is a
    cumulative sum of each base's score. Base scores are calculated from
    their quality values.

    More about the trimming algorithm:
    http://www.phrap.org/phredphrap/phred.html
    http://www.clcbio.com/manual/genomics/Quality_abif_trimming.html
    """

    start = False   # flag for starting position of trimmed seq_recorduence
    segment = 20    # minimum seq_recorduence length
    trim_start = 0  # init start index
    cutoff = 0.05   # default cutoff value for calculating base score

    if len(seq_record) <= segment:
        raise ValueError('Sequence not trimmed because length is shorter than \
                         minimum required length (20).')
    else:
        # calculate base score
        score_list = [cutoff - (10 ** (qual/-10.0)) for qual in
                      seq_record.letter_annotations['phred_quality']]

        # calculate cummulative score
        # if cummulative value < 0, set it to 0
        # first value is set to 0, because of the assumption that
        # the first base will always be trimmed out
        cummul_score = [0]
        for i in xrange(1, len(score_list)):
            score = cummul_score[-1] + score_list[i]
            if score < 0:
                cummul_score.append(0)
            else:
                cummul_score.append(score)
                if not start:
                    # trim_start = value when cummulative score is first > 0
                    trim_start = i
                    start = True
        
        # trim_finish = index of highest cummulative score,
        # marking the end of sequence segment with highest cummulative score
        trim_finish = cummul_score.index(max(cummul_score))
                         
        return seq_record[trim_start:trim_finish]

if __name__ == '__main__':
    print "Self-testing AbifIO..."
