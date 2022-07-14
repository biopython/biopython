from Bio.SeqIO.Interfaces import SequenceWriter
from Bio.SeqIO.FastaIO import FastaIterator, FastaWriter
from Bio.SeqIO.QualityIO import FastqPhredIterator, FastqPhredWriter

import os

def importPyleon():
    '''
    ensure pyleon is imported correctly
    '''

    try:
        from pyleon import PyLeon
    except ImportError:
        print('Error: pyleon could not be imported!')
        print('1)   git clone https://github.com/jdvne/pyleon')
        print('2)   follow installation instructions in README.md')

    return PyLeon


def LeonIterator(handle):
    if type(handle) != str:
        print('Error: source must be of type str')
        raise TypeError

    pyleon = importPyleon()(silent=True)
    filepath = pyleon.decompress(handle)

    # return iterator based upon decompressed form
    if 'fasta' in os.path.basename(filepath):
        return FastaIterator(filepath)
    if 'fastq' in os.path.basename(filepath):
        return FastqPhredIterator(filepath)

    print(f'Error: LEON failed to decompress {handle}')
    raise FileNotFoundError


class LeonWriter(SequenceWriter):
    def __init__(self, handle, kmer_size=None, abundance=None, nb_cores=None, 
        lossless=True, seq_only=False, noheader=False, noqual=False):
        if type(handle) != str:
            print('Error: source must be of type str')
            raise TypeError

        if type(handle) != str:
            print('Error: source must be of type str')
            raise TypeError

        self.handle = handle
        self.kwargs = {
            'kmer_size': kmer_size,
            'abundance': abundance,
            'nb_cores': nb_cores,
            'lossless': lossless,
            'seq_only': seq_only,
            'noheader': noheader,
            'noqual': noqual
        }


    def write_file(self, sequences):
        pyleon = importPyleon()(silent=True)

        # attempt to write as fastq
        try:
            filepath = self.handle.replace('.leon', '')
            if '.fastq' not in filepath: filepath += '.fastq'
            
            fastq_phred_writer = FastqPhredWriter(filepath)
            fastq_phred_writer._header_written = True
            fastq_phred_writer._footer_written = False

            count = fastq_phred_writer.write_file(sequences)

        # if unable to write as fastq, write as fasta
        except ValueError:
            # handle failed fastq file
            os.remove(filepath)

            filepath = self.handle.replace('.leon', '')
            if '.fasta' not in filepath: filepath += '.fasta'

            count = FastaWriter(filepath).write_file(sequences)
        
        pyleon.compress(filepath, **self.kwargs)
        os.remove(filepath)
        return count