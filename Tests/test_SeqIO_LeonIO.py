
import unittest
import os

from Bio import SeqIO

class LeonParseWrite(unittest.TestCase):
    def test(self):
        for file in os.listdir('./Leon'):
            if '.leon' not in file:
                continue

            sequences = list(SeqIO.parse(f'./Leon/{file}', 'leon'))
            SeqIO.write(sequences, f'./Leon/temp_{file}', 'leon')
            os.remove(f'./Leon/{file.replace(".leon", ".d")}')
            os.remove(f'./Leon/temp_{file}')