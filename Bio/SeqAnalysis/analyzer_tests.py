#Work in progress, not ready for reviewing


from Bio.SeqAnalysis.SeqAnalyzer.py import
import unittest
from Bio import SeqIO
#cant create objects
records = {
    'czlowiek': 'MGKKRAPQPIVEKLISNFSQK',
    'pies': 'MSVKPSKKKSKRSKVKKKISFDFSDDDDSEIGVSFR',
}
analyzer = Bio.SeqAnalysis.SeqAnalyzer(records)
reults = analyzer.results()

SeqAnalyzer = SeqAnalyzer(input_data)

print(SeqAnalyzer.results()[0], '\n', SeqAnalyzer.results()[1])
