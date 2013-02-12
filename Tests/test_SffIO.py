import re
import unittest
from Bio import SeqIO

# sffinfo E3MFGYR02_random_10_reads.sff | sed -n '/>\|Run Prefix\|Region\|XY/p'
test_data = """
>E3MFGYR02JWQ7T
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3946_2103
>E3MFGYR02JA6IL
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3700_3115
>E3MFGYR02JHD4H
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  3771_2095
>E3MFGYR02GFKUC
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2520_2738
>E3MFGYR02FTGED
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2268_2739
>E3MFGYR02FR9G7
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2255_0361
>E3MFGYR02GAZMS
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2468_1618
>E3MFGYR02HHZ8O
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2958_1574
>E3MFGYR02GPGB1
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2633_0607
>E3MFGYR02F7Z7G
  Run Prefix:   R_2008_01_09_16_16_00_
  Region #:     2
  XY Location:  2434_1658"""


class TestUAN(unittest.TestCase):
    def setUp(self):
        self.records = [record for record in SeqIO.parse('Roche/E3MFGYR02_random_10_reads.sff', 'sff')]
        self.test_annotations = {}
        for line in test_data.splitlines():
            fields = re.split(r"\s+", line.strip())
            if '>' in line:
                current_name = fields[0].lstrip('>')
                self.test_annotations[current_name] = {}
            elif 'Prefix' in line:
                time_list = map(int, fields[2].split('_')[1:-1])
                self.test_annotations[current_name]["time"] = time_list
            elif 'Region' in line:
                region = int(fields[-1])
                self.test_annotations[current_name]["region"] = region
            elif 'XY' in line:
                x, y = map(int, fields[-1].split('_'))
                self.test_annotations[current_name]["coords"] = (x, y)

    def test_time(self):
        for record in self.records:
            self.assertEqual(record.annotations["time"], self.test_annotations[record.name]["time"])

    def test_region(self):
        for record in self.records:
            self.assertEqual(record.annotations["region"], self.test_annotations[record.name]["region"])

    def test_coords(self):
        for record in self.records:
            self.assertEqual(record.annotations["coords"], self.test_annotations[record.name]["coords"])

if __name__ == '__main__':
    unittest.main()
