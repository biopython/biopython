# Copyright 2014 Marco Galardini.  All rights reserved.
# Adapted from test_Mymodule.py by Jeff Chang
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

try:
    import numpy
except ImportError:
    from Bio import MissingExternalDependencyError
    raise MissingExternalDependencyError(
        "Install NumPy if you want to use Bio.Phenotype.")

import os
import json
import unittest

from Bio._py3k import StringIO

from Bio import Phenotype

# Example plate files
JSON_PLATE = 'Phenotype/Plate.json'
JSON_PLATE_2 = 'Phenotype/Plate_2.json'
JSON_PLATE_3 = 'Phenotype/Plate_3.json'
JSON_BAD = 'Phenotype/BadPlate.json'

CSV_PLATES = 'Phenotype/Plates.csv'

class TestPhenoMicro(unittest.TestCase):
    
    def test_PhenMicroIO(self):
        '''Test basic functionalities of Phenotype IO methods'''
        self.assertRaises(ValueError, Phenotype.read, CSV_PLATES, 'pm-csv')
        self.assertRaises(ValueError, Phenotype.read, CSV_PLATES, 'pm-json')
        self.assertRaises(ValueError, Phenotype.read, CSV_PLATES, 'pm-noformat')
        self.assertRaises(ValueError, Phenotype.read, CSV_PLATES, 'PM-CSV')
        self.assertRaises(TypeError, Phenotype.read, CSV_PLATES, 1)
        self.assertRaises(KeyError, Phenotype.read, JSON_BAD, 'pm-json')
        
        p1 = Phenotype.read(JSON_PLATE_3, 'pm-json')
        p2 = next(Phenotype.parse(CSV_PLATES, 'pm-csv'))
        
        handle = StringIO()
        
        c = Phenotype.write([p1, p2], handle, 'pm-json')
        self.assertEqual(c, 2)
        
        handle.flush()
        handle.seek(0)
        #Now ready to read back from the handle...
        try:
            records = list(Phenotype.parse(handle, 'pm-json'))
        except ValueError as e:
            #This is BAD.  We can't read our own output.
            #I want to see the output when called from the test harness,
            #run_tests.py (which can be funny about new lines on Windows)
            handle.seek(0)
            raise ValueError("%s\n\n%s\n\n%s"
                              % (str(e), repr(handle.read()), repr(records)))
         
        self.assertEqual(p1, records[0])

        handle.close()
        handle = StringIO()
        self.assertRaises(TypeError, Phenotype.write, p1, handle, 1)
        self.assertRaises(ValueError, Phenotype.write, p1, handle, 'PM-JSON')
        self.assertRaises(ValueError, Phenotype.write, p1, handle, 'pm-csv')
        handle.close()
    
    def test_PlateRecord(self):
        '''Test basic functionalities of PlateRecord objects'''
        self.assertRaises(ValueError,
                Phenotype.PhenMicro.PlateRecord, 'test', [1,2,3])
        self.assertRaises(TypeError, 
                Phenotype.PhenMicro.PlateRecord, 'test', 1)
        
        handle = open(JSON_PLATE)
        j = json.load(handle)
        handle.close()
        
        p = Phenotype.PhenMicro.PlateRecord(j['csv_data']['Plate Type'])

        times = j['measurements']['Hour']
        for k in j['measurements']:
            if k == 'Hour':continue
            p[k] = Phenotype.PhenMicro.WellRecord(k, signals=
                        {times[i]:j['measurements'][k][i]
                        for i in range(len(times))})
        
        del j['measurements']
        p.qualifiers = j
            
        self.assertEqual(p.id, 'PM01')
        self.assertEqual(len(p), 96)
        self.assertEqual(p.qualifiers, j)
        self.assertRaises(ValueError, p._is_well, 'a')
        self.assertEqual(p['A01'].id, 'A01')
        self.assertRaises(KeyError, p.__getitem__, 'test')
        self.assertEqual(len(p[1]), 12)
        self.assertEqual(len(p[1:4:2]), 24)
        self.assertEqual(p[1,2], p['B03'])
        self.assertEqual(len(p[:,1]), 8)
        self.assertEqual(len(p[:,1:4:2]), 16)
        self.assertRaises(TypeError, p.__getitem__, 1, 2, 3)
        self.assertRaises(IndexError, p.__getitem__, 13)        
        self.assertRaises(ValueError, p.__setitem__, 'A02', p['A01'])
        self.assertRaises(ValueError, p.__setitem__, 'A02', 'a')
        p['A02'] = p['A02']
        for w in p:
            pass
        self.assertEqual('A01' in p, True)
        self.assertEqual('test' in p, False)
        self.assertRaises(ValueError, next, p.get_row('test'))
        self.assertEqual(next(p.get_row('A')), p['A01'])
        self.assertRaises(ValueError, next, p.get_column('test'))
        self.assertEqual(next(p.get_column('12')), p['A12'])
        self.assertEqual(next(p.get_column('1')), p['A01'])
        self.assertRaises(ValueError, p.subtract_control, 'A121')
        self.assertRaises(ValueError, p.subtract_control, wells=['A121'])
        p2 = p.subtract_control()
        self.assertEqual(p2.id, p.id)
        self.assertEqual(p2['A02'], p['A02'] - p['A01'])
        self.assertEqual(repr(p), "PlateRecord('WellRecord['A01'], WellRecord"+
        "['A02'], WellRecord['A03'], WellRecord['A04']...WellRecord['H12']')")
        self.assertEqual(str(p), "Plate ID: PM01\nWell: 96\nRows: 8\nColumns: "+
        "12\nPlateRecord('WellRecord['A01'], WellRecord['A02'], WellRecord"+
        "['A03'], WellRecord['A04']...WellRecord['H12']')")
        
        handle = open(JSON_PLATE_2)
        j = json.load(handle)
        handle.close()
        
        p1 = Phenotype.PhenMicro.PlateRecord(j['csv_data']['Plate Type'])

        times = j['measurements']['Hour']
        for k in j['measurements']:
            if k == 'Hour':continue
            p1[k] = Phenotype.PhenMicro.WellRecord(k, signals=
                        {times[i]:j['measurements'][k][i]
                        for i in range(len(times))})              
        
        del j['measurements']
        p1.qualifiers = j
        
        self.assertRaises(TypeError, p.__add__, 'a')
        self.assertRaises(TypeError, p.__sub__, 'a')
        
        p3 = p + p1
        self.assertEqual(p3['A02'], p['A02'] + p1['A02'])
        
        p3 = p - p1
        self.assertEqual(p3['A02'], p['A02'] - p1['A02'])
        
        del p['A02']
        self.assertRaises(ValueError, p.__add__, p1)
        self.assertRaises(ValueError, p.__sub__, p1)
        
    def test_WellRecord(self):
        '''Test basic functionalities of WellRecord objects'''
        handle = open(JSON_PLATE)
        p = json.load(handle)
        handle.close()
        
        times = p['measurements']['Hour']
        w = Phenotype.PhenMicro.WellRecord('A10', signals=
                        {times[i]:p['measurements']['A10'][i]
                            for i in range(len(times))})
        
        w1 = Phenotype.PhenMicro.WellRecord('H12', signals=
                        {times[i]:p['measurements']['H12'][i]
                            for i in range(len(times))})
                            
        self.assertIsInstance(w.plate,
                    Phenotype.PhenMicro.PlateRecord)    
        self.assertEqual(w.id, 'A10')
        self.assertEqual(len(w), len(times))
        self.assertEqual(len(w), 384)
        self.assertEqual(max(w), (95.75, 217.0))
        self.assertEqual(min(w), (0.0, 37.0))
        self.assertEqual(max(w, key=lambda x: x[1]),
                (16.75, 313.0))
        self.assertEqual(min(w, key=lambda x: x[1]),
                (0.25, 29.0))
        self.assertEqual(len(w[:]), 96)
        self.assertEqual(w[1], 29.)
        self.assertEqual(w[12], 272.)
        self.assertEqual(w[1:5], [29.,35.,39.,43.])
        self.assertRaises(ValueError, w.__getitem__, 'a')
        self.assertAlmostEqual(w[1:2:.25][0], 29.)
        self.assertAlmostEqual(w[1.3567], 33.7196)
        self.assertEqual(w.get_raw()[0], (0.0, 37.0))
        self.assertEqual(w.get_raw()[-1], (95.75, 217.0))
        self.assertEqual(w.get_times()[0], 0.0)
        self.assertEqual(w.get_times()[-1], 95.75)
        self.assertEqual(w.get_signals()[0], 37.0)
        self.assertEqual(w.get_signals()[-1], 217.0)
        self.assertEqual(repr(w),
                        "WellRecord('(0.0, 37.0), (0.25, 29.0), (0.5, 32.0),"+
                        " (0.75, 30.0), (1.0, 29.0)...(95.75, 217.0)')")
        self.assertEqual(str(w),
                        "Well ID: A10\nTime points: 384\nMinum signal 0.25 at "+
                        "time 29.00\nMaximum signal 16.75 at time "+
                        "313.00\nWellRecord('(0.0, 37.0), (0.25, 29.0), "+
                        "(0.5, 32.0), (0.75, 30.0), "+
                        "(1.0, 29.0)...(95.75, 217.0)')")
                        
        w.fit()
        try:
            import scipy
            self.assertAlmostEqual(w.area, 20879.5)
            self.assertAlmostEqual(w.model, 'gompertz')
            self.assertAlmostEqual(w.lag, 6.0425868725090357)
            self.assertAlmostEqual(w.plateau, 188.51404344898586)
            self.assertAlmostEqual(w.slope, 48.190618284831132)
            self.assertAlmostEqual(w.v, 0.10000000000000001)
            self.assertAlmostEqual(w.y0, 45.879770069807989)
        except ImportError:
            self.assertAlmostEqual(w.area, None)
            self.assertAlmostEqual(w.model, None)
            self.assertAlmostEqual(w.lag, None)
            self.assertAlmostEqual(w.plateau, None)
            self.assertAlmostEqual(w.slope, None)
            self.assertAlmostEqual(w.v, None)
            self.assertAlmostEqual(w.y0, None)
        self.assertEqual(w.max, 313.0)
        self.assertEqual(w.min, 29.0)
        self.assertEqual(w.average_height, 217.82552083333334)

        self.assertRaises(TypeError, w.__add__, 'a')

        w2 = w + w1
        self.assertEqual(w2.id, 'A10')
        self.assertEqual(len(w2), len(times))
        self.assertEqual(len(w2), 384)
        self.assertEqual(max(w2), (95.75, 327.0))
        self.assertEqual(min(w2), (0.0, 63.0))
        self.assertEqual(max(w2, key=lambda x: x[1]),
                (18.25, 357.0))
        self.assertEqual(min(w2, key=lambda x: x[1]),
                (0.25, 55.0))
        self.assertEqual(w2[1], 71.)
        self.assertEqual(w2[12], 316.)
        self.assertEqual(w2[1:5], [71.0, 88.0, 94.0, 94.0])
        self.assertAlmostEqual(w2[1:2:.25][0], 71.0)
        self.assertAlmostEqual(w2[1.3567], 77.7196)
        self.assertEqual(w2.get_raw()[0], (0.0, 63.0))
        self.assertEqual(w2.get_raw()[-1], (95.75, 327.0))
        self.assertEqual(w2.get_times()[0], 0.0)
        self.assertEqual(w2.get_times()[-1], 95.75)
        self.assertEqual(w2.get_signals()[0], 63.0)
        self.assertEqual(w2.get_signals()[-1], 327.0)

        self.assertRaises(TypeError, w.__sub__, 'a')

        w2 = w - w1
        self.assertEqual(w2.id, 'A10')
        self.assertEqual(len(w2), len(times))
        self.assertEqual(len(w2), 384)
        self.assertEqual(max(w2), (95.75, 107.0))
        self.assertEqual(min(w2), (0.0, 11.0))
        self.assertEqual(max(w2, key=lambda x: x[1]),
                (15.75, 274.0))
        self.assertEqual(min(w2, key=lambda x: x[1]),
                (3.25, -20.0))
        self.assertEqual(w2[1], -13.)
        self.assertEqual(w2[12], 228.)
        self.assertEqual(w2[1:5], [-13.0, -18.0, -16.0, -8.0])
        self.assertAlmostEqual(w2[1:2:.25][0], -13.0)
        self.assertAlmostEqual(w2[1.3567], -10.2804)
        self.assertEqual(w2.get_raw()[0], (0.0, 11.0))
        self.assertEqual(w2.get_raw()[-1], (95.75, 107.0))
        self.assertEqual(w2.get_times()[0], 0.0)
        self.assertEqual(w2.get_times()[-1], 95.75)
        self.assertEqual(w2.get_signals()[0], 11.0)
        self.assertEqual(w2.get_signals()[-1], 107.0)
        
        w[1] = 1
            
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
