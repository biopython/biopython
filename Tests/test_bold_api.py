# -*- coding: utf-8 -*-
# Copyright 2014 by Carlos Pena.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
import unittest
import warnings

from Bio import BiopythonWarning
from Bio._py3k import HTTPError
from Bio import MissingExternalDependencyError

from Bio import bold
from Bio.bold import api


class TestApi(unittest.TestCase):

    def setUp(self):
        warnings.simplefilter('ignore', BiopythonWarning)

    def test_call_id(self):
        seq = "TTTTTGGTATTTGAGCAGGAATAGTAGGAACTTCTCTCAGTTTAATTATTCGAATAGAATTAGGTAATCCAGGTTTCTTAATTGGAGATGATCAAATTTATAATACTATTGTAACAGCCCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTGTAATTGGAGGATTTGGAAATTGACTAGTTCCCCTAATATTAGGTGCACCTGATATAGCTTTCCCTCGTATAAATAATATAAGATATTGACTACTTCCACCATCTTTAATATTATTAATTTCAAGTAGTATTGTAGAAAATGGAGCTGGAACAGGTTGAACAGTTTACCCCCCTCTTTCCTCTAATATTGCTCATAGAGGAACCTCAGTAGACTTAGCAATTTTTTCTCTTCATTTAGCTGGTATTTCTTCTATTTTAGGAGCTATTAATTTTATTACTACAATTATTAATATACGAGTTAATGGAATATCCTATGATCAAATACCTTTATTTGTTTGAGCTGTTGGAATTACAGCTCTTCTTTTACTTCTTTCTTTACCTGTTTTAGCAGGAGCTATCACAATACTTCTTACAGATCGAAATTTAAATACATCATTTTTTGATCCTGCAGGAGGAGGTGATCCAATTTTATACCAACATTTATTTTGATTTTTTGGTCACCC"
        db = "COX1_SPECIES_PUBLIC"
        res = bold.call_id(seq, db)
        for item in res.items:
            if item['similarity'] == 1:
                self.assertEqual('Euptychia ordinata', item['taxonomic_identification'])

    def test_call_taxon_search(self):
        taxonomic_identification = 'Euptychia ordinata'
        expected = 302603
        res = bold.call_taxon_search(taxonomic_identification, fuzzy=False)
        item = res.items[0]
        self.assertEqual(expected, item['tax_id'])

        taxonomic_identification = 'Fabaceae'
        res = bold.call_taxon_search(taxonomic_identification, fuzzy=False)
        item = res.items[0]
        self.assertEqual('Plants', item['tax_division'])
        self.assertEqual(187, item['parent_id'])
        self.assertEqual('Fabales', item['parent_name'])
        self.assertEqual('Fabaceae', item['taxon_rep'])

        taxonomic_identification = 'Diplura'
        res = bold.call_taxon_search(taxonomic_identification, fuzzy=False)
        self.assertEqual(2, len(res.items))

    def test_call_taxon_search_returns_empty(self):
        taxonomic_identification = 'Fake species name'
        self.assertRaises(ValueError, bold.call_taxon_search, taxonomic_identification, fuzzy=False)

    def test_call_taxon_search_fuzzy_true(self):
        taxonomic_identification = 'Fabaceae'
        res = bold.call_taxon_search(taxonomic_identification, fuzzy=True)
        self.assertEqual(5, len(res.items))

    def test_call_taxon_search_fuzzy_error(self):
        self.assertRaises(ValueError, bold.call_taxon_search, 'Fabaceae', 'true')

    def test_call_specimen_data(self):
        taxon = 'Euptychia'
        res = bold.call_specimen_data(taxon)
        item = res.items[0]
        self.assertEqual('Nymphalidae', item['taxonomy_family_taxon_name'])

    def test_call_specimen_data_several_taxa(self):
        taxon = 'Euptychia|Mycalesis'
        res = bold.call_specimen_data(taxon)
        self.assertTrue('Mycalesis' in [item['taxonomy_genus_taxon_name'] for item in res.items])

    def test_call_specimen_data_bin(self):
        bin = 'BOLD:AAE2777'
        res = bold.call_specimen_data(bin=bin)
        taxonomy_identifications = []
        append = taxonomy_identifications.append
        for item in res.items:
            if 'taxonomy_identification_provided_by' in item:
                append(item['taxonomy_identification_provided_by'])
        self.assertTrue('Jose Montero' in taxonomy_identifications)

    def test_call_specimen_data_container(self):
        container = 'ACRJP'
        try:
            res = bold.call_specimen_data(container=container)
        except HTTPError:
            # e.g. due to timeout
            raise MissingExternalDependencyError("internet connection failed")

        taxonomy_identifications = []
        append = taxonomy_identifications.append
        for item in res.items:
            if 'taxonomy_identification_provided_by' in item:
                append(item['taxonomy_identification_provided_by'])
        self.assertTrue('Jacques L. Pierre' in taxonomy_identifications)

    def test_call_specimen_data_institutions(self):
        institutions = 'University of Turku'
        res = bold.call_specimen_data(institutions=institutions)
        taxonomy_identifications = []
        append = taxonomy_identifications.append
        for item in res.items:
            if 'taxonomy_identification_provided_by' in item:
                append(item['taxonomy_identification_provided_by'])
        self.assertTrue('Meri Lindqvist' in taxonomy_identifications)

    def test_call_specimen_data_researchers(self):
        researchers = 'Thibaud Decaens'
        res = bold.call_specimen_data(researchers=researchers)
        collection_event_countries = []
        append = collection_event_countries.append
        for item in res.items:
            if 'collection_event_country' in item:
                append(item['collection_event_country'])
        self.assertTrue('Peru' in collection_event_countries)

    def test_call_specimen_data_geo(self):
        geo = 'Iceland'
        res = bold.call_specimen_data(geo=geo)
        collection_event_countries = []
        append = collection_event_countries.append
        for item in res.items:
            if 'collection_event_country' in item:
                append(item['collection_event_country'])
        self.assertTrue('Iceland' in collection_event_countries)

    def test_call_specimen_data_format_tsv(self):
        geo = 'Iceland'
        res = bold.call_specimen_data(geo=geo, format='tsv')
        self.assertTrue('Iceland' in res.items)

    def test_call_specimen_data_wrong_format(self):
        geo = 'Iceland'
        self.assertRaises(ValueError, bold.call_specimen_data, geo=geo, format='csv')

    def test_call_specimen_data_return_empty(self):
        geo = 'Fake country name'
        self.assertRaises(ValueError, bold.call_specimen_data, geo=geo)

    def test_call_taxon_data_basic(self):
        tax_id = 302603
        # using default datatype='basic'
        res = bold.call_taxon_data(tax_id, data_type='basic')
        item = res.items[0]
        self.assertEqual(7044, item['parent_id'])

    def test_call_taxon_data_basic_empty(self):
        tax_id = 302603
        res = bold.call_taxon_data(tax_id)
        item = res.items[0]
        self.assertEqual(7044, item['parent_id'])

    def test_call_taxon_data_includetree_false(self):
        tax_id = 302603
        # using default datatype='basic'
        res = bold.call_taxon_data(tax_id, data_type='basic', include_tree=False)
        item = res.items[0]
        self.assertEqual(7044, item['parent_id'])

    def test_call_taxon_data_includetree_true(self):
        tax_id = 302603
        # using default datatype='basic'
        res = bold.call_taxon_data(tax_id, data_type='basic', include_tree=True)
        self.assertEqual(7, len(res.items))

    def test_call_taxon_data_includetree_error(self):
        tax_id = 302603
        # using default datatype='basic'
        self.assertRaises(ValueError, bold.call_taxon_data, (tax_id, 'basic', 'true'))

    def test_call_sequence_data(self):
        taxon = 'Hermeuptychia'
        geo = 'Peru'
        res = bold.call_sequence_data(taxon=taxon, geo=geo)
        items = res.items
        seq_record_ids = [item.id for item in items]
        self.assertTrue('GBLN4477-14|Hermeuptychia' in seq_record_ids)

    def test_call_sequence_data_returns_empty(self):
        taxon = 'Fake taxon'
        geo = 'Fake country'
        self.assertRaises(ValueError, bold.call_sequence_data, taxon, geo)

    def test_call_full_data(self):
        taxon = 'Hermeuptychia'
        geo = 'Peru'
        res = bold.call_full_data(taxon=taxon, geo=geo)
        genbank_accession_numbers = [item['specimen_identifiers_sample_id'] for item in res.items]
        self.assertTrue('KF466142' in genbank_accession_numbers)

    def test_call_full_data_invalid(self):
        geo = 'Peru'
        format = 'csv'
        self.assertRaises(ValueError, bold.call_full_data, geo=geo, format=format)

    def test_call_trace_files(self):
        taxon = 'Euptychia mollis'
        institutions = 'York University'
        res = bold.call_trace_files(taxon=taxon,
                                    institutions=institutions)
        self.assertNotEqual(res.file_contents, None)

    def test_parse_json(self):
        res = api.Response()

        # call_taxon_search
        json_string = '{"302603":{"taxid":302603,"taxon":"Euptychia ordinata","tax_rank":"species","tax_division":"Animals","parentid":7044,"parentname":"Euptychia"}}'
        res._parse_json(json_string)
        item = res.items[0]
        self.assertEqual(302603, item['tax_id'])
        self.assertEqual(7044, item['parent_id'])

        # data_type = basic
        json_string = '{"taxid":891,"taxon":"Fabaceae","tax_rank":"family","tax_division":"Plants","parentid":187,"parentname":"Fabales","taxonrep":"Fabaceae"}'
        res._parse_json(json_string)
        item = res.items[0]
        self.assertEqual('Fabaceae', item['taxon'])
        self.assertEqual('Plants', item['tax_division'])

        # data_type = images
        json_string = '{"images":[{"copyright_institution":"Smithsonian Tropical Research Institute","specimenid":2616716,"copyright":"Matthew J. MIller","imagequality":4,"photographer":"Oscar Lopez","image":"BSPBB\/MJM_7364_IMG_2240_d+1345758620.JPG","fieldnum":"MJM 7364","sampleid":"MJM 7364","mam_uri":null,"copyright_license":"CreativeCommons - Attribution Non-Commercial","meta":"Dorsal","copyright_holder":"Matthew J. MIller","catalognum":"","copyright_contact":"millerm@si.edu","copyright_year":"2012","taxonrep":"Momotus momota","aspectratio":1.608,"original":true,"external":null}]}'
        res._parse_json(json_string)
        item = res.items[0]
        self.assertEqual('Oscar Lopez', item['images'][0]['photographer'])

        # data_type = geo
        json_string = '{"country":{"Brazil":3,"Mexico":2,"Panama":10,"Guatemala":1,"Peru":13,"Bolivia":6,"Ecuador":2},"sitemap":"http:\/\/www.boldsystems.org\/index.php\/TaxBrowser_Maps_CollectionSites?taxid=88899"}'
        res._parse_json(json_string)
        item = res.items[0]
        self.assertTrue('Brazil' in item['country'].keys())

        # data_type = stats
        json_string = '{"stats":{"publicspecies":2,"publicbins":3,"publicmarkersequences":{"COI-5P":6},"publicrecords":6,"specimenrecords":"45","sequencedspecimens":"25","barcodespecimens":"22","species":"3","barcodespecies":"3"}}'
        res._parse_json(json_string)
        item = res.items[0]
        self.assertTrue('publicspecies' in item['stats'].keys())

        # data_type = sequencinlabs
        json_string = '{"sequencinglabs":{"Smithsonian Tropical Research Institute":7,"Biodiversity Institute of Ontario":13,"Universidade Federal de Minas Gerais":1,"Mined from GenBank":2,"Royal Ontario Museum":2}}'
        res._parse_json(json_string)
        item = res.items[0]
        self.assertTrue('Royal Ontario Museum' in item['sequencinglabs'].keys())

        # data_type = thirdparty
        json_string = r'{"taxid": 88899, "taxon": "Momotus", "tax_rank": "genus", "tax_division": "Animals", "parentid": 88898, "parentname": "Momotidae", "wikipedia_summary": "Momotus</b></i> is a small genus of the motmots, a family of near passerine birds found in forest and woodland of the Neotropics. They have a colourful plumage, which is green on the back becoming blue on the flight feathers and the long tails. The barbs near the ends of the two longest central tail feathers fall off, leaving a length of bare shaft so that tails appear racket-shaped. \n\nMomotus</i> species, like other motmots, eat small prey such as insects and lizards, and will also take fruit. They nest in tunnels in banks, laying about four white eggs.", "wikipedia_link": "http://en.wikipedia.org/wiki/Momotus", "gbif_map": "http://data.gbif.org/species/2475289/overviewMap.png"}'
        res._parse_json(json_string)
        item = res.items[0]
        self.assertTrue('wikipedia_summary' in item.keys())

    def test_parse_data_empty(self):
        result_string = ''
        response = api.Response()
        self.assertRaises(ValueError, response._parse_data, 'call_id', result_string)

    def tearDown(self):
        pass


if __name__ == '__main__':
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
