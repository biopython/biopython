#!/usr/bin/env python
# Copyright 2010 by Andrea Pierleoni
# Revisions copyright 2010-2013 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Test for the Uniprot parser on Uniprot XML files.
"""
import os
import unittest

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

#Left as None if the import within UniProtIO fails
if SeqIO.UniprotIO.ElementTree is None:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError("No ElementTree module was found. "
                            "Use Python 2.5+, lxml or elementtree if you "
                            "want to use Bio.SeqIO.UniprotIO.")

from seq_tests_common import compare_reference, compare_record


class TestUniprot(unittest.TestCase):

    def test_uni001(self):
        "Parsing Uniprot file uni001"
        filename = 'uni001'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        with open(datafile) as test_handle:
            seq_record = SeqIO.read(test_handle, "uniprot-xml")

        self.assertTrue(isinstance(seq_record, SeqRecord))

        # test a couple of things on the record -- this is not exhaustive
        self.assertEqual(seq_record.id, "Q91G55")
        self.assertEqual(seq_record.name, "043L_IIV6")
        self.assertEqual(seq_record.description, "Uncharacterized protein 043L")
        self.assertEqual(repr(seq_record.seq), "Seq('MDLINNKLNIEIQKFCLDLEKKYNINYNNLIDLWFNKESTERLIKCEVNLENKI...IPI', ProteinAlphabet())")

        # self.assertEqual(seq_record.accessions, ['Q91G55']) #seq_record.accessions does not exist
        # self.assertEqual(seq_record.organism_classification, ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Mammalia', 'Eutheria', 'Primates', 'Catarrhini', 'Hominidae', 'Homo'])
        # self.assertEqual(record.seqinfo, (348, 39676, '75818910'))

        self.assertEqual(len(seq_record.features), 1)
        self.assertEqual(repr(seq_record.features[0]), "SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(116)), type='chain', id='PRO_0000377969')")

        self.assertEqual(len(seq_record.annotations['references']), 2)
        self.assertEqual(seq_record.annotations['references'][0].authors, 'Jakob N.J., Mueller K., Bahr U., Darai G.')
        self.assertEqual(seq_record.annotations['references'][0].title, 'Analysis of the first complete DNA sequence of an invertebrate iridovirus: coding strategy of the genome of Chilo iridescent virus.')
        self.assertEqual(seq_record.annotations['references'][0].journal, 'Virology 286:182-196(2001)')
        self.assertEqual(seq_record.annotations['references'][0].comment, 'journal article | 2001 | Scope: NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA] | ')

        self.assertEqual(len(seq_record.dbxrefs), 11)
        self.assertEqual(seq_record.dbxrefs[0], 'DOI:10.1006/viro.2001.0963')

        self.assertEqual(seq_record.annotations['sequence_length'], 116)
        self.assertEqual(seq_record.annotations['sequence_checksum'], '4A29B35FB716523C')
        self.assertEqual(seq_record.annotations['modified'], '2009-07-07')
        self.assertEqual(seq_record.annotations['accessions'], ['Q91G55'])
        self.assertEqual(seq_record.annotations['taxonomy'], ['Viruses', 'dsDNA viruses, no RNA stage', 'Iridoviridae', 'Iridovirus'])
        self.assertEqual(seq_record.annotations['sequence_mass'], 13673)
        self.assertEqual(seq_record.annotations['dataset'], 'Swiss-Prot')
        self.assertEqual(seq_record.annotations['gene_name_ORF'], ['IIV6-043L'])
        self.assertEqual(seq_record.annotations['version'], 21)
        self.assertEqual(seq_record.annotations['sequence_modified'], '2001-12-01')
        self.assertEqual(seq_record.annotations['keywords'], ['Complete proteome', 'Virus reference strain'])
        self.assertEqual(seq_record.annotations['organism_host'], ['Acheta domesticus', 'House cricket', 'Chilo suppressalis', 'striped riceborer', 'Gryllus bimaculatus', 'Two-spotted cricket', 'Gryllus campestris', 'Spodoptera frugiperda', 'Fall armyworm'])
        self.assertEqual(seq_record.annotations['created'], '2009-06-16')
        self.assertEqual(seq_record.annotations['organism_name'], ['Chilo iridescent virus'])
        self.assertEqual(seq_record.annotations['organism'], 'Invertebrate iridescent virus 6 (IIV-6)')
        self.assertEqual(seq_record.annotations['recommendedName_fullName'], ['Uncharacterized protein 043L'])
        self.assertEqual(seq_record.annotations['sequence_version'], 1)
        self.assertEqual(seq_record.annotations['proteinExistence'], ['Predicted'])

    def test_uni003(self):
        "Parsing Uniprot file uni003"
        filename = 'uni003'
        # test the record parser

        datafile = os.path.join('SwissProt', filename)

        test_handle = open(datafile)
        seq_record = SeqIO.read(test_handle, "uniprot-xml")
        test_handle.close()

        self.assertTrue(isinstance(seq_record, SeqRecord))

        # test general record entries
        self.assertEqual(seq_record.id, "O44185")
        self.assertEqual(seq_record.name, "FLP13_CAEEL")
        self.assertEqual(seq_record.description,
                "FMRFamide-like neuropeptides 13")
        self.assertEqual(repr(seq_record.seq),
                "Seq('MMTSLLTISMFVVAIQAFDSSEIRMLDEQYDTKNPFFQFLENSKRSDRPTRAMD...GRK', ProteinAlphabet())")

        self.assertEqual(len(seq_record.annotations['references']), 7)
        self.assertEqual(seq_record.annotations['references'][5].authors,
                'Kim K., Li C.')
        self.assertEqual(seq_record.annotations['references'][5].title,
                'Expression and regulation of an FMRFamide-related '
                'neuropeptide gene family in Caenorhabditis elegans.')
        self.assertEqual(seq_record.annotations['references'][5].journal,
                'J. Comp. Neurol. 475:540-550(2004)')
        self.assertEqual(seq_record.annotations['references'][5].comment,
                'journal article | 2004 | Scope: TISSUE SPECIFICITY, '
                'DEVELOPMENTAL STAGE | ')

        self.assertEqual(seq_record.annotations["accessions"], ['O44185'])
        self.assertEqual(seq_record.annotations["created"], "2004-05-10")
        self.assertEqual(seq_record.annotations["dataset"], "Swiss-Prot")
        self.assertEqual(seq_record.annotations["gene_name_ORF"], ['F33D4.3'])
        self.assertEqual(seq_record.annotations["gene_name_primary"], "flp-13")
        self.assertEqual(seq_record.annotations["keywords"],
                ['Amidation', 'Cleavage on pair of basic residues',
                'Complete proteome', 'Direct protein sequencing',
                'Neuropeptide', 'Reference proteome', 'Repeat',
                'Secreted', 'Signal'])
        self.assertEqual(seq_record.annotations["modified"], "2012-11-28")
        self.assertEqual(seq_record.annotations["organism"],
                "Caenorhabditis elegans")
        self.assertEqual(seq_record.annotations["proteinExistence"],
                ['evidence at protein level'])
        self.assertEqual(seq_record.annotations["recommendedName_fullName"],
                ['FMRFamide-like neuropeptides 13'])
        self.assertEqual(seq_record.annotations["sequence_length"], 160)
        self.assertEqual(seq_record.annotations["sequence_checksum"],
                "BE4C24E9B85FCD11")
        self.assertEqual(seq_record.annotations["sequence_mass"], 17736)
        self.assertEqual(seq_record.annotations["sequence_modified"], "1998-06-01")
        self.assertEqual(seq_record.annotations["sequence_precursor"], "true")
        self.assertEqual(seq_record.annotations["sequence_version"], 1)
        self.assertEqual(seq_record.annotations["taxonomy"],
                ['Eukaryota', 'Metazoa', 'Ecdysozoa', 'Nematoda',
                'Chromadorea', 'Rhabditida', 'Rhabditoidea', 'Rhabditidae',
                'Peloderinae', 'Caenorhabditis'])
        self.assertEqual(seq_record.annotations["type"],
                ['ECO:0000006', 'ECO:0000001'])
        self.assertEqual(seq_record.annotations["version"], 74)

        # test comment entries
        self.assertEqual(seq_record.annotations["comment_allergen"],
                ['Causes an allergic reaction in human.'])
        self.assertEqual(seq_record.annotations["comment_alternativeproducts_isoform"],
                ['Q8W1X2-1', 'Q8W1X2-2'])
        self.assertEqual(seq_record.annotations["comment_biotechnology"],
                ['Green fluorescent protein has been engineered to produce a '
                'vast number of variously colored mutants, fusion proteins, '
                'and biosensors. Fluorescent proteins and its mutated allelic '
                'forms, blue, cyan and yellow have become a useful and '
                'ubiquitous tool for making chimeric proteins, where they '
                'function as a fluorescent protein tag. Typically they '
                'tolerate N- and C-terminal fusion to a broad variety of '
                'proteins. They have been expressed in most known cell types '
                'and are used as a noninvasive fluorescent marker in living '
                'cells and organisms. They enable a wide range of applications '
                'where they have functioned as a cell lineage tracer, reporter '
                'of gene expression, or as a measure of protein-protein '
                'interactions.', 'Can also be used as a molecular thermometer, '
                'allowing accurate temperature measurements in fluids. The '
                'measurement process relies on the detection of the blinking '
                'of GFP using fluorescence correlation spectroscopy.'])
        self.assertEqual(seq_record.annotations["comment_catalyticactivity"],
                ['ATP + acetyl-CoA + HCO(3)(-) = ADP + phosphate + malonyl-CoA.',
                'ATP + biotin-[carboxyl-carrier-protein] + CO(2) = ADP + '
                'phosphate + carboxy-biotin-[carboxyl-carrier-protein].'])
        self.assertEqual(seq_record.annotations["comment_caution"],
                ['Could be the product of a pseudogene. The existence of a '
                'transcript at this locus is supported by only one sequence '
                'submission (PubMed:2174397).'])
        self.assertEqual(seq_record.annotations["comment_cofactor"],
                ['Biotin (By similarity).', 'Binds 2 manganese ions per '
                'subunit (By similarity).'])
        self.assertEqual(seq_record.annotations["comment_developmentalstage"],
                ['Expressed from the comma stage of embryogenesis, during all '
                'larval stages, and in low levels in adults.'])
        self.assertEqual(seq_record.annotations["comment_disease"],
                ['Defects in MC2R are the cause of glucocorticoid deficiency '
                'type 1 (GCCD1) [MIM:202200]; also known as familial '
                'glucocorticoid deficiency type 1 (FGD1). GCCD1 is an '
                'autosomal recessive disorder due to congenital '
                'insensitivity or resistance to adrenocorticotropin (ACTH). '
                'It is characterized by progressive primary adrenal '
                'insufficiency, without mineralocorticoid deficiency.'])
        self.assertEqual(seq_record.annotations["comment_disruptionphenotype"],
                ['Mice display impaired B-cell development which does not '
                'progress pass the progenitor stage.'])
        self.assertEqual(seq_record.annotations["comment_domain"],
                ['Two regions, an N-terminal (aa 96-107) and a C-terminal '
                '(aa 274-311) are required for binding FGF2.'])
        self.assertEqual(seq_record.annotations["comment_enzymeregulation"],
                ['By phosphorylation. The catalytic activity is inhibited by '
                'soraphen A, a polyketide isolated from the myxobacterium '
                'Sorangium cellulosum and a potent inhibitor of fungal growth.'])
        self.assertEqual(seq_record.annotations["comment_function"],
                ['FMRFamides and FMRFamide-like peptides are neuropeptides. '
                'AADGAPLIRF-amide and APEASPFIRF-amide inhibit muscle tension '
                'in somatic muscle. APEASPFIRF-amide is a potent inhibitor of '
                'the activity of dissected pharyngeal myogenic muscle system.'])
        self.assertEqual(seq_record.annotations["comment_induction"],
                ['Repressed in presence of fatty acids. Repressed 3-fold by '
                'lipid precursors, inositol and choline, and also controlled '
                'by regulatory factors INO2, INO4 and OPI1.'])
        self.assertEqual(seq_record.annotations["comment_interaction_intactId"],
                ['EBI-356720', 'EBI-746969', 'EBI-720116'])
        self.assertEqual(seq_record.annotations["comment_massspectrometry"],
                ['88..98:1032|MALDI', '100..110:1133.7|MALDI'])
        self.assertEqual(seq_record.annotations["comment_miscellaneous"],
                ['Present with 20200 molecules/cell in log phase SD medium.'])
        self.assertEqual(seq_record.annotations["comment_onlineinformation"],
                ['NIEHS-SNPs@http://egp.gs.washington.edu/data/api5/'])
        self.assertEqual(seq_record.annotations["comment_pathway"],
                ['Lipid metabolism; malonyl-CoA biosynthesis; malonyl-CoA '
                'from acetyl-CoA: step 1/1.'])
        self.assertEqual(seq_record.annotations["comment_RNAediting"],
                ['Partially edited. RNA editing generates receptor isoforms '
                'that differ in their ability to interact with the '
                'phospholipase C signaling cascade in a transfected cell '
                'line, suggesting that this RNA processing event may '
                'contribute to the modulation of serotonergic '
                'neurotransmission in the central nervous system.'])
        self.assertEqual(seq_record.annotations["comment_PTM"],
                ['Acetylation at Lys-251 impairs antiapoptotic function.'])
        self.assertEqual(seq_record.annotations["comment_pharmaceutical"],
                ['Could be used as a possible therapeutic agent for treating '
                'rheumatoid arthritis.'])
        self.assertEqual(seq_record.annotations["comment_polymorphism"],
                ['Position 23 is polymorphic; the frequencies in unrelated '
                'Caucasians are 0.87 for Cys and 0.13 for Ser.'])
        self.assertEqual(seq_record.annotations["comment_similarity"],
                ['Belongs to the FARP (FMRFamide related peptide) family.'])
        self.assertEqual(seq_record.annotations["comment_subcellularlocation_location"],
                ['Secreted'])
        self.assertEqual(seq_record.annotations["comment_subunit"],
                ['Homodimer.'])
        self.assertEqual(seq_record.annotations["comment_tissuespecificity"],
                ['Each flp gene is expressed in a distinct set of neurons. '
                'Flp-13 is expressed in the ASE sensory neurons, the DD motor '
                'neurons, the 15, M3 and M5 cholinergic pharyngeal '
                'motoneurons, and the ASG, ASK and BAG neurons.'])
        self.assertEqual(seq_record.annotations["comment_toxicdose"],
                ['LD(50) is 50 ug/kg in mouse by intracerebroventricular '
                'injection and 600 ng/g in Blatella germanica.'])

    def compare_txt_xml(self, old, new):
        self.assertEqual(old.id, new.id)
        self.assertEqual(old.name, new.name)
        self.assertEqual(len(old), len(new))
        self.assertEqual(str(old.seq), str(new.seq))
        for key in set(old.annotations).intersection(new.annotations):
            if key == "references":
                self.assertEqual(len(old.annotations[key]),
                                 len(new.annotations[key]))
                for r1, r2 in zip(old.annotations[key], new.annotations[key]):
                    #Tweak for line breaks in plain text SwissProt
                    r1.title = r1.title.replace("- ", "-")
                    r2.title = r2.title.replace("- ", "-")
                    r1.journal = r1.journal.rstrip(".")  # Should parser do this?
                    r1.medline_id = ""  # Missing in UniPort XML? TODO - check
                    #Lots of extra comments in UniProt XML
                    r1.comment = ""
                    r2.comment = ""
                    if not r2.journal:
                        r1.journal = ""
                    compare_reference(r1, r2)
            elif old.annotations[key] == new.annotations[key]:
                pass
            elif key in ["date"]:
                #TODO - Why is this a list vs str?
                pass
            elif type(old.annotations[key]) != type(new.annotations[key]):
                raise TypeError("%s gives %s vs %s" %
                                 (key, old.annotations[key], new.annotations[key]))
            elif key in ["organism"]:
                if old.annotations[key] == new.annotations[key]:
                    pass
                elif old.annotations[key].startswith(new.annotations[key]+" "):
                    pass
                else:
                    raise ValueError(key)
            elif isinstance(old.annotations[key], list) \
            and sorted(old.annotations[key]) == sorted(new.annotations[key]):
                pass
            else:
                raise ValueError("%s gives %s vs %s" %
                                 (key, old.annotations[key], new.annotations[key]))
        self.assertEqual(len(old.features), len(new.features),
                         "Features in %s, %i vs %i" %
                         (old.id, len(old.features), len(new.features)))
        for f1, f2 in zip(old.features, new.features):
            """
            self.assertEqual(f1.location.nofuzzy_start, f2.location.nofuzzy_start,
                             "%s %s vs %s %s" %
                             (f1.location, f1.type, f2.location, f2.type))
            self.assertEqual(f1.location.nofuzzy_end, f2.location.nofuzzy_end,
                             "%s %s vs %s %s" %
                             (f1.location, f1.type, f2.location, f2.type))
            """
            self.assertEqual(repr(f1.location), repr(f2.location),
                            "%s %s vs %s %s" %
                            (f1.location, f1.type, f2.location, f2.type))

    def test_Q13639(self):
        """Compare SwissProt text and uniprot XML versions of Q13639."""
        old = SeqIO.read("SwissProt/Q13639.txt", "swiss")
        new = SeqIO.read("SwissProt/Q13639.xml", "uniprot-xml")
        self.compare_txt_xml(old, new)

    def test_multi_ex(self):
        """Compare SwissProt text and uniprot XML versions of several examples."""
        txt_list = list(SeqIO.parse("SwissProt/multi_ex.txt", "swiss"))
        xml_list = list(SeqIO.parse("SwissProt/multi_ex.xml", "uniprot-xml"))
        fas_list = list(SeqIO.parse("SwissProt/multi_ex.fasta", "fasta"))
        with open("SwissProt/multi_ex.list") as handle:
            ids = [x.strip() for x in handle]
        self.assertEqual(len(txt_list), len(ids))
        self.assertEqual(len(txt_list), len(fas_list))
        self.assertEqual(len(txt_list), len(xml_list))
        for txt, xml, fas, id in zip(txt_list, xml_list, fas_list, ids):
            self.assertEqual(txt.id, id)
            self.assertTrue(txt.id in fas.id.split("|"))
            self.assertEqual(str(txt.seq), str(fas.seq))
            self.compare_txt_xml(txt, xml)

    def test_multi_ex_index(self):
        """Index SwissProt text and uniprot XML versions of several examples."""
        txt_list = list(SeqIO.parse("SwissProt/multi_ex.txt", "swiss"))
        xml_list = list(SeqIO.parse("SwissProt/multi_ex.xml", "uniprot-xml"))
        with open("SwissProt/multi_ex.list") as handle:
            ids = [x.strip() for x in handle]
        txt_index = SeqIO.index("SwissProt/multi_ex.txt", "swiss")
        xml_index = SeqIO.index("SwissProt/multi_ex.xml", "uniprot-xml")
        self.assertEqual(sorted(txt_index), sorted(ids))
        self.assertEqual(sorted(xml_index), sorted(ids))
        #Check SeqIO.parse() versus SeqIO.index() for plain text "swiss"
        for old in txt_list:
            new = txt_index[old.id]
            compare_record(old, new)
        #Check SeqIO.parse() versus SeqIO.index() for XML "uniprot-xml"
        for old in xml_list:
            new = xml_index[old.id]
            compare_record(old, new)
        txt_index.close()
        xml_index.close()

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
