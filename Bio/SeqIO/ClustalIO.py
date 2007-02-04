# Copyright 2006, 2007 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Clustal parsing based on earlier code by Thomas Sicheritz-Ponten,
# copyright 2001.  See Bio/SeqIO/generic.py

#For reading alignments:
from Bio.Alphabet import generic_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#For writing alignments:
from Bio.SeqIO.Interfaces import SequenceWriter
from Bio.Clustalw import ClustalAlignment

#This is a generator function!
def ClustalIterator(handle, alphabet = generic_alphabet) :
    """Reads a Clustalw file returning a SeqRecord object iterator

    The entire file is loaded at once, but the SeqRecord objects
    are only created "on request".

    For more information on the file format, please see:
    http://www.bioperl.org/wiki/ClustalW_multiple_alignment_format

    You might like to look at Bio.Clustalw which has an interface
    to the command line tool clustalw, and can also clustal alignment
    files into Bio.Clustalw.ClustalAlignment objects.
         
    We call this the "clustal" format which is consistent with EMBOSS.
    Sadly BioPerl calls it the "clustalw" format, so we can't match
    them both.
    """
    line = handle.readline()
    if not line: return
    if not line[:7] == 'CLUSTAL':
        raise SyntaxError("Did not find CLUSTAL header")

    #If the alignment contains entries with the same sequence
    #identifier (not a good idea - but seems possible), then this
    #dictionary based parser will merge their sequences.  Fix this?
    seqs = {}
    ids = []
    while True:
        line = handle.readline()
        if not line: break
        if line[0] == ' ': continue
        fields = line.rstrip().split()
        if not len(fields): continue
        
        #We expect there to be two fields, but on older files
        #there may be a third entry containing a letter count.
        if len(fields) < 2 or len(fields) > 3:
            raise SyntaxError("Could not parse line:\n%s" % line)

        name, seq = fields[0], fields[1]
        if not name in ids: ids.append(name)
        seqs.setdefault(name, '')
        seqs[name] += seq.upper()

        if len(fields) == 3 :
            #This MAY be an old style file with a letter count...
            try :
                oddity = int(fields[2])
            except ValueError :
                raise SyntaxError("Could not parse line, odd third field:\n%s" % line)
            #Check this equals the number of letters (excluding gaps)
            #so far for this sequence?
            if len(seqs[name].replace("-","")) <> oddity :
                raise SyntaxError("Could not parse line, odd third field:\n%s" % line)

    for id in ids :
        yield SeqRecord(Seq(seqs[id], alphabet), id=id)
    
class ClustalWriter(SequenceWriter):
    """Write Clustal sequence alignments"""
    def __init__(self, handle, truncate=10):
        """Creates the writer object

        Use the method write_file() to actually record your sequence records."""
        self.handle = handle
        self.truncate = truncate
    
    def write_file(self, records) :
        """Use this to write an entire file containing the given records.

        records - a SeqRecord iterator, or list of SeqRecords

        Right now this code uses Bio.Clustalw.ClustalAlignment to do
        the hard work - this may change in the future.
        """
        # ToDo - decide if using Bio.Clustalw.ClustalAlignment is
        # actually the best way to handle this.
        #
        # Copying that thirty lines of code (with slight tweaks)
        # would be much simpler, and would probably run quicker and
        # use less memory as it doesn't build a ClustalAlignment
        # object.
        #
        # The downside is code duplication.
        alignment_length = None
        alignment = ClustalAlignment()
        for record in records :
            if alignment_length is None :
                alignment_length = len(record.seq)
            elif alignment_length <> len(record.seq) :
                raise ValueError, "Sequences of different lengths"
            
            #ToDo, check alphabet for this sequence matches that
            #specified for the alignment.  Not sure how the
            #alphabet.contains() method is intended to be used,
            #but it doesn't make sense to me right now.

            #Doing this works, but ClustalAlignment currently uses
            #the record.descrption when outputing the records.
            #alignment._records.append(record)

            #Make sure we don't get any spaces in the record
            #identifier when output in the file by replacing
            #them with underscores:
            alignment.add_sequence(record.id.replace(" ","_"),
                                   record.seq.tostring())

        self.handle.write(str(alignment))
        #Don't close the handle.  Doing so would prevent this code
        #from writing concatenated Clustal files which might be used
        #in phylogenetic bootstrapping (very common with phylip).
        #self.handle.close()


if __name__ == "__main__" :
    # Run a quick self-test

    #This is a truncated version of the example in Tests/cw02.aln
    aln_example1 = \
"""CLUSTAL W (1.81) multiple sequence alignment


gi|4959044|gb|AAD34209.1|AF069      MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN 50
gi|671626|emb|CAA85685.1|           ---------MSPQTETKASVGFKAGVKEYKLTYYTPEYETKDTDILAAFR 41
                                              * *: ::    :.   :*  :  :. : . :*  ::   .

gi|4959044|gb|AAD34209.1|AF069      LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW 100
gi|671626|emb|CAA85685.1|           VTPQPG-----------------VPPEEAGAAVAAESSTGT--------- 65
                                    :   **                  **:...   *.*** ..         

gi|4959044|gb|AAD34209.1|AF069      LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT 150
gi|671626|emb|CAA85685.1|           WTTVWTDGLTSLDRYKG-----RCYHIEPVPG------------------ 92
                                     .:*   * *: .* :*        : :* .*                  

gi|4959044|gb|AAD34209.1|AF069      SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE 200
gi|671626|emb|CAA85685.1|           -EKDQCICYVAYPLDLFEEGSVTNMFTSIVGNVFGFKALRALRLEDLRIP 141
                                     *::.  .    .:: :*..*  :* .*   .. .  :    .  :    

gi|4959044|gb|AAD34209.1|AF069      VPTTRAQRRA 210
gi|671626|emb|CAA85685.1|           VAYVKTFQGP 151
                                    *. .:: : .
                                     
"""                 

    #This example is a truncated version of the dataset used here:
    #http://virgil.ruc.dk/kurser/Sekvens/Treedraw.htm
    aln_example2 = \
"""CLUSTAL X (1.83) multiple sequence alignment


V_Harveyi_PATH                 --MKNWIKVAVAAIA--LSAA------------------TVQAATEVKVG
B_subtilis_YXEM                MKMKKWTVLVVAALLAVLSACG------------NGNSSSKEDDNVLHVG
B_subtilis_GlnH_homo_YCKK      MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVG
YA80_HAEIN                     MKKLLFTTALLTGAIAFSTF-----------SHAGEIADRVEKTKTLLVG
FLIY_ECOLI                     MKLAHLGRQALMGVMAVALVAG---MSVKSFADEG-LLNKVKERGTLLVG
E_coli_GlnH                    --MKSVLKVSLAALTLAFAVS------------------SHAADKKLVVA
Deinococcus_radiodurans        -MKKSLLSLKLSGLLVPSVLALS--------LSACSSPSSTLNQGTLKIA
HISJ_E_COLI                    MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG
                                         : .                                 : :.

V_Harveyi_PATH                 MSGRYFPFTFVKQ--DKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGL
B_subtilis_YXEM                ATGQSYPFAYKEN--GKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGE
B_subtilis_GlnH_homo_YCKK      TEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAG
YA80_HAEIN                     TEGTYAPFTFHDK-SGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAG
FLIY_ECOLI                     LEGTYPPFSFQGD-DGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLAS
E_coli_GlnH                    TDTAFVPFEFKQG--DKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPA
Deinococcus_radiodurans        MEGTYPPFTSKNE-QGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAG
HISJ_E_COLI                    TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS
                                     **       .:  *::::.   : :.   .        ..:   

V_Harveyi_PATH                 LETGRIDTISNQITMTDARKAKYLFADPYVVDG-AQI
B_subtilis_YXEM                LQTGKLDTISNQVAVTDERKETYNFTKPYAYAG-TQI
B_subtilis_GlnH_homo_YCKK      LNSKRFDVVANQVG-KTDREDKYDFSDKYTTSR-AVV
YA80_HAEIN                     LNAKRFDVIANQTNPSPERLKKYSFTTPYNYSG-GVI
FLIY_ECOLI                     LDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQAL
E_coli_GlnH                    LQTKNVDLALAGITITDERKKAIDFSDGYYKSG-LLV
Deinococcus_radiodurans        LQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEII
HISJ_E_COLI                    LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLV
                               *.: . *        .  *     *:          :

"""

    from StringIO import StringIO

    records = list(ClustalIterator(StringIO(aln_example1)))
    assert 2 == len(records)
    assert records[0].id == "gi|4959044|gb|AAD34209.1|AF069"
    assert records[1].id == "gi|671626|emb|CAA85685.1|"
    assert records[0].seq.tostring() == \
          "MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN" + \
          "LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW" + \
          "LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT" + \
          "SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE" + \
          "VPTTRAQRRA"

    records = list(ClustalIterator(StringIO(aln_example2)))
    assert 8 == len(records)
    assert records[-1].id == "HISJ_E_COLI"
    assert records[-1].seq.tostring() == \
          "MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG" + \
          "TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS" + \
          "LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLV"
