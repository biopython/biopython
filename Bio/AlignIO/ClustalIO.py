# Copyright 2006-2010 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""
Bio.AlignIO support for the "clustal" output from CLUSTAL W and other tools.

You are expected to use this module via the Bio.AlignIO functions (or the
Bio.SeqIO functions if you want to work directly with the gapped sequences).
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Interfaces import AlignmentIterator, SequentialAlignmentWriter

class ClustalWriter(SequentialAlignmentWriter):
    """Clustalw alignment writer."""
    def write_alignment(self, alignment):
        """Use this to write (another) single alignment to an open file."""

        if len(alignment) == 0:
            raise ValueError("Must have at least one sequence")
        if alignment.get_alignment_length() == 0:
            #This doubles as a check for an alignment object    
            raise ValueError("Non-empty sequences are required")

        #Old versions of the parser in Bio.Clustalw used a ._version property,
        try:
            version = str(alignment._version)
        except AttributeError:
            version = ""
        if not version:
            version = '1.81'
        if version.startswith("2."):
            #e.g. 2.0.x
            output = "CLUSTAL %s multiple sequence alignment\n\n\n" % version
        else:
            #e.g. 1.81 or 1.83
            output = "CLUSTAL X (%s) multiple sequence alignment\n\n\n" % version
        
        cur_char = 0
        max_length = len(alignment[0])

        if max_length <= 0:
            raise ValueError("Non-empty sequences are required")

        # keep displaying sequences until we reach the end
        while cur_char != max_length:
            # calculate the number of sequences to show, which will
            # be less if we are at the end of the sequence
            if (cur_char + 50) > max_length:
                show_num = max_length - cur_char
            else:
                show_num = 50

            # go through all of the records and print out the sequences
            # when we output, we do a nice 80 column output, although this
            # may result in truncation of the ids.
            for record in alignment:
                #Make sure we don't get any spaces in the record
                #identifier when output in the file by replacing
                #them with underscores:
                line = record.id[0:30].replace(" ","_").ljust(36)
                line += record.seq[cur_char:(cur_char + show_num)].tostring()
                output += line + "\n"

            # now we need to print out the star info, if we've got it
            # This was stored by Bio.Clustalw using a ._star_info property.
            if hasattr(alignment, "_star_info") and alignment._star_info != '':
                output += (" " * 36) + \
                     alignment._star_info[cur_char:(cur_char + show_num)] + "\n"

            output += "\n"
            cur_char += show_num

        # Want a trailing blank new line in case the output is concatenated
        self.handle.write(output + "\n")

class ClustalIterator(AlignmentIterator):
    """Clustalw alignment iterator."""

    def next(self):
        handle = self.handle
        try:
            #Header we saved from when we were parsing
            #the previous alignment.
            line = self._header
            del self._header
        except AttributeError:      
            line = handle.readline()
        if not line:
            raise StopIteration

        #Whitelisted headers we know about
        known_headers = ['CLUSTAL', 'PROBCONS', 'MUSCLE']
        if line.strip().split()[0] not in known_headers:
            raise ValueError("%s is not a known CLUSTAL header: %s" % \
                             (line.strip().split()[0],
                              ", ".join(known_headers)))

        # find the clustal version in the header line
        version = None
        for word in line.split():
            if word[0]=='(' and word[-1]==')':
                word = word[1:-1]
            if word[0] in '0123456789':
                version = word
                break

        #There should be two blank lines after the header line
        line = handle.readline()
        while line.strip() == "":
            line = handle.readline()

        #If the alignment contains entries with the same sequence
        #identifier (not a good idea - but seems possible), then this
        #dictionary based parser will merge their sequences.  Fix this?
        ids = []
        seqs = []
        consensus = ""
        seq_cols = None #: Used to extract the consensus

        #Use the first block to get the sequence identifiers
        while True:
            if line[0] != " " and line.strip() != "":
                #Sequences identifier...
                fields = line.rstrip().split()

                #We expect there to be two fields, there can be an optional
                #"sequence number" field containing the letter count.
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError("Could not parse line:\n%s" % line)

                ids.append(fields[0])
                seqs.append(fields[1])

                #Record the sequence position to get the consensus
                if seq_cols is None:
                    start = len(fields[0]) + line[len(fields[0]):].find(fields[1])
                    end = start + len(fields[1])
                    seq_cols = slice(start, end)
                    del start, end
                assert fields[1] == line[seq_cols]

                if len(fields) == 3:
                    #This MAY be an old style file with a letter count...
                    try:
                        letters = int(fields[2])
                    except ValueError:
                        raise ValueError("Could not parse line, bad sequence number:\n%s" % line)
                    if len(fields[1].replace("-","")) != letters:
                        raise ValueError("Could not parse line, invalid sequence number:\n%s" % line)
            elif line[0] == " ":
                #Sequence consensus line...
                assert len(ids) == len(seqs)
                assert len(ids) > 0
                assert seq_cols is not None
                consensus = line[seq_cols]
                assert not line[:seq_cols.start].strip()
                assert not line[seq_cols.stop:].strip()
                #Check for blank line (or end of file)
                line = handle.readline()
                assert line.strip() == ""
                break
            else:
                #No consensus
                break
            line = handle.readline()
            if not line : break #end of file

        assert line.strip() == ""
        assert seq_cols is not None

        #Confirm all same length
        for s in seqs:
            assert len(s) == len(seqs[0])
        if consensus:
            assert len(consensus) == len(seqs[0])

        #Loop over any remaining blocks...
        done = False
        while not done:
            #There should be a blank line between each block.
            #Also want to ignore any consensus line from the
            #previous block.
            while (not line) or line.strip() == "":
                line = handle.readline()
                if not line : break # end of file
            if not line : break # end of file

            if line.split(None,1)[0] in known_headers:
                #Found concatenated alignment.
                done = True
                self._header = line
                break

            for i in range(len(ids)):
                assert line[0] != " ", "Unexpected line:\n%s" % repr(line)
                fields = line.rstrip().split()
                
                #We expect there to be two fields, there can be an optional
                #"sequence number" field containing the letter count.
                if len(fields) < 2 or len(fields) > 3:
                    raise ValueError("Could not parse line:\n%s" % repr(line))

                if fields[0] != ids[i]:
                    raise ValueError("Identifiers out of order? Got '%s' but expected '%s'" \
                                      % (fields[0], ids[i]))

                if fields[1] != line[seq_cols]:
                    start = len(fields[0]) + line[len(fields[0]):].find(fields[1])
                    assert start == seq_cols.start, 'Old location %s -> %i:XX' % (seq_cols, start)
                    end = start + len(fields[1])
                    seq_cols = slice(start, end)
                    del start, end

                #Append the sequence
                seqs[i] += fields[1]
                assert len(seqs[i]) == len(seqs[0])

                if len(fields) == 3:
                    #This MAY be an old style file with a letter count...
                    try:
                        letters = int(fields[2])
                    except ValueError:
                        raise ValueError("Could not parse line, bad sequence number:\n%s" % line)
                    if len(seqs[i].replace("-","")) != letters:
                        raise ValueError("Could not parse line, invalid sequence number:\n%s" % line)

                #Read in the next line
                line = handle.readline()
            #There should now be a consensus line
            if consensus:
                assert line[0] == " "
                assert seq_cols is not None
                consensus += line[seq_cols]
                assert len(consensus) == len(seqs[0])
                assert not line[:seq_cols.start].strip()
                assert not line[seq_cols.stop:].strip()
                #Read in the next line
                line = handle.readline()
            

        assert len(ids) == len(seqs)
        if len(seqs) == 0 or len(seqs[0]) == 0:
            raise StopIteration

        if self.records_per_alignment is not None \
        and self.records_per_alignment != len(ids):
            raise ValueError("Found %i records in this alignment, told to expect %i" \
                             % (len(ids), self.records_per_alignment))

        records = (SeqRecord(Seq(s, self.alphabet), id=i, description=i) \
                   for (i,s) in zip(ids, seqs)) 
        alignment = MultipleSeqAlignment(records, self.alphabet)
        #TODO - Handle alignment annotation better, for now
        #mimic the old parser in Bio.Clustalw
        if version:
            alignment._version = version
        if consensus:
            alignment_length = len(seqs[0])
            assert len(consensus) == alignment_length, \
                   "Alignment length is %i, consensus length is %i, '%s'" \
                   % (alignment_length, len(consensus), consensus)
            alignment._star_info = consensus
        return alignment
    
if __name__ == "__main__":
    print "Running a quick self-test"

    #This is a truncated version of the example in Tests/cw02.aln
    #Notice the inclusion of sequence numbers (right hand side)
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
    #with the last record repeated twice (deliberate toture test)
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
HISJ_E_COLI                    LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLV
                               *.: . *        .  *     *:          :

"""

    from StringIO import StringIO

    alignments = list(ClustalIterator(StringIO(aln_example1)))
    assert 1 == len(alignments)
    assert alignments[0]._version == "1.81"
    alignment = alignments[0]
    assert 2 == len(alignment)
    assert alignment[0].id == "gi|4959044|gb|AAD34209.1|AF069"
    assert alignment[1].id == "gi|671626|emb|CAA85685.1|"
    assert alignment[0].seq.tostring() == \
          "MENSDSNDKGSDQSAAQRRSQMDRLDREEAFYQFVNNLSEEDYRLMRDNN" + \
          "LLGTPGESTEEELLRRLQQIKEGPPPQSPDENRAGESSDDVTNSDSIIDW" + \
          "LNSVRQTGNTTRSRQRGNQSWRAVSRTNPNSGDFRFSLEINVNRNNGSQT" + \
          "SENESEPSTRRLSVENMESSSQRQMENSASESASARPSRAERNSTEAVTE" + \
          "VPTTRAQRRA"

    alignments = list(ClustalIterator(StringIO(aln_example2)))
    assert 1 == len(alignments)
    assert alignments[0]._version == "1.83"
    alignment = alignments[0]
    assert 9 == len(alignment)
    assert alignment[-1].id == "HISJ_E_COLI"
    assert alignment[-1].seq.tostring() == \
          "MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG" + \
          "TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS" + \
          "LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLV"

    for alignment in ClustalIterator(StringIO(aln_example2 + aln_example1)):
        print "Alignment with %i records of length %i" \
              % (len(alignment),
                 alignment.get_alignment_length())

    print "Checking empty file..."
    assert 0 == len(list(ClustalIterator(StringIO(""))))

    print "Checking write/read..."
    alignments = list(ClustalIterator(StringIO(aln_example1))) \
               + list(ClustalIterator(StringIO(aln_example2)))*2
    handle = StringIO()
    ClustalWriter(handle).write_file(alignments)
    handle.seek(0)
    for i,a in enumerate(ClustalIterator(handle)):
        assert a.get_alignment_length() == alignments[i].get_alignment_length()
    handle.seek(0)

    print "Testing write/read when there is only one sequence..."
    alignment = alignment[0:1]
    handle = StringIO()
    ClustalWriter(handle).write_file([alignment])
    handle.seek(0)
    for i,a in enumerate(ClustalIterator(handle)):
        assert a.get_alignment_length() == alignment.get_alignment_length()
        assert len(a) == 1

    aln_example3 = \
"""CLUSTAL 2.0.9 multiple sequence alignment


Test1seq             ------------------------------------------------------------
AT3G20900.1-SEQ      ATGAACAAAGTAGCGAGGAAGAACAAAACATCAGGTGAACAAAAAAAAAACTCAATCCAC
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             -----AGTTACAATAACTGACGAAGCTAAGTAGGCTACTAATTAACGTCATCAACCTAAT
AT3G20900.1-SEQ      ATCAAAGTTACAATAACTGACGAAGCTAAGTAGGCTAGAAATTAAAGTCATCAACCTAAT
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             ACATAGCACTTAGAAAAAAGTGAAGTAAGAAAATATAAAATAATAAAAGGGTGGGTTATC
AT3G20900.1-SEQ      ACATAGCACTTAGAAAAAAGTGAAGCAAGAAAATATAAAATAATAAAAGGGTGGGTTATC
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             AATTGATAGTGTAAATCATCGTATTCCGGTGATATACCCTACCACAAAAACTCAAACCGA
AT3G20900.1-SEQ      AATTGATAGTGTAAATCATAGTTGATTTTTGATATACCCTACCACAAAAACTCAAACCGA
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             CTTGATTCAAATCATCTCAATAAATTAGCGCCAAAATAATGAAAAAAATAATAACAAACA
AT3G20900.1-SEQ      CTTGATTCAAATCATCTCAAAAAACAAGCGCCAAAATAATGAAAAAAATAATAACAAAAA
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             AAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATT
AT3G20900.1-SEQ      CAAACAAACCAAAATAAGAAAAAACATTACGCAAAACATAATAATTTACTCTTCGTTATT
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             GTATTAACAAATCAAAGAGCTGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGAT
AT3G20900.1-SEQ      GTATTAACAAATCAAAGAGATGAATTTTGATCACCTGCTAATACTACTTTCTGTATTGAT
AT3G20900.1-CDS      ------------------------------------------------------------
                                                                                 

Test1seq             CCTATATCAACGTAAACAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGT
AT3G20900.1-SEQ      CCTATATCAAAAAAAAAAAAGATACTAATAATTAACTAAAAGTACGTTCATCGATCGTGT
AT3G20900.1-CDS      ------------------------------------------------------ATGAAC
                                                                             *   

Test1seq             TCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGT
AT3G20900.1-SEQ      GCGTTGACGAAGAAGAGCTCTATCTCCGGCGGAGCAAAGAAAACGATCTGTCTCCGTCGT
AT3G20900.1-CDS      AAAGTAGCGAGGAAGAACAAAACATC------AGCAAAGAAAACGATCTGTCTCCGTCGT
                         *  *** ***** *   *  **      ****************************

Test1seq             AACACACGGTCGCTAGAGAAACTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
AT3G20900.1-SEQ      AACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
AT3G20900.1-CDS      AACACACAGTTTTTCGAGACCCTTTGCTTCTTCGGCGCCGGTGGACACGTCAGCATCTCC
                     ******* **   * ****  ***************************************

Test1seq             GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCGTGGTGACGTCAGCACCGCT
AT3G20900.1-SEQ      GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCT
AT3G20900.1-CDS      GGTATCCTAGACTTCTTGGCTTTCGGGGTACAACAACCGCCTGGTGACGTCAGCACCGCT
                     **************************************** *******************

Test1seq             GCTGGGGATGGAGAGGGAACAGAGTT-
AT3G20900.1-SEQ      GCTGGGGATGGAGAGGGAACAGAGTAG
AT3G20900.1-CDS      GCTGGGGATGGAGAGGGAACAGAGTAG
                     *************************  
"""
    alignments = list(ClustalIterator(StringIO(aln_example3)))
    assert 1 == len(alignments)
    assert alignments[0]._version == "2.0.9"

    print "The End"
