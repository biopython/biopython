# Copyright 2006, 2007 by Peter Cock.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Bio.SeqIO support for the "phylip" (PHYLIP) file format.

You were expected to use this module via the Bio.SeqIO functions.
This module has now been replaced by Bio.AlignIO.PhylipIO, and is
deprecated."""

import warnings
warnings.warn("Bio.SeqIO.PhylipIO is deprecated.  You can continue to read" \
              + " and write 'clustal' files with Bio.SeqIO, but this is now" \
              + " handled via Bio.AlignIO internally.",
              DeprecationWarning)

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Interfaces import SequenceWriter
from sets import Set

#This is a generator function!
#TODO - Should the default be Gapped(single_letter_alphabet) instead?
def PhylipIterator(handle, alphabet = single_letter_alphabet) :
    """Reads a Phylip alignment file returning a SeqRecord object iterator.

    Record identifiers are limited to at most 10 characters.

    It only copes with interlaced phylip files!  Sequential files won't work
    where the sequences are split over multiple lines.

    For more information on the file format, please see:
    http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles
    """
    line = handle.readline()
    if not line: return
    line = line.strip()
    parts = filter(None, line.split())
    if len(parts)<>2 :
        raise ValueError("First line should have two integers")
    try :
        number_of_seqs = int(parts[0])
        length_of_seqs = int(parts[1])
    except ValueError:
        raise ValueError("First line should have two integers")

    ids = []
    seqs = []

    #Expects STRICT truncation/padding to 10 characters
    #Does not require any white space between name and seq.
    for i in range(0,number_of_seqs) :
        line = handle.readline().rstrip()
        ids.append(line[:10].strip()) #first ten characters
        seqs.append([line[10:].strip().replace(" ","")])

    line=""
    while True :
        #Skip any blank lines between blocks...
        while ""==line.strip():
            line = handle.readline()
            if not line : break #end of file
        if not line : break
        #print "New block..."
        for i in range(0,number_of_seqs) :
            seqs[i].append(line.strip().replace(" ",""))
            line = handle.readline()
            if (not line) and i+1 < number_of_seqs :
                raise ValueError("End of file mid-block")
        if not line : break #end of file

    for i in range(0,number_of_seqs) :
        seq = "".join(seqs[i])
        if len(seq)<>length_of_seqs :
            raise ValueError("Sequence %i length %i, expected length %i" \
                              % (i+1, len(seq), length_of_seqs))
        yield SeqRecord(Seq(seq, alphabet), id=ids[i], name=ids[i], description="")

class PhylipWriter(SequenceWriter):
    """Write interlaced Phylip sequence alignments.

    For more information on the file format, please see:
    http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    http://evolution.genetics.washington.edu/phylip/doc/main.html#inputfiles

    All sequences must be the same length."""
    def __init__(self, handle, truncate=10):
        """Creates the writer object

        Use the method write_file() to actually record your sequence records."""
        self.handle = handle
        self.truncate = truncate
    
    def write_file(self, records) :
        """Use this to write an entire file containing the given records.

        If records is an iterator that does not support len(records) or
        records[index] then it is converted into a list.
        """
        #Need length, and multiple passes - and iterator will not do.
        records = list(records)

        if len(records)==0 :
            raise ValueError("Must have at least one sequence")
        length_of_sequences = len(records[0].seq)
        for record in records :
            if length_of_sequences <> len(record.seq) :
                raise ValueError("Sequences must all be the same length")
        if length_of_sequences <= 0 :
            raise ValueError("Non-empty sequences are required")
        
        if len(records) > len(Set([r.id[:self.truncate] for r in records])) :
            raise ValueError("Repeated identifier, possibly due to truncation")

        handle = self.handle

        # From experimentation, the use of tabs is not understood by the
        # EMBOSS suite.  The nature of the expected white space is not
        # defined, simply "These are in free format, separated by blanks"
        handle.write(" %i %s\n" % (len(records), length_of_sequences))
        block=0
        while True :
            for record in records :
                if block==0 :
                    #Write name (truncated/padded to 10 characters)
                    """
                    Quoting the PHYLIP version 3.6 documentation:
                    
                    The name should be ten characters in length, filled out to
                    the full ten characters by blanks if shorter. Any printable
                    ASCII/ISO character is allowed in the name, except for
                    parentheses ("(" and ")"), square brackets ("[" and "]"),
                    colon (":"), semicolon (";") and comma (","). If you forget
                    to extend the names to ten characters in length by blanks,
                    the program [i.e. PHYLIP] will get out of synchronization
                    with the contents of the data file, and an error message
                    will result.

                    Note that Tab characters count as only one character in the
                    species names. Their inclusion can cause trouble.
                    """
                    name = record.id.strip()
                    #Either remove the banned characters, or map them to something
                    #else like an underscore "_" or pipe "|" character...
                    for char in "[]()," :
                        name = name.replace(char,"")
                    for char in ":;" :
                        name = name.replace(char,"|")

                    #Now truncate and right pad to expected length.
                    handle.write(name[:self.truncate].ljust(self.truncate))
                else :
                    #write 10 space indent
                    handle.write(" "*self.truncate)
                #Write five chunks of ten letters per line...
                for chunk in range(0,5) :
                    i = block*50 + chunk*10
                    seq_segment = record.seq.tostring()[i:i+10]
                    #TODO - Force any gaps to be '-' character?  Look at the alphabet...
                    #TODO - How to cope with '?' or '.' in the sequence?
                    handle.write(" %s" % seq_segment)
                    if i+10 > length_of_sequences : break
                handle.write("\n")
            block=block+1
            if block*50 > length_of_sequences : break
            handle.write("\n")

        #Don't close the handle.  Doing so would prevent this code
        #from writing concatenated phylip files which are used
        #in phylogenetic bootstrapping
        #handle.close()
if __name__=="__main__" :
    print "Testing"

    phylip_text="""     8    286
V_Harveyi_ --MKNWIKVA VAAIA--LSA A--------- ---------T VQAATEVKVG 
B_subtilis MKMKKWTVLV VAALLAVLSA CG-------- ----NGNSSS KEDDNVLHVG 
B_subtilis MKKALLALFM VVSIAALAAC GAGNDNQSKD NAKDGDLWAS IKKKGVLTVG 
YA80_HAEIN MKKLLFTTAL LTGAIAFSTF ---------- -SHAGEIADR VEKTKTLLVG 
FLIY_ECOLI MKLAHLGRQA LMGVMAVALV AG---MSVKS FADEG-LLNK VKERGTLLVG 
E_coli_Gln --MKSVLKVS LAALTLAFAV S--------- ---------S HAADKKLVVA 
Deinococcu -MKKSLLSLK LSGLLVPSVL ALS------- -LSACSSPSS TLNQGTLKIA 
HISJ_E_COL MKKLVLSLSL VLAFSSATAA F--------- ---------- AAIPQNIRIG 

           MSGRYFPFTF VKQ--DKLQG FEVDMWDEIG KRNDYKIEYV TANFSGLFGL 
           ATGQSYPFAY KEN--GKLTG FDVEVMEAVA KKIDMKLDWK LLEFSGLMGE 
           TEGTYEPFTY HDKDTDKLTG YDVEVITEVA KRLGLKVDFK ETQWGSMFAG 
           TEGTYAPFTF HDK-SGKLTG FDVEVIRKVA EKLGLKVEFK ETQWDAMYAG 
           LEGTYPPFSF QGD-DGKLTG FEVEFAQQLA KHLGVEASLK PTKWDGMLAS 
           TDTAFVPFEF KQG--DKYVG FDVDLWAAIA KELKLDYELK PMDFSGIIPA 
           MEGTYPPFTS KNE-QGELVG FDVDIAKAVA QKLNLKPEFV LTEWSGILAG 
           TDPTYAPFES KNS-QGELVG FDIDLAKELC KRINTQCTFV ENPLDALIPS 

           LETGRIDTIS NQITMTDARK AKYLFADPYV VDG-AQITVR KGNDSIQGVE 
           LQTGKLDTIS NQVAVTDERK ETYNFTKPYA YAG-TQIVVK KDNTDIKSVD 
           LNSKRFDVVA NQVG-KTDRE DKYDFSDKYT TSR-AVVVTK KDNNDIKSEA 
           LNAKRFDVIA NQTNPSPERL KKYSFTTPYN YSG-GVIVTK SSDNSIKSFE 
           LDSKRIDVVI NQVTISDERK KKYDFSTPYT ISGIQALVKK GNEGTIKTAD 
           LQTKNVDLAL AGITITDERK KAIDFSDGYY KSG-LLVMVK ANNNDVKSVK 
           LQANKYDVIV NQVGITPERQ NSIGFSQPYA YSRPEIIVAK NNTFNPQSLA 
           LKAKKIDAIM SSLSITEKRQ QEIAFTDKLY AADSRLVVAK NSDIQP-TVE 

           DLAGKTVAVN LGSNFEQLLR DYDKDGKINI KTYDT--GIE HDVALGRADA 
           DLKGKTVAAV LGSNHAKNLE SKDPDKKINI KTYETQEGTL KDVAYGRVDA 
           DVKGKTSAQS LTSNYNKLAT N----AGAKV EGVEGMAQAL QMIQQARVDM 
           DLKGRKSAQS ATSNWGKDAK A----AGAQI LVVDGLAQSL ELIKQGRAEA 
           DLKGKKVGVG LGTNYEEWLR QNV--QGVDV RTYDDDPTKY QDLRVGRIDA 
           DLDGKVVAVK SGTGSVDYAK AN--IKTKDL RQFPNIDNAY MELGTNRADA 
           DLKGKRVGST LGSNYEKQLI DTG---DIKI VTYPGAPEIL ADLVAGRIDA 
           SLKGKRVGVL QGTTQETFGN EHWAPKGIEI VSYQGQDNIY SDLTAGRIDA 

           FIMDRLSALE -LIKKT-GLP LQLAGEPFET I-----QNAW PFVDNEKGRK 
           YVNSRTVLIA -QIKKT-GLP LKLAGDPIVY E-----QVAF PFAKDDAHDK 
           TYNDKLAVLN -YLKTSGNKN VKIAFETGEP Q-----STYF TFRKGS--GE 
           TINDKLAVLD -YFKQHPNSG LKIAYDRGDK T-----PTAF AFLQGE--DA 
           ILVDRLAALD -LVKKT-NDT LAVTGEAFSR Q-----ESGV ALRKGN--ED 
           VLHDTPNILY -FIKTAGNGQ FKAVGDSLEA Q-----QYGI AFPKGS--DE 
           AYNDRLVVNY -IINDQ-KLP VRGAGQIGDA A-----PVGI ALKKGN--SA 
           AFQDEVAASE GFLKQPVGKD YKFGGPSVKD EKLFGVGTGM GLRKED--NE 

           LQAEVNKALA EMRADGTVEK ISVKWFGADI TK----
           LRKKVNKALD ELRKDGTLKK LSEKYFNEDI TVEQKH
           VVDQVNKALK EMKEDGTLSK ISKKWFGEDV SK----
           LITKFNQVLE ALRQDGTLKQ ISIEWFGYDI TQ----
           LLKAVNDAIA EMQKDGTLQA LSEKWFGADV TK----
           LRDKVNGALK TLRENGTYNE IYKKWFGTEP K-----
           LKDQIDKALT EMRSDGTFEK ISQKWFGQDV GQP---
           LREALNKAFA EMRADGTYEK LAKKYFDFDV YGG---
"""

    from cStringIO import StringIO
    handle = StringIO(phylip_text)
    count=0
    for record in PhylipIterator(handle) :
        count=count+1
        print record.id
        #print record.seq.tostring()
    assert count == 8

    expected="""mkklvlslsl vlafssataa faaipqniri gtdptyapfe sknsqgelvg
    fdidlakelc krintqctfv enpldalips lkakkidaim sslsitekrq qeiaftdkly
    aadsrlvvak nsdiqptves lkgkrvgvlq gttqetfgne hwapkgieiv syqgqdniys
    dltagridaafqdevaaseg flkqpvgkdy kfggpsvkde klfgvgtgmg lrkednelre
    alnkafaemradgtyeklak kyfdfdvygg""".replace(" ","").replace("\n","").upper()
    assert record.seq.tostring().replace("-","") == expected

    #From here:
    #http://atgc.lirmm.fr/phyml/usersguide.html
    phylip_text2="""5 60
Tax1        CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAG
Tax2        CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGG
Tax3        CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGG
Tax4        TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGG
Tax5        CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGG

GAAATGGTCAATATTACAAGGT
GAAATGGTCAACATTAAAAGAT
GAAATCGTCAATATTAAAAGGT
GAAATGGTCAATCTTAAAAGGT
GAAATGGTCAATATTAAAAGGT"""

    phylip_text3="""5 60
Tax1        CCATCTCACGGTCGGTACGATACACCTGCTTTTGGCAGGAAATGGTCAATATTACAAGGT
Tax2        CCATCTCACGGTCAGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAACATTAAAAGAT
Tax3        CCATCTCCCGCTCAGTAAGATACCCCTGCTGTTGGCGGGAAATCGTCAATATTAAAAGGT
Tax4        TCATCTCATGGTCAATAAGATACTCCTGCTTTTGGCGGGAAATGGTCAATCTTAAAAGGT
Tax5        CCATCTCACGGTCGGTAAGATACACCTGCTTTTGGCGGGAAATGGTCAATATTAAAAGGT"""

    handle = StringIO(phylip_text2)
    list2 = list(PhylipIterator(handle))
    handle.close()
    assert len(list2)==5

    handle = StringIO(phylip_text3)
    list3 = list(PhylipIterator(handle))
    handle.close()
    assert len(list3)==5

    for i in range(0,5) :
        list2[i].id == list3[i].id
        list2[i].seq.tostring() == list3[i].seq.tostring()

    #From here:
    #http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    #Note the lack of any white space between names 2 and 3 and their seqs.
    phylip_text4="""  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT
Salmo gairAAGCCTTGGC AGTGCAGGGT
H. SapiensACCGGTTGGC CGTTCAGGGT
Chimp     AAACCCTTGC CGTTACGCTT
Gorilla   AAACCCTTGC CGGTACGCTT

GAGCCCGGGC AATACAGGGT AT
GAGCCGTGGC CGGGCACGGT AT
ACAGGTTGGC CGTTCAGGGT AA
AAACCGAGGC CGGGACACTC AT
AAACCATTGC CGGTACGCTT AA"""

    #From here:
    #http://evolution.genetics.washington.edu/phylip/doc/sequence.html
    phylip_text5="""  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT
GAGCCCGGGC AATACAGGGT AT
Salmo gairAAGCCTTGGC AGTGCAGGGT
GAGCCGTGGC CGGGCACGGT AT
H. SapiensACCGGTTGGC CGTTCAGGGT
ACAGGTTGGC CGTTCAGGGT AA
Chimp     AAACCCTTGC CGTTACGCTT
AAACCGAGGC CGGGACACTC AT
Gorilla   AAACCCTTGC CGGTACGCTT
AAACCATTGC CGGTACGCTT AA"""

    phylip_text5a="""  5    42
Turkey    AAGCTNGGGC ATTTCAGGGT GAGCCCGGGC AATACAGGGT AT
Salmo gairAAGCCTTGGC AGTGCAGGGT GAGCCGTGGC CGGGCACGGT AT
H. SapiensACCGGTTGGC CGTTCAGGGT ACAGGTTGGC CGTTCAGGGT AA
Chimp     AAACCCTTGC CGTTACGCTT AAACCGAGGC CGGGACACTC AT
Gorilla   AAACCCTTGC CGGTACGCTT AAACCATTGC CGGTACGCTT AA"""

    handle = StringIO(phylip_text4)
    list4 = list(PhylipIterator(handle))
    handle.close()
    assert len(list4)==5

    handle = StringIO(phylip_text5)
    try :
        list5 = list(PhylipIterator(handle))
        assert len(list5)==5
        print "That should have failed..."
    except ValueError :
        print "Evil multiline non-interlaced example failed as expected"
    handle.close()

    handle = StringIO(phylip_text5a)
    list5 = list(PhylipIterator(handle))
    handle.close()
    assert len(list5)==5

    for i in range(0,5) :
        list4[i].id == list5[i].id
        list4[i].seq.tostring() == list5[i].seq.tostring()

    

    """
    handle = StringIO(phylip_text)
    out_handle=open("/tmp/test.phy","w")
    writer = PhylipWriter(out_handle)
    writer.write_file(PhylipIterator(handle))
    out_handle.close()

    print "---------------------"

    print open("/tmp/test.phy").read()
    """
