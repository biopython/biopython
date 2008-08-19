# Copyright 2008 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Multiple sequence alignment input/output as Alignment objects.

The Bio.AlignIO interface is deliberately very similar to Bio.SeqIO, and in
fact the two are connected internally.  Both modules use the same set of file
format names (lower case strings).  From the user's perspective, you can read
in a PHYLIP file containing one or more alignments using Bio.AlignIO, or you
can read in the sequences within these alignmenta using Bio.SeqIO.

Input
=====
For the typical special case when your file or handle contains one and only
one alignment, use the function Bio.AlignIO.read().  This takes an input file
handle, format string and optional number of sequences per alignment.  It will
return a single Alignment object (or raise an exception if there isn't just
one alignment):

    from Bio import AlignIO
    handle = open("example.aln", "rU")
    align = AlignIO.read(handle, "clustal")
    handle.close()
    print align

For the general case, when the handle could contain any number of alignments,
use the function Bio.AlignIO.parse(...) which takes the same arguments, but
returns an iterator giving Alignment objects.  For example, using the output
from the EMBOSS water or needle pairwise alignment prorams:

    from Bio import AlignIO
    handle = open("example.txt", "rU")
    for alignment in AlignIO.parse(handle, "emboss") :
        print alignment

If you want random access to the alignments by number, turn this into a list:

    from Bio import AlignIO
    handle = open("example.aln", "rU")
    alignments = list(AlignIO.parse(handle, "clustal"))
    print alignments[0]

Most alignment file formats can be concatenated so as to hold as many
different multiple sequence alignments as possible.  One common example
is the output of the tool seqboot in the PHLYIP suite.  Sometimes there
can be a file header and footer, as seen in the EMBOSS alignment output.

There is an optional argument for the number of sequences per alignment which
is usually only needed with the alignments stored in the FASTA format.
Without this information, there is no clear way to tell if you have say a
single alignment of 20 sequences, or four alignments of 5 sequences.  e.g.

    from Bio import AlignIO
    handle = open("example.faa", "rU")
    for alignment in AlignIO.parse(handle, "fasta", seq_count=5) :
        print alignment

The above code would split up the FASTA files, and try and batch every five
sequences into an alignment.

Output
======
Use the function Bio.AlignIO.write(...), which takes a complete set of
Alignment objects (either as a list, or an iterator), an output file handle
and of course the file format.

    from Bio import AlignIO
    alignments = ...
    handle = open("example.faa", "w")
    alignment = SeqIO.write(alignments, handle, "fasta")
    handle.close()

In general, you are expected to call this function once (with all your
alignments) and then close the file handle.  However, for file formats
like PHYLIP where multiple alignments are stored sequentially (with no file
header and footer), then multiple calls to the write function should work as
expected.

File Formats
============
When specifying the file format, use lowercase strings.  The same format
names are also used in Bio.SeqIO and include the following:

clustal   - Ouput from Clustal W or X, see also the module Bio.Clustalw
            which can be used to run the command line tool from Biopython.
emboss    - The "pairs" and "simple" alignment format from the EMBOSS tools.
fasta     - The generic sequence file format where each record starts with a
            identifer line starting with a ">" character, followed by lines
            of sequence.
fasta-m10 - For the pairswise alignments output by Bill Pearson's FASTA
            tools when used with the -m 10 command line option for machine
            readable output.
ig        - The IntelliGenetics file format, apparently the same as the
            MASE alignment format.
nexus     - Output from NEXUS, see also the module Bio.Nexus which can also
            read any phylogenetic trees in these files.
phylip    - Used by the PHLIP tools.
stockholm - A richly annotated alignment file format used by PFAM.

Note that while Bio.AlignIO can read all the above file formats, it cannot
write to all of them.

You can also use any file format supported by Bio.SeqIO, such as "fasta" or
"ig" (which are listed above), PROVIDED the sequences in your file are all the
same length.

Further Information
===================
See the wiki page http://biopython.org/wiki/AlignIO and also the Bio.AlignIO
chapter in the Biopython Tutorial and Cookbook which is also available online:

http://biopython.org/DIST/docs/tutorial/Tutorial.html
http://biopython.org/DIST/docs/tutorial/Tutorial.pdf
"""

#TODO
# - define policy on reading aligned sequences with gaps in
#   (e.g. - and . characters) including how the alphabet interacts
#
# - Can we build the to_alignment(...) functionality
#   into the generic Alignment class instead?
#
# - How best to handle unique/non unique record.id when writing.
#   For most file formats reading such files is fine; The stockholm
#   parser would fail.
#
# - MSF multiple alignment format, aka GCG, aka PileUp format (*.msf)
#   http://www.bioperl.org/wiki/MSF_multiple_alignment_format 

import os
#from cStringIO import StringIO
from StringIO import StringIO
from Bio.Alphabet import generic_alphabet, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Generic import Alignment

import StockholmIO
import ClustalIO
import NexusIO
import PhylipIO
import EmbossIO
import FastaIO

#Convention for format names is "mainname-subtype" in lower case.
#Please use the same names as BioPerl and EMBOSS where possible.

_FormatToIterator ={#"fasta" is done via Bio.SeqIO
                    "clustal" : ClustalIO.ClustalIterator,
                    "emboss" : EmbossIO.EmbossIterator,
                    "fasta-m10" : FastaIO.FastaM10Iterator,
                    "nexus" : NexusIO.NexusIterator,
                    "phylip" : PhylipIO.PhylipIterator,
                    "stockholm" : StockholmIO.StockholmIterator,
                    }

_FormatToWriter ={#"fasta" is done via Bio.SeqIO
                  #"emboss" : EmbossIO.EmbossWriter, (unfinished)
                  "nexus" : NexusIO.NexusWriter,
                  "phylip" : PhylipIO.PhylipWriter,
                  "stockholm" : StockholmIO.StockholmWriter,
                  "clustal" : ClustalIO.ClustalWriter,
                  }

#This is a generator function
def _multiple_alignment_record_iterator(alignments) :
    """Private function to loop over all the records in many alignments."""
    for alignment in alignments :
        for record in alignment :
            yield record

def write(alignments, handle, format) :
    """Write complete set of alignments to a file.

    sequences - A list (or iterator) of Alignment objects
    handle    - File handle object to write to
    format    - lower case string describing the file format to write.

    You should close the handle after calling this function.

    There is no return value.
    """
    from Bio import SeqIO

    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError("Need a file handle, not a string (i.e. not a filename)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format <> format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)
    if isinstance(alignments, Alignment) :
        raise TypeError("Need an Alignment list/iterator, not just a single Alignment")

    #Map the file format to a writer class
    if format in _FormatToIterator :
        writer_class = _FormatToWriter[format]
        writer_class(handle).write_file(alignments)
    elif format in SeqIO._FormatToWriter :
        #Exploit the existing SeqIO parser to the dirty work!
        #TODO - Once we drop support for Python 2.3, this helper function can be
        #replaced with a generator expression.
        SeqIO.write(_multiple_alignment_record_iterator(alignments), handle, format)
    elif format in _FormatToIterator or format in SeqIO._FormatToIterator :
        raise ValueError("Reading format '%s' is supported, but not writing" % format)
    else :
        raise ValueError("Unknown format '%s'" % format)

    return

#This is a generator function!
def _SeqIO_to_alignment_iterator(handle, format, seq_count=None) :
    """Private function, uses Bio.SeqIO to create an Alignment iterator.

    handle   - handle to the file.
    format   - string describing the file format.
    seq_count- Optional integer, number of sequences expected in
               each alignment.  Recommended for fasta format files.

    If count is omitted (default) then all the sequences in
    the file are combined into a single Alignment.
    """
    from Bio import SeqIO

    assert format in SeqIO._FormatToIterator

    if seq_count :
        #Use the count to split the records into batches.
        seq_record_iterator = SeqIO.parse(handle, format)

        records = []
        for record in seq_record_iterator :
            records.append(record)
            if len(records) == seq_count :
                yield SeqIO.to_alignment(records)
                records = []
        if len(records) > 0 :
            raise ValueError("Check seq_count argument, not enough sequences?")
    else :
        #Must assume that there is a single alignment using all
        #the SeqRecord objects:
        records = list(SeqIO.parse(handle, format))
        if records :
            yield SeqIO.to_alignment(records)
        else :
            #No alignment found!
            pass
    
def parse(handle, format, seq_count=None) :
    """Turns a sequence file into an iterator returning Alignment objects.

    handle   - handle to the file.
    format   - string describing the file format.
    seq_count- Optional integer, number of sequences expected in
               each alignment.  Recommended for fasta format files.

    If you have the file name in a string 'filename', use:

    from Bio import AlignIO
    my_iterator = AlignIO.parse(open(filename,"rU"), format)

    If you have a string 'data' containing the file contents, use:

    from Bio import AlignIO
    from StringIO import StringIO
    my_iterator = AlignIO.parse(StringIO(data), format)

    Use the Bio.AlignIO.read(handle, format[, seq_count]) function when
    you expect a single record only.
    """
    from Bio import SeqIO

    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError("Need a file handle, not a string (i.e. not a filename)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format <> format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)

    #Map the file format to a sequence iterator:
    if format in _FormatToIterator :
        iterator_generator = _FormatToIterator[format]
        return iterator_generator(handle, seq_count)
    elif format in SeqIO._FormatToIterator :
        #Exploit the existing SeqIO parser to the dirty work!
        return _SeqIO_to_alignment_iterator(handle, format, seq_count)
    else :
        raise ValueError("Unknown format '%s'" % format)

def read(handle, format, seq_count=None) :
    """Turns an alignment file into a single Alignment object.

    handle   - handle to the file.
    format   - string describing the file format.
    seq_count- Optional interger, number of sequences expected in
               the alignment to check you got what you expected.

    If the handle contains no alignments, or more than one alignment,
    an exception is raised.  For example, using a PFAM/Stockholm file
    containing one alignment:

    from Bio import AlignIO
    align = AlignIO.read(open("example.sth"), "stockholm")

    If however you want the first alignment from a file containing
    multiple alignments this function would raise an exception.
    Instead use:

    from Bio import AlignIO
    align = AlignIO.parse(open("example.sth"), "stockholm").next()

    Use the Bio.AlignIO.parse() function if you want to read multiple
    records from the handle.
    """
    iterator = parse(handle, format, seq_count)
    try :
        first = iterator.next()
    except StopIteration :
        first = None
    if first is None :
        raise ValueError, "No records found in handle"
    try :
        second = iterator.next()
    except StopIteration :
        second = None
    if second is not None :
        raise ValueError, "More than one record found in handle"
    if seq_count :
        assert len(first.get_all_seqs())==seq_count
    return first

           
if __name__ == "__main__" :
    #Run some tests...
    from Bio.Alphabet import generic_nucleotide
    from sets import Set

    for format in _FormatToIterator :
        print "parse(handle to empty file)"
        iterator = parse(StringIO(""), format=format)
        assert len(list(iterator))==0
        iterator = parse(StringIO(""), format=format, seq_count = 42)
        assert len(list(iterator))==0
    print

    def align_cmp(align1, align2) :
        if align1.get_alignment_length() <> align2.get_alignment_length() :
            return False
        recs1 = align1.get_all_seqs()
        recs2 = align2.get_all_seqs()
        if len(recs1) <> len(recs2) :
            return False
        for r1, r2 in zip(recs1, recs2) :
            if r1.seq.tostring() <> r2.seq.tostring() :
                return False
        return True
    
    # Fasta file with unusual layout, from here:
    # http://virgil.ruc.dk/kurser/Sekvens/Treedraw.htm
    faa_example = \
""">V_Harveyi_PATH
mknwikvava aialsaatvq aatevkvgms gryfpftfvk qdklqgfevd mwdeigkrnd
ykieyvtanf sglfglletg ridtisnqit mtdarkakyl fadpyvvdga qitvrkgnds
iqgvedlagk tvavnlgsnf eqllrdydkd gkiniktydt giehdvalgr adafimdrls
alelikktgl plqlagepfe tiqnawpfvd nekgrklqae vnkalaemra dgtvekisvk
wfgaditk
>B_subtilis_YXEM
mkmkkwtvlv vaallavlsa cgngnssske ddnvlhvgat gqsypfayke ngkltgfdve
vmeavakkid mkldwkllef sglmgelqtg kldtisnqva vtderketyn ftkpyayagt
qivvkkdntd iksvddlkgk tvaavlgsnh aknleskdpd kkiniktyet qegtlkdvay
grvdayvnsr tvliaqikkt glplklagdp ivyeqvafpf akddahdklr kkvnkaldel
rkdgtlkkls ekyfneditv eqkh
>FLIY_ECOLI
mklahlgrqa lmgvmavalv agmsvksfad egllnkvker gtllvglegt yppfsfqgdd
gkltgfevef aqqlakhlgv easlkptkwd gmlasldskr idvvinqvti sderkkkydf
stpytisgiq alvkkgnegt iktaddlkgk kvgvglgtny eewlrqnvqg vdvrtydddp
tkyqdlrvgr idailvdrla aldlvkktnd tlavtgeafs rqesgvalrk gnedllkavn
daiaemqkdg tlqalsekwf gadvtk
>Deinococcus_radiodurans
mkksllslkl sgllvpsvla lslsacssps stlnqgtlki amegtyppft skneqgelvg
fdvdiakava qklnlkpefv ltewsgilag lqankydviv nqvgitperq nsigfsqpya
ysrpeiivak nntfnpqsla dlkgkrvgst lgsnyekqli dtgdikivty pgapeiladl
vagridaayn drlvvnyiin dqklpvrgag qigdaapvgi alkkgnsalk dqidkaltem
rsdgtfekis qkwfgqdvgq p
>B_subtilis_GlnH_homo_YCKK
mkkallalfm vvsiaalaac gagndnqskd nakdgdlwas ikkkgvltvg tegtyepfty
hdkdtdkltg ydveviteva krlglkvdfk etqwgsmfag lnskrfdvva nqvgktdred
kydfsdkytt sravvvtkkd nndikseadv kgktsaqslt snynklatna gakvegvegm
aqalqmiqqa rvdmtyndkl avlnylktsg nknvkiafet gepqstyftf rkgsgevvdq
vnkalkemke dgtlskiskk wfgedvsk
>YA80_HAEIN
mkkllfttal ltgaiafstf shageiadrv ektktllvgt egtyapftfh dksgkltgfd
vevirkvaek lglkvefket qwdamyagln akrfdvianq tnpsperlkk ysfttpynys
ggvivtkssd nsiksfedlk grksaqsats nwgkdakaag aqilvvdgla qslelikqgr
aeatindkla vldyfkqhpn sglkiaydrg dktptafafl qgedalitkf nqvlealrqd
gtlkqisiew fgyditq
>E_coli_GlnH
mksvlkvsla altlafavss haadkklvva tdtafvpfef kqgdkyvgfd vdlwaaiake
lkldyelkpm dfsgiipalq tknvdlalag ititderkka idfsdgyyks gllvmvkann
ndvksvkdld gkvvavksgt gsvdyakani ktkdlrqfpn idnaymelgt nradavlhdt
pnilyfikta gngqfkavgd sleaqqygia fpkgsdelrd kvngalktlr engtyneiyk
kwfgtepk
>HISJ_E_COLI
mkklvlslsl vlafssataa faaipqniri gtdptyapfe sknsqgelvg fdidlakelc
krintqctfv enpldalips lkakkidaim sslsitekrq qeiaftdkly aadsrlvvak
nsdiqptves lkgkrvgvlq gttqetfgne hwapkgieiv syqgqdniys dltagridaa
fqdevaaseg flkqpvgkdy kfggpsvkde klfgvgtgmg lrkednelre alnkafaemr
adgtyeklak kyfdfdvygg"""

    # This alignment was created from the fasta example given above
    aln_example = \
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

V_Harveyi_PATH                 LETGRIDTISNQITMTDARKAKYLFADPYVVDG-AQITVRKGNDSIQGVE
B_subtilis_YXEM                LQTGKLDTISNQVAVTDERKETYNFTKPYAYAG-TQIVVKKDNTDIKSVD
B_subtilis_GlnH_homo_YCKK      LNSKRFDVVANQVG-KTDREDKYDFSDKYTTSR-AVVVTKKDNNDIKSEA
YA80_HAEIN                     LNAKRFDVIANQTNPSPERLKKYSFTTPYNYSG-GVIVTKSSDNSIKSFE
FLIY_ECOLI                     LDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTAD
E_coli_GlnH                    LQTKNVDLALAGITITDERKKAIDFSDGYYKSG-LLVMVKANNNDVKSVK
Deinococcus_radiodurans        LQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLA
HISJ_E_COLI                    LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQP-TVE
                               *.: . *        .  *     *:          :  : .        

V_Harveyi_PATH                 DLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDT--GIEHDVALGRADA
B_subtilis_YXEM                DLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEGTLKDVAYGRVDA
B_subtilis_GlnH_homo_YCKK      DVKGKTSAQSLTSNYNKLATN----AGAKVEGVEGMAQALQMIQQARVDM
YA80_HAEIN                     DLKGRKSAQSATSNWGKDAKA----AGAQILVVDGLAQSLELIKQGRAEA
FLIY_ECOLI                     DLKGKKVGVGLGTNYEEWLRQNV--QGVDVRTYDDDPTKYQDLRVGRIDA
E_coli_GlnH                    DLDGKVVAVKSGTGSVDYAKAN--IKTKDLRQFPNIDNAYMELGTNRADA
Deinococcus_radiodurans        DLKGKRVGSTLGSNYEKQLIDTG---DIKIVTYPGAPEILADLVAGRIDA
HISJ_E_COLI                    SLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDA
                               .: *:  .    :               .:            :   * : 

V_Harveyi_PATH                 FIMDRLSALE-LIKKT-GLPLQLAGEPFETI-----QNAWPFVDNEKGRK
B_subtilis_YXEM                YVNSRTVLIA-QIKKT-GLPLKLAGDPIVYE-----QVAFPFAKDDAHDK
B_subtilis_GlnH_homo_YCKK      TYNDKLAVLN-YLKTSGNKNVKIAFETGEPQ-----STYFTFRKGS--GE
YA80_HAEIN                     TINDKLAVLD-YFKQHPNSGLKIAYDRGDKT-----PTAFAFLQGE--DA
FLIY_ECOLI                     ILVDRLAALD-LVKKT-NDTLAVTGEAFSRQ-----ESGVALRKGN--ED
E_coli_GlnH                    VLHDTPNILY-FIKTAGNGQFKAVGDSLEAQ-----QYGIAFPKGS--DE
Deinococcus_radiodurans        AYNDRLVVNY-IINDQ-KLPVRGAGQIGDAA-----PVGIALKKGN--SA
HISJ_E_COLI                    AFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKED--NE
                                  .        .:                           : . .    

V_Harveyi_PATH                 LQAEVNKALAEMRADGTVEKISVKWFGADITK----
B_subtilis_YXEM                LRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH
B_subtilis_GlnH_homo_YCKK      VVDQVNKALKEMKEDGTLSKISKKWFGEDVSK----
YA80_HAEIN                     LITKFNQVLEALRQDGTLKQISIEWFGYDITQ----
FLIY_ECOLI                     LLKAVNDAIAEMQKDGTLQALSEKWFGADVTK----
E_coli_GlnH                    LRDKVNGALKTLRENGTYNEIYKKWFGTEPK-----
Deinococcus_radiodurans        LKDQIDKALTEMRSDGTFEKISQKWFGQDVGQP---
HISJ_E_COLI                    LREALNKAFAEMRADGTYEKLAKKYFDFDVYGG---
                               :   .: .:  :: :** . :  ::*. :       
"""

    # This is the clustal example (above) but output in phylip format,
    # with truncated names.  Note there is an ambiguity here: two
    # different sequences both called "B_subtilis", originally
    # "B_subtilis_YXEM" and "B_subtilis_GlnH_homo_YCKK"
    phy_example = \
"""     8    286
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
    # This is the clustal example (above) but output in phylip format,
    nxs_example = \
"""#NEXUS
BEGIN DATA;
dimensions ntax=8 nchar=286;
format missing=?
symbols="ABCDEFGHIKLMNPQRSTUVWXYZ"
interleave datatype=PROTEIN gap= -;

matrix
V_Harveyi_PATH             --MKNWIKVAVAAIA--LSAA------------------TVQAATEVKVG
B_subtilis_YXEM            MKMKKWTVLVVAALLAVLSACG------------NGNSSSKEDDNVLHVG
B_subtilis_GlnH_homo_YCKK  MKKALLALFMVVSIAALAACGAGNDNQSKDNAKDGDLWASIKKKGVLTVG
YA80_HAEIN                 MKKLLFTTALLTGAIAFSTF-----------SHAGEIADRVEKTKTLLVG
FLIY_ECOLI                 MKLAHLGRQALMGVMAVALVAG---MSVKSFADEG-LLNKVKERGTLLVG
E_coli_GlnH                --MKSVLKVSLAALTLAFAVS------------------SHAADKKLVVA
Deinococcus_radiodurans    -MKKSLLSLKLSGLLVPSVLALS--------LSACSSPSSTLNQGTLKIA
HISJ_E_COLI                MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG

V_Harveyi_PATH             MSGRYFPFTFVKQ--DKLQGFEVDMWDEIGKRNDYKIEYVTANFSGLFGL
B_subtilis_YXEM            ATGQSYPFAYKEN--GKLTGFDVEVMEAVAKKIDMKLDWKLLEFSGLMGE
B_subtilis_GlnH_homo_YCKK  TEGTYEPFTYHDKDTDKLTGYDVEVITEVAKRLGLKVDFKETQWGSMFAG
YA80_HAEIN                 TEGTYAPFTFHDK-SGKLTGFDVEVIRKVAEKLGLKVEFKETQWDAMYAG
FLIY_ECOLI                 LEGTYPPFSFQGD-DGKLTGFEVEFAQQLAKHLGVEASLKPTKWDGMLAS
E_coli_GlnH                TDTAFVPFEFKQG--DKYVGFDVDLWAAIAKELKLDYELKPMDFSGIIPA
Deinococcus_radiodurans    MEGTYPPFTSKNE-QGELVGFDVDIAKAVAQKLNLKPEFVLTEWSGILAG
HISJ_E_COLI                TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS

V_Harveyi_PATH             LETGRIDTISNQITMTDARKAKYLFADPYVVDG-AQITVRKGNDSIQGVE
B_subtilis_YXEM            LQTGKLDTISNQVAVTDERKETYNFTKPYAYAG-TQIVVKKDNTDIKSVD
B_subtilis_GlnH_homo_YCKK  LNSKRFDVVANQVG-KTDREDKYDFSDKYTTSR-AVVVTKKDNNDIKSEA
YA80_HAEIN                 LNAKRFDVIANQTNPSPERLKKYSFTTPYNYSG-GVIVTKSSDNSIKSFE
FLIY_ECOLI                 LDSKRIDVVINQVTISDERKKKYDFSTPYTISGIQALVKKGNEGTIKTAD
E_coli_GlnH                LQTKNVDLALAGITITDERKKAIDFSDGYYKSG-LLVMVKANNNDVKSVK
Deinococcus_radiodurans    LQANKYDVIVNQVGITPERQNSIGFSQPYAYSRPEIIVAKNNTFNPQSLA
HISJ_E_COLI                LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQP-TVE

V_Harveyi_PATH             DLAGKTVAVNLGSNFEQLLRDYDKDGKINIKTYDT--GIEHDVALGRADA
B_subtilis_YXEM            DLKGKTVAAVLGSNHAKNLESKDPDKKINIKTYETQEGTLKDVAYGRVDA
B_subtilis_GlnH_homo_YCKK  DVKGKTSAQSLTSNYNKLATN----AGAKVEGVEGMAQALQMIQQARVDM
YA80_HAEIN                 DLKGRKSAQSATSNWGKDAKA----AGAQILVVDGLAQSLELIKQGRAEA
FLIY_ECOLI                 DLKGKKVGVGLGTNYEEWLRQNV--QGVDVRTYDDDPTKYQDLRVGRIDA
E_coli_GlnH                DLDGKVVAVKSGTGSVDYAKAN--IKTKDLRQFPNIDNAYMELGTNRADA
Deinococcus_radiodurans    DLKGKRVGSTLGSNYEKQLIDTG---DIKIVTYPGAPEILADLVAGRIDA
HISJ_E_COLI                SLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDA

V_Harveyi_PATH             FIMDRLSALE-LIKKT-GLPLQLAGEPFETI-----QNAWPFVDNEKGRK
B_subtilis_YXEM            YVNSRTVLIA-QIKKT-GLPLKLAGDPIVYE-----QVAFPFAKDDAHDK
B_subtilis_GlnH_homo_YCKK  TYNDKLAVLN-YLKTSGNKNVKIAFETGEPQ-----STYFTFRKGS--GE
YA80_HAEIN                 TINDKLAVLD-YFKQHPNSGLKIAYDRGDKT-----PTAFAFLQGE--DA
FLIY_ECOLI                 ILVDRLAALD-LVKKT-NDTLAVTGEAFSRQ-----ESGVALRKGN--ED
E_coli_GlnH                VLHDTPNILY-FIKTAGNGQFKAVGDSLEAQ-----QYGIAFPKGS--DE
Deinococcus_radiodurans    AYNDRLVVNY-IINDQ-KLPVRGAGQIGDAA-----PVGIALKKGN--SA
HISJ_E_COLI                AFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKED--NE

V_Harveyi_PATH             LQAEVNKALAEMRADGTVEKISVKWFGADITK----
B_subtilis_YXEM            LRKKVNKALDELRKDGTLKKLSEKYFNEDITVEQKH
B_subtilis_GlnH_homo_YCKK  VVDQVNKALKEMKEDGTLSKISKKWFGEDVSK----
YA80_HAEIN                 LITKFNQVLEALRQDGTLKQISIEWFGYDITQ----
FLIY_ECOLI                 LLKAVNDAIAEMQKDGTLQALSEKWFGADVTK----
E_coli_GlnH                LRDKVNGALKTLRENGTYNEIYKKWFGTEPK-----
Deinococcus_radiodurans    LKDQIDKALTEMRSDGTFEKISQKWFGQDVGQP---
HISJ_E_COLI                LREALNKAFAEMRADGTYEKLAKKYFDFDVYGG---
;
end;
"""

    # This example uses DNA, from here:
    # http://www.molecularevolution.org/resources/fileformats/
    nxs_example2 = \
"""#NEXUS 

Begin data;
        Dimensions ntax=10 nchar=705;
        Format datatype=dna interleave=yes gap=- missing=?;
        Matrix
Cow     ATGGCATATCCCATACAACTAGGATTCCAAGATGCAACATCACCAATCATAGAAGAACTA
Carp    ATGGCACACCCAACGCAACTAGGTTTCAAGGACGCGGCCATACCCGTTATAGAGGAACTT
Chicken ATGGCCAACCACTCCCAACTAGGCTTTCAAGACGCCTCATCCCCCATCATAGAAGAGCTC
Human   ATGGCACATGCAGCGCAAGTAGGTCTACAAGACGCTACTTCCCCTATCATAGAAGAGCTT
Loach   ATGGCACATCCCACACAATTAGGATTCCAAGACGCGGCCTCACCCGTAATAGAAGAACTT
Mouse   ATGGCCTACCCATTCCAACTTGGTCTACAAGACGCCACATCCCCTATTATAGAAGAGCTA
Rat     ATGGCTTACCCATTTCAACTTGGCTTACAAGACGCTACATCACCTATCATAGAAGAACTT
Seal    ATGGCATACCCCCTACAAATAGGCCTACAAGATGCAACCTCTCCCATTATAGAGGAGTTA
Whale   ATGGCATATCCATTCCAACTAGGTTTCCAAGATGCAGCATCACCCATCATAGAAGAGCTC
Frog    ATGGCACACCCATCACAATTAGGTTTTCAAGACGCAGCCTCTCCAATTATAGAAGAATTA

Cow     CTTCACTTTCATGACCACACGCTAATAATTGTCTTCTTAATTAGCTCATTAGTACTTTAC
Carp    CTTCACTTCCACGACCACGCATTAATAATTGTGCTCCTAATTAGCACTTTAGTTTTATAT
Chicken GTTGAATTCCACGACCACGCCCTGATAGTCGCACTAGCAATTTGCAGCTTAGTACTCTAC
Human   ATCACCTTTCATGATCACGCCCTCATAATCATTTTCCTTATCTGCTTCCTAGTCCTGTAT
Loach   CTTCACTTCCATGACCATGCCCTAATAATTGTATTTTTGATTAGCGCCCTAGTACTTTAT
Mouse   ATAAATTTCCATGATCACACACTAATAATTGTTTTCCTAATTAGCTCCTTAGTCCTCTAT
Rat     ACAAACTTTCATGACCACACCCTAATAATTGTATTCCTCATCAGCTCCCTAGTACTTTAT
Seal    CTACACTTCCATGACCACACATTAATAATTGTGTTCCTAATTAGCTCATTAGTACTCTAC
Whale   CTACACTTTCACGATCATACACTAATAATCGTTTTTCTAATTAGCTCTTTAGTTCTCTAC
Frog    CTTCACTTCCACGACCATACCCTCATAGCCGTTTTTCTTATTAGTACGCTAGTTCTTTAC

Cow     ATTATTTCACTAATACTAACGACAAAGCTGACCCATACAAGCACGATAGATGCACAAGAA
Carp    ATTATTACTGCAATGGTATCAACTAAACTTACTAATAAATATATTCTAGACTCCCAAGAA
Chicken CTTCTAACTCTTATACTTATAGAAAAACTATCA---TCAAACACCGTAGATGCCCAAGAA
Human   GCCCTTTTCCTAACACTCACAACAAAACTAACTAATACTAACATCTCAGACGCTCAGGAA
Loach   GTTATTATTACAACCGTCTCAACAAAACTCACTAACATATATATTTTGGACTCACAAGAA
Mouse   ATCATCTCGCTAATATTAACAACAAAACTAACACATACAAGCACAATAGATGCACAAGAA
Rat     ATTATTTCACTAATACTAACAACAAAACTAACACACACAAGCACAATAGACGCCCAAGAA
Seal    ATTATCTCACTTATACTAACCACGAAACTCACCCACACAAGTACAATAGACGCACAAGAA
Whale   ATTATTACCCTAATGCTTACAACCAAATTAACACATACTAGTACAATAGACGCCCAAGAA
Frog    ATTATTACTATTATAATAACTACTAAACTAACTAATACAAACCTAATGGACGCACAAGAG

Cow     GTAGAGACAATCTGAACCATTCTGCCCGCCATCATCTTAATTCTAATTGCTCTTCCTTCT
Carp    ATCGAAATCGTATGAACCATTCTACCAGCCGTCATTTTAGTACTAATCGCCCTGCCCTCC
Chicken GTTGAACTAATCTGAACCATCCTACCCGCTATTGTCCTAGTCCTGCTTGCCCTCCCCTCC
Human   ATAGAAACCGTCTGAACTATCCTGCCCGCCATCATCCTAGTCCTCATCGCCCTCCCATCC
Loach   ATTGAAATCGTATGAACTGTGCTCCCTGCCCTAATCCTCATTTTAATCGCCCTCCCCTCA
Mouse   GTTGAAACCATTTGAACTATTCTACCAGCTGTAATCCTTATCATAATTGCTCTCCCCTCT
Rat     GTAGAAACAATTTGAACAATTCTCCCAGCTGTCATTCTTATTCTAATTGCCCTTCCCTCC
Seal    GTGGAAACGGTGTGAACGATCCTACCCGCTATCATTTTAATTCTCATTGCCCTACCATCA
Whale   GTAGAAACTGTCTGAACTATCCTCCCAGCCATTATCTTAATTTTAATTGCCTTGCCTTCA
Frog    ATCGAAATAGTGTGAACTATTATACCAGCTATTAGCCTCATCATAATTGCCCTTCCATCC

Cow     TTACGAATTCTATACATAATAGATGAAATCAATAACCCATCTCTTACAGTAAAAACCATA
Carp    CTACGCATCCTGTACCTTATAGACGAAATTAACGACCCTCACCTGACAATTAAAGCAATA
Chicken CTCCAAATCCTCTACATAATAGACGAAATCGACGAACCTGATCTCACCCTAAAAGCCATC
Human   CTACGCATCCTTTACATAACAGACGAGGTCAACGATCCCTCCCTTACCATCAAATCAATT
Loach   CTACGAATTCTATATCTTATAGACGAGATTAATGACCCCCACCTAACAATTAAGGCCATG
Mouse   CTACGCATTCTATATATAATAGACGAAATCAACAACCCCGTATTAACCGTTAAAACCATA
Rat     CTACGAATTCTATACATAATAGACGAGATTAATAACCCAGTTCTAACAGTAAAAACTATA
Seal    TTACGAATCCTCTACATAATGGACGAGATCAATAACCCTTCCTTGACCGTAAAAACTATA
Whale   TTACGGATCCTTTACATAATAGACGAAGTCAATAACCCCTCCCTCACTGTAAAAACAATA
Frog    CTTCGTATCCTATATTTAATAGATGAAGTTAATGATCCACACTTAACAATTAAAGCAATC

Cow     GGACATCAGTGATACTGAAGCTATGAGTATACAGATTATGAGGACTTAAGCTTCGACTCC
Carp    GGACACCAATGATACTGAAGTTACGAGTATACAGACTATGAAAATCTAGGATTCGACTCC
Chicken GGACACCAATGATACTGAACCTATGAATACACAGACTTCAAGGACCTCTCATTTGACTCC
Human   GGCCACCAATGGTACTGAACCTACGAGTACACCGACTACGGCGGACTAATCTTCAACTCC
Loach   GGGCACCAATGATACTGAAGCTACGAGTATACTGATTATGAAAACTTAAGTTTTGACTCC
Mouse   GGGCACCAATGATACTGAAGCTACGAATATACTGACTATGAAGACCTATGCTTTGATTCA
Rat     GGACACCAATGATACTGAAGCTATGAATATACTGACTATGAAGACCTATGCTTTGACTCC
Seal    GGACATCAGTGATACTGAAGCTATGAGTACACAGACTACGAAGACCTGAACTTTGACTCA
Whale   GGTCACCAATGATATTGAAGCTATGAGTATACCGACTACGAAGACCTAAGCTTCGACTCC
Frog    GGCCACCAATGATACTGAAGCTACGAATATACTAACTATGAGGATCTCTCATTTGACTCT

Cow     TACATAATTCCAACATCAGAATTAAAGCCAGGGGAGCTACGACTATTAGAAGTCGATAAT
Carp    TATATAGTACCAACCCAAGACCTTGCCCCCGGACAATTCCGACTTCTGGAAACAGACCAC
Chicken TACATAACCCCAACAACAGACCTCCCCCTAGGCCACTTCCGCCTACTAGAAGTCGACCAT
Human   TACATACTTCCCCCATTATTCCTAGAACCAGGCGACCTGCGACTCCTTGACGTTGACAAT
Loach   TACATAATCCCCACCCAGGACCTAACCCCTGGACAATTCCGGCTACTAGAGACAGACCAC
Mouse   TATATAATCCCAACAAACGACCTAAAACCTGGTGAACTACGACTGCTAGAAGTTGATAAC
Rat     TACATAATCCCAACCAATGACCTAAAACCAGGTGAACTTCGTCTATTAGAAGTTGATAAT
Seal    TATATGATCCCCACACAAGAACTAAAGCCCGGAGAACTACGACTGCTAGAAGTAGACAAT
Whale   TATATAATCCCAACATCAGACCTAAAGCCAGGAGAACTACGATTATTAGAAGTAGATAAC
Frog    TATATAATTCCAACTAATGACCTTACCCCTGGACAATTCCGGCTGCTAGAAGTTGATAAT

Cow     CGAGTTGTACTACCAATAGAAATAACAATCCGAATGTTAGTCTCCTCTGAAGACGTATTA
Carp    CGAATAGTTGTTCCAATAGAATCCCCAGTCCGTGTCCTAGTATCTGCTGAAGACGTGCTA
Chicken CGCATTGTAATCCCCATAGAATCCCCCATTCGAGTAATCATCACCGCTGATGACGTCCTC
Human   CGAGTAGTACTCCCGATTGAAGCCCCCATTCGTATAATAATTACATCACAAGACGTCTTG
Loach   CGAATGGTTGTTCCCATAGAATCCCCTATTCGCATTCTTGTTTCCGCCGAAGATGTACTA
Mouse   CGAGTCGTTCTGCCAATAGAACTTCCAATCCGTATATTAATTTCATCTGAAGACGTCCTC
Rat     CGGGTAGTCTTACCAATAGAACTTCCAATTCGTATACTAATCTCATCCGAAGACGTCCTG
Seal    CGAGTAGTCCTCCCAATAGAAATAACAATCCGCATACTAATCTCATCAGAAGATGTACTC
Whale   CGAGTTGTCTTACCTATAGAAATAACAATCCGAATATTAGTCTCATCAGAAGACGTACTC
Frog    CGAATAGTAGTCCCAATAGAATCTCCAACCCGACTTTTAGTTACAGCCGAAGACGTCCTC

Cow     CACTCATGAGCTGTGCCCTCTCTAGGACTAAAAACAGACGCAATCCCAGGCCGTCTAAAC
Carp    CATTCTTGAGCTGTTCCATCCCTTGGCGTAAAAATGGACGCAGTCCCAGGACGACTAAAT
Chicken CACTCATGAGCCGTACCCGCCCTCGGGGTAAAAACAGACGCAATCCCTGGACGACTAAAT
Human   CACTCATGAGCTGTCCCCACATTAGGCTTAAAAACAGATGCAATTCCCGGACGTCTAAAC
Loach   CACTCCTGGGCCCTTCCAGCCATGGGGGTAAAGATAGACGCGGTCCCAGGACGCCTTAAC
Mouse   CACTCATGAGCAGTCCCCTCCCTAGGACTTAAAACTGATGCCATCCCAGGCCGACTAAAT
Rat     CACTCATGAGCCATCCCTTCACTAGGGTTAAAAACCGACGCAATCCCCGGCCGCCTAAAC
Seal    CACTCATGAGCCGTACCGTCCCTAGGACTAAAAACTGATGCTATCCCAGGACGACTAAAC
Whale   CACTCATGGGCCGTACCCTCCTTGGGCCTAAAAACAGATGCAATCCCAGGACGCCTAAAC
Frog    CACTCGTGAGCTGTACCCTCCTTGGGTGTCAAAACAGATGCAATCCCAGGACGACTTCAT

Cow     CAAACAACCCTTATATCGTCCCGTCCAGGCTTATATTACGGTCAATGCTCAGAAATTTGC
Carp    CAAGCCGCCTTTATTGCCTCACGCCCAGGGGTCTTTTACGGACAATGCTCTGAAATTTGT
Chicken CAAACCTCCTTCATCACCACTCGACCAGGAGTGTTTTACGGACAATGCTCAGAAATCTGC
Human   CAAACCACTTTCACCGCTACACGACCGGGGGTATACTACGGTCAATGCTCTGAAATCTGT
Loach   CAAACCGCCTTTATTGCCTCCCGCCCCGGGGTATTCTATGGGCAATGCTCAGAAATCTGT
Mouse   CAAGCAACAGTAACATCAAACCGACCAGGGTTATTCTATGGCCAATGCTCTGAAATTTGT
Rat     CAAGCTACAGTCACATCAAACCGACCAGGTCTATTCTATGGCCAATGCTCTGAAATTTGC
Seal    CAAACAACCCTAATAACCATACGACCAGGACTGTACTACGGTCAATGCTCAGAAATCTGT
Whale   CAAACAACCTTAATATCAACACGACCAGGCCTATTTTATGGACAATGCTCAGAGATCTGC
Frog    CAAACATCATTTATTGCTACTCGTCCGGGAGTATTTTACGGACAATGTTCAGAAATTTGC

Cow     GGGTCAAACCACAGTTTCATACCCATTGTCCTTGAGTTAGTCCCACTAAAGTACTTTGAA
Carp    GGAGCTAATCACAGCTTTATACCAATTGTAGTTGAAGCAGTACCTCTCGAACACTTCGAA
Chicken GGAGCTAACCACAGCTACATACCCATTGTAGTAGAGTCTACCCCCCTAAAACACTTTGAA
Human   GGAGCAAACCACAGTTTCATGCCCATCGTCCTAGAATTAATTCCCCTAAAAATCTTTGAA
Loach   GGAGCAAACCACAGCTTTATACCCATCGTAGTAGAAGCGGTCCCACTATCTCACTTCGAA
Mouse   GGATCTAACCATAGCTTTATGCCCATTGTCCTAGAAATGGTTCCACTAAAATATTTCGAA
Rat     GGCTCAAATCACAGCTTCATACCCATTGTACTAGAAATAGTGCCTCTAAAATATTTCGAA
Seal    GGTTCAAACCACAGCTTCATACCTATTGTCCTCGAATTGGTCCCACTATCCCACTTCGAG
Whale   GGCTCAAACCACAGTTTCATACCAATTGTCCTAGAACTAGTACCCCTAGAAGTCTTTGAA
Frog    GGAGCAAACCACAGCTTTATACCAATTGTAGTTGAAGCAGTACCGCTAACCGACTTTGAA

Cow     AAATGATCTGCGTCAATATTA---------------------TAA
Carp    AACTGATCCTCATTAATACTAGAAGACGCCTCGCTAGGAAGCTAA
Chicken GCCTGATCCTCACTA------------------CTGTCATCTTAA
Human   ATA---------------------GGGCCCGTATTTACCCTATAG
Loach   AACTGGTCCACCCTTATACTAAAAGACGCCTCACTAGGAAGCTAA
Mouse   AACTGATCTGCTTCAATAATT---------------------TAA
Rat     AACTGATCAGCTTCTATAATT---------------------TAA
Seal    AAATGATCTACCTCAATGCTT---------------------TAA
Whale   AAATGATCTGTATCAATACTA---------------------TAA
Frog    AACTGATCTTCATCAATACTA---GAAGCATCACTA------AGA
        ;
End;
"""

    # This example uses amino acids, from here:
    # http://www.molecularevolution.org/resources/fileformats/
    nxs_example3 = \
"""#NEXUS 

Begin data;
        Dimensions ntax=10 nchar=234;
        Format datatype=protein gap=- interleave;
        Matrix
Cow     MAYPMQLGFQDATSPIMEELLHFHDHTLMIVFLISSLVLYIISLMLTTKLTHTSTMDAQE
Carp    MAHPTQLGFKDAAMPVMEELLHFHDHALMIVLLISTLVLYIITAMVSTKLTNKYILDSQE
Chicken MANHSQLGFQDASSPIMEELVEFHDHALMVALAICSLVLYLLTLMLMEKLS-SNTVDAQE
Human   MAHAAQVGLQDATSPIMEELITFHDHALMIIFLICFLVLYALFLTLTTKLTNTNISDAQE
Loach   MAHPTQLGFQDAASPVMEELLHFHDHALMIVFLISALVLYVIITTVSTKLTNMYILDSQE
Mouse   MAYPFQLGLQDATSPIMEELMNFHDHTLMIVFLISSLVLYIISLMLTTKLTHTSTMDAQE
Rat     MAYPFQLGLQDATSPIMEELTNFHDHTLMIVFLISSLVLYIISLMLTTKLTHTSTMDAQE
Seal    MAYPLQMGLQDATSPIMEELLHFHDHTLMIVFLISSLVLYIISLMLTTKLTHTSTMDAQE
Whale   MAYPFQLGFQDAASPIMEELLHFHDHTLMIVFLISSLVLYIITLMLTTKLTHTSTMDAQE
Frog    MAHPSQLGFQDAASPIMEELLHFHDHTLMAVFLISTLVLYIITIMMTTKLTNTNLMDAQE

Cow     VETIWTILPAIILILIALPSLRILYMMDEINNPSLTVKTMGHQWYWSYEYTDYEDLSFDS
Carp    IEIVWTILPAVILVLIALPSLRILYLMDEINDPHLTIKAMGHQWYWSYEYTDYENLGFDS
Chicken VELIWTILPAIVLVLLALPSLQILYMMDEIDEPDLTLKAIGHQWYWTYEYTDFKDLSFDS
Human   METVWTILPAIILVLIALPSLRILYMTDEVNDPSLTIKSIGHQWYWTYEYTDYGGLIFNS
Loach   IEIVWTVLPALILILIALPSLRILYLMDEINDPHLTIKAMGHQWYWSYEYTDYENLSFDS
Mouse   VETIWTILPAVILIMIALPSLRILYMMDEINNPVLTVKTMGHQWYWSYEYTDYEDLCFDS
Rat     VETIWTILPAVILILIALPSLRILYMMDEINNPVLTVKTMGHQWYWSYEYTDYEDLCFDS
Seal    VETVWTILPAIILILIALPSLRILYMMDEINNPSLTVKTMGHQWYWSYEYTDYEDLNFDS
Whale   VETVWTILPAIILILIALPSLRILYMMDEVNNPSLTVKTMGHQWYWSYEYTDYEDLSFDS
Frog    IEMVWTIMPAISLIMIALPSLRILYLMDEVNDPHLTIKAIGHQWYWSYEYTNYEDLSFDS

Cow     YMIPTSELKPGELRLLEVDNRVVLPMEMTIRMLVSSEDVLHSWAVPSLGLKTDAIPGRLN
Carp    YMVPTQDLAPGQFRLLETDHRMVVPMESPVRVLVSAEDVLHSWAVPSLGVKMDAVPGRLN
Chicken YMTPTTDLPLGHFRLLEVDHRIVIPMESPIRVIITADDVLHSWAVPALGVKTDAIPGRLN
Human   YMLPPLFLEPGDLRLLDVDNRVVLPIEAPIRMMITSQDVLHSWAVPTLGLKTDAIPGRLN
Loach   YMIPTQDLTPGQFRLLETDHRMVVPMESPIRILVSAEDVLHSWALPAMGVKMDAVPGRLN
Mouse   YMIPTNDLKPGELRLLEVDNRVVLPMELPIRMLISSEDVLHSWAVPSLGLKTDAIPGRLN
Rat     YMIPTNDLKPGELRLLEVDNRVVLPMELPIRMLISSEDVLHSWAIPSLGLKTDAIPGRLN
Seal    YMIPTQELKPGELRLLEVDNRVVLPMEMTIRMLISSEDVLHSWAVPSLGLKTDAIPGRLN
Whale   YMIPTSDLKPGELRLLEVDNRVVLPMEMTIRMLVSSEDVLHSWAVPSLGLKTDAIPGRLN
Frog    YMIPTNDLTPGQFRLLEVDNRMVVPMESPTRLLVTAEDVLHSWAVPSLGVKTDAIPGRLH

Cow     QTTLMSSRPGLYYGQCSEICGSNHSFMPIVLELVPLKYFEKWSASML-------
Carp    QAAFIASRPGVFYGQCSEICGANHSFMPIVVEAVPLEHFENWSSLMLEDASLGS
Chicken QTSFITTRPGVFYGQCSEICGANHSYMPIVVESTPLKHFEAWSSL------LSS
Human   QTTFTATRPGVYYGQCSEICGANHSFMPIVLELIPLKIFEM-------GPVFTL
Loach   QTAFIASRPGVFYGQCSEICGANHSFMPIVVEAVPLSHFENWSTLMLKDASLGS
Mouse   QATVTSNRPGLFYGQCSEICGSNHSFMPIVLEMVPLKYFENWSASMI-------
Rat     QATVTSNRPGLFYGQCSEICGSNHSFMPIVLEMVPLKYFENWSASMI-------
Seal    QTTLMTMRPGLYYGQCSEICGSNHSFMPIVLELVPLSHFEKWSTSML-------
Whale   QTTLMSTRPGLFYGQCSEICGSNHSFMPIVLELVPLEVFEKWSVSML-------
Frog    QTSFIATRPGVFYGQCSEICGANHSFMPIVVEAVPLTDFENWSSSML-EASL--
        ;
End;
"""
    
    # This example with its slightly odd (partial) annotation is from here:
    # http://www.cgb.ki.se/cgb/groups/sonnhammer/Stockholm.html
    sth_example = \
"""# STOCKHOLM 1.0
#=GF ID CBS
#=GF AC PF00571
#=GF DE CBS domain
#=GF AU Bateman A
#=GF CC CBS domains are small intracellular modules mostly found  
#=GF CC in 2 or four copies within a protein. 
#=GF SQ 67
#=GS O31698/18-71 AC O31698
#=GS O83071/192-246 AC O83071
#=GS O83071/259-312 AC O83071
#=GS O31698/88-139 AC O31698
#=GS O31698/88-139 OS Bacillus subtilis
O83071/192-246          MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRVPVYERS
#=GR O83071/192-246 SA  999887756453524252..55152525....36463774777
O83071/259-312          MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVAIVLDEY
#=GR O83071/259-312 SS  CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEE
O31698/18-71            MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAIPVLDPS
#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH
O31698/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31698/88-139 SS   CCCCCCCHHHHHHHHHHH..HEEEEEEE....EEEEEEEEEEH
#=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH
O31699/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31699/88-139 AS   ________________*__________________________
#=GR_O31699/88-139_IN   ____________1______________2__________0____
//
"""

    # Interlaced example from BioPerl documentation.  Also note the blank line.
    # http://www.bioperl.org/wiki/Stockholm_multiple_alignment_format
    sth_example2 = \
"""# STOCKHOLM 1.0
#=GC SS_cons       .................<<<<<<<<...<<<<<<<........>>>>>>>..
AP001509.1         UUAAUCGAGCUCAACACUCUUCGUAUAUCCUC-UCAAUAUGG-GAUGAGGGU
#=GR AP001509.1 SS -----------------<<<<<<<<---..<<-<<-------->>->>..--
AE007476.1         AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGU
#=GR AE007476.1 SS -----------------<<<<<<<<-----<<.<<-------->>.>>----

#=GC SS_cons       ......<<<<<<<.......>>>>>>>..>>>>>>>>...............
AP001509.1         CUCUAC-AGGUA-CCGUAAA-UACCUAGCUACGAAAAGAAUGCAGUUAAUGU
#=GR AP001509.1 SS -------<<<<<--------->>>>>--->>>>>>>>---------------
AE007476.1         UUCUACAAGGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU
#=GR AE007476.1 SS ------.<<<<<--------->>>>>.-->>>>>>>>---------------
//"""

    # Sample GenBank record from here:
    # http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
    gbk_example = \
"""LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999
DEFINITION  Saccharomyces cerevisiae TCP1-beta gene, partial cds, and Axl2p
            (AXL2) and Rev7p (REV7) genes, complete cds.
ACCESSION   U49845
VERSION     U49845.1  GI:1293613
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae (baker's yeast)
  ORGANISM  Saccharomyces cerevisiae
            Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;
            Saccharomycetales; Saccharomycetaceae; Saccharomyces.
REFERENCE   1  (bases 1 to 5028)
  AUTHORS   Torpey,L.E., Gibbs,P.E., Nelson,J. and Lawrence,C.W.
  TITLE     Cloning and sequence of REV7, a gene whose function is required for
            DNA damage-induced mutagenesis in Saccharomyces cerevisiae
  JOURNAL   Yeast 10 (11), 1503-1509 (1994)
  PUBMED    7871890
REFERENCE   2  (bases 1 to 5028)
  AUTHORS   Roemer,T., Madden,K., Chang,J. and Snyder,M.
  TITLE     Selection of axial growth sites in yeast requires Axl2p, a novel
            plasma membrane glycoprotein
  JOURNAL   Genes Dev. 10 (7), 777-793 (1996)
  PUBMED    8846915
REFERENCE   3  (bases 1 to 5028)
  AUTHORS   Roemer,T.
  TITLE     Direct Submission
  JOURNAL   Submitted (22-FEB-1996) Terry Roemer, Biology, Yale University, New
            Haven, CT, USA
FEATURES             Location/Qualifiers
     source          1..5028
                     /organism="Saccharomyces cerevisiae"
                     /db_xref="taxon:4932"
                     /chromosome="IX"
                     /map="9"
     CDS             <1..206
                     /codon_start=3
                     /product="TCP1-beta"
                     /protein_id="AAA98665.1"
                     /db_xref="GI:1293614"
                     /translation="SSIYNGISTSGLDLNNGTIADMRQLGIVESYKLKRAVVSSASEA
                     AEVLLRVDNIIRARPRTANRQHM"
     gene            687..3158
                     /gene="AXL2"
     CDS             687..3158
                     /gene="AXL2"
                     /note="plasma membrane glycoprotein"
                     /codon_start=1
                     /function="required for axial budding pattern of S.
                     cerevisiae"
                     /product="Axl2p"
                     /protein_id="AAA98666.1"
                     /db_xref="GI:1293615"
                     /translation="MTQLQISLLLTATISLLHLVVATPYEAYPIGKQYPPVARVNESF
                     TFQISNDTYKSSVDKTAQITYNCFDLPSWLSFDSSSRTFSGEPSSDLLSDANTTLYFN
                     VILEGTDSADSTSLNNTYQFVVTNRPSISLSSDFNLLALLKNYGYTNGKNALKLDPNE
                     VFNVTFDRSMFTNEESIVSYYGRSQLYNAPLPNWLFFDSGELKFTGTAPVINSAIAPE
                     TSYSFVIIATDIEGFSAVEVEFELVIGAHQLTTSIQNSLIINVTDTGNVSYDLPLNYV
                     YLDDDPISSDKLGSINLLDAPDWVALDNATISGSVPDELLGKNSNPANFSVSIYDTYG
                     DVIYFNFEVVSTTDLFAISSLPNINATRGEWFSYYFLPSQFTDYVNTNVSLEFTNSSQ
                     DHDWVKFQSSNLTLAGEVPKNFDKLSLGLKANQGSQSQELYFNIIGMDSKITHSNHSA
                     NATSTRSSHHSTSTSSYTSSTYTAKISSTSAAATSSAPAALPAANKTSSHNKKAVAIA
                     CGVAIPLGVILVALICFLIFWRRRRENPDDENLPHAISGPDLNNPANKPNQENATPLN
                     NPFDDDASSYDDTSIARRLAALNTLKLDNHSATESDISSVDEKRDSLSGMNTYNDQFQ
                     SQSKEELLAKPPVQPPESPFFDPQNRSSSVYMDSEPAVNKSWRYTGNLSPVSDIVRDS
                     YGSQKTVDTEKLFDLEAPEKEKRTSRDVTMSSLDPWNSNISPSPVRKSVTPSPYNVTK
                     HRNRHLQNIQDSQSGKNGITPTTMSTSSSDDFVPVKDGENFCWVHSMEPDRRPSKKRL
                     VDFSNKSNVNVGQVKDIHGRIPEML"
     gene            complement(3300..4037)
                     /gene="REV7"
     CDS             complement(3300..4037)
                     /gene="REV7"
                     /codon_start=1
                     /product="Rev7p"
                     /protein_id="AAA98667.1"
                     /db_xref="GI:1293616"
                     /translation="MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYTTYQSFNLPQ
                     FVPINRHPALIDYIEELILDVLSKLTHVYRFSICIINKKNDLCIEKYVLDFSELQHVD
                     KDDQIITETEVFDEFRSSLNSLIMHLEKLPKVNDDTITFEAVINAIELELGHKLDRNR
                     RVDSLEEKAEIERDSNWVKCQEDENLPDNNGFQPPKIKLTSLVGSDVGPLIIHQFSEK
                     LISGDDKILNGVYSQYEEGESIFGSLF"
ORIGIN
        1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg
       61 ccgacatgag acagttaggt atcgtcgaga gttacaagct aaaacgagca gtagtcagct
      121 ctgcatctga agccgctgaa gttctactaa gggtggataa catcatccgt gcaagaccaa
      181 gaaccgccaa tagacaacat atgtaacata tttaggatat acctcgaaaa taataaaccg
      241 ccacactgtc attattataa ttagaaacag aacgcaaaaa ttatccacta tataattcaa
      301 agacgcgaaa aaaaaagaac aacgcgtcat agaacttttg gcaattcgcg tcacaaataa
      361 attttggcaa cttatgtttc ctcttcgagc agtactcgag ccctgtctca agaatgtaat
      421 aatacccatc gtaggtatgg ttaaagatag catctccaca acctcaaagc tccttgccga
      481 gagtcgccct cctttgtcga gtaattttca cttttcatat gagaacttat tttcttattc
      541 tttactctca catcctgtag tgattgacac tgcaacagcc accatcacta gaagaacaga
      601 acaattactt aatagaaaaa ttatatcttc ctcgaaacga tttcctgctt ccaacatcta
      661 cgtatatcaa gaagcattca cttaccatga cacagcttca gatttcatta ttgctgacag
      721 ctactatatc actactccat ctagtagtgg ccacgcccta tgaggcatat cctatcggaa
      781 aacaataccc cccagtggca agagtcaatg aatcgtttac atttcaaatt tccaatgata
      841 cctataaatc gtctgtagac aagacagctc aaataacata caattgcttc gacttaccga
      901 gctggctttc gtttgactct agttctagaa cgttctcagg tgaaccttct tctgacttac
      961 tatctgatgc gaacaccacg ttgtatttca atgtaatact cgagggtacg gactctgccg
     1021 acagcacgtc tttgaacaat acataccaat ttgttgttac aaaccgtcca tccatctcgc
     1081 tatcgtcaga tttcaatcta ttggcgttgt taaaaaacta tggttatact aacggcaaaa
     1141 acgctctgaa actagatcct aatgaagtct tcaacgtgac ttttgaccgt tcaatgttca
     1201 ctaacgaaga atccattgtg tcgtattacg gacgttctca gttgtataat gcgccgttac
     1261 ccaattggct gttcttcgat tctggcgagt tgaagtttac tgggacggca ccggtgataa
     1321 actcggcgat tgctccagaa acaagctaca gttttgtcat catcgctaca gacattgaag
     1381 gattttctgc cgttgaggta gaattcgaat tagtcatcgg ggctcaccag ttaactacct
     1441 ctattcaaaa tagtttgata atcaacgtta ctgacacagg taacgtttca tatgacttac
     1501 ctctaaacta tgtttatctc gatgacgatc ctatttcttc tgataaattg ggttctataa
     1561 acttattgga tgctccagac tgggtggcat tagataatgc taccatttcc gggtctgtcc
     1621 cagatgaatt actcggtaag aactccaatc ctgccaattt ttctgtgtcc atttatgata
     1681 cttatggtga tgtgatttat ttcaacttcg aagttgtctc cacaacggat ttgtttgcca
     1741 ttagttctct tcccaatatt aacgctacaa ggggtgaatg gttctcctac tattttttgc
     1801 cttctcagtt tacagactac gtgaatacaa acgtttcatt agagtttact aattcaagcc
     1861 aagaccatga ctgggtgaaa ttccaatcat ctaatttaac attagctgga gaagtgccca
     1921 agaatttcga caagctttca ttaggtttga aagcgaacca aggttcacaa tctcaagagc
     1981 tatattttaa catcattggc atggattcaa agataactca ctcaaaccac agtgcgaatg
     2041 caacgtccac aagaagttct caccactcca cctcaacaag ttcttacaca tcttctactt
     2101 acactgcaaa aatttcttct acctccgctg ctgctacttc ttctgctcca gcagcgctgc
     2161 cagcagccaa taaaacttca tctcacaata aaaaagcagt agcaattgcg tgcggtgttg
     2221 ctatcccatt aggcgttatc ctagtagctc tcatttgctt cctaatattc tggagacgca
     2281 gaagggaaaa tccagacgat gaaaacttac cgcatgctat tagtggacct gatttgaata
     2341 atcctgcaaa taaaccaaat caagaaaacg ctacaccttt gaacaacccc tttgatgatg
     2401 atgcttcctc gtacgatgat acttcaatag caagaagatt ggctgctttg aacactttga
     2461 aattggataa ccactctgcc actgaatctg atatttccag cgtggatgaa aagagagatt
     2521 ctctatcagg tatgaataca tacaatgatc agttccaatc ccaaagtaaa gaagaattat
     2581 tagcaaaacc cccagtacag cctccagaga gcccgttctt tgacccacag aataggtctt
     2641 cttctgtgta tatggatagt gaaccagcag taaataaatc ctggcgatat actggcaacc
     2701 tgtcaccagt ctctgatatt gtcagagaca gttacggatc acaaaaaact gttgatacag
     2761 aaaaactttt cgatttagaa gcaccagaga aggaaaaacg tacgtcaagg gatgtcacta
     2821 tgtcttcact ggacccttgg aacagcaata ttagcccttc tcccgtaaga aaatcagtaa
     2881 caccatcacc atataacgta acgaagcatc gtaaccgcca cttacaaaat attcaagact
     2941 ctcaaagcgg taaaaacgga atcactccca caacaatgtc aacttcatct tctgacgatt
     3001 ttgttccggt taaagatggt gaaaattttt gctgggtcca tagcatggaa ccagacagaa
     3061 gaccaagtaa gaaaaggtta gtagattttt caaataagag taatgtcaat gttggtcaag
     3121 ttaaggacat tcacggacgc atcccagaaa tgctgtgatt atacgcaacg atattttgct
     3181 taattttatt ttcctgtttt attttttatt agtggtttac agatacccta tattttattt
     3241 agtttttata cttagagaca tttaatttta attccattct tcaaatttca tttttgcact
     3301 taaaacaaag atccaaaaat gctctcgccc tcttcatatt gagaatacac tccattcaaa
     3361 attttgtcgt caccgctgat taatttttca ctaaactgat gaataatcaa aggccccacg
     3421 tcagaaccga ctaaagaagt gagttttatt ttaggaggtt gaaaaccatt attgtctggt
     3481 aaattttcat cttcttgaca tttaacccag tttgaatccc tttcaatttc tgctttttcc
     3541 tccaaactat cgaccctcct gtttctgtcc aacttatgtc ctagttccaa ttcgatcgca
     3601 ttaataactg cttcaaatgt tattgtgtca tcgttgactt taggtaattt ctccaaatgc
     3661 ataatcaaac tatttaagga agatcggaat tcgtcgaaca cttcagtttc cgtaatgatc
     3721 tgatcgtctt tatccacatg ttgtaattca ctaaaatcta aaacgtattt ttcaatgcat
     3781 aaatcgttct ttttattaat aatgcagatg gaaaatctgt aaacgtgcgt taatttagaa
     3841 agaacatcca gtataagttc ttctatatag tcaattaaag caggatgcct attaatggga
     3901 acgaactgcg gcaagttgaa tgactggtaa gtagtgtagt cgaatgactg aggtgggtat
     3961 acatttctat aaaataaaat caaattaatg tagcatttta agtataccct cagccacttc
     4021 tctacccatc tattcataaa gctgacgcaa cgattactat tttttttttc ttcttggatc
     4081 tcagtcgtcg caaaaacgta taccttcttt ttccgacctt ttttttagct ttctggaaaa
     4141 gtttatatta gttaaacagg gtctagtctt agtgtgaaag ctagtggttt cgattgactg
     4201 atattaagaa agtggaaatt aaattagtag tgtagacgta tatgcatatg tatttctcgc
     4261 ctgtttatgt ttctacgtac ttttgattta tagcaagggg aaaagaaata catactattt
     4321 tttggtaaag gtgaaagcat aatgtaaaag ctagaataaa atggacgaaa taaagagagg
     4381 cttagttcat cttttttcca aaaagcaccc aatgataata actaaaatga aaaggatttg
     4441 ccatctgtca gcaacatcag ttgtgtgagc aataataaaa tcatcacctc cgttgccttt
     4501 agcgcgtttg tcgtttgtat cttccgtaat tttagtctta tcaatgggaa tcataaattt
     4561 tccaatgaat tagcaatttc gtccaattct ttttgagctt cttcatattt gctttggaat
     4621 tcttcgcact tcttttccca ttcatctctt tcttcttcca aagcaacgat ccttctaccc
     4681 atttgctcag agttcaaatc ggcctctttc agtttatcca ttgcttcctt cagtttggct
     4741 tcactgtctt ctagctgttg ttctagatcc tggtttttct tggtgtagtt ctcattatta
     4801 gatctcaagt tattggagtc ttcagccaat tgctttgtat cagacaattg actctctaac
     4861 ttctccactt cactgtcgag ttgctcgttt ttagcggaca aagatttaat ctcgttttct
     4921 ttttcagtgt tagattgctc taattctttg agctgttctc tcagctcctc atatttttct
     4981 tgccatgact cagattctaa ttttaagcta ttcaatttct ctttgatc
//"""

    # GenBank format protein (aka GenPept) file from:
    # http://www.molecularevolution.org/resources/fileformats/
    gbk_example2 = \
"""LOCUS       AAD51968                 143 aa            linear   BCT 21-AUG-2001
DEFINITION  transcriptional regulator RovA [Yersinia enterocolitica].
ACCESSION   AAD51968
VERSION     AAD51968.1  GI:5805369
DBSOURCE    locus AF171097 accession AF171097.1
KEYWORDS    .
SOURCE      Yersinia enterocolitica
  ORGANISM  Yersinia enterocolitica
            Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales;
            Enterobacteriaceae; Yersinia.
REFERENCE   1  (residues 1 to 143)
  AUTHORS   Revell,P.A. and Miller,V.L.
  TITLE     A chromosomally encoded regulator is required for expression of the
            Yersinia enterocolitica inv gene and for virulence
  JOURNAL   Mol. Microbiol. 35 (3), 677-685 (2000)
  MEDLINE   20138369
   PUBMED   10672189
REFERENCE   2  (residues 1 to 143)
  AUTHORS   Revell,P.A. and Miller,V.L.
  TITLE     Direct Submission
  JOURNAL   Submitted (22-JUL-1999) Molecular Microbiology, Washington
            University School of Medicine, Campus Box 8230, 660 South Euclid,
            St. Louis, MO 63110, USA
COMMENT     Method: conceptual translation.
FEATURES             Location/Qualifiers
     source          1..143
                     /organism="Yersinia enterocolitica"
                     /mol_type="unassigned DNA"
                     /strain="JB580v"
                     /serotype="O:8"
                     /db_xref="taxon:630"
     Protein         1..143
                     /product="transcriptional regulator RovA"
                     /name="regulates inv expression"
     CDS             1..143
                     /gene="rovA"
                     /coded_by="AF171097.1:380..811"
                     /note="regulator of virulence"
                     /transl_table=11
ORIGIN      
        1 mestlgsdla rlvrvwrali dhrlkplelt qthwvtlhni nrlppeqsqi qlakaigieq
       61 pslvrtldql eekglitrht candrrakri klteqsspii eqvdgvicst rkeilggisp
      121 deiellsgli dklerniiql qsk
//"""



    print "#########################################################"
    print "# Sequence Input Tests                                  #"
    print "#########################################################"

    #ToDo - Check alphabet, or at least DNA/amino acid, for those
    #       filetype that specify it (e.g. Nexus, GenBank)
    tests = [
         (aln_example,  "clustal",   8, "HISJ_E_COLI",
          "MKKLVLSLSLVLAFSSATAAF-------------------AAIPQNIRIG" + \
          "TDPTYAPFESKNS-QGELVGFDIDLAKELCKRINTQCTFVENPLDALIPS" + \
          "LKAKKIDAIMSSLSITEKRQQEIAFTDKLYAADSRLVVAKNSDIQP-TVE" + \
          "SLKGKRVGVLQGTTQETFGNEHWAPKGIEIVSYQGQDNIYSDLTAGRIDA" + \
          "AFQDEVAASEGFLKQPVGKDYKFGGPSVKDEKLFGVGTGMGLRKED--NE" + \
          "LREALNKAFAEMRADGTYEKLAKKYFDFDVYGG---", True),
         (phy_example,  "phylip",    8, "HISJ_E_COL", None, False),
         (nxs_example,  "nexus",     8, "HISJ_E_COLI", None, True),
         (nxs_example2, "nexus",    10, "Frog",
          "ATGGCACACCCATCACAATTAGGTTTTCAAGACGCAGCCTCTCCAATTATAGAAGAATTA" + \
          "CTTCACTTCCACGACCATACCCTCATAGCCGTTTTTCTTATTAGTACGCTAGTTCTTTAC" + \
          "ATTATTACTATTATAATAACTACTAAACTAACTAATACAAACCTAATGGACGCACAAGAG" + \
          "ATCGAAATAGTGTGAACTATTATACCAGCTATTAGCCTCATCATAATTGCCCTTCCATCC" + \
          "CTTCGTATCCTATATTTAATAGATGAAGTTAATGATCCACACTTAACAATTAAAGCAATC" + \
          "GGCCACCAATGATACTGAAGCTACGAATATACTAACTATGAGGATCTCTCATTTGACTCT" + \
          "TATATAATTCCAACTAATGACCTTACCCCTGGACAATTCCGGCTGCTAGAAGTTGATAAT" + \
          "CGAATAGTAGTCCCAATAGAATCTCCAACCCGACTTTTAGTTACAGCCGAAGACGTCCTC" + \
          "CACTCGTGAGCTGTACCCTCCTTGGGTGTCAAAACAGATGCAATCCCAGGACGACTTCAT" + \
          "CAAACATCATTTATTGCTACTCGTCCGGGAGTATTTTACGGACAATGTTCAGAAATTTGC" + \
          "GGAGCAAACCACAGCTTTATACCAATTGTAGTTGAAGCAGTACCGCTAACCGACTTTGAA" + \
          "AACTGATCTTCATCAATACTA---GAAGCATCACTA------AGA", True),
         (nxs_example3, "nexus",    10, "Frog",
          'MAHPSQLGFQDAASPIMEELLHFHDHTLMAVFLISTLVLYIITIMMTTKLTNTNLMDAQE' + \
          'IEMVWTIMPAISLIMIALPSLRILYLMDEVNDPHLTIKAIGHQWYWSYEYTNYEDLSFDS' + \
          'YMIPTNDLTPGQFRLLEVDNRMVVPMESPTRLLVTAEDVLHSWAVPSLGVKTDAIPGRLH' + \
          'QTSFIATRPGVFYGQCSEICGANHSFMPIVVEAVPLTDFENWSSSML-EASL--', True),
         (sth_example,  "stockholm", 5, "O31699/88-139",
          'EVMLTDIPRLHINDPIMK--GFGMVINN------GFVCVENDE', True),
         (sth_example2, "stockholm", 2, "AE007476.1",
          'AAAAUUGAAUAUCGUUUUACUUGUUUAU-GUCGUGAAU-UGG-CACGA-CGU' + \
          'UUCUACAAGGUG-CCGG-AA-CACCUAACAAUAAGUAAGUCAGCAGUGAGAU', True),
         (gbk_example, "genbank", 1, "U49845.1", None, True),
         (gbk_example2,"genbank", 1, 'AAD51968.1',
          "MESTLGSDLARLVRVWRALIDHRLKPLELTQTHWVTLHNINRLPPEQSQIQLAKAIGIEQ" + \
          "PSLVRTLDQLEEKGLITRHTCANDRRAKRIKLTEQSSPIIEQVDGVICSTRKEILGGISP" + \
          "DEIELLSGLIDKLERNIIQLQSK", True),
    ]
    
    for (data, format, rec_count, last_id, last_seq, dict_check) in tests:
        #Note all the examples have one and only one alignment
        print "%s file with %i records" % (format, rec_count)

        print "Bio.AlignIO.read(handle, format)"
        alignment = read(StringIO(data), format)
        assert len(alignment.get_all_seqs()) == rec_count

        print "Bio.AlignIO.read(handle, format, seq_count)"
        alignment = read(StringIO(data), format, rec_count)
        assert len(alignment.get_all_seqs()) == rec_count
        
        print "Bio.AlignIO.parse(handle, format)"
        #Basic check, turning the iterator into a list...
        #This uses "for x in iterator" interally.
        iterator = parse(StringIO(data), format=format)
        as_list = list(iterator)

        assert len(as_list) == 1
        assert len(as_list[0].get_all_seqs()) == rec_count, \
            "Expected %i records, found %i" \
            % (rec_count, len(as_list[0].get_all_seqs()))
        assert as_list[0].get_all_seqs()[-1].id == last_id, \
            "Expected '%s' as last record ID, found '%s'" \
            % (last_id, as_list[0].get_all_seqs()[-1].id)
        if last_seq :
            assert as_list[0].get_all_seqs()[-1].seq.tostring() == last_seq

        print "Bio.AlignIO.parse(handle, format, seq_count)"
        as_list2 = list(parse(StringIO(data), format=format, seq_count=rec_count))
        assert len(as_list2) == len(as_list)
        for a1, a2 in zip(as_list, as_list2) :
            assert align_cmp(a1, a2)

        half = rec_count / 2
        if half*2 == rec_count :
            #Even... try splitting the alignments in two.
            #This should work for things parsed by Bio.SeqIO like fasta,
            #but fail for things parsed by Bio.AlignIO itself, like phylip,
            #clustal, stockholm, ...
            try :
                list(parse(StringIO(data), format=format, seq_count=half))
                assert format not in _FormatToIterator
            except ValueError, e :
                assert format in _FormatToIterator, \
                       "Format %s, %s" % (format, str(e))
        del half
            

        print "Iteration using .next()"
        #Test iteration including use of the next() method and "for x in iterator"
        iterator = parse(StringIO(data), format=format)
        count = 1
        alignment = iterator.next()
        assert alignment is not None
        assert str(alignment.__class__) == "Bio.Align.Generic.Alignment"
        #print record
        for alignment in iterator :
            assert len(alignment.get_all_seqs()) == len(as_list[0].get_all_seqs())
            count = count + 1
        assert count == len(as_list)

        #Test iteration using just next() method
        iterator = parse(StringIO(data), format=format)
        count = 0
        while True :
            try :
                alignment = iterator.next()
            except StopIteration :
                break
            if alignment is None : break
            assert len(alignment.get_all_seqs()) == len(as_list[0].get_all_seqs())
            count = count + 1
        assert count == len(as_list)

        print "parse(handle)"
        iterator = parse(StringIO(data), format=format)
        for (i, alignment) in enumerate(iterator) :
            pass
        assert i+1 == len(as_list)

        if format not in ["nexus"] :
            print "Triple copy of data"
            #Basic check, turning the iterator into a list...
            #This uses "for x in iterator" interally.
            iterator = parse(StringIO(data + "\n" + data + "\n" + data), format=format)
            triple_list = list(iterator)
            if format in _FormatToIterator :
                #This format should be understood, so three alignments
                assert len(triple_list) == 3
                for a in triple_list :
                    assert len(a.get_all_seqs()) == rec_count
            else :
                #We have forced this format into a single alignment
                assert len(triple_list) == 1
                assert len(triple_list[0].get_all_seqs()) == 3 * rec_count

            #Try with explicit count argument
            assert 3==len(list(parse(StringIO(data + "\n" + data + "\n" + data), format, rec_count)))

            try :
                alignment = read(StringIO(data + "\n" + data + "\n" + data), format, rec_count)
                assert False, "Should have failed"
            except ValueError :
                #Failed, good
                pass

        for out_format in _FormatToWriter :
            print "writing to %s" % out_format
            #Going to write to a handle...
            handle = StringIO()
            
            try :
                write(as_list, handle=handle, format=out_format)
            except ValueError, e :
                #This is often expected to happen, for example when we try and
                #write sequences of different lengths to an alignment file.
                print "Failed: %s" % str(e)
                #Carry on to the next format:
                continue
        
        print

    print "#########################################################"
    print "# AlignIO Tests finished                                #"
    print "#########################################################"
