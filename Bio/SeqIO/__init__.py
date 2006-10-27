# Copyright 2006 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
#New Wiki page:
# http://biopython.org/wiki/SeqIO
#
#Nice link:
# http://www.ebi.ac.uk/help/formats_frame.html

"""Sequence input/output designed to look similar to the bioperl design.

Input
=====
There are four helper functions which all take a filename/handle/data, and
optional format

File2SequenceIterator - SeqRecord iterator (low memory, forward access only)
File2SequenceList     - List of SeqRecord objects
File2SequenceDict     - Dictionary of SeqRecord objects by record ID
File2Alignment        - Alignment from a multiple sequence alignment file

For non-interlaced files (e.g. Fasta, GenBank, EMBL) with multiple records using
a sequence iterator can save you a lot of memory (RAM).  The saving for interlaced
file formats (e.g. most multiple alignment file formats).  However, you will only
be able to access the records one by one.

These will all invoke the relevant parser with default settings.  You may want
more control, in which case you need to create a sequence iterator directly.

Output
======
There is a single helper function which takes a complete set of records (either
a list, or an iterator) and a filename:

Sequences2File        - Creates a file, writes provided sequences and close file.

For some writers (e.g. phylip and interlaced file formats) if you provide a
SeqRecord iterator, it will be converted into a list.

If you are using a sequential file format, you may want to write out the records
one at a time.  To do this, you will need to create a sequence writer directly.

File Formats
============
When specifying formats, use lowercase strings.

Old Files
=========
The modules Bio.SeqIO.FASTA and Bio.SeqIO.generic are going to be marked depreciated
"""

#TODO
# - define policy on reading aligned sequences with gaps in
#   (e.g. - and . characters) including how the alphabet interacts
#
# - How best to handle unique/non unique record.id when writing.
#   For most file formats reading such files is fine; The stockholm
#   parser would fail.
#
# - EMBL sequence format, ideally combined with GenBank nicely.
#   http://www.bioperl.org/wiki/EMBL_sequence_format
#   Possibly do this in Bio.GenBank?
#   See http://bugzilla.open-bio.org/show_bug.cgi?id=2059#c11
#
# - MSF multiple alignment format, aka GCG, aka PileUp format (*.msf)
#   http://www.bioperl.org/wiki/MSF_multiple_alignment_format 
#
# - Writing NEXUS multiple alignment format (*.nxs)
#   http://www.bioperl.org/wiki/NEXUS_multiple_alignment_format
#   Can be simply offload to Bio.Nexus for this?

"""
FAO BioPython Developers
========================
The way I envision this SeqIO system working as that for any sequence file format
we have an iterator that returns SeqRecord objects.

This also applies to interlaced fileformats (like clustal) where the file cannot
be read record by record.  You should still return an iterator!

These file format specific sequence iterators may be implemented as:
* Classes which take a handle for __init__ and provide the __iter__ method
* Functions that take a handle, and return an iterator object
* Generator functions that take a handle, and yeild SeqRecord objects

It is then trivial to turn this iterator into a list of SeqRecord objects, an in
memory dictionary, or a multiple sequence alignment object.

For building the dictionary by default the id propery of each SeqRecord is used
as the key.  You should always populate the id property, and it should be unique.
For some file formats the accession number is a good choice.

When adding a new file format, please use the same lower case format name as
BioPerl, or if they have not defined one, try the names used by EMBOSS.
"""

import os
#from cStringIO import StringIO
from StringIO import StringIO
from Bio.Alphabet import generic_alphabet, generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Generic import Alignment

from Interfaces import *
import FastaIO
#import EmblGenBankIO
import GenBankIO
import StockholmIO
import ClustalIO
import PhylipIO
import NexusIO

#Convention for format names is "mainname-subtype" in lower case.
#Please use the same names as BioPerl where possible.
#
#Note that this simple system copes with defining
#multiple possible iterators for a given format/extension
#with the -subtype suffix
_ext2format = {"fasta"     : "fasta",
               "faa"       : "fasta", #Used by the NCBI
               "fna"       : "fasta", #Used by the NCBI
               "fnn"       : "fasta", #Used by the NCBI
               "mfasta"    : "fasta", #Used for multiple alignments in Fasta format
               "fa"        : "fasta",
               "genbank"   : "genbank",
               "gbk"       : "genbank", #Used by the NCBI
               "gb"        : "genbank",
               "aln"       : "clustal", #aln is almost always clustal format
               "phy"       : "phylip", #phy is used by clustal
               "phylip"    : "phylip",
               "stockholm" : "stockholm",
               "sth"       : "stockholm", #Used by PFAM (Sanger Inst)
               "pfam"      : "stockholm",
               "nexus"     : "nexus",
               "nxs"       : "nexus", #nxs is used by clustal
               "nex"       : "nexus", #nxs is used by clustal
               }

_format2iterator = {"fasta" : FastaIO.FastaIterator,
                    "genbank" : GenBankIO.GenBankIterator,
                    #See http://bugzilla.open-bio.org/show_bug.cgi?id=2059#c11
                    #"genbank" : EmblGenBankIO.GenBankIterator,
                    #"genbank-cds" : EmblGenBankIO.GenBankCdsFeatureIterator,
                    #"embl" : EmblGenBankIO.EmblIterator, #Not written yet
                    #"embl-cds" : EmblGenBankIO.EmblCdsFeatureIterator,
                    "clustal" : ClustalIO.ClustalIterator,
                    "phylip" : PhylipIO.PhylipIterator,
                    "nexus" : NexusIO.NexusIterator,
                    "stockholm" : StockholmIO.StockholmIterator,
                    }

_format2writer = {"fasta" : FastaIO.FastaWriter,
                  "phylip" : PhylipIO.PhylipWriter,
                  "stockholm" : StockholmIO.StockholmWriter,
                  }

def _filename2format(filename) :
    """Helper function to guess file format based on extension"""
    parts = os.path.basename(filename).split(os.path.extsep)
    assert len(parts) > 1, "No filename extension"
    extension = parts[-1]
    try :
        format = _ext2format[extension.lower()]
    except KeyError :
        assert False, "Unknown extension, " + extension
    #This shouldn't be needed:
    assert format == format.lower().strip()
    return format


def Sequences2File(sequences, filename=None, format=None, handle=None) :
    """Write sequences to a file (and closes the file).

    sequences - A list (or iterator) of SeqRecord objects

    You should also supply at least one of:
    filename - Where to write the records
    handle   - File handle object to write to (which will get closed)

    You are strongly recommended to supply:
    format   - What format to use.  If ommitted, then the filename
               extension will be used to try guess.
    """
    assert filename is not None or handle is not None, \
        "filename or handle is required"
    if not format :
        assert filename is not None, \
           "A filename and/or file format must be supplied"
        #This may raise an exception...
        format =  _filename2format(filename)
    try :
        writer_class = _format2writer[format]
    except KeyError :
        assert False, "Unknown format, " + format

    if handle is None:
        handle = open(filename,"w")
    writer_class(handle).write_file(sequences)
    handle.close() #just in case the writer object forgot
    
def File2SequenceIterator(filename=None, format=None, handle=None, contents=None) :
    """Turns a sequence file into a iterator returning SeqRecords

    You must supply the data using one of the following arguments
    (which are used in preference order: contents, handle, filename)
    
    filename - Path to a local file containing the sequences.
    handle   - Handle to the file.
    contents - String containing full contents of the file.

    format   - String describing the file format.  If omitted,
               then then filename (if given) will be used to
               guess the format.

    Note that file will be parsed with default settings,
    which may result in a generic alphabet or other non-ideal
    settings.  For more control, use the format specific
    iterator directly."""

    assert (filename is not None) \
    or (handle is not None) or (contents is not None), \
    "A filename, file contents or file handle must be supplied"

    if filename is not None :
        assert os.path.isfile(filename), "File not found: " + filename

    if not format :
        assert filename is not None, \
           "A filename and/or file format must be supplied"
        #This may raise an exception...
        format =  _filename2format(filename)
    try :
        iterator_generator = _format2iterator[format]
    except KeyError :
        assert False, "Unknown format, " + format

    if contents is not None :
        #This StringIO object will go out of scope and thus get closed.
        return iterator_generator(StringIO(contents))
    elif handle is not None :
        #Its up to the caller to close this handle - they opened it.
        return iterator_generator(handle)
    else :
        #TODO - I can't see a nice way to explicitly close this handle...
        return iterator_generator(open(filename,"rU"))

def File2SequenceList(filename=None, format=None, handle=None, contents=None) :
    """Turns a sequence file into a list of SeqRecords

    See File2SequenceIterator for details on the arguments."""
    iterator = File2SequenceIterator(filename=filename, format=format,
                                     handle=handle, contents=contents)
    return SequenceList(iterator)

def File2SequenceDict(filename=None, format=None, handle=None, contents=None, record2key=None) :
    """Turns a sequence file into a dictionary of SeqRecords

    If no function record2key is provided, then each record's
    id is used as its key.

    See File2SequenceIterator for details on the other four arguments."""
    iterator = File2SequenceIterator(filename=filename, format=format,
                                     handle=handle, contents=contents)
    return SequenceDict(iterator, record2key)

def Iter2Alignment(SeqIterator, alphabet=generic_alphabet, strict=True) :
    """Returns a multiple sequence alignment

    alphabet - Optional alphabet.  Stongly recommended.
    strict   - Optional, defaults to True.  Should error checking
               be done?
    """
    alignment_length = None
    alignment = Alignment(alphabet)
    for record in iterator :
        if strict :
            if alignment_length is None :
                alignment_length = len(record.seq)
            elif alignment_length <> len(record.seq) :
                raise ValueError, "Sequences of different lengths"
            
            #ToDo, check alphabet for this sequence matches that
            #specified for the alignment.  Not sure how the
            #alphabet.contains() method is intended to be used,
            #but it doesn't make sense to me right now.
            
        #This is abusing the "private" records list,
        #we should really have a method like add_sequence
        #but which takes SeqRecord objects.
        alignment._records.append(record)
    return alignment

def File2Alignment(filename=None, format=None, handle=None, contents=None,
                   alphabet=generic_alphabet, strict=True) :
    """Returns a multiple sequence alignment

    alphabet - Optional alphabet.  Stongly recommended.

    See File2SequenceIterator for details on the other four arguments."""
    iterator = File2SequenceIterator(filename=filename, format=format,
                                     handle=handle, contents=contents)
    return Iter2Alignment(iterator, alphabet=alphabet)
           
if __name__ == "__main__" :
    #Run some tests...
    from Bio.Alphabet import generic_nucleotide
    from sets import Set
    
    # Fasta file with unusual lay out, from here:
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

    #ToDo - Check alphabet, or at least DNA/amino acid, or those
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
         (faa_example,  "fasta",     8, "HISJ_E_COLI", 
          'mkklvlslslvlafssataafaaipqnirigtdptyapfesknsqgelvgfdidlakelc' + \
          'krintqctfvenpldalipslkakkidaimsslsitekrqqeiaftdklyaadsrlvvak' + \
          'nsdiqptveslkgkrvgvlqgttqetfgnehwapkgieivsyqgqdniysdltagridaa' + \
          'fqdevaasegflkqpvgkdykfggpsvkdeklfgvgtgmglrkednelrealnkafaemr' + \
          'adgtyeklakkyfdfdvygg', True),
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
        
        print "%s file with %i records" % (format, rec_count)
        
        print "File2SequenceIterator(contents of file)"

        #Basic check, turning the iterator into a list...
        #This uses "for x in iterator" interally.
        iterator = File2SequenceIterator(contents=data, format=format)
        as_list = list(iterator)
        assert len(as_list) == rec_count
        assert as_list[-1].id == last_id
        if last_seq :
            assert as_list[-1].seq.tostring() == last_seq

        #Test iteration including use of the next() method and "for x in iterator"
        iterator = File2SequenceIterator(contents=data, format=format)
        count = 1
        record = iterator.next()
        assert record is not None
        assert str(record.__class__) == "Bio.SeqRecord.SeqRecord"
        #print record
        for record in iterator :
            assert record.id == as_list[count].id
            assert record.seq.tostring() == as_list[count].seq.tostring()
            count = count + 1
        assert count == rec_count
        assert record is not None
        assert record.id == last_id

        #Test iteration using just next() method
        iterator = File2SequenceIterator(contents=data, format=format)
        count = 0
        while True :
            try :
                record = iterator.next()
            except StopIteration :
                break
            if record is None : break
            assert record.id == as_list[count].id
            assert record.seq.tostring() == as_list[count].seq.tostring()
            count=count+1
        assert count == rec_count

        print "File2SequenceIterator(handle to file)"
        iterator = File2SequenceIterator(handle=StringIO(data), format=format)
        for (i, record) in enumerate(iterator) :
            assert record.id == as_list[i].id
            assert record.seq.tostring() == as_list[i].seq.tostring()            
        assert i+1 == rec_count

        print "File2SequenceIterator(handle to empty file)"
        iterator = File2SequenceIterator(handle=StringIO(""), format=format)
        assert len(list(iterator))==0

        print "File2SequenceIterator(empty contents)"
        iterator = File2SequenceIterator(contents="", format=format)
        assert len(list(iterator))==0

        print "File2SequenceList(contents of file)"
        seq_list = File2SequenceList(contents=data, format=format)
        assert len(seq_list) == len(as_list)
        assert [r.id for r in seq_list] == [r.id for r in as_list]
        for i in range(0, rec_count) :
            seq_list[i].seq.tostring() == as_list[i].seq.tostring()

        if dict_check :
            print "File2SequenceDict(contents of file)"
            seq_dict = File2SequenceDict(contents=data, format=format)
            assert Set(seq_dict.keys()) == Set([r.id for r in as_list])
            assert last_id in seq_dict
            assert seq_dict[last_id].seq.tostring() == as_list[-1].seq.tostring()

        print
        
    print "Checking phy <-> aln examples agree using File2SequenceList"
    #Only compare the first 10 characters of the record.id as they
    #are truncated in the phylip file.  Cannot use File2SequenceDict
    #on the phylip file as there is a repeared id.
    aln_list = File2SequenceList(contents=aln_example, format="clustal")
    phy_list = File2SequenceList(contents=phy_example, format="phylip")
    assert len(aln_list) == len(phy_list)
    assert Set([r.id[0:10] for r in aln_list]) == Set([r.id for r in phy_list])
    for i in range(0, len(aln_list)) :
        assert aln_list[i].id[0:10] == phy_list[i].id
        assert aln_list[i].seq.tostring() == phy_list[i].seq.tostring()
        
    print "Checking nxs <-> aln examples agree using File2SequenceIterator"
    #Only compare the first 10 characters of the record.id as they
    #are truncated in the phylip file.  Cannot use File2SequenceDict
    #on the phylip file as there is a repeared id.
    aln_iter = File2SequenceIterator(contents=aln_example, format="clustal")
    nxs_iter = File2SequenceIterator(contents=nxs_example, format="nexus")
    while True :
        try :
            aln_record = aln_iter.next()
        except StopIteration :
            aln_record = None
        try :
            nxs_record = nxs_iter.next()
        except StopIteration :
            nxs_record = None
        if aln_record is None or nxs_record is None :
            assert aln_record is None
            assert nxs_record is None
            break
        assert aln_record.id == nxs_record.id
        assert aln_record.seq.tostring() == nxs_record.seq.tostring()
    
    print "Checking faa <-> aln examples agree using File2SequenceDict"
    #In my examples, aln_example is an alignment of faa_example
    aln_dict = File2SequenceDict(contents=aln_example, format="clustal")
    faa_dict = File2SequenceDict(contents=faa_example, format="fasta")

    ids = Set(aln_dict.keys())
    assert ids == Set(faa_dict.keys())

    for id in ids :
        #The aln file contains gaps as "-", and this fasta file does not
        assert aln_dict[id].seq.tostring().upper().replace("-","") == \
               faa_dict[id].seq.tostring().upper()

    print
    print "#########################################################"
    print "# Sequence Output Tests                                 #"
    print "#########################################################"
    print

    general_output_formats = ["fasta"]
    alignment_formats = ["phylip","stockholm"]
    for (in_data, in_format, rec_count, last_id, last_seq, unique_ids) in tests:
        if unique_ids :
            in_list =  File2SequenceList(contents=in_data, format=in_format)
            seq_lengths = [len(r.seq) for r in in_list]
            output_formats = general_output_formats[:]
            if min(seq_lengths)==max(seq_lengths) :
                output_formats.extend(alignment_formats)
                print "Checking conversion from %s (including to alignment formats)" % in_format
            else :
                print "Checking conversion from %s (excluding alignment formats)" % in_format
            for out_format in output_formats :
                print "Converting %s iterator -> %s" % (in_format, out_format)
                output = open("temp.txt","w")
                iterator = File2SequenceIterator(contents=in_data, format=in_format)
                #I am using an iterator here deliberately, as some format
                #writers (e.g. phylip and stockholm) will have to cope with
                #this and get the record count).
                Sequences2File(iterator, handle=output, format=out_format)
                output.close()

                print "Checking %s <-> %s" % (in_format, out_format)
                out_list = File2SequenceList(filename="temp.txt", format=out_format)

                assert rec_count == len(out_list)
                if last_seq :
                    assert last_seq == out_list[-1].seq.tostring()
                if out_format=="phylip" :
                    assert last_id[0:10] == out_list[-1].id
                else :
                    assert last_id == out_list[-1].id

                for i in range(0, rec_count) :
                    assert in_list[-1].seq.tostring() == out_list[-1].seq.tostring()
                    if out_format=="phylip" :
                        assert in_list[i].id[0:10] == out_list[i].id
                    else :
                        assert in_list[i].id == out_list[i].id
            print

    print "#########################################################"
    print "# SeqIO Tests finished                                  #"
    print "#########################################################"


    in_list = File2SequenceList(r"c:\temp\nexus_etc\example.nexus")
