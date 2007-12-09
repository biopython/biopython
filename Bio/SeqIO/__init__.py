# Copyright 2006-2007 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
#Nice link:
# http://www.ebi.ac.uk/help/formats_frame.html

"""Sequence input/output designed to look similar to the bioperl design.

The Bio.SeqIO module is also documented by a whole chapter in the Biopython
tutorial, and by this wiki on the webage, http://biopython.org/wiki/SeqIO

Input
=====
The main function is Bio.SeqIO.parse(...) which takes an input file handle,
and format string.  This returns an iterator giving SeqRecord objects.

    from Bio import SeqIO
    handle = open("example.fasta", "rU")
    for record in SeqIO.parse(handle, "fasta") :
        print record
    handle.close()

Note that the parse() function will all invoke the relevant parser for
the format with its default settings.  You may want more control, in which case
you need to create a format specific sequence iterator directly.

For non-interlaced files (e.g. Fasta, GenBank, EMBL) with multiple records
using a sequence iterator can save you a lot of memory (RAM).  There is less
benefit for interlaced file formats (e.g. most multiple alignment file formats).
However, an iterator only lets you access the records one by one.

If you want random access to the records by number, turn this into a list:

    from Bio import SeqIO
    handle = open("example.fasta", "rU")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    print records[0]

If you want random access to the records by a key such as the record id, turn
the iterator into a dictionary:

    from Bio import SeqIO
    handle = open("example.fasta", "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    print record["gi:12345678"]

If you expect your file to contain one-and-only-one record, then we provide
the following 'helper' function which will return a single SeqRecord, or raise
an exception if there are no records or more than one record:

    from Bio import SeqIO
    handle = open("example.fasta", "rU")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    print record

This style is useful when you expect a single record only (and would consider
multiple records an error).  For example, when dealing with GenBank files for
bacterial genomes or chromosomes, there is normally only a single record.
Alternatively, use this with a handle when download a single record from the
internet.

However, if you just want the first record from a file containing multiple
record, use the iterator's next() method:

    from Bio import SeqIO
    handle = open("example.fasta", "rU")
    record = SeqIO.parse(handle, "fasta").next()
    handle.close()
    print record

The above code will work as long as the file contains at least one record.
Note that if there is more than one record, the remaining records will be
silently ignored.

Input - Alignments
==================
Currently an alignment class cannot be created from SeqRecord objects.
Instead, use the to_alignment(...) function, like so:

    from Bio import SeqIO
    handle = open("example.aln", "rU")
    alignment = SeqIO.to_alignment(SeqIO.parse(handle, "clustal"))
    handle.close()

This function may be removed in future once alignments can be created
directly from SeqRecord objects.

Output
======
Use the function Bio.SeqIO.write(...), which takes a complete set of SeqRecord
objects (either as a list, or an iterator), an output file handle and of course
the file format.

    from Bio import SeqIO
    records = ...
    handle = open("example.faa", "w")
    SeqIO.write(records, handle, "fasta")
    handle.close()

In general, you are expected to call this function once (with all your records)
and then close the file handle.

Output - Advanced
=================
The effect of calling write() multiple times on a single file will vary
depending on the file format, and is best avoided unless you have a strong reason
to do so.

Trying this for certain alignment formats (e.g. phylip, clustal, stockholm) would
have the effect of concatenating several multiple sequence alignments together.
Such files are created by the PHYLIP suite of programs for bootstrap analysis.

For sequential files formats (e.g. fasta, genbank) each "record block" holds a
single sequence.  For these files it would probably be safe to call write()
multiple times.

File Formats
============
When specifying formats, use lowercase strings.

Old Files
=========
The modules Bio.SeqIO.FASTA and Bio.SeqIO.generic are depreciated and may be
removed.
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

import FastaIO
import InsdcIO #EMBL and GenBank
import StockholmIO
import ClustalIO
import PhylipIO
import NexusIO
import SwissIO

#Convention for format names is "mainname-subtype" in lower case.
#Please use the same names as BioPerl where possible.
#
#Note that this simple system copes with defining
#multiple possible iterators for a given format/extension
#with the -subtype suffix

_FormatToIterator ={"fasta" : FastaIO.FastaIterator,
                    "genbank" : InsdcIO.GenBankIterator,
                    "genbank-cds" : InsdcIO.GenBankCdsFeatureIterator,
                    "embl" : InsdcIO.EmblIterator,
                    "embl-cds" : InsdcIO.EmblCdsFeatureIterator,
                    "clustal" : ClustalIO.ClustalIterator,
                    "phylip" : PhylipIO.PhylipIterator,
                    "nexus" : NexusIO.NexusIterator,
                    "stockholm" : StockholmIO.StockholmIterator,
                    "swiss" : SwissIO.SwissIterator,
                    }

_FormatToWriter ={"fasta" : FastaIO.FastaWriter,
                  "phylip" : PhylipIO.PhylipWriter,
                  "stockholm" : StockholmIO.StockholmWriter,
                  "clustal" : ClustalIO.ClustalWriter,
                  }

def write(sequences, handle, format) :
    """Write complete set of sequences to a file

    sequences - A list (or iterator) of SeqRecord objects
    handle    - File handle object to write to
    format    - What format to use.

    You should close the handle after calling this function.

    There is no return value.
    """

    #Try and give helpful error messages:
    if isinstance(handle, basestring) :
        raise TypeError("Need a file handle, not a string (i.e. not a filename)")
    if not isinstance(format, basestring) :
        raise TypeError("Need a string for the file format (lower case)")
    if not format :
        raise ValueError("Format required (lower case string)")
    if format <> format.lower() :
        raise ValueError("Format string '%s' should be lower case" % format)
    if isinstance(sequences,SeqRecord):
        raise ValueError("Use a SeqRecord list/iterator, not just a single SeqRecord")

    #Map the file format to a writer class
    try :
        writer_class = _FormatToWriter[format]
    except KeyError :
        raise ValueError("Unknown format '%s'" % format)

    writer_class(handle).write_file(sequences)
    #Don't close the file, as that would prevent things like
    #creating concatenated phylip files for bootstrapping.
    #handle.close()
    return
    
def parse(handle, format) :
    """Turns a sequence file into a iterator returning SeqRecords.

    handle   - handle to the file.
    format   - string describing the file format.

    If you have the file name in a string 'filename', use:

    from Bio import SeqIO
    my_iterator = SeqIO.parse(open(filename,"rU"), format)

    If you have a string 'data' containing the file contents, use:

    from Bio import SeqIO
    from StringIO import StringIO
    my_iterator = SeqIO.parse(StringIO(data), format)

    Note that file will be parsed with default settings,
    which may result in a generic alphabet or other non-ideal
    settings.  For more control, you must use the format specific
    iterator directly...

    Use the Bio.SeqIO.read(handle, format) function when you expect
    a single record only.
    """

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
    try :
        iterator_generator = _FormatToIterator[format]
    except KeyError :
        raise ValueError("Unknown format '%s'" % format)

    #Its up to the caller to close this handle - they opened it.
    return iterator_generator(handle)

def read(handle, format) :
    """Turns a sequence file into a single SeqRecord.

    handle   - handle to the file.
    format   - string describing the file format.

    If the handle contains no records, or more than one record,
    an exception is raised.  For example, using a GenBank file
    containing one record:

    from Bio import SeqIO
    record = SeqIO.read(open("example.gbk"), "genbank")

    If however you want the first record from a file containing,
    multiple records this function would raise an exception.
    Instead use:

    from Bio import SeqIO
    record = SeqIO.parse(open("example.gbk"), "genbank").next()

    Use the Bio.SeqIO.parse(handle, format) function if you want
    to read multiple records from the handle.
    """
    iterator = parse(handle, format)
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
    return first

def to_dict(sequences, key_function=None) :
    """Turns a sequence iterator or list into a dictionary.

    sequences  - An iterator that returns SeqRecord objects,
                 or simply a list of SeqRecord objects.
    key_function - Optional function which when given a SeqRecord
                   returns a unique string for the dictionary key.

    e.g. key_function = lambda rec : rec.name
    or,  key_function = lambda rec : rec.description.split()[0]

    If key_function is ommitted then record.id is used, on the
    assumption that the records objects returned are SeqRecords
    with a unique id field.

    If there are duplicate keys, an error is raised.

    Example usage:

    from Bio import SeqIO
    filename = "example.fasta"
    d = SeqIO.to_dict(SeqIO.parse(open(faa_filename, "rU")),
        key_function = lambda rec : rec.description.split()[0])
    print len(d)
    print d.keys()[0:10]
    key = d.keys()[0]
    print d[key]
    """    
    if key_function is None :
        key_function = lambda rec : rec.id

    d = dict()
    for record in sequences :
        key = key_function(record)
        if key in d :
            raise ValueError("Duplicate key '%s'" % key)
        d[key] = record
    return d

def to_alignment(sequences, alphabet=generic_alphabet, strict=True) :
    """Returns a multiple sequence alignment.

    sequences -An iterator that returns SeqRecord objects,
               or simply a list of SeqRecord objects.
               All the record sequences must be the same length.
    alphabet - Optional alphabet.  Stongly recommended.
    strict   - Optional, defaults to True.  Should error checking
               be done?
    """
    #TODO - Move this functionality into the Alignment class instead?
    alignment_length = None
    alignment = Alignment(alphabet)
    for record in sequences :
        if strict :
            if alignment_length is None :
                alignment_length = len(record.seq)
            elif alignment_length <> len(record.seq) :
                raise ValueError("Sequences of different lengths")
            
            if not isinstance(record.seq.alphabet, alphabet.__class__) :
                raise ValueError("Incompatible sequence alphabet")
            
            #ToDo, additional checks on the specified alignment...
            #Should we look at the alphabet.contains() method?
            
        #This is abusing the "private" records list,
        #we should really have a method like add_sequence
        #but which takes SeqRecord objects.  See also Bug 1944
        alignment._records.append(record)
    return alignment
           
if __name__ == "__main__" :
    #Run some tests...
    from Bio.Alphabet import generic_nucleotide
    from sets import Set
    
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


    swiss_example = \
"""ID   104K_THEAN              Reviewed;         893 AA.
AC   Q4U9M9;
DT   18-APR-2006, integrated into UniProtKB/Swiss-Prot.
DT   05-JUL-2005, sequence version 1.
DT   31-OCT-2006, entry version 8.
DE   104 kDa microneme-rhoptry antigen precursor (p104).
GN   ORFNames=TA08425;
OS   Theileria annulata.
OC   Eukaryota; Alveolata; Apicomplexa; Piroplasmida; Theileriidae;
OC   Theileria.
OX   NCBI_TaxID=5874;
RN   [1]
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RC   STRAIN=Ankara;
RX   PubMed=15994557; DOI=10.1126/science.1110418;
RA   Pain A., Renauld H., Berriman M., Murphy L., Yeats C.A., Weir W.,
RA   Kerhornou A., Aslett M., Bishop R., Bouchier C., Cochet M.,
RA   Coulson R.M.R., Cronin A., de Villiers E.P., Fraser A., Fosker N.,
RA   Gardner M., Goble A., Griffiths-Jones S., Harris D.E., Katzer F.,
RA   Larke N., Lord A., Maser P., McKellar S., Mooney P., Morton F.,
RA   Nene V., O'Neil S., Price C., Quail M.A., Rabbinowitsch E.,
RA   Rawlings N.D., Rutter S., Saunders D., Seeger K., Shah T., Squares R.,
RA   Squares S., Tivey A., Walker A.R., Woodward J., Dobbelaere D.A.E.,
RA   Langsley G., Rajandream M.A., McKeever D., Shiels B., Tait A.,
RA   Barrell B.G., Hall N.;
RT   "Genome of the host-cell transforming parasite Theileria annulata
RT   compared with T. parva.";
RL   Science 309:131-133(2005).
CC   -!- SUBCELLULAR LOCATION: Cell membrane; lipid-anchor; GPI-anchor
CC       (Potential). In microneme/rhoptry complexes (By similarity).
DR   EMBL; CR940353; CAI76474.1; -; Genomic_DNA.
DR   InterPro; IPR007480; DUF529.
DR   Pfam; PF04385; FAINT; 4.
KW   Complete proteome; GPI-anchor; Lipoprotein; Membrane; Repeat; Signal;
KW   Sporozoite.
FT   SIGNAL        1     19       Potential.
FT   CHAIN        20    873       104 kDa microneme-rhoptry antigen.
FT                                /FTId=PRO_0000232680.
FT   PROPEP      874    893       Removed in mature form (Potential).
FT                                /FTId=PRO_0000232681.
FT   COMPBIAS    215    220       Poly-Leu.
FT   COMPBIAS    486    683       Lys-rich.
FT   COMPBIAS    854    859       Poly-Arg.
FT   LIPID       873    873       GPI-anchor amidated aspartate
FT                                (Potential).
SQ   SEQUENCE   893 AA;  101921 MW;  2F67CEB3B02E7AC1 CRC64;
     MKFLVLLFNI LCLFPILGAD ELVMSPIPTT DVQPKVTFDI NSEVSSGPLY LNPVEMAGVK
     YLQLQRQPGV QVHKVVEGDI VIWENEEMPL YTCAIVTQNE VPYMAYVELL EDPDLIFFLK
     EGDQWAPIPE DQYLARLQQL RQQIHTESFF SLNLSFQHEN YKYEMVSSFQ HSIKMVVFTP
     KNGHICKMVY DKNIRIFKAL YNEYVTSVIG FFRGLKLLLL NIFVIDDRGM IGNKYFQLLD
     DKYAPISVQG YVATIPKLKD FAEPYHPIIL DISDIDYVNF YLGDATYHDP GFKIVPKTPQ
     CITKVVDGNE VIYESSNPSV ECVYKVTYYD KKNESMLRLD LNHSPPSYTS YYAKREGVWV
     TSTYIDLEEK IEELQDHRST ELDVMFMSDK DLNVVPLTNG NLEYFMVTPK PHRDIIIVFD
     GSEVLWYYEG LENHLVCTWI YVTEGAPRLV HLRVKDRIPQ NTDIYMVKFG EYWVRISKTQ
     YTQEIKKLIK KSKKKLPSIE EEDSDKHGGP PKGPEPPTGP GHSSSESKEH EDSKESKEPK
     EHGSPKETKE GEVTKKPGPA KEHKPSKIPV YTKRPEFPKK SKSPKRPESP KSPKRPVSPQ
     RPVSPKSPKR PESLDIPKSP KRPESPKSPK RPVSPQRPVS PRRPESPKSP KSPKSPKSPK
     VPFDPKFKEK LYDSYLDKAA KTKETVTLPP VLPTDESFTH TPIGEPTAEQ PDDIEPIEES
     VFIKETGILT EEVKTEDIHS ETGEPEEPKR PDSPTKHSPK PTGTHPSMPK KRRRSDGLAL
     STTDLESEAG RILRDPTGKI VTMKRSKSFD DLTTVREKEH MGAEIRKIVV DDDGTEADDE
     DTHPSKEKHL STVRRRRPRP KKSSKSSKPR KPDSAFVPSI IFIFLVSLIV GIL
//
ID   104K_THEPA              Reviewed;         924 AA.
AC   P15711; Q4N2B5;
DT   01-APR-1990, integrated into UniProtKB/Swiss-Prot.
DT   01-APR-1990, sequence version 1.
DT   31-OCT-2006, entry version 31.
DE   104 kDa microneme-rhoptry antigen precursor (p104).
GN   OrderedLocusNames=TP04_0437;
OS   Theileria parva.
OC   Eukaryota; Alveolata; Apicomplexa; Piroplasmida; Theileriidae;
OC   Theileria.
OX   NCBI_TaxID=5875;
RN   [1]
RP   NUCLEOTIDE SEQUENCE [GENOMIC DNA].
RC   STRAIN=Muguga;
RX   MEDLINE=90158697; PubMed=1689460; DOI=10.1016/0166-6851(90)90007-9;
RA   Iams K.P., Young J.R., Nene V., Desai J., Webster P., Ole-Moiyoi O.K.,
RA   Musoke A.J.;
RT   "Characterisation of the gene encoding a 104-kilodalton microneme-
RT   rhoptry protein of Theileria parva.";
RL   Mol. Biochem. Parasitol. 39:47-60(1990).
RN   [2]
RP   NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA].
RC   STRAIN=Muguga;
RX   PubMed=15994558; DOI=10.1126/science.1110439;
RA   Gardner M.J., Bishop R., Shah T., de Villiers E.P., Carlton J.M.,
RA   Hall N., Ren Q., Paulsen I.T., Pain A., Berriman M., Wilson R.J.M.,
RA   Sato S., Ralph S.A., Mann D.J., Xiong Z., Shallom S.J., Weidman J.,
RA   Jiang L., Lynn J., Weaver B., Shoaibi A., Domingo A.R., Wasawo D.,
RA   Crabtree J., Wortman J.R., Haas B., Angiuoli S.V., Creasy T.H., Lu C.,
RA   Suh B., Silva J.C., Utterback T.R., Feldblyum T.V., Pertea M.,
RA   Allen J., Nierman W.C., Taracha E.L.N., Salzberg S.L., White O.R.,
RA   Fitzhugh H.A., Morzaria S., Venter J.C., Fraser C.M., Nene V.;
RT   "Genome sequence of Theileria parva, a bovine pathogen that transforms
RT   lymphocytes.";
RL   Science 309:134-137(2005).
CC   -!- SUBCELLULAR LOCATION: Cell membrane; lipid-anchor; GPI-anchor
CC       (Potential). In microneme/rhoptry complexes.
CC   -!- DEVELOPMENTAL STAGE: Sporozoite antigen.
DR   EMBL; M29954; AAA18217.1; -; Unassigned_DNA.
DR   EMBL; AAGK01000004; EAN31789.1; -; Genomic_DNA.
DR   PIR; A44945; A44945.
DR   InterPro; IPR007480; DUF529.
DR   Pfam; PF04385; FAINT; 4.
KW   Complete proteome; GPI-anchor; Lipoprotein; Membrane; Repeat; Signal;
KW   Sporozoite.
FT   SIGNAL        1     19       Potential.
FT   CHAIN        20    904       104 kDa microneme-rhoptry antigen.
FT                                /FTId=PRO_0000046081.
FT   PROPEP      905    924       Removed in mature form (Potential).
FT                                /FTId=PRO_0000232679.
FT   COMPBIAS    508    753       Pro-rich.
FT   COMPBIAS    880    883       Poly-Arg.
FT   LIPID       904    904       GPI-anchor amidated aspartate
FT                                (Potential).
SQ   SEQUENCE   924 AA;  103626 MW;  289B4B554A61870E CRC64;
     MKFLILLFNI LCLFPVLAAD NHGVGPQGAS GVDPITFDIN SNQTGPAFLT AVEMAGVKYL
     QVQHGSNVNI HRLVEGNVVI WENASTPLYT GAIVTNNDGP YMAYVEVLGD PNLQFFIKSG
     DAWVTLSEHE YLAKLQEIRQ AVHIESVFSL NMAFQLENNK YEVETHAKNG ANMVTFIPRN
     GHICKMVYHK NVRIYKATGN DTVTSVVGFF RGLRLLLINV FSIDDNGMMS NRYFQHVDDK
     YVPISQKNYE TGIVKLKDYK HAYHPVDLDI KDIDYTMFHL ADATYHEPCF KIIPNTGFCI
     TKLFDGDQVL YESFNPLIHC INEVHIYDRN NGSIICLHLN YSPPSYKAYL VLKDTGWEAT
     THPLLEEKIE ELQDQRACEL DVNFISDKDL YVAALTNADL NYTMVTPRPH RDVIRVSDGS
     EVLWYYEGLD NFLVCAWIYV SDGVASLVHL RIKDRIPANN DIYVLKGDLY WTRITKIQFT
     QEIKRLVKKS KKKLAPITEE DSDKHDEPPE GPGASGLPPK APGDKEGSEG HKGPSKGSDS
     SKEGKKPGSG KKPGPAREHK PSKIPTLSKK PSGPKDPKHP RDPKEPRKSK SPRTASPTRR
     PSPKLPQLSK LPKSTSPRSP PPPTRPSSPE RPEGTKIIKT SKPPSPKPPF DPSFKEKFYD
     DYSKAASRSK ETKTTVVLDE SFESILKETL PETPGTPFTT PRPVPPKRPR TPESPFEPPK
     DPDSPSTSPS EFFTPPESKR TRFHETPADT PLPDVTAELF KEPDVTAETK SPDEAMKRPR
     SPSEYEDTSP GDYPSLPMKR HRLERLRLTT TEMETDPGRM AKDASGKPVK LKRSKSFDDL
     TTVELAPEPK ASRIVVDDEG TEADDEETHP PEERQKTEVR RRRPPKKPSK SPRPSKPKKP
     KKPDSAYIPS ILAILVVSLI VGIL
//
ID   108_SOLLC               Reviewed;         102 AA.
AC   Q43495;
DT   15-JUL-1999, integrated into UniProtKB/Swiss-Prot.
DT   01-NOV-1996, sequence version 1.
DT   31-OCT-2006, entry version 37.
DE   Protein 108 precursor.
OS   Solanum lycopersicum (Tomato) (Lycopersicon esculentum).
OC   Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
OC   Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons;
OC   asterids; lamiids; Solanales; Solanaceae; Solanum; Lycopersicon.
OX   NCBI_TaxID=4081;
RN   [1]
RP   NUCLEOTIDE SEQUENCE [MRNA].
RC   STRAIN=cv. VF36; TISSUE=Anther;
RX   MEDLINE=94143497; PubMed=8310077; DOI=10.1104/pp.101.4.1413;
RA   Chen R., Smith A.G.;
RT   "Nucleotide sequence of a stamen- and tapetum-specific gene from
RT   Lycopersicon esculentum.";
RL   Plant Physiol. 101:1413-1413(1993).
CC   -!- TISSUE SPECIFICITY: Stamen- and tapetum-specific.
CC   -!- SIMILARITY: Belongs to the A9/FIL1 family.
DR   EMBL; Z14088; CAA78466.1; -; mRNA.
DR   PIR; S26409; S26409.
DR   InterPro; IPR013770; LPT_helical.
DR   InterPro; IPR003612; LTP/seed_store/tryp_amyl_inhib.
DR   Pfam; PF00234; Tryp_alpha_amyl; 1.
DR   SMART; SM00499; AAI; 1.
KW   Signal.
FT   SIGNAL        1     30       Potential.
FT   CHAIN        31    102       Protein 108.
FT                                /FTId=PRO_0000000238.
FT   DISULFID     41     77       By similarity.
FT   DISULFID     51     66       By similarity.
FT   DISULFID     67     92       By similarity.
FT   DISULFID     79     99       By similarity.
SQ   SEQUENCE   102 AA;  10576 MW;  CFBAA1231C3A5E92 CRC64;
     MASVKSSSSS SSSSFISLLL LILLVIVLQS QVIECQPQQS CTASLTGLNV CAPFLVPGSP
     TASTECCNAV QSINHDCMCN TMRIAAQIPA QCNLPPLSCS AN
//
"""

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
         (gbk_example, "genbank-cds", 3, "AAA98667.1",
          'MNRWVEKWLRVYLKCYINLILFYRNVYPPQSFDYTTYQSFNLPQFVPINRHPALIDYIEE' + \
          'LILDVLSKLTHVYRFSICIINKKNDLCIEKYVLDFSELQHVDKDDQIITETEVFDEFRSS' + \
          'LNSLIMHLEKLPKVNDDTITFEAVINAIELELGHKLDRNRRVDSLEEKAEIERDSNWVKC' + \
          'QEDENLPDNNGFQPPKIKLTSLVGSDVGPLIIHQFSEKLISGDDKILNGVYSQYEEGESI' + \
          'FGSLF', True),
          (swiss_example,"swiss", 3, "Q43495",
          "MASVKSSSSSSSSSFISLLLLILLVIVLQSQVIECQPQQSCTASLTGLNVCAPFLVPGSP" + \
          "TASTECCNAVQSINHDCMCNTMRIAAQIPAQCNLPPLSCSAN", True),
    ]
    
    for (data, format, rec_count, last_id, last_seq, dict_check) in tests:
        
        print "%s file with %i records" % (format, rec_count)
        
        print "Bio.SeqIO.parse(handle)"

        #Basic check, turning the iterator into a list...
        #This uses "for x in iterator" interally.
        iterator = parse(StringIO(data), format=format)
        as_list = list(iterator)
        assert len(as_list) == rec_count, \
            "Expected %i records, found %i" \
            % (rec_count, len(as_list))
        assert as_list[-1].id == last_id, \
            "Expected '%s' as last record ID, found '%s'" \
            % (last_id, as_list[-1].id)
        if last_seq :
            assert as_list[-1].seq.tostring() == last_seq

        #Test iteration including use of the next() method and "for x in iterator"
        iterator = parse(StringIO(data), format=format)
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
        iterator = parse(StringIO(data), format=format)
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

        print "parse(...)"
        iterator = parse(StringIO(data), format=format)
        for (i, record) in enumerate(iterator) :
            assert record.id == as_list[i].id
            assert record.seq.tostring() == as_list[i].seq.tostring()            
        assert i+1 == rec_count

        print "parse(handle to empty file)"
        iterator = parse(StringIO(""), format=format)
        assert len(list(iterator))==0

        if dict_check :
            print "to_dict(parse(...))"
            seq_dict = to_dict(parse(StringIO(data), format=format))
            assert Set(seq_dict.keys()) == Set([r.id for r in as_list])
            assert last_id in seq_dict
            assert seq_dict[last_id].seq.tostring() == as_list[-1].seq.tostring()

        if len(Set([len(r.seq) for r in as_list]))==1 :
            #All the sequences in the example are the same length,
            #so it make sense to try turning this file into an alignment.
            print "to_alignment(parse(handle))"
            alignment = to_alignment(parse(handle = StringIO(data), format=format))
            assert len(alignment._records)==rec_count
            assert alignment.get_alignment_length() == len(as_list[0].seq)
            for i in range(0, rec_count) :
                assert as_list[i].id == alignment._records[i].id
                assert as_list[i].id == alignment.get_all_seqs()[i].id
                assert as_list[i].seq.tostring() == alignment._records[i].seq.tostring()
                assert as_list[i].seq.tostring() == alignment.get_all_seqs()[i].seq.tostring()

        print "read(...)"
        if rec_count == 1 :
            record = read(StringIO(data), format)
            assert isinstance(record, SeqRecord)
        else :
            try :
                record = read(StringIO(data), format)
                assert False, "Should have failed"
            except ValueError :
                #Expected to fail
                pass

        print
        
    print "Checking phy <-> aln examples agree using list(parse(...))"
    #Only compare the first 10 characters of the record.id as they
    #are truncated in the phylip file.  Cannot use to_dict(parse(...))
    #on the phylip file as there is a repeared id.
    aln_list = list(parse(StringIO(aln_example), format="clustal"))
    phy_list = list(parse(StringIO(phy_example), format="phylip"))
    assert len(aln_list) == len(phy_list)
    assert Set([r.id[0:10] for r in aln_list]) == Set([r.id for r in phy_list])
    for i in range(0, len(aln_list)) :
        assert aln_list[i].id[0:10] == phy_list[i].id
        assert aln_list[i].seq.tostring() == phy_list[i].seq.tostring()
        
    print "Checking nxs <-> aln examples agree using parse"
    #Only compare the first 10 characters of the record.id as they
    #are truncated in the phylip file.  Cannot use to_dict(parse(...))
    #on the phylip file as there is a repeared id.
    aln_iter = parse(StringIO(aln_example), format="clustal")
    nxs_iter = parse(StringIO(nxs_example), format="nexus")
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
    
    print "Checking faa <-> aln examples agree using to_dict(parse(...)"
    #In my examples, aln_example is an alignment of faa_example
    aln_dict = to_dict(parse(StringIO(aln_example), format="clustal"))
    faa_dict = to_dict(parse(StringIO(faa_example), format="fasta"))

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
    alignment_formats = ["phylip","stockholm","clustal"]
    for (in_data, in_format, rec_count, last_id, last_seq, unique_ids) in tests:
        if unique_ids :
            in_list =  list(parse(StringIO(in_data), format=in_format))
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
                iterator = parse(StringIO(in_data), format=in_format)
                #I am using an iterator here deliberately, as some format
                #writers (e.g. phylip and stockholm) will have to cope with
                #this and get the record count.

                try :
                    write(iterator, output, out_format)
                except ValueError, e:
                    print "FAILED: %s" % str(e)
                    #Try next format instead...
                    continue

                output.close()

                print "Checking %s <-> %s" % (in_format, out_format)
                out_list = list(parse(open("temp.txt","rU"), format=out_format))

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
