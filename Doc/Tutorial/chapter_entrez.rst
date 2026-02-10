.. _`chapter:entrez`:

Accessing NCBI’s Entrez databases
=================================

Entrez (https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html) is a data
retrieval system that provides users access to NCBI’s databases such as
PubMed, GenBank, GEO, and many others. You can access Entrez from a web
browser to manually enter queries, or you can use Biopython’s
``Bio.Entrez`` module for programmatic access to Entrez. The latter
allows you for example to search PubMed or download GenBank records from
within a Python script.

The ``Bio.Entrez`` module makes use of the Entrez Programming Utilities
(also known as EUtils), consisting of eight tools that are described in
detail on NCBI’s page at https://www.ncbi.nlm.nih.gov/books/NBK25501/.
Each of these tools corresponds to one Python function in the
``Bio.Entrez`` module, as described in the sections below. This module
makes sure that the correct URL is used for the queries, and that NCBI’s
guidelines for responsible data access are being followed.

The output returned by the Entrez Programming Utilities is typically in
XML format. To parse such output, you have several options:

#. Use ``Bio.Entrez``\ ’s parser to parse the XML output into a Python
   object;

#. Use one of the XML parsers available in Python’s standard library;

#. Read the XML output as raw text, and parse it by string searching and
   manipulation.

See the Python documentation for a description of the XML parsers in
Python’s standard library. Here, we discuss the parser in Biopython’s
``Bio.Entrez`` module. This parser can be used to parse data as provided
through ``Bio.Entrez``\ ’s programmatic access functions to Entrez, but
can also be used to parse XML data from NCBI Entrez that are stored in a
file. In the latter case, the XML file should be opened in binary mode
(e.g. ``open("myfile.xml", "rb")``) for the XML parser in ``Bio.Entrez``
to work correctly. Alternatively, you can pass the file name or path to
the XML file, and let ``Bio.Entrez`` take care of opening and closing
the file.

NCBI uses DTD (Document Type Definition) files to describe the structure
of the information contained in XML files. Most of the DTD files used by
NCBI are included in the Biopython distribution. The ``Bio.Entrez``
parser makes use of the DTD files when parsing an XML file returned by
NCBI Entrez.

Occasionally, you may find that the DTD file associated with a specific
XML file is missing in the Biopython distribution. In particular, this
may happen when NCBI updates its DTD files. If this happens,
``Entrez.read`` will show a warning message with the name and URL of the
missing DTD file. The parser will proceed to access the missing DTD file
through the internet, allowing the parsing of the XML file to continue.
However, the parser is much faster if the DTD file is available locally.
For this purpose, please download the DTD file from the URL in the
warning message and place it in the directory
``...site-packages/Bio/Entrez/DTDs``, containing the other DTD files. If
you don’t have write access to this directory, you can also place the
DTD file in ``~/.biopython/Bio/Entrez/DTDs``, where ``~`` represents
your home directory. Since this directory is read before the directory
``...site-packages/Bio/Entrez/DTDs``, you can also put newer versions of
DTD files there if the ones in ``...site-packages/Bio/Entrez/DTDs``
become outdated. Alternatively, if you installed Biopython from source,
you can add the DTD file to the source code’s ``Bio/Entrez/DTDs``
directory, and reinstall Biopython. This will install the new DTD file
in the correct location together with the other DTD files.

The Entrez Programming Utilities can also generate output in other
formats, such as the Fasta or GenBank file formats for sequence
databases, or the MedLine format for the literature database, discussed
in Section :ref:`sec:entrez-specialized-parsers`.

The functions in ``Bio.Entrez`` for programmatic access to Entrez return
data either in binary format or in text format, depending on the type of
data requested. In most cases, these functions return data in text
format by decoding the data obtained from NCBI Entrez to Python strings
under the assumption that the encoding is UTF-8. However, XML data are
returned in binary format. The reason for this is that the encoding is
specified in the XML document itself, which means that we won’t know the
correct encoding to use until we start parsing the file.
``Bio.Entrez``\ ’s parser therefore accepts data in binary format,
extracts the encoding from the XML, and uses it to decode all text in
the XML document to Python strings, ensuring that all text (in
particular in languages other than English) are interpreted correctly.
This is also the reason why you should open an XML file a binary mode
when you want to use ``Bio.Entrez``\ ’s parser to parse the file.

.. _`sec:entrez-guidelines`:

Entrez Guidelines
-----------------

Before using Biopython to access the NCBI’s online resources (via
``Bio.Entrez`` or some of the other modules), please read the `NCBI’s
Entrez User
Requirements <https://www.ncbi.nlm.nih.gov/books/NBK25497/>`__. If the
NCBI finds you are abusing their systems, they can and will ban your
access!

To paraphrase:

-  For any series of more than 100 requests, do this at weekends or
   outside USA peak times. This is up to you to obey.

-  Use the https://eutils.ncbi.nlm.nih.gov address, not the standard
   NCBI Web address. Biopython uses this web address.

-  If you are using a API key, you can make at most 10 queries per
   second, otherwise at most 3 queries per second. This is automatically
   enforced by Biopython. Include ``api_key="MyAPIkey"`` in the argument
   list or set it as a module level variable:

   .. code:: pycon

      >>> from Bio import Entrez
      >>> Entrez.api_key = "MyAPIkey"

-  Use the optional email parameter so the NCBI can contact you if there
   is a problem. You can either explicitly set this as a parameter with
   each call to Entrez (e.g. include ``email="A.N.Other@example.com"``
   in the argument list), or you can set a global email address:

   .. doctest

   .. code:: pycon

      >>> from Bio import Entrez
      >>> Entrez.email = "A.N.Other@example.com"

   ``Bio.Entrez`` will then use this email address with each call to
   Entrez. The ``example.com`` address is a reserved domain name
   specifically for documentation (RFC 2606). Please DO NOT use a random
   email – it’s better not to give an email at all. The email parameter
   has been mandatory since June 1, 2010. In case of excessive usage,
   NCBI will attempt to contact a user at the e-mail address provided
   prior to blocking access to the E-utilities.

-  If you are using Biopython within some larger software suite, use the
   tool parameter to specify this. You can either explicitly set the
   tool name as a parameter with each call to Entrez (e.g. include
   ``tool="MyLocalScript"`` in the argument list), or you can set a
   global tool name:

   .. doctest

   .. code:: pycon

      >>> from Bio import Entrez
      >>> Entrez.tool = "MyLocalScript"

   The tool parameter will default to Biopython.

-  For large queries, the NCBI also recommend using their session
   history feature (the WebEnv session cookie string, see
   Section :ref:`sec:entrez-webenv`). This is only slightly more
   complicated.

In conclusion, be sensible with your usage levels. If you plan to
download lots of data, consider other options. For example, if you want
easy access to all the human genes, consider fetching each chromosome by
FTP as a GenBank file, and importing these into your own BioSQL database
(see Section :ref:`sec:BioSQL`).

.. _`sec:entrez-einfo`:

EInfo: Obtaining information about the Entrez databases
-------------------------------------------------------

EInfo provides field index term counts, last update, and available links
for each of NCBI’s databases. In addition, you can use EInfo to obtain a
list of all database names accessible through the Entrez utilities:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.einfo()
   >>> result = stream.read()
   >>> stream.close()

The variable ``result`` now contains a list of databases in XML format:

.. code:: pycon

   >>> print(result)
   <?xml version="1.0"?>
   <!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD eInfoResult, 11 May 2002//EN"
    "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/eInfo_020511.dtd">
   <eInfoResult>
   <DbList>
           <DbName>pubmed</DbName>
           <DbName>protein</DbName>
           <DbName>nucleotide</DbName>
           <DbName>nuccore</DbName>
           <DbName>nucgss</DbName>
           <DbName>nucest</DbName>
           <DbName>structure</DbName>
           <DbName>genome</DbName>
           <DbName>books</DbName>
           <DbName>cancerchromosomes</DbName>
           <DbName>cdd</DbName>
           <DbName>gap</DbName>
           <DbName>domains</DbName>
           <DbName>gene</DbName>
           <DbName>genomeprj</DbName>
           <DbName>gensat</DbName>
           <DbName>geo</DbName>
           <DbName>gds</DbName>
           <DbName>homologene</DbName>
           <DbName>journals</DbName>
           <DbName>mesh</DbName>
           <DbName>ncbisearch</DbName>
           <DbName>nlmcatalog</DbName>
           <DbName>omia</DbName>
           <DbName>omim</DbName>
           <DbName>pmc</DbName>
           <DbName>popset</DbName>
           <DbName>probe</DbName>
           <DbName>proteinclusters</DbName>
           <DbName>pcassay</DbName>
           <DbName>pccompound</DbName>
           <DbName>pcsubstance</DbName>
           <DbName>snp</DbName>
           <DbName>taxonomy</DbName>
           <DbName>toolkit</DbName>
           <DbName>unigene</DbName>
           <DbName>unists</DbName>
   </DbList>
   </eInfoResult>

Since this is a fairly simple XML file, we could extract the information
it contains simply by string searching. Using ``Bio.Entrez``\ ’s parser
instead, we can directly parse this XML file into a Python object:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> stream = Entrez.einfo()
   >>> record = Entrez.read(stream)

Now ``record`` is a dictionary with three keys:

.. cont-doctest

.. code:: pycon

   >>> sorted(record.keys())
   ['DbInfo', 'DbList', 'ERROR']

The values stored in the ``'DbList`` key is the list of database
names shown in the XML above:

.. code:: pycon

   >>> record["DbList"]
   ['pubmed', 'protein', 'nucleotide', 'nuccore', 'nucgss', 'nucest',
    'structure', 'genome', 'books', 'cancerchromosomes', 'cdd', 'gap',
    'domains', 'gene', 'genomeprj', 'gensat', 'geo', 'gds', 'homologene',
    'journals', 'mesh', 'ncbisearch', 'nlmcatalog', 'omia', 'omim', 'pmc',
    'popset', 'probe', 'proteinclusters', 'pcassay', 'pccompound',
    'pcsubstance', 'snp', 'taxonomy', 'toolkit', 'unigene', 'unists']

For each of these databases, we can use EInfo again to obtain more
information:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.einfo(db="pubmed")
   >>> record = Entrez.read(stream)
   >>> record["DbInfo"][0]["Description"]
   'PubMed bibliographic record'

.. code:: pycon

   >>> record["DbInfo"][0]["Count"]
   '17989604'
   >>> record["DbInfo"][0]["LastUpdate"]
   '2008/05/24 06:45'

Try ``record["DbInfo"].keys()`` for other information stored in this
record. One of the most useful is a list of possible search fields for
use with ESearch:

.. code:: pycon

   >>> for field in record["DbInfo"]["FieldList"]:
   ...     print("%(Name)s, %(FullName)s, %(Description)s" % field)
   ...
   ALL, All Fields, All terms from all searchable fields
   UID, UID, Unique number assigned to publication
   FILT, Filter, Limits the records
   TITL, Title, Words in title of publication
   WORD, Text Word, Free text associated with publication
   MESH, MeSH Terms, Medical Subject Headings assigned to publication
   MAJR, MeSH Major Topic, MeSH terms of major importance to publication
   AUTH, Author, Author(s) of publication
   JOUR, Journal, Journal abbreviation of publication
   AFFL, Affiliation, Author's institutional affiliation and address
   ...

That’s a long list, but indirectly this tells you that for the PubMed
database, you can do things like ``Jones[AUTH]`` to search the author
field, or ``Sanger[AFFL]`` to restrict to authors at the Sanger Centre.
This can be very handy - especially if you are not so familiar with a
particular database.

.. _`sec:entrez-esearch`:

ESearch: Searching the Entrez databases
---------------------------------------

To search any of these databases, we use ``Bio.Entrez.esearch()``. For
example, let’s search in PubMed for publications that include Biopython
in their title:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esearch(db="pubmed", term="biopython[title]", retmax="40")
   >>> record = Entrez.read(stream)
   >>> "19304878" in record["IdList"]
   True

.. code:: pycon

   >>> print(record["IdList"])
   ['22909249', '19304878']

In this output, you see PubMed IDs (including 19304878 which is the PMID
for the Biopython application note), which can be retrieved by EFetch
(see section :ref:`sec:efetch`).

You can also use ESearch to search GenBank. Here we’ll do a quick search
for the *matK* gene in *Cypripedioideae* orchids (see
Section :ref:`sec:entrez-einfo` about EInfo for one way to find out
which fields you can search in each Entrez database):

.. code:: pycon

   >>> stream = Entrez.esearch(
   ...     db="nucleotide", term="Cypripedioideae[Orgn] AND matK[Gene]", idtype="acc"
   ... )
   >>> record = Entrez.read(stream)
   >>> record["Count"]
   '348'
   >>> record["IdList"]
   ['JQ660909.1', 'JQ660908.1', 'JQ660907.1', 'JQ660906.1', ..., 'JQ660890.1']

Each of the IDs (JQ660909.1, JQ660908.1, JQ660907.1, …) is a GenBank
identifier (Accession number). See section :ref:`sec:efetch` for
information on how to actually download these GenBank records.

Note that instead of a species name like ``Cypripedioideae[Orgn]``, you
can restrict the search using an NCBI taxon identifier, here this would
be ``txid158330[Orgn]``. This isn’t currently documented on the ESearch
help page - the NCBI explained this in reply to an email query. You can
often deduce the search term formatting by playing with the Entrez web
interface. For example, including ``complete[prop]`` in a genome search
restricts to just completed genomes.

As a final example, let’s get a list of computational journal titles:

.. code:: pycon

   >>> stream = Entrez.esearch(db="nlmcatalog", term="computational[Journal]", retmax="20")
   >>> record = Entrez.read(stream)
   >>> print("{} computational journals found".format(record["Count"]))
   117 computational Journals found
   >>> print("The first 20 are\n{}".format(record["IdList"]))
   ['101660833', '101664671', '101661657', '101659814', '101657941',
    '101653734', '101669877', '101649614', '101647835', '101639023',
    '101627224', '101647801', '101589678', '101585369', '101645372',
    '101586429', '101582229', '101574747', '101564639', '101671907']

Again, we could use EFetch to obtain more information for each of these
journal IDs.

ESearch has many useful options — see the `ESearch help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch>`__
for more information.

EPost: Uploading a list of identifiers
--------------------------------------

EPost uploads a list of UIs for use in subsequent search strategies; see
the `EPost help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EPost>`__
for more information. It is available from Biopython through the
``Bio.Entrez.epost()`` function.

To give an example of when this is useful, suppose you have a long list
of IDs you want to download using EFetch (maybe sequences, maybe
citations – anything). When you make a request with EFetch your list of
IDs, the database etc, are all turned into a long URL sent to the
server. If your list of IDs is long, this URL gets long, and long URLs
can break (e.g. some proxies don’t cope well).

Instead, you can break this up into two steps, first uploading the list
of IDs using EPost (this uses an “HTML post” internally, rather than an
“HTML get”, getting round the long URL problem). With the history
support, you can then refer to this long list of IDs, and download the
associated data with EFetch.

Let’s look at a simple example to see how EPost works – uploading some
PubMed identifiers:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> id_list = ["19304878", "18606172", "16403221", "16377612", "14871861", "14630660"]
   >>> print(Entrez.epost("pubmed", id=",".join(id_list)).read())
   <?xml version="1.0"?>
   <!DOCTYPE ePostResult PUBLIC "-//NLM//DTD ePostResult, 11 May 2002//EN"
    "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/ePost_020511.dtd">
   <ePostResult>
       <QueryKey>1</QueryKey>
       <WebEnv>NCID_01_206841095_130.14.22.101_9001_1242061629</WebEnv>
   </ePostResult>

The returned XML includes two important strings, ``QueryKey`` and
``WebEnv`` which together define your history session. You would extract
these values for use with another Entrez call such as EFetch:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> id_list = ["19304878", "18606172", "16403221", "16377612", "14871861", "14630660"]
   >>> search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(id_list)))
   >>> webenv = search_results["WebEnv"]
   >>> query_key = search_results["QueryKey"]

Section :ref:`sec:entrez-webenv` shows how to use the history
feature.

ESummary: Retrieving summaries from primary IDs
-----------------------------------------------

ESummary retrieves document summaries from a list of primary IDs (see
the `ESummary help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESummary>`__
for more information). In Biopython, ESummary is available as
``Bio.Entrez.esummary()``. Using the search result above, we can for
example find out more about the journal with ID 30367:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esummary(db="nlmcatalog", id="101660833")
   >>> record = Entrez.read(stream)
   >>> info = record[0]["TitleMainList"][0]
   >>> print("Journal info\nid: {}\nTitle: {}".format(record[0]["Id"], info["Title"]))
   Journal info
   id: 101660833
   Title: IEEE transactions on computational imaging.

.. _`sec:efetch`:

EFetch: Downloading full records from Entrez
--------------------------------------------

EFetch is what you use when you want to retrieve a full record from
Entrez. This covers several possible databases, as described on the main
`EFetch Help page <https://www.ncbi.nlm.nih.gov/books/NBK3837/>`__.

For most of their databases, the NCBI support several different file
formats. Requesting a specific file format from Entrez using
``Bio.Entrez.efetch()`` requires specifying the ``rettype`` and/or
``retmode`` optional arguments. The different combinations are described
for each database type on the pages linked to on `NCBI efetch
webpage <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`__.

One common usage is downloading sequences in the FASTA or
GenBank/GenPept plain text formats (which can then be parsed with
``Bio.SeqIO``, see
Sections :ref:`sec:SeqIO_GenBank_Online`
and :ref:`sec:efetch`). From the *Cypripedioideae* example above, we
can download GenBank record EU490707 using ``Bio.Entrez.efetch``:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
   >>> print(stream.read())
   LOCUS       EU490707                1302 bp    DNA     linear   PLN 26-JUL-2016
   DEFINITION  Selenipedium aequinoctiale maturase K (matK) gene, partial cds;
               chloroplast.
   ACCESSION   EU490707
   VERSION     EU490707.1
   KEYWORDS    .
   SOURCE      chloroplast Selenipedium aequinoctiale
     ORGANISM  Selenipedium aequinoctiale
               Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
               Spermatophyta; Magnoliopsida; Liliopsida; Asparagales; Orchidaceae;
               Cypripedioideae; Selenipedium.
   REFERENCE   1  (bases 1 to 1302)
     AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A., Endara,L.,
               Williams,N.H. and Moore,M.
     TITLE     Phylogenetic utility of ycf1 in orchids: a plastid gene more
               variable than matK
     JOURNAL   Plant Syst. Evol. 277 (1-2), 75-84 (2009)
   REFERENCE   2  (bases 1 to 1302)
     AUTHORS   Neubig,K.M., Whitten,W.M., Carlsward,B.S., Blanco,M.A.,
               Endara,C.L., Williams,N.H. and Moore,M.J.
     TITLE     Direct Submission
     JOURNAL   Submitted (14-FEB-2008) Department of Botany, University of
               Florida, 220 Bartram Hall, Gainesville, FL 32611-8526, USA
   FEATURES             Location/Qualifiers
        source          1..1302
                        /organism="Selenipedium aequinoctiale"
                        /organelle="plastid:chloroplast"
                        /mol_type="genomic DNA"
                        /specimen_voucher="FLAS:Blanco 2475"
                        /db_xref="taxon:256374"
        gene            <1..>1302
                        /gene="matK"
        CDS             <1..>1302
                        /gene="matK"
                        /codon_start=1
                        /transl_table=11
                        /product="maturase K"
                        /protein_id="ACC99456.1"
                        /translation="IFYEPVEIFGYDNKSSLVLVKRLITRMYQQNFLISSVNDSNQKG
                        FWGHKHFFSSHFSSQMVSEGFGVILEIPFSSQLVSSLEEKKIPKYQNLRSIHSIFPFL
                        EDKFLHLNYVSDLLIPHPIHLEILVQILQCRIKDVPSLHLLRLLFHEYHNLNSLITSK
                        KFIYAFSKRKKRFLWLLYNSYVYECEYLFQFLRKQSSYLRSTSSGVFLERTHLYVKIE
                        HLLVVCCNSFQRILCFLKDPFMHYVRYQGKAILASKGTLILMKKWKFHLVNFWQSYFH
                        FWSQPYRIHIKQLSNYSFSFLGYFSSVLENHLVVRNQMLENSFIINLLTKKFDTIAPV
                        ISLIGSLSKAQFCTVLGHPISKPIWTDFSDSDILDRFCRICRNLCRYHSGSSKKQVLY
                        RIKYILRLSCARTLARKHKSTVRTFMRRLGSGLLEEFFMEEE"
   ORIGIN      
           1 attttttacg aacctgtgga aatttttggt tatgacaata aatctagttt agtacttgtg
          61 aaacgtttaa ttactcgaat gtatcaacag aattttttga tttcttcggt taatgattct
         121 aaccaaaaag gattttgggg gcacaagcat tttttttctt ctcatttttc ttctcaaatg
         181 gtatcagaag gttttggagt cattctggaa attccattct cgtcgcaatt agtatcttct
         241 cttgaagaaa aaaaaatacc aaaatatcag aatttacgat ctattcattc aatatttccc
         301 tttttagaag acaaattttt acatttgaat tatgtgtcag atctactaat accccatccc
         361 atccatctgg aaatcttggt tcaaatcctt caatgccgga tcaaggatgt tccttctttg
         421 catttattgc gattgctttt ccacgaatat cataatttga atagtctcat tacttcaaag
         481 aaattcattt acgccttttc aaaaagaaag aaaagattcc tttggttact atataattct
         541 tatgtatatg aatgcgaata tctattccag tttcttcgta aacagtcttc ttatttacga
         601 tcaacatctt ctggagtctt tcttgagcga acacatttat atgtaaaaat agaacatctt
         661 ctagtagtgt gttgtaattc ttttcagagg atcctatgct ttctcaagga tcctttcatg
         721 cattatgttc gatatcaagg aaaagcaatt ctggcttcaa agggaactct tattctgatg
         781 aagaaatgga aatttcatct tgtgaatttt tggcaatctt attttcactt ttggtctcaa
         841 ccgtatagga ttcatataaa gcaattatcc aactattcct tctcttttct ggggtatttt
         901 tcaagtgtac tagaaaatca tttggtagta agaaatcaaa tgctagagaa ttcatttata
         961 ataaatcttc tgactaagaa attcgatacc atagccccag ttatttctct tattggatca
        1021 ttgtcgaaag ctcaattttg tactgtattg ggtcatccta ttagtaaacc gatctggacc
        1081 gatttctcgg attctgatat tcttgatcga ttttgccgga tatgtagaaa tctttgtcgt
        1141 tatcacagcg gatcctcaaa aaaacaggtt ttgtatcgta taaaatatat acttcgactt
        1201 tcgtgtgcta gaactttggc acggaaacat aaaagtacag tacgcacttt tatgcgaaga
        1261 ttaggttcgg gattattaga agaattcttt atggaagaag aa
   //
   <BLANKLINE>
   <BLANKLINE>

Please be aware that as of October 2016 GI identifiers are discontinued
in favor of accession numbers. You can still fetch sequences based on
their GI, but new sequences are no longer given this identifier. You
should instead refer to them by the “Accession number” as done in the
example.

The arguments ``rettype="gb"`` and ``retmode="text"`` let us download
this record in the GenBank format.

Note that until Easter 2009, the Entrez EFetch API let you use “genbank”
as the return type, however the NCBI now insist on using the official
return types of “gb” or “gbwithparts” (or “gp” for proteins) as
described on online. Also note that until Feb 2012, the Entrez EFetch
API would default to returning plain text files, but now defaults to
XML.

Alternatively, you could for example use ``rettype="fasta"`` to get the
Fasta-format; see the `EFetch Sequences Help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`__
for other options. Remember – the available formats depend on which
database you are downloading from - see the main `EFetch Help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`__.

If you fetch the record in one of the formats accepted by ``Bio.SeqIO``
(see Chapter :ref:`chapter:seqio`), you could directly
parse it into a ``SeqRecord``:

.. doctest . internet

.. code:: pycon

   >>> from Bio import SeqIO
   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
   >>> record = SeqIO.read(stream, "genbank")
   >>> stream.close()
   >>> print(record.id)
   EU490707.1
   >>> print(record.name)
   EU490707
   >>> print(record.description)
   Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast
   >>> print(len(record.features))
   3
   >>> record.seq
   Seq('ATTTTTTACGAACCTGTGGAAATTTTTGGTTATGACAATAAATCTAGTTTAGTA...GAA')

Note that a more typical use would be to save the sequence data to a
local file, and *then* parse it with ``Bio.SeqIO``. This can save you
having to re-download the same file repeatedly while working on your
script, and places less load on the NCBI’s servers. For example:

.. code:: python

   import os
   from Bio import SeqIO
   from Bio import Entrez

   Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   filename = "EU490707.gbk"
   if not os.path.isfile(filename):
       # Downloading...
       stream = Entrez.efetch(db="nucleotide", id="EU490707", rettype="gb", retmode="text")
       output = open(filename, "w")
       output.write(streame.read())
       output.close()
       stream.close()
       print("Saved")

   print("Parsing...")
   record = SeqIO.read(filename, "genbank")
   print(record)

To get the output in XML format, which you can parse using the
``Bio.Entrez.read()`` function, use ``retmode="xml"``:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.efetch(db="nucleotide", id="EU490707", retmode="xml")
   >>> record = Entrez.read(stream)
   >>> stream.close()
   >>> record[0]["GBSeq_definition"]
   'Selenipedium aequinoctiale maturase K (matK) gene, partial cds; chloroplast'
   >>> record[0]["GBSeq_source"]
   'chloroplast Selenipedium aequinoctiale'

So, that dealt with sequences. For examples of parsing file formats
specific to the other databases (e.g. the ``MEDLINE`` format used in
PubMed), see Section :ref:`sec:entrez-specialized-parsers`.

If you want to perform a search with ``Bio.Entrez.esearch()``, and then
download the records with ``Bio.Entrez.efetch()``, you should use the
WebEnv history feature – see Section :ref:`sec:entrez-webenv`.

.. _`sec:elink`:

ELink: Searching for related items in NCBI Entrez
-------------------------------------------------

ELink, available from Biopython as ``Bio.Entrez.elink()``, can be used
to find related items in the NCBI Entrez databases. For example, you can
us this to find nucleotide entries for an entry in the gene database,
and other cool stuff.

Let’s use ELink to find articles related to the Biopython application
note published in *Bioinformatics* in 2009. The PubMed ID of this
article is 19304878:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> pmid = "19304878"
   >>> record = Entrez.read(Entrez.elink(dbfrom="pubmed", id=pmid))

The ``record`` variable consists of a Python list, one for each database
in which we searched. Since we specified only one PubMed ID to search
for, ``record`` contains only one item. This item is a dictionary
containing information about our search term, as well as all the related
items that were found:

.. cont-doctest

.. code:: pycon

   >>> record[0]["DbFrom"]
   'pubmed'
   >>> record[0]["IdList"]
   ['19304878']

The ``"LinkSetDb"`` key contains the search results, stored as a list
consisting of one item for each target database. In our search results,
we only find hits in the PubMed database (although sub-divided into
categories):

.. cont-doctest

.. code:: pycon

   >>> len(record[0]["LinkSetDb"])
   8

The exact numbers should increase over time:

.. code:: pycon

   >>> for linksetdb in record[0]["LinkSetDb"]:
   ...     print(linksetdb["DbTo"], linksetdb["LinkName"], len(linksetdb["Link"]))
   ...
   pubmed pubmed_pubmed 284
   pubmed pubmed_pubmed_alsoviewed 7
   pubmed pubmed_pubmed_citedin 926
   pubmed pubmed_pubmed_combined 6
   pubmed pubmed_pubmed_five 6
   pubmed pubmed_pubmed_refs 17
   pubmed pubmed_pubmed_reviews 12
   pubmed pubmed_pubmed_reviews_five 6

The actual search results are stored as under the ``"Link"`` key.

Let’s now at the first search result:

.. code:: pycon

   >>> record[0]["LinkSetDb"][0]["Link"][0]
   {'Id': '19304878'}

This is the article we searched for, which doesn’t help us much, so
let’s look at the second search result:

.. code:: pycon

   >>> record[0]["LinkSetDb"][0]["Link"][1]
   {'Id': '14630660'}

This paper, with PubMed ID 14630660, is about the Biopython PDB parser.

We can use a loop to print out all PubMed IDs:

.. code:: pycon

   >>> for link in record[0]["LinkSetDb"][0]["Link"]:
   ...     print(link["Id"])
   ...
   19304878
   14630660
   18689808
   17121776
   16377612
   12368254
   ......

Now that was nice, but personally I am often more interested to find out
if a paper has been cited. Well, ELink can do that too – at least for
journals in Pubmed Central (see
Section :ref:`sec:elink-citations`).

For help on ELink, see the `ELink help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ELink>`__.
There is an entire sub-page just for the `link
names <https://eutils.ncbi.nlm.nih.gov/corehtml/query/static/entrezlinks.html>`__,
describing how different databases can be cross referenced.

EGQuery: Global Query - counts for search terms
-----------------------------------------------

EGQuery provides counts for a search term in each of the Entrez
databases (i.e. a global query). This is particularly useful to find out
how many items your search terms would find in each database without
actually performing lots of separate searches with ESearch (see the
example in :ref:`sec:entrez_example_genbank` below).

In this example, we use ``Bio.Entrez.egquery()`` to obtain the counts
for “Biopython”:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.egquery(term="biopython")
   >>> record = Entrez.read(stream)
   >>> for row in record["eGQueryResult"]:
   ...     print(row["DbName"], row["Count"])
   ...
   pubmed 6
   pmc 62
   journals 0
   ...

See the `EGQuery help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EGQuery>`__
for more information.

ESpell: Obtaining spelling suggestions
--------------------------------------

ESpell retrieves spelling suggestions. In this example, we use
``Bio.Entrez.espell()`` to obtain the correct spelling of Biopython:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.espell(term="biopythooon")
   >>> record = Entrez.read(stream)
   >>> record["Query"]
   'biopythooon'
   >>> record["CorrectedQuery"]
   'biopython'

See the `ESpell help
page <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESpell>`__
for more information. The main use of this is for GUI tools to provide
automatic suggestions for search terms.

Parsing huge Entrez XML files
-----------------------------

The ``Entrez.read`` function reads the entire XML file returned by
Entrez into a single Python object, which is kept in memory. To parse
Entrez XML files too large to fit in memory, you can use the function
``Entrez.parse``. This is a generator function that reads records in the
XML file one by one. This function is only useful if the XML file
reflects a Python list object (in other words, if ``Entrez.read`` on a
computer with infinite memory resources would return a Python list).

For example, you can download the entire Entrez Gene database for a
given organism as a file from NCBI’s ftp site. These files can be very
large. As an example, on September 4, 2009, the file
``Homo_sapiens.ags.gz``, containing the Entrez Gene database for human,
had a size of 116576 kB. This file, which is in the ``ASN`` format, can
be converted into an XML file using NCBI’s ``gene2xml`` program (see
NCBI’s ftp site for more information):

.. code:: console

   $ gene2xml -b T -i Homo_sapiens.ags -o Homo_sapiens.xml

The resulting XML file has a size of 6.1 GB. Attempting ``Entrez.read``
on this file will result in a ``MemoryError`` on many computers.

The XML file ``Homo_sapiens.xml`` consists of a list of Entrez gene
records, each corresponding to one Entrez gene in human.
``Entrez.parse`` retrieves these gene records one by one. You can then
print out or store the relevant information in each record by iterating
over the records. For example, this script iterates over the Entrez gene
records and prints out the gene numbers and names for all current genes:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = open("Homo_sapiens.xml", "rb")
   >>> records = Entrez.parse(stream)

Alternatively, you can use

.. code:: pycon

   >>> records = Entrez.parse("Homo_sapiens.xml")

and let ``Bio.Entrez`` take care of opening and closing the file. This
is safer, as the file will then automatically be closed after parsing
it, or if an error occurs.

.. code:: pycon

   >>> for record in records:
   ...     status = record["Entrezgene_track-info"]["Gene-track"]["Gene-track_status"]
   ...     if status.attributes["value"] == "discontinued":
   ...         continue
   ...     geneid = record["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]
   ...     genename = record["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
   ...     print(geneid, genename)
   ...
   1 A1BG
   2 A2M
   3 A2MP
   8 AA
   9 NAT1
   10 NAT2
   11 AACP
   12 SERPINA3
   13 AADAC
   14 AAMP
   15 AANAT
   16 AARS
   17 AAVS1
   ...

HTML escape characters
----------------------

Pubmed records may contain HTML tags to indicate e.g. subscripts,
superscripts, or italic text, as well as mathematical symbols via
MathML. By default, the ``Bio.Entrez`` parser treats all text as plain
text without markup; for example, the fragment “:math:`P < 0.05`” in the
abstract of a Pubmed record, which is encoded as

.. code:: text

   <i>P</i> &lt; 0.05

in the XML returned by Entrez, is converted to the Python string

.. code:: text

   '<i>P</i> < 0.05'

by the ``Bio.Entrez`` parser. While this is more human-readable, it is
not valid HTML due to the less-than sign, and makes further processing
of the text e.g. by an HTML parser impractical. To ensure that all
strings returned by the parser are valid HTML, call ``Entrez.read`` or
``Entrez.parse`` with the ``escape`` argument set to ``True``:

.. code:: pycon

   >>> record = Entrez.read(stream, escape=True)

The parser will then replace all characters disallowed in HTML by their
HTML-escaped equivalent; in the example above, the parser will generate

.. code:: text

   '<i>P</i> &lt; 0.05'

which is a valid HTML fragment. By default, ``escape`` is ``False``.

Handling errors
---------------

The file is not an XML file
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For example, this error occurs if you try to parse a Fasta file as if it
were an XML file:

.. doctest ../Tests/GenBank

.. code:: pycon

   >>> from Bio import Entrez
   >>> stream = open("NC_005816.fna", "rb")  # a Fasta file
   >>> record = Entrez.read(stream)
   Traceback (most recent call last):
     ...
   Bio.Entrez.Parser.NotXMLError: Failed to parse the XML data (syntax error: line 1, column 0). Please make sure that the input data are in XML format.

Here, the parser didn’t find the ``<?xml ...`` tag with which an XML
file is supposed to start, and therefore decides (correctly) that the
file is not an XML file.

The file ends prematurely or is otherwise corrupted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When your file is in the XML format but is corrupted (for example, by
ending prematurely), the parser will raise a CorruptedXMLError.

Here is an example of an XML file that ends prematurely:

.. code:: text

   <?xml version="1.0"?>
   <!DOCTYPE eInfoResult PUBLIC "-//NLM//DTD eInfoResult, 11 May 2002//EN" "https://www.ncbi.nlm.nih.gov/entrez/query/DTD/eInfo_020511.dtd">
   <eInfoResult>
   <DbList>
           <DbName>pubmed</DbName>
           <DbName>protein</DbName>
           <DbName>nucleotide</DbName>
           <DbName>nuccore</DbName>
           <DbName>nucgss</DbName>
           <DbName>nucest</DbName>
           <DbName>structure</DbName>
           <DbName>genome</DbName>
           <DbName>books</DbName>
           <DbName>cancerchromosomes</DbName>
           <DbName>cdd</DbName>

which will generate the following traceback:

.. code:: pycon

   >>> Entrez.read(stream)
   Traceback (most recent call last):
     ...
   Bio.Entrez.Parser.CorruptedXMLError: Failed to parse the XML data (no element found: line 16, column 0). Please make sure that the input data are not corrupted.

Note that the error message tells you at what point in the XML file the
error was detected.

The file contains items that are missing from the associated DTD
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some extraordinary cases, the XML file is not fully consistent with the
associated DTD file, and contains some tags that are not defined in the
DTD file. By default, the parser will stop and raise a ``ValidationError``
if it cannot find some tag in the DTD. You can instruct the parser to skip
such tags instead of raising a ``ValidationError`` by calling ``Entrez.read``
or ``Entrez.parse`` with the argument ``validate`` equal to False:

.. code:: pycon

   >>> from Bio import Entrez
   >>> stream = open("myxmlfile.xml", "rb")
   >>> record = Entrez.read(stream, validate=False)
   >>> stream.close()

Of course, the information contained in the XML tags that are not in the
DTD will not be present in the record returned by ``Entrez.read``.

The file contains an error message
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This may occur, for example, when you attempt to access a PubMed record
for a nonexistent PubMed ID. By default, this will raise a
``RuntimeError``:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esummary(db="pubmed", id="99999999")
   >>> record = Entrez.read(stream)
   Traceback (most recent call last):
   ...
   RuntimeError: UID=99999999: cannot get document summary

If you are accessing multiple PubMed records, the ``RuntimeError`` would
prevent you from receiving results for any of the PubMed records if one
of the PubMed IDs is incorrect. To circumvent this, you can set the
``ignore_errors`` argument to ``True``. This will return the requested
results for the valid PubMed IDs, and an ``ErrorElement`` for the
incorrect ID:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esummary(db="pubmed", id="19304878,99999999,31278684")
   >>> record = Entrez.read(stream, ignore_errors=True)
   >>> len(record)
   3
   >>> record[0].tag
   'DocSum'
   >>> record[0]["Title"]
   'Biopython: freely available Python tools for computational molecular biology and bioinformatics.'
   >>> record[1].tag
   'ERROR'
   >>> record[1]
   ErrorElement('UID=99999999: cannot get document summary')
   >>> record[2].tag
   'DocSum'
   >>> record[2]["Title"]
   'Sharing Programming Resources Between Bio* Projects.'

.. _`sec:entrez-specialized-parsers`:

Specialized parsers
-------------------

The ``Bio.Entrez.read()`` function can parse most (if not all) XML
output returned by Entrez. Entrez typically allows you to retrieve
records in other formats, which may have some advantages compared to the
XML format in terms of readability (or download size).

To request a specific file format from Entrez using
``Bio.Entrez.efetch()`` requires specifying the ``rettype`` and/or
``retmode`` optional arguments. The different combinations are described
for each database type on the `NCBI efetch
webpage <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`__.

One obvious case is you may prefer to download sequences in the FASTA or
GenBank/GenPept plain text formats (which can then be parsed with
``Bio.SeqIO``, see
Sections :ref:`sec:SeqIO_GenBank_Online`
and :ref:`sec:efetch`). For the literature databases, Biopython
contains a parser for the ``MEDLINE`` format used in PubMed.

.. _`sec:entrez-and-medline`:

Parsing Medline records
~~~~~~~~~~~~~~~~~~~~~~~

You can find the Medline parser in ``Bio.Medline``. Suppose we want to
parse the file ``pubmed_result1.txt``, containing one Medline record.
You can find this file in Biopython’s ``Tests\Medline`` directory. The
file looks like this:

.. code:: text

   PMID- 12230038
   OWN - NLM
   STAT- MEDLINE
   DA  - 20020916
   DCOM- 20030606
   LR  - 20041117
   PUBM- Print
   IS  - 1467-5463 (Print)
   VI  - 3
   IP  - 3
   DP  - 2002 Sep
   TI  - The Bio* toolkits--a brief overview.
   PG  - 296-302
   AB  - Bioinformatics research is often difficult to do with commercial software. The
         Open Source BioPerl, BioPython and Biojava projects provide toolkits with
   ...

We first open the file and then parse it:

.. doctest ../Tests/Medline

.. code:: pycon

   >>> from Bio import Medline
   >>> with open("pubmed_result1.txt") as stream:
   ...     record = Medline.read(stream)
   ...

The ``record`` now contains the Medline record as a Python dictionary:

.. cont-doctest

.. code:: pycon

   >>> record["PMID"]
   '12230038'

.. code:: pycon

   >>> record["AB"]
   'Bioinformatics research is often difficult to do with commercial software.
   The Open Source BioPerl, BioPython and Biojava projects provide toolkits with
   multiple functionality that make it easier to create customized pipelines or
   analysis. This review briefly compares the quirks of the underlying languages
   and the functionality, documentation, utility and relative advantages of the
   Bio counterparts, particularly from the point of view of the beginning
   biologist programmer.'

The key names used in a Medline record can be rather obscure; use

.. code:: pycon

   >>> help(record)

for a brief summary.

To parse a file containing multiple Medline records, you can use the
``parse`` function instead:

.. doctest ../Tests/Medline

.. code:: pycon

   >>> from Bio import Medline
   >>> with open("pubmed_result2.txt") as stream:
   ...     for record in Medline.parse(stream):
   ...         print(record["TI"])
   ...
   A high level interface to SCOP and ASTRAL implemented in python.
   GenomeDiagram: a python package for the visualization of large-scale genomic data.
   Open source clustering software.
   PDB file parser and structure class implemented in Python.

Instead of parsing Medline records stored in files, you can also parse
Medline records downloaded by ``Bio.Entrez.efetch``. For example, let’s
look at all Medline records in PubMed related to Biopython:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esearch(db="pubmed", term="biopython")
   >>> record = Entrez.read(stream)
   >>> record["IdList"]
   ['19304878', '18606172', '16403221', '16377612', '14871861', '14630660', '12230038']

We now use ``Bio.Entrez.efetch`` to download these Medline records:

.. code:: pycon

   >>> idlist = record["IdList"]
   >>> stream = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")

Here, we specify ``rettype="medline", retmode="text"`` to obtain the
Medline records in plain-text Medline format. Now we use ``Bio.Medline``
to parse these records:

.. code:: pycon

   >>> from Bio import Medline
   >>> records = Medline.parse(stream)
   >>> for record in records:
   ...     print(record["AU"])
   ...
   ['Cock PJ', 'Antao T', 'Chang JT', 'Chapman BA', 'Cox CJ', 'Dalke A', ..., 'de Hoon MJ']
   ['Munteanu CR', 'Gonzalez-Diaz H', 'Magalhaes AL']
   ['Casbon JA', 'Crooks GE', 'Saqi MA']
   ['Pritchard L', 'White JA', 'Birch PR', 'Toth IK']
   ['de Hoon MJ', 'Imoto S', 'Nolan J', 'Miyano S']
   ['Hamelryck T', 'Manderick B']
   ['Mangalam H']

For comparison, here we show an example using the XML format:

.. code:: pycon

   >>> stream = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="xml")
   >>> records = Entrez.read(stream)
   >>> for record in records["PubmedArticle"]:
   ...     print(record["MedlineCitation"]["Article"]["ArticleTitle"])
   ...
   Biopython: freely available Python tools for computational molecular biology and
    bioinformatics.
   Enzymes/non-enzymes classification model complexity based on composition, sequence,
    3D and topological indices.
   A high level interface to SCOP and ASTRAL implemented in python.
   GenomeDiagram: a python package for the visualization of large-scale genomic data.
   Open source clustering software.
   PDB file parser and structure class implemented in Python.
   The Bio* toolkits--a brief overview.

Note that in both of these examples, for simplicity we have naively
combined ESearch and EFetch. In this situation, the NCBI would expect
you to use their history feature, as illustrated in
Section :ref:`sec:entrez-webenv`.

Parsing GEO records
~~~~~~~~~~~~~~~~~~~

GEO (`Gene Expression Omnibus <https://www.ncbi.nlm.nih.gov/geo/>`__) is
a data repository of high-throughput gene expression and hybridization
array data. The ``Bio.Geo`` module can be used to parse GEO-formatted
data.

The following code fragment shows how to parse the example GEO file
``GSE16.txt`` into a record and print the record:

.. code:: pycon

   >>> from Bio import Geo
   >>> stream = open("GSE16.txt")
   >>> records = Geo.parse(stream)
   >>> for record in records:
   ...     print(record)
   ...

You can search the “gds” database (GEO datasets) with ESearch:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esearch(db="gds", term="GSE16")
   >>> record = Entrez.read(stream)
   >>> stream.close()
   >>> record["Count"]
   '27'

.. code:: pycon

   >>> record["IdList"]
   ['200000016', '100000028', ...]

From the Entrez website, UID “200000016” is GDS16 while the other hit
“100000028” is for the associated platform, GPL28. Unfortunately, at the
time of writing the NCBI don’t seem to support downloading GEO files
using Entrez (not as XML, nor in the *Simple Omnibus Format in Text*
(SOFT) format).

However, it is actually pretty straight forward to download the GEO
files by FTP from ftp://ftp.ncbi.nih.gov/pub/geo/ instead. In this case
you might want
ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SOFT/by_series/GSE16/GSE16_family.soft.gz
(a compressed file, see the Python module gzip).

Parsing UniGene records
~~~~~~~~~~~~~~~~~~~~~~~

UniGene is an NCBI database of the transcriptome, with each UniGene
record showing the set of transcripts that are associated with a
particular gene in a specific organism. A typical UniGene record looks
like this:

.. code:: text

   ID          Hs.2
   TITLE       N-acetyltransferase 2 (arylamine N-acetyltransferase)
   GENE        NAT2
   CYTOBAND    8p22
   GENE_ID     10
   LOCUSLINK   10
   HOMOL       YES
   EXPRESS      bone| connective tissue| intestine| liver| liver tumor| normal| soft tissue/muscle tissue tumor| adult
   RESTR_EXPR   adult
   CHROMOSOME  8
   STS         ACC=PMC310725P3 UNISTS=272646
   STS         ACC=WIAF-2120 UNISTS=44576
   STS         ACC=G59899 UNISTS=137181
   ...
   STS         ACC=GDB:187676 UNISTS=155563
   PROTSIM     ORG=10090; PROTGI=6754794; PROTID=NP_035004.1; PCT=76.55; ALN=288
   PROTSIM     ORG=9796; PROTGI=149742490; PROTID=XP_001487907.1; PCT=79.66; ALN=288
   PROTSIM     ORG=9986; PROTGI=126722851; PROTID=NP_001075655.1; PCT=76.90; ALN=288
   ...
   PROTSIM     ORG=9598; PROTGI=114619004; PROTID=XP_519631.2; PCT=98.28; ALN=288

   SCOUNT      38
   SEQUENCE    ACC=BC067218.1; NID=g45501306; PID=g45501307; SEQTYPE=mRNA
   SEQUENCE    ACC=NM_000015.2; NID=g116295259; PID=g116295260; SEQTYPE=mRNA
   SEQUENCE    ACC=D90042.1; NID=g219415; PID=g219416; SEQTYPE=mRNA
   SEQUENCE    ACC=D90040.1; NID=g219411; PID=g219412; SEQTYPE=mRNA
   SEQUENCE    ACC=BC015878.1; NID=g16198419; PID=g16198420; SEQTYPE=mRNA
   SEQUENCE    ACC=CR407631.1; NID=g47115198; PID=g47115199; SEQTYPE=mRNA
   SEQUENCE    ACC=BG569293.1; NID=g13576946; CLONE=IMAGE:4722596; END=5'; LID=6989; SEQTYPE=EST; TRACE=44157214
   ...
   SEQUENCE    ACC=AU099534.1; NID=g13550663; CLONE=HSI08034; END=5'; LID=8800; SEQTYPE=EST
   //

This particular record shows the set of transcripts (shown in the
``SEQUENCE`` lines) that originate from the human gene NAT2, encoding en
N-acetyltransferase. The ``PROTSIM`` lines show proteins with
significant similarity to NAT2, whereas the ``STS`` lines show the
corresponding sequence-tagged sites in the genome.

To parse UniGene files, use the ``Bio.UniGene`` module:

.. code:: pycon

   >>> from Bio import UniGene
   >>> input = open("myunigenefile.data")
   >>> record = UniGene.read(input)

The ``record`` returned by ``UniGene.read`` is a Python object with
attributes corresponding to the fields in the UniGene record. For
example,

.. code:: pycon

   >>> record.ID
   "Hs.2"
   >>> record.title
   "N-acetyltransferase 2 (arylamine N-acetyltransferase)"

The ``EXPRESS`` and ``RESTR_EXPR`` lines are stored as Python lists of
strings:

.. code:: python

   [
       "bone",
       "connective tissue",
       "intestine",
       "liver",
       "liver tumor",
       "normal",
       "soft tissue/muscle tissue tumor",
       "adult",
   ]

Specialized objects are returned for the ``STS``, ``PROTSIM``, and
``SEQUENCE`` lines, storing the keys shown in each line as attributes:

.. code:: pycon

   >>> record.sts[0].acc
   'PMC310725P3'
   >>> record.sts[0].unists
   '272646'

and similarly for the ``PROTSIM`` and ``SEQUENCE`` lines.

To parse a file containing more than one UniGene record, use the
``parse`` function in ``Bio.UniGene``:

.. code:: pycon

   >>> from Bio import UniGene
   >>> input = open("unigenerecords.data")
   >>> records = UniGene.parse(input)
   >>> for record in records:
   ...     print(record.ID)
   ...

Using a proxy
-------------

Normally you won’t have to worry about using a proxy, but if this is an
issue on your network here is how to deal with it. Internally,
``Bio.Entrez`` uses the standard Python library ``urllib`` for accessing
the NCBI servers. This will check an environment variable called
``http_proxy`` to configure any simple proxy automatically.
Unfortunately this module does not support the use of proxies which
require authentication.

You may choose to set the ``http_proxy`` environment variable once (how
you do this will depend on your operating system). Alternatively you can
set this within Python at the start of your script, for example:

.. code:: python

   import os

   os.environ["http_proxy"] = "http://proxyhost.example.com:8080"

See the `urllib
documentation <https://docs.python.org/2/library/urllib.html>`__ for
more details.

.. _`sec:entrez_examples`:

Examples
--------

.. _`sec:pub_med`:

PubMed and Medline
~~~~~~~~~~~~~~~~~~

If you are in the medical field or interested in human issues (and many
times even if you are not!), PubMed
(https://www.ncbi.nlm.nih.gov/PubMed/) is an excellent source of all
kinds of goodies. So like other things, we’d like to be able to grab
information from it and use it in Python scripts.

In this example, we will query PubMed for all articles having to do with
orchids (see section :ref:`sec:orchids` for our
motivation). We first check how many of such articles there are:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.egquery(term="orchid")
   >>> record = Entrez.read(stream)
   >>> for row in record["eGQueryResult"]:
   ...     if row["DbName"] == "pubmed":
   ...         print(row["Count"])
   ...
   463

Now we use the ``Bio.Entrez.efetch`` function to download the PubMed IDs
of these 463 articles:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esearch(db="pubmed", term="orchid", retmax=463)
   >>> record = Entrez.read(stream)
   >>> stream.close()
   >>> idlist = record["IdList"]

This returns a Python list containing all of the PubMed IDs of articles
related to orchids:

.. code:: pycon

   >>> print(idlist)
   ['18680603', '18665331', '18661158', '18627489', '18627452', '18612381',
   '18594007', '18591784', '18589523', '18579475', '18575811', '18575690',
   ...

Now that we’ve got them, we obviously want to get the corresponding
Medline records and extract the information from them. Here, we’ll
download the Medline records in the Medline flat-file format, and use
the ``Bio.Medline`` module to parse them:

.. cont-doctest

.. code:: pycon

   >>> from Bio import Medline
   >>> stream = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
   >>> records = Medline.parse(stream)

NOTE - We’ve just done a separate search and fetch here, the NCBI much
prefer you to take advantage of their history support in this situation.
See Section :ref:`sec:entrez-webenv`.

Keep in mind that ``records`` is an iterator, so you can iterate through
the records only once. If you want to save the records, you can convert
them to a list:

.. cont-doctest

.. code:: pycon

   >>> records = list(records)

Let’s now iterate over the records to print out some information about
each record:

.. code:: pycon

   >>> for record in records:
   ...     print("title:", record.get("TI", "?"))
   ...     print("authors:", record.get("AU", "?"))
   ...     print("source:", record.get("SO", "?"))
   ...     print("")
   ...

The output for this looks like:

.. code:: text

   title: Sex pheromone mimicry in the early spider orchid (ophrys sphegodes):
   patterns of hydrocarbons as the key mechanism for pollination by sexual
   deception [In Process Citation]
   authors: ['Schiestl FP', 'Ayasse M', 'Paulus HF', 'Lofstedt C', 'Hansson BS',
   'Ibarra F', 'Francke W']
   source: J Comp Physiol [A] 2000 Jun;186(6):567-74

Especially interesting to note is the list of authors, which is returned
as a standard Python list. This makes it easy to manipulate and search
using standard Python tools. For instance, we could loop through a whole
bunch of entries searching for a particular author with code like the
following:

.. code:: pycon

   >>> search_author = "Waits T"
   >>> for record in records:
   ...     if not "AU" in record:
   ...         continue
   ...     if search_author in record["AU"]:
   ...         print("Author %s found: %s" % (search_author, record["SO"]))
   ...

Hopefully this section gave you an idea of the power and flexibility of
the Entrez and Medline interfaces and how they can be used together.

.. _`sec:entrez_example_genbank`:

Searching, downloading, and parsing Entrez Nucleotide records
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we’ll show a simple example of performing a remote Entrez query. In
section :ref:`sec:orchids` of the parsing examples, we
talked about using NCBI’s Entrez website to search the NCBI nucleotide
databases for info on Cypripedioideae, our friends the lady slipper
orchids. Now, we’ll look at how to automate that process using a Python
script. In this example, we’ll just show how to connect, get the
results, and parse them, with the Entrez module doing all of the work.

First, we use EGQuery to find out the number of results we will get
before actually downloading them. EGQuery will tell us how many search
results were found in each of the databases, but for this example we are
only interested in nucleotides:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.egquery(term="Cypripedioideae")
   >>> record = Entrez.read(stream)
   >>> for row in record["eGQueryResult"]:
   ...     if row["DbName"] == "nuccore":
   ...         print(row["Count"])
   ...
   4457

So, we expect to find 4457 Entrez Nucleotide records (this increased
from 814 records in 2008; it is likely to continue to increase in the
future). If you find some ridiculously high number of hits, you may want
to reconsider if you really want to download all of them, which is our
next step. Let’s use the ``retmax`` argument to restrict the maximum
number of records retrieved to the number available in 2008:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esearch(
   ...     db="nucleotide", term="Cypripedioideae", retmax=814, idtype="acc"
   ... )
   >>> record = Entrez.read(stream)
   >>> stream.close()

Here, ``record`` is a Python dictionary containing the search results
and some auxiliary information. Just for information, let’s look at what
is stored in this dictionary:

.. code:: pycon

   >>> print(record.keys())
   ['Count', 'RetMax', 'IdList', 'TranslationSet', 'RetStart', 'QueryTranslation']

First, let’s check how many results were found:

.. code:: pycon

   >>> print(record["Count"])
   '4457'

You might have expected this to be 814, the maximum number of records we
asked to retrieve. However, ``Count`` represents the total number of
records available for that search, not how many were retrieved. The
retrieved records are stored in ``record['IdList']``, which should
contain the total number we asked for:

.. code:: pycon

   >>> len(record["IdList"])
   814

Let’s look at the first five results:

.. code:: pycon

   >>> record["IdList"][:5]
   ['KX265015.1', 'KX265014.1', 'KX265013.1', 'KX265012.1', 'KX265011.1']

We can download these records using ``efetch``. While you could download
these records one by one, to reduce the load on NCBI’s servers, it is
better to fetch a bunch of records at the same time, shown below.
However, in this situation you should ideally be using the history
feature described later in Section :ref:`sec:entrez-webenv`.

.. code:: pycon

   >>> idlist = ",".join(record["IdList"][:5])
   >>> print(idlist)
   KX265015.1, KX265014.1, KX265013.1, KX265012.1, KX265011.1]
   >>> stream = Entrez.efetch(db="nucleotide", id=idlist, retmode="xml")
   >>> records = Entrez.read(stream)
   >>> len(records)
   5

Each of these records corresponds to one GenBank record.

.. code:: pycon

   >>> print(records[0].keys())
   ['GBSeq_moltype', 'GBSeq_source', 'GBSeq_sequence',
    'GBSeq_primary-accession', 'GBSeq_definition', 'GBSeq_accession-version',
    'GBSeq_topology', 'GBSeq_length', 'GBSeq_feature-table',
    'GBSeq_create-date', 'GBSeq_other-seqids', 'GBSeq_division',
    'GBSeq_taxonomy', 'GBSeq_references', 'GBSeq_update-date',
    'GBSeq_organism', 'GBSeq_locus', 'GBSeq_strandedness']

   >>> print(records[0]["GBSeq_primary-accession"])
   DQ110336

   >>> print(records[0]["GBSeq_other-seqids"])
   ['gb|DQ110336.1|', 'gi|187237168']

   >>> print(records[0]["GBSeq_definition"])
   Cypripedium calceolus voucher Davis 03-03 A maturase (matR) gene, partial cds;
   mitochondrial

   >>> print(records[0]["GBSeq_organism"])
   Cypripedium calceolus

You could use this to quickly set up searches – but for heavy usage, see
Section :ref:`sec:entrez-webenv`.

.. _`sec:entrez-search-fetch-genbank`:

Searching, downloading, and parsing GenBank records
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GenBank record format is a very popular method of holding
information about sequences, sequence features, and other associated
sequence information. The format is a good way to get information from
the NCBI databases at https://www.ncbi.nlm.nih.gov/.

In this example we’ll show how to query the NCBI databases,to retrieve
the records from the query, and then parse them using ``Bio.SeqIO`` -
something touched on in
Section :ref:`sec:SeqIO_GenBank_Online`. For
simplicity, this example *does not* take advantage of the WebEnv history
feature – see Section :ref:`sec:entrez-webenv` for this.

First, we want to make a query and find out the ids of the records to
retrieve. Here we’ll do a quick search for one of our favorite
organisms, *Opuntia* (prickly-pear cacti). We can do quick search and
get back the GIs (GenBank identifiers) for all of the corresponding
records. First we check how many records there are:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.egquery(term="Opuntia AND rpl16")
   >>> record = Entrez.read(stream)
   >>> for row in record["eGQueryResult"]:
   ...     if row["DbName"] == "nuccore":
   ...         print(row["Count"])
   ...
   9

Now we download the list of GenBank identifiers:

.. code:: pycon

   >>> stream = Entrez.esearch(db="nuccore", term="Opuntia AND rpl16")
   >>> record = Entrez.read(stream)
   >>> gi_list = record["IdList"]
   >>> gi_list
   ['57240072', '57240071', '6273287', '6273291', '6273290', '6273289', '6273286',
   '6273285', '6273284']

Now we use these GIs to download the GenBank records - note that with
older versions of Biopython you had to supply a comma separated list of
GI numbers to Entrez, as of Biopython 1.59 you can pass a list and this
is converted for you:

.. code:: pycon

   >>> gi_str = ",".join(gi_list)
   >>> stream = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")

If you want to look at the raw GenBank files, you can read from this
stream and print out the result:

.. code:: pycon

   >>> text = stream.read()
   >>> print(text)
   LOCUS       AY851612                 892 bp    DNA     linear   PLN 10-APR-2007
   DEFINITION  Opuntia subulata rpl16 gene, intron; chloroplast.
   ACCESSION   AY851612
   VERSION     AY851612.1  GI:57240072
   KEYWORDS    .
   SOURCE      chloroplast Austrocylindropuntia subulata
     ORGANISM  Austrocylindropuntia subulata
               Eukaryota; Viridiplantae; Streptophyta; Embryophyta; Tracheophyta;
               Spermatophyta; Magnoliophyta; eudicotyledons; core eudicotyledons;
               Caryophyllales; Cactaceae; Opuntioideae; Austrocylindropuntia.
   REFERENCE   1  (bases 1 to 892)
     AUTHORS   Butterworth,C.A. and Wallace,R.S.
   ...

In this case, we are just getting the raw records. To get the records in
a more Python-friendly form, we can use ``Bio.SeqIO`` to parse the
GenBank data into ``SeqRecord`` objects, including ``SeqFeature``
objects (see Chapter :ref:`chapter:seqio`):

.. code:: pycon

   >>> from Bio import SeqIO
   >>> stream = Entrez.efetch(db="nuccore", id=gi_str, rettype="gb", retmode="text")
   >>> records = SeqIO.parse(stream, "gb")

We can now step through the records and look at the information we are
interested in:

.. code:: pycon

   >>> for record in records:
   ...     print(f"{record.name}, length {len(record)}, with {len(record.features)} features")
   ...
   AY851612, length 892, with 3 features
   AY851611, length 881, with 3 features
   AF191661, length 895, with 3 features
   AF191665, length 902, with 3 features
   AF191664, length 899, with 3 features
   AF191663, length 899, with 3 features
   AF191660, length 893, with 3 features
   AF191659, length 894, with 3 features
   AF191658, length 896, with 3 features

Using these automated query retrieval functionality is a big plus over
doing things by hand. Although the module should obey the NCBI’s max
three queries per second rule, the NCBI have other recommendations like
avoiding peak hours. See Section :ref:`sec:entrez-guidelines`. In
particular, please note that for simplicity, this example does not use
the WebEnv history feature. You should use this for any non-trivial
search and download work, see Section :ref:`sec:entrez-webenv`.

Finally, if plan to repeat your analysis, rather than downloading the
files from the NCBI and parsing them immediately (as shown in this
example), you should just download the records *once* and save them to
your hard disk, and then parse the local file.

Finding the lineage of an organism
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Staying with a plant example, let’s now find the lineage of the
Cypripedioideae orchid family. First, we search the Taxonomy database
for Cypripedioideae, which yields exactly one NCBI taxonomy identifier:

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esearch(db="Taxonomy", term="Cypripedioideae")
   >>> record = Entrez.read(stream)
   >>> record["IdList"]
   ['158330']
   >>> record["IdList"][0]
   '158330'

Now, we use ``efetch`` to download this entry in the Taxonomy database,
and then parse it:

.. cont-doctest

.. code:: pycon

   >>> stream = Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
   >>> records = Entrez.read(stream)

Again, this record stores lots of information:

.. code:: pycon

   >>> records[0].keys()
   ['Lineage', 'Division', 'ParentTaxId', 'PubDate', 'LineageEx',
    'CreateDate', 'TaxId', 'Rank', 'GeneticCode', 'ScientificName',
    'MitoGeneticCode', 'UpdateDate']

We can get the lineage directly from this record:

.. code:: pycon

   >>> records[0]["Lineage"]
   'cellular organisms; Eukaryota; Viridiplantae; Streptophyta; Streptophytina;
    Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida;
    Liliopsida; Asparagales; Orchidaceae'

The record data contains much more than just the information shown here
- for example look under ``"LineageEx"`` instead of ``"Lineage"`` and
you’ll get the NCBI taxon identifiers of the lineage entries too.

.. _`sec:entrez-webenv`:

Using the history and WebEnv
----------------------------

Often you will want to make a series of linked queries. Most typically,
running a search, perhaps refining the search, and then retrieving
detailed search results. You *can* do this by making a series of
separate calls to Entrez. However, the NCBI prefer you to take advantage
of their history support - for example combining ESearch and EFetch.

Another typical use of the history support would be to combine EPost and
EFetch. You use EPost to upload a list of identifiers, which starts a
new history session. You then download the records with EFetch by
referring to the session (instead of the identifiers).

Searching for and downloading sequences using the history
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Suppose we want to search and download all the *Opuntia* rpl16
nucleotide sequences, and store them in a FASTA file. As shown in
Section :ref:`sec:entrez-search-fetch-genbank`, we can naively
combine ``Bio.Entrez.esearch()`` to get a list of Accession numbers, and
then call ``Bio.Entrez.efetch()`` to download them all.

However, the approved approach is to run the search with the history
feature. Then, we can fetch the results by reference to the search
results - which the NCBI can anticipate and cache.

To do this, call ``Bio.Entrez.esearch()`` as normal, but with the
additional argument of ``usehistory="y"``,

.. doctest . internet

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "history.user@example.com"  # Always tell NCBI who you are
   >>> stream = Entrez.esearch(
   ...     db="nucleotide", term="Opuntia[orgn] and rpl16", usehistory="y", idtype="acc"
   ... )
   >>> search_results = Entrez.read(stream)
   >>> stream.close()

As before (see Section :ref:`sec:entrez_example_genbank`), the
XML output includes the first ``retmax`` search results, with ``retmax``
defaulting to 20:

.. cont-doctest

.. code:: pycon

   >>> acc_list = search_results["IdList"]
   >>> count = int(search_results["Count"])
   >>> len(acc_list)
   20

.. code:: pycon

   >>> count
   28

You also get given two additional pieces of information, the ``WebEnv``
session cookie, and the ``QueryKey``:

.. cont-doctest

.. code:: pycon

   >>> webenv = search_results["WebEnv"]
   >>> query_key = search_results["QueryKey"]

Having stored these values in variables ``session_cookie`` and
``query_key`` we can use them as parameters to ``Bio.Entrez.efetch()``
instead of giving the GI numbers as identifiers.

While for small searches you might be OK downloading everything at once,
it is better to download in batches. You use the ``retstart`` and
``retmax`` parameters to specify which range of search results you want
returned (starting entry using zero-based counting, and maximum number
of results to return). Note that if Biopython encounters a transient
failure like a HTTP 500 response when communicating with NCBI, it will
automatically try again a couple of times. For example,

.. code:: python

   # This assumes you have already run a search as shown above,
   # and set the variables count, webenv, query_key

   batch_size = 3
   output = open("orchid_rpl16.fasta", "w")
   for start in range(0, count, batch_size):
       end = min(count, start + batch_size)
       print("Going to download record %i to %i" % (start + 1, end))
       stream = Entrez.efetch(
           db="nucleotide",
           rettype="fasta",
           retmode="text",
           retstart=start,
           retmax=batch_size,
           webenv=webenv,
           query_key=query_key,
           idtype="acc",
       )
       data = stream.read()
       stream.close()
       output.write(data)
   output.close()

For illustrative purposes, this example downloaded the FASTA records in
batches of three. Unless you are downloading genomes or chromosomes, you
would normally pick a larger batch size.

Searching for and downloading abstracts using the history
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is another history example, searching for papers published in the
last year about the *Opuntia*, and then downloading them into a file in
MedLine format:

.. code:: python

   from Bio import Entrez

   Entrez.email = "history.user@example.com"
   search_results = Entrez.read(
       Entrez.esearch(
           db="pubmed", term="Opuntia[ORGN]", reldate=365, datetype="pdat", usehistory="y"
       )
   )
   count = int(search_results["Count"])
   print("Found %i results" % count)

   batch_size = 10
   output = open("recent_orchid_papers.txt", "w")
   for start in range(0, count, batch_size):
       end = min(count, start + batch_size)
       print("Going to download record %i to %i" % (start + 1, end))
       stream = Entrez.efetch(
           db="pubmed",
           rettype="medline",
           retmode="text",
           retstart=start,
           retmax=batch_size,
           webenv=search_results["WebEnv"],
           query_key=search_results["QueryKey"],
       )
       data = stream.read()
       stream.close()
       output.write(data)
   output.close()

At the time of writing, this gave 28 matches - but because this is a
date dependent search, this will of course vary. As described in
Section :ref:`sec:entrez-and-medline` above, you can then use
``Bio.Medline`` to parse the saved records.

.. _`sec:elink-citations`:

Searching for citations
~~~~~~~~~~~~~~~~~~~~~~~

Back in Section :ref:`sec:elink` we mentioned ELink can be used to
search for citations of a given paper. Unfortunately this only covers
journals indexed for PubMed Central (doing it for all the journals in
PubMed would mean a lot more work for the NIH). Let’s try this for the
Biopython PDB parser paper, PubMed ID 14630660:

.. code:: pycon

   >>> from Bio import Entrez
   >>> Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are
   >>> pmid = "14630660"
   >>> results = Entrez.read(
   ...     Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", id=pmid)
   ... )
   >>> pmc_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
   >>> pmc_ids
   ['2744707', '2705363', '2682512', ..., '1190160']

Great - eleven articles. But why hasn’t the Biopython application note
been found (PubMed ID 19304878)? Well, as you might have guessed from
the variable names, there are not actually PubMed IDs, but PubMed
Central IDs. Our application note is the third citing paper in that
list, PMCID 2682512.

So, what if (like me) you’d rather get back a list of PubMed IDs? Well
we can call ELink again to translate them. This becomes a two step
process, so by now you should expect to use the history feature to
accomplish it (Section :ref:`sec:entrez-webenv`).

But first, taking the more straightforward approach of making a second
(separate) call to ELink:

.. code:: pycon

   >>> results2 = Entrez.read(
   ...     Entrez.elink(dbfrom="pmc", db="pubmed", LinkName="pmc_pubmed", id=",".join(pmc_ids))
   ... )
   >>> pubmed_ids = [link["Id"] for link in results2[0]["LinkSetDb"][0]["Link"]]
   >>> pubmed_ids
   ['19698094', '19450287', '19304878', ..., '15985178']

This time you can immediately spot the Biopython application note as the
third hit (PubMed ID 19304878).

Now, let’s do that all again but with the history … *TODO*.

And finally, don’t forget to include your *own* email address in the
Entrez calls.
