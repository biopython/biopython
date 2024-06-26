.. _`chapter:uniprot`:

Swiss-Prot and ExPASy
=====================

Parsing Swiss-Prot files
------------------------

Swiss-Prot (https://web.expasy.org/docs/swiss-prot_guideline.html) is a
hand-curated database of protein sequences. Biopython can parse the
“plain text” Swiss-Prot file format, which is still used for the UniProt
Knowledgebase which combined Swiss-Prot, TrEMBL and PIR-PSD.

Although in the following we focus on the older human readable plain
text format, ``Bio.SeqIO`` can read both this and the newer UniProt XML
file format for annotated protein sequences.

Parsing Swiss-Prot records
~~~~~~~~~~~~~~~~~~~~~~~~~~

In
Section :ref:`sec:SeqIO_ExPASy_and_SwissProt`,
we described how to extract the sequence of a Swiss-Prot record as a
``SeqRecord`` object. Alternatively, you can store the Swiss-Prot record
in a ``Bio.SwissProt.Record`` object, which in fact stores the complete
information contained in the Swiss-Prot record. In this section, we
describe how to extract ``Bio.SwissProt.Record`` objects from a
Swiss-Prot file.

To parse a Swiss-Prot record, we first get a handle to a Swiss-Prot
record. There are several ways to do so, depending on where and how the
Swiss-Prot record is stored:

-  Open a Swiss-Prot file locally:

   .. doctest ../Tests

   .. code:: pycon

      >>> handle = open("SwissProt/F2CXE6.txt")

-  Open a gzipped Swiss-Prot file:

   .. code:: pycon

      >>> import gzip
      >>> handle = gzip.open("myswissprotfile.dat.gz", "rt")

-  Open a Swiss-Prot file over the internet:

   .. code:: pycon

      >>> from urllib.request import urlopen
      >>> url = "https://raw.githubusercontent.com/biopython/biopython/master/Tests/SwissProt/F2CXE6.txt"
      >>> handle = urlopen(url)

   to open the file stored on the Internet before calling ``read``.

-  Open a Swiss-Prot file over the internet from the ExPASy database
   (see section :ref:`sec:expasy_swissprot`):

   .. code:: pycon

      >>> from Bio import ExPASy
      >>> handle = ExPASy.get_sprot_raw("F2CXE6")

The key point is that for the parser, it doesn’t matter how the handle
was created, as long as it points to data in the Swiss-Prot format. The
parser will automatically decode the data as ASCII (the encoding used by
Swiss-Prot) if the handle was opened in binary mode.

We can use ``Bio.SeqIO`` as described in
Section :ref:`sec:SeqIO_ExPASy_and_SwissProt`
to get file format agnostic ``SeqRecord`` objects. Alternatively, we can
use ``Bio.SwissProt`` get ``Bio.SwissProt.Record`` objects, which are a
much closer match to the underlying file format.

To read one Swiss-Prot record from the handle, we use the function
``read()``:

.. cont-doctest

.. code:: pycon

   >>> from Bio import SwissProt
   >>> record = SwissProt.read(handle)

This function should be used if the handle points to exactly one
Swiss-Prot record. It raises a ``ValueError`` if no Swiss-Prot record
was found, and also if more than one record was found.

We can now print out some information about this record:

.. cont-doctest

.. code:: pycon

   >>> print(record.description)
   SubName: Full=Plasma membrane intrinsic protein {ECO:0000313|EMBL:BAN04711.1}; SubName: Full=Predicted protein {ECO:0000313|EMBL:BAJ87517.1};
   >>> for ref in record.references:
   ...     print("authors:", ref.authors)
   ...     print("title:", ref.title)
   ...
   authors: Matsumoto T., Tanaka T., Sakai H., Amano N., Kanamori H., Kurita K., Kikuta A., Kamiya K., Yamamoto M., Ikawa H., Fujii N., Hori K., Itoh T., Sato K.
   title: Comprehensive sequence analysis of 24,783 barley full-length cDNAs derived from 12 clone libraries.
   authors: Shibasaka M., Sasano S., Utsugi S., Katsuhara M.
   title: Functional characterization of a novel plasma membrane intrinsic protein2 in barley.
   authors: Shibasaka M., Katsuhara M., Sasano S.
   title: 
   >>> print(record.organism_classification)
   ['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Tracheophyta', 'Spermatophyta', 'Magnoliophyta', 'Liliopsida', 'Poales', 'Poaceae', 'BEP clade', 'Pooideae', 'Triticeae', 'Hordeum']

To parse a file that contains more than one Swiss-Prot record, we use
the ``parse`` function instead. This function allows us to iterate over
the records in the file.

For example, let’s parse the full Swiss-Prot database and collect all
the descriptions. You can download this from the `ExPASy FTP
site <ftp://ftp.expasy.org/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz>`__
as a single gzipped-file ``uniprot_sprot.dat.gz`` (about 300MB). This is
a compressed file containing a single file, ``uniprot_sprot.dat`` (over
1.5GB).

As described at the start of this section, you can use the Python
library ``gzip`` to open and uncompress a ``.gz`` file, like this:

.. code:: pycon

   >>> import gzip
   >>> handle = gzip.open("uniprot_sprot.dat.gz", "rt")

However, uncompressing a large file takes time, and each time you open
the file for reading in this way, it has to be decompressed on the fly.
So, if you can spare the disk space you’ll save time in the long run if
you first decompress the file to disk, to get the ``uniprot_sprot.dat``
file inside. Then you can open the file for reading as usual:

.. code:: pycon

   >>> handle = open("uniprot_sprot.dat")

As of June 2009, the full Swiss-Prot database downloaded from ExPASy
contained 468851 Swiss-Prot records. One concise way to build up a list
of the record descriptions is with a list comprehension:

.. code:: pycon

   >>> from Bio import SwissProt
   >>> handle = open("uniprot_sprot.dat")
   >>> descriptions = [record.description for record in SwissProt.parse(handle)]
   >>> len(descriptions)
   468851
   >>> descriptions[:5]
   ['RecName: Full=Protein MGF 100-1R;',
    'RecName: Full=Protein MGF 100-1R;',
    'RecName: Full=Protein MGF 100-1R;',
    'RecName: Full=Protein MGF 100-1R;',
    'RecName: Full=Protein MGF 100-2L;']

Or, using a for loop over the record iterator:

.. code:: pycon

   >>> from Bio import SwissProt
   >>> descriptions = []
   >>> handle = open("uniprot_sprot.dat")
   >>> for record in SwissProt.parse(handle):
   ...     descriptions.append(record.description)
   ...
   >>> len(descriptions)
   468851

Because this is such a large input file, either way takes about eleven
minutes on my new desktop computer (using the uncompressed
``uniprot_sprot.dat`` file as input).

It is equally easy to extract any kind of information you’d like from
Swiss-Prot records. To see the members of a Swiss-Prot record, use

.. code:: pycon

   >>> dir(record)
   ['__doc__', '__init__', '__module__', 'accessions', 'annotation_update',
   'comments', 'created', 'cross_references', 'data_class', 'description',
   'entry_name', 'features', 'gene_name', 'host_organism', 'keywords',
   'molecule_type', 'organelle', 'organism', 'organism_classification',
   'references', 'seqinfo', 'sequence', 'sequence_length',
   'sequence_update', 'taxonomy_id']

Parsing the Swiss-Prot keyword and category list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Swiss-Prot also distributes a file ``keywlist.txt``, which lists the
keywords and categories used in Swiss-Prot. The file contains entries in
the following form:

.. code:: text

   ID   2Fe-2S.
   AC   KW-0001
   DE   Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron
   DE   atoms complexed to 2 inorganic sulfides and 4 sulfur atoms of
   DE   cysteines from the protein.
   SY   Fe2S2; [2Fe-2S] cluster; [Fe2S2] cluster; Fe2/S2 (inorganic) cluster;
   SY   Di-mu-sulfido-diiron; 2 iron, 2 sulfur cluster binding.
   GO   GO:0051537; 2 iron, 2 sulfur cluster binding
   HI   Ligand: Iron; Iron-sulfur; 2Fe-2S.
   HI   Ligand: Metal-binding; 2Fe-2S.
   CA   Ligand.
   //
   ID   3D-structure.
   AC   KW-0002
   DE   Protein, or part of a protein, whose three-dimensional structure has
   DE   been resolved experimentally (for example by X-ray crystallography or
   DE   NMR spectroscopy) and whose coordinates are available in the PDB
   DE   database. Can also be used for theoretical models.
   HI   Technical term: 3D-structure.
   CA   Technical term.
   //
   ID   3Fe-4S.
   ...

The entries in this file can be parsed by the ``parse`` function in the
``Bio.SwissProt.KeyWList`` module. Each entry is then stored as a
``Bio.SwissProt.KeyWList.Record``, which is a Python dictionary.

.. code:: pycon

   >>> from Bio.SwissProt import KeyWList
   >>> handle = open("keywlist.txt")
   >>> records = KeyWList.parse(handle)
   >>> for record in records:
   ...     print(record["ID"])
   ...     print(record["DE"])
   ...

This prints

.. code:: text

   2Fe-2S.
   Protein which contains at least one 2Fe-2S iron-sulfur cluster: 2 iron atoms
   complexed to 2 inorganic sulfides and 4 sulfur atoms of cysteines from the
   protein.
   ...

Parsing Prosite records
-----------------------

Prosite is a database containing protein domains, protein families,
functional sites, as well as the patterns and profiles to recognize
them. Prosite was developed in parallel with Swiss-Prot. In Biopython, a
Prosite record is represented by the ``Bio.ExPASy.Prosite.Record``
class, whose members correspond to the different fields in a Prosite
record.

In general, a Prosite file can contain more than one Prosite records.
For example, the full set of Prosite records, which can be downloaded as
a single file (``prosite.dat``) from the `ExPASy FTP
site <ftp://ftp.expasy.org/databases/prosite/prosite.dat>`__, contains
2073 records (version 20.24 released on 4 December 2007). To parse such
a file, we again make use of an iterator:

.. code:: pycon

   >>> from Bio.ExPASy import Prosite
   >>> handle = open("myprositefile.dat")
   >>> records = Prosite.parse(handle)

We can now take the records one at a time and print out some
information. For example, using the file containing the complete Prosite
database, we’d find

.. code:: pycon

   >>> from Bio.ExPASy import Prosite
   >>> handle = open("prosite.dat")
   >>> records = Prosite.parse(handle)
   >>> record = next(records)
   >>> record.accession
   'PS00001'
   >>> record.name
   'ASN_GLYCOSYLATION'
   >>> record.pdoc
   'PDOC00001'
   >>> record = next(records)
   >>> record.accession
   'PS00004'
   >>> record.name
   'CAMP_PHOSPHO_SITE'
   >>> record.pdoc
   'PDOC00004'
   >>> record = next(records)
   >>> record.accession
   'PS00005'
   >>> record.name
   'PKC_PHOSPHO_SITE'
   >>> record.pdoc
   'PDOC00005'

and so on. If you’re interested in how many Prosite records there are,
you could use

.. code:: pycon

   >>> from Bio.ExPASy import Prosite
   >>> handle = open("prosite.dat")
   >>> records = Prosite.parse(handle)
   >>> n = 0
   >>> for record in records:
   ...     n += 1
   ...
   >>> n
   2073

To read exactly one Prosite from the handle, you can use the ``read``
function:

.. code:: pycon

   >>> from Bio.ExPASy import Prosite
   >>> handle = open("mysingleprositerecord.dat")
   >>> record = Prosite.read(handle)

This function raises a ValueError if no Prosite record is found, and
also if more than one Prosite record is found.

Parsing Prosite documentation records
-------------------------------------

In the Prosite example above, the ``record.pdoc`` accession numbers
``'PDOC00001'``, ``'PDOC00004'``, ``'PDOC00005'`` and so on refer to
Prosite documentation. The Prosite documentation records are available
from ExPASy as individual files, and as one file (``prosite.doc``)
containing all Prosite documentation records.

We use the parser in ``Bio.ExPASy.Prodoc`` to parse Prosite
documentation records. For example, to create a list of all accession
numbers of Prosite documentation record, you can use

.. code:: pycon

   >>> from Bio.ExPASy import Prodoc
   >>> handle = open("prosite.doc")
   >>> records = Prodoc.parse(handle)
   >>> accessions = [record.accession for record in records]

Again a ``read()`` function is provided to read exactly one Prosite
documentation record from the handle.

Parsing Enzyme records
----------------------

ExPASy’s Enzyme database is a repository of information on enzyme
nomenclature. A typical Enzyme record looks as follows:

.. code:: text

   ID   3.1.1.34
   DE   Lipoprotein lipase.
   AN   Clearing factor lipase.
   AN   Diacylglycerol lipase.
   AN   Diglyceride lipase.
   CA   Triacylglycerol + H(2)O = diacylglycerol + a carboxylate.
   CC   -!- Hydrolyzes triacylglycerols in chylomicrons and very low-density
   CC       lipoproteins (VLDL).
   CC   -!- Also hydrolyzes diacylglycerol.
   PR   PROSITE; PDOC00110;
   DR   P11151, LIPL_BOVIN ;  P11153, LIPL_CAVPO ;  P11602, LIPL_CHICK ;
   DR   P55031, LIPL_FELCA ;  P06858, LIPL_HUMAN ;  P11152, LIPL_MOUSE ;
   DR   O46647, LIPL_MUSVI ;  P49060, LIPL_PAPAN ;  P49923, LIPL_PIG   ;
   DR   Q06000, LIPL_RAT   ;  Q29524, LIPL_SHEEP ;
   //

In this example, the first line shows the EC (Enzyme Commission) number
of lipoprotein lipase (second line). Alternative names of lipoprotein
lipase are "clearing factor lipase", "diacylglycerol lipase", and
"diglyceride lipase" (lines 3 through 5). The line starting with "CA"
shows the catalytic activity of this enzyme. Comment lines start with
"CC". The "PR" line shows references to the Prosite Documentation
records, and the "DR" lines show references to Swiss-Prot records. Not
of these entries are necessarily present in an Enzyme record.

In Biopython, an Enzyme record is represented by the
``Bio.ExPASy.Enzyme.Record`` class. This record derives from a Python
dictionary and has keys corresponding to the two-letter codes used in
Enzyme files. To read an Enzyme file containing one Enzyme record, use
the ``read`` function in ``Bio.ExPASy.Enzyme``:

.. doctest ../Tests/Enzymes

.. code:: pycon

   >>> from Bio.ExPASy import Enzyme
   >>> with open("lipoprotein.txt") as handle:
   ...     record = Enzyme.read(handle)
   ...
   >>> record["ID"]
   '3.1.1.34'
   >>> record["DE"]
   'Lipoprotein lipase.'
   >>> record["AN"]
   ['Clearing factor lipase.', 'Diacylglycerol lipase.', 'Diglyceride lipase.']
   >>> record["CA"]
   'Triacylglycerol + H(2)O = diacylglycerol + a carboxylate.'
   >>> record["PR"]
   ['PDOC00110']

.. code:: pycon

   >>> record["CC"]
   ['Hydrolyzes triacylglycerols in chylomicrons and very low-density lipoproteins
   (VLDL).', 'Also hydrolyzes diacylglycerol.']
   >>> record["DR"]
   [['P11151', 'LIPL_BOVIN'], ['P11153', 'LIPL_CAVPO'], ['P11602', 'LIPL_CHICK'],
   ['P55031', 'LIPL_FELCA'], ['P06858', 'LIPL_HUMAN'], ['P11152', 'LIPL_MOUSE'],
   ['O46647', 'LIPL_MUSVI'], ['P49060', 'LIPL_PAPAN'], ['P49923', 'LIPL_PIG'],
   ['Q06000', 'LIPL_RAT'], ['Q29524', 'LIPL_SHEEP']]

The ``read`` function raises a ValueError if no Enzyme record is found,
and also if more than one Enzyme record is found.

The full set of Enzyme records can be downloaded as a single file
(``enzyme.dat``) from the `ExPASy FTP
site <ftp://ftp.expasy.org/databases/enzyme/enzyme.dat>`__, containing
4877 records (release of 3 March 2009). To parse such a file containing
multiple Enzyme records, use the ``parse`` function in
``Bio.ExPASy.Enzyme`` to obtain an iterator:

.. code:: pycon

   >>> from Bio.ExPASy import Enzyme
   >>> handle = open("enzyme.dat")
   >>> records = Enzyme.parse(handle)

We can now iterate over the records one at a time. For example, we can
make a list of all EC numbers for which an Enzyme record is available:

.. code:: pycon

   >>> ecnumbers = [record["ID"] for record in records]

Accessing the ExPASy server
---------------------------

Swiss-Prot, Prosite, and Prosite documentation records can be downloaded
from the ExPASy web server at https://www.expasy.org. Four kinds of
queries are available from ExPASy:

get_prodoc_entry
   To download a Prosite documentation record in HTML format

get_prosite_entry
   To download a Prosite record in HTML format

get_prosite_raw
   To download a Prosite or Prosite documentation record in raw format

get_sprot_raw
   To download a Swiss-Prot record in raw format

To access this web server from a Python script, we use the
``Bio.ExPASy`` module.

.. _`sec:expasy_swissprot`:

Retrieving a Swiss-Prot record
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let’s say we are looking at chalcone synthases for Orchids (see
section :ref:`sec:orchids` for some justification for
looking for interesting things about orchids). Chalcone synthase is
involved in flavanoid biosynthesis in plants, and flavanoids make lots
of cool things like pigment colors and UV protectants.

If you do a search on Swiss-Prot, you can find three orchid proteins for
Chalcone Synthase, id numbers O23729, O23730, O23731. Now, let’s write a
script which grabs these, and parses out some interesting information.

First, we grab the records, using the ``get_sprot_raw()`` function of
``Bio.ExPASy``. This function is very nice since you can feed it an id
and get back a handle to a raw text record (no HTML to mess with!). We
can the use ``Bio.SwissProt.read`` to pull out the Swiss-Prot record, or
``Bio.SeqIO.read`` to get a SeqRecord. The following code accomplishes
what I just wrote:

.. code:: pycon

   >>> from Bio import ExPASy
   >>> from Bio import SwissProt

   >>> accessions = ["O23729", "O23730", "O23731"]
   >>> records = []

   >>> for accession in accessions:
   ...     handle = ExPASy.get_sprot_raw(accession)
   ...     record = SwissProt.read(handle)
   ...     records.append(record)
   ...

If the accession number you provided to ``ExPASy.get_sprot_raw`` does
not exist, then ``SwissProt.read(handle)`` will raise a ``ValueError``.
You can catch ``ValueException`` exceptions to detect invalid accession
numbers:

.. code:: pycon

   >>> for accession in accessions:
   ...     handle = ExPASy.get_sprot_raw(accession)
   ...     try:
   ...         record = SwissProt.read(handle)
   ...     except ValueException:
   ...         print("WARNING: Accession %s not found" % accession)
   ...     records.append(record)
   ...

Searching with UniProt
~~~~~~~~~~~~~~~~~~~~~~

Now, you may remark that I knew the records’ accession numbers
beforehand. Indeed, ``get_sprot_raw()`` needs either the entry name or
an accession number. When you don't have them handy, you can use
`UniProt <https://www.uniprot.org/>`_ to search for them.
You can also use the UniProt package to programmatically search for proteins.

For example, let's search for proteins from a specific organism (organism ID: 2697049)
that have been reviewed. We can do this with the following code:

.. code:: pycon

   >>> from Bio import UniProt

   >>> query = "(organism_id:2697049) AND (reviewed:true)"
   >>> results = list(UniProt.search(query))

The ``UniProt.search`` method returns an iterator over the search results.
The iterator returns one result at a time, fetching more results from UniProt as needed until all results are returned.
We can efficiently create a list from this iterator for this specific query
because this query only returns a few results (17 at the time of writing).

Let's try a search that returns more results.
At the time of writing, there are 5,147 results for the query "Insulin AND (reviewed:true)".
We can use slicing to get a list of the first 50 results.

.. code:: pycon

   >>> from Bio import UniProt
   >>> from itertools import islice

   >>> query = "Insulin AND (reviewed:true)"
   >>> results = UniProt.search(query, batch_size=50)[:50]

You can get the total number of search results (regardless of the batch size)
with the ``len`` method:

.. code:: pycon

   >>> from Bio import UniProt

   >>> query = "Insulin AND (reviewed:true)"
   >>> result_iterator = UniProt.search(query, batch_size=0)
   >>> len(result_iterator)
   5147

Retrieving Prosite and Prosite documentation records
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prosite and Prosite documentation records can be retrieved either in
HTML format, or in raw format. To parse Prosite and Prosite
documentation records with Biopython, you should retrieve the records in
raw format. For other purposes, however, you may be interested in these
records in HTML format.

To retrieve a Prosite or Prosite documentation record in raw format, use
``get_prosite_raw()``. For example, to download a Prosite record and
print it out in raw text format, use

.. code:: pycon

   >>> from Bio import ExPASy
   >>> handle = ExPASy.get_prosite_raw("PS00001")
   >>> text = handle.read()
   >>> print(text)

To retrieve a Prosite record and parse it into a ``Bio.Prosite.Record``
object, use

.. code:: pycon

   >>> from Bio import ExPASy
   >>> from Bio import Prosite
   >>> handle = ExPASy.get_prosite_raw("PS00001")
   >>> record = Prosite.read(handle)

The same function can be used to retrieve a Prosite documentation record
and parse it into a ``Bio.ExPASy.Prodoc.Record`` object:

.. code:: pycon

   >>> from Bio import ExPASy
   >>> from Bio.ExPASy import Prodoc
   >>> handle = ExPASy.get_prosite_raw("PDOC00001")
   >>> record = Prodoc.read(handle)

For non-existing accession numbers, ``ExPASy.get_prosite_raw`` returns a
handle to an empty string. When faced with an empty string,
``Prosite.read`` and ``Prodoc.read`` will raise a ValueError. You can
catch these exceptions to detect invalid accession numbers.

The functions ``get_prosite_entry()`` and ``get_prodoc_entry()`` are
used to download Prosite and Prosite documentation records in HTML
format. To create a web page showing one Prosite record, you can use

.. code:: pycon

   >>> from Bio import ExPASy
   >>> handle = ExPASy.get_prosite_entry("PS00001")
   >>> html = handle.read()
   >>> with open("myprositerecord.html", "w") as out_handle:
   ...     out_handle.write(html)
   ...

and similarly for a Prosite documentation record:

.. code:: pycon

   >>> from Bio import ExPASy
   >>> handle = ExPASy.get_prodoc_entry("PDOC00001")
   >>> html = handle.read()
   >>> with open("myprodocrecord.html", "w") as out_handle:
   ...     out_handle.write(html)
   ...

For these functions, an invalid accession number returns an error
message in HTML format.

Scanning the Prosite database
-----------------------------

`ScanProsite <https://prosite.expasy.org/prosite.html>`__ allows you to
scan protein sequences online against the Prosite database by providing
a UniProt or PDB sequence identifier or the sequence itself. For more
information about ScanProsite, please see the `ScanProsite
documentation <https://prosite.expasy.org/prosite_doc.html>`__ as well
as the `documentation for programmatic access of
ScanProsite <https://prosite.expasy.org/scanprosite/scanprosite_doc.html#rest>`__.

You can use Biopython’s ``Bio.ExPASy.ScanProsite`` module to scan the
Prosite database from Python. This module both helps you to access
ScanProsite programmatically, and to parse the results returned by
ScanProsite. To scan for Prosite patterns in the following protein
sequence:

.. code:: text

   MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFT
   CRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN

you can use the following code:

.. doctest . internet

.. code:: pycon

   >>> sequence = (
   ...     "MEHKEVVLLLLLFLKSGQGEPLDDYVNTQGASLFSVTKKQLGAGSIEECAAKCEEDEEFT"
   ...     "CRAFQYHSKEQQCVIMAENRKSSIIIRMRDVVLFEKKVYLSECKTGNGKNYRGTMSKTKN"
   ... )
   >>> from Bio.ExPASy import ScanProsite
   >>> handle = ScanProsite.scan(seq=sequence)

By executing ``handle.read()``, you can obtain the search results in raw
XML format. Instead, let’s use ``Bio.ExPASy.ScanProsite.read`` to parse
the raw XML into a Python object:

.. cont-doctest

.. code:: pycon

   >>> result = ScanProsite.read(handle)
   >>> type(result)
   <class 'Bio.ExPASy.ScanProsite.Record'>

A ``Bio.ExPASy.ScanProsite.Record`` object is derived from a list, with
each element in the list storing one ScanProsite hit. This object also
stores the number of hits, as well as the number of search sequences, as
returned by ScanProsite. This ScanProsite search resulted in six hits:

.. cont-doctest

.. code:: pycon

   >>> result.n_seq
   1
   >>> result.n_match
   1
   >>> len(result)
   1
   >>> result[0]
   {'sequence_ac': 'USERSEQ1', 'start': 16, 'stop': 98, 'signature_ac': 'PS50948', 'score': '8.873', 'level': '0'}

Other ScanProsite parameters can be passed as keyword arguments; see the
`documentation for programmatic access of
ScanProsite <https://prosite.expasy.org/scanprosite/scanprosite_doc.html#rest>`__
for more information. As an example, passing ``lowscore=1`` to include
matches with low level scores lets use find one additional hit:

.. cont-doctest

.. code:: pycon

   >>> handle = ScanProsite.scan(seq=sequence, lowscore=1)
   >>> result = ScanProsite.read(handle)
   >>> result.n_match
   2
