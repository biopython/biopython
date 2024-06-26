.. _`chapter:searchio`:

BLAST and other sequence search tools
=====================================

Biological sequence identification is an integral part of
bioinformatics. Several tools are available for this, each with their
own algorithms and approaches, such as BLAST (arguably the most
popular), FASTA, HMMER, and many more. In general, these tools usually
use your sequence to search a database of potential matches. With the
growing number of known sequences (hence the growing number of potential
matches), interpreting the results becomes increasingly hard as there
could be hundreds or even thousands of potential matches. Naturally,
manual interpretation of these searches’ results is out of the question.
Moreover, you often need to work with several sequence search tools,
each with its own statistics, conventions, and output format. Imagine
how daunting it would be when you need to work with multiple sequences
using multiple search tools.

We know this too well ourselves, which is why we created the
``Bio.SearchIO`` submodule in Biopython. ``Bio.SearchIO`` allows you to
extract information from your search results in a convenient way, while
also dealing with the different standards and conventions used by
different search tools. The name ``SearchIO`` is a homage to BioPerl’s
module of the same name.

In this chapter, we’ll go through the main features of ``Bio.SearchIO``
to show what it can do for you. We’ll use two popular search tools along
the way: BLAST and BLAT. They are used merely for illustrative purposes,
and you should be able to adapt the workflow to any other search tools
supported by ``Bio.SearchIO`` in a breeze. You’re very welcome to follow
along with the search output files we’ll be using. The BLAST output file
can be downloaded
`here <https://github.com/biopython/biopython/blob/master/Doc/examples/my_blast.xml>`__,
and the BLAT output file
`here <https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/my_blat.psl>`__
or are included with the Biopython source code under the
``Doc/examples/`` folder. Both output files were generated using this
sequence:

.. code:: text

   >mystery_seq
   CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

The BLAST result is an XML file generated using ``blastn`` against the
NCBI ``refseq_rna`` database. For BLAT, the sequence database was the
February 2009 ``hg19`` human genome draft and the output format is PSL.

We’ll start from an introduction to the ``Bio.SearchIO`` object model.
The model is the representation of your search results, thus it is core
to ``Bio.SearchIO`` itself. After that, we’ll check out the main
functions in ``Bio.SearchIO`` that you may often use.

Now that we’re all set, let’s go to the first step: introducing the core
object model.

.. _`sec:searchio-model`:

The SearchIO object model
-------------------------

Despite the wildly differing output styles among many sequence search
tools, it turns out that their underlying concept is similar:

-  The output file may contain results from one or more search queries.

-  In each search query, you will see one or more hits from the given
   search database.

-  In each database hit, you will see one or more regions containing the
   actual sequence alignment between your query sequence and the
   database sequence.

-  Some programs like BLAT or Exonerate may further split these regions
   into several alignment fragments (or blocks in BLAT and possibly
   exons in exonerate). This is not something you always see, as
   programs like BLAST and HMMER do not do this.

Realizing this generality, we decided use it as base for creating the
``Bio.SearchIO`` object model. The object model consists of a nested
hierarchy of Python objects, each one representing one concept outlined
above. These objects are:

-  ``QueryResult``, to represent a single search query.

-  ``Hit``, to represent a single database hit. ``Hit`` objects are
   contained within ``QueryResult`` and in each ``QueryResult`` there is
   zero or more ``Hit`` objects.

-  ``HSP`` (short for high-scoring pair), to represent region(s) of
   significant alignments between query and hit sequences. ``HSP``
   objects are contained within ``Hit`` objects and each ``Hit`` has one
   or more ``HSP`` objects.

-  ``HSPFragment``, to represent a single contiguous alignment between
   query and hit sequences. ``HSPFragment`` objects are contained within
   ``HSP`` objects. Most sequence search tools like BLAST and HMMER
   unify ``HSP`` and ``HSPFragment`` objects as each ``HSP`` will only
   have a single ``HSPFragment``. However there are tools like BLAT and
   Exonerate that produce ``HSP`` containing multiple ``HSPFragment``.
   Don’t worry if this seems a tad confusing now, we’ll elaborate more
   on these two objects later on.

These four objects are the ones you will interact with when you use
``Bio.SearchIO``. They are created using one of the main
``Bio.SearchIO`` methods: ``read``, ``parse``, ``index``, or
``index_db``. The details of these methods are provided in later
sections. For this section, we’ll only be using read and parse. These
functions behave similarly to their ``Bio.SeqIO`` and ``Bio.AlignIO``
counterparts:

-  ``read`` is used for search output files with a single query and
   returns a ``QueryResult`` object

-  ``parse`` is used for search output files with multiple queries and
   returns a generator that yields ``QueryResult`` objects

With that settled, let’s start probing each ``Bio.SearchIO`` object,
beginning with ``QueryResult``.

.. _`sec:searchio-qresult`:

QueryResult
~~~~~~~~~~~

The QueryResult object represents a single search query and contains
zero or more Hit objects. Let’s see what it looks like using the BLAST
file we have:

.. doctest examples

.. code:: pycon

   >>> from Bio import SearchIO
   >>> blast_qresult = SearchIO.read("my_blast.xml", "blast-xml")
   >>> print(blast_qresult)
   Program: blastn (2.2.27+)
     Query: 42291 (61)
            mystery_seq
    Target: refseq_rna
      Hits: ----  -----  ----------------------------------------------------------
               #  # HSP  ID + description
            ----  -----  ----------------------------------------------------------
               0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
               1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
               2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...
               3      2  gi|301171322|ref|NR_035857.1|  Pan troglodytes microRNA...
               4      1  gi|301171267|ref|NR_035851.1|  Pan troglodytes microRNA...
               5      2  gi|262205330|ref|NR_030198.1|  Homo sapiens microRNA 52...
               6      1  gi|262205302|ref|NR_030191.1|  Homo sapiens microRNA 51...
               7      1  gi|301171259|ref|NR_035850.1|  Pan troglodytes microRNA...
               8      1  gi|262205451|ref|NR_030222.1|  Homo sapiens microRNA 51...
               9      2  gi|301171447|ref|NR_035871.1|  Pan troglodytes microRNA...
              10      1  gi|301171276|ref|NR_035852.1|  Pan troglodytes microRNA...
              11      1  gi|262205290|ref|NR_030188.1|  Homo sapiens microRNA 51...
              12      1  gi|301171354|ref|NR_035860.1|  Pan troglodytes microRNA...
              13      1  gi|262205281|ref|NR_030186.1|  Homo sapiens microRNA 52...
              14      2  gi|262205298|ref|NR_030190.1|  Homo sapiens microRNA 52...
              15      1  gi|301171394|ref|NR_035865.1|  Pan troglodytes microRNA...
              16      1  gi|262205429|ref|NR_030218.1|  Homo sapiens microRNA 51...
              17      1  gi|262205423|ref|NR_030217.1|  Homo sapiens microRNA 52...
              18      1  gi|301171401|ref|NR_035866.1|  Pan troglodytes microRNA...
              19      1  gi|270133247|ref|NR_032574.1|  Macaca mulatta microRNA ...
              20      1  gi|262205309|ref|NR_030193.1|  Homo sapiens microRNA 52...
              21      2  gi|270132717|ref|NR_032716.1|  Macaca mulatta microRNA ...
              22      2  gi|301171437|ref|NR_035870.1|  Pan troglodytes microRNA...
              23      2  gi|270133306|ref|NR_032587.1|  Macaca mulatta microRNA ...
              24      2  gi|301171428|ref|NR_035869.1|  Pan troglodytes microRNA...
              25      1  gi|301171211|ref|NR_035845.1|  Pan troglodytes microRNA...
              26      2  gi|301171153|ref|NR_035838.1|  Pan troglodytes microRNA...
              27      2  gi|301171146|ref|NR_035837.1|  Pan troglodytes microRNA...
              28      2  gi|270133254|ref|NR_032575.1|  Macaca mulatta microRNA ...
              29      2  gi|262205445|ref|NR_030221.1|  Homo sapiens microRNA 51...
              ~~~
              97      1  gi|356517317|ref|XM_003527287.1|  PREDICTED: Glycine ma...
              98      1  gi|297814701|ref|XM_002875188.1|  Arabidopsis lyrata su...
              99      1  gi|397513516|ref|XM_003827011.1|  PREDICTED: Pan panisc...

We’ve just begun to scratch the surface of the object model, but you can
see that there’s already some useful information. By invoking ``print``
on the ``QueryResult`` object, you can see:

-  The program name and version (blastn version 2.2.27+)

-  The query ID, description, and its sequence length (ID is 42291,
   description is ‘mystery_seq’, and it is 61 nucleotides long)

-  The target database to search against (refseq_rna)

-  A quick overview of the resulting hits. For our query sequence, there
   are 100 potential hits (numbered 0–99 in the table). For each hit, we
   can also see how many HSPs it contains, its ID, and a snippet of its
   description. Notice here that ``Bio.SearchIO`` truncates the hit
   table overview, by showing only hits numbered 0–29, and then 97–99.

Now let’s check our BLAT results using the same procedure as above:

.. cont-doctest

.. code:: pycon

   >>> blat_qresult = SearchIO.read("my_blat.psl", "blat-psl")
   >>> print(blat_qresult)
   Program: blat (<unknown version>)
     Query: mystery_seq (61)
            <unknown description>
    Target: <unknown target>
      Hits: ----  -----  ----------------------------------------------------------
               #  # HSP  ID + description
            ----  -----  ----------------------------------------------------------
               0     17  chr19  <unknown description>

You’ll immediately notice that there are some differences. Some of these
are caused by the way PSL format stores its details, as you’ll see. The
rest are caused by the genuine program and target database differences
between our BLAST and BLAT searches:

-  The program name and version. ``Bio.SearchIO`` knows that the program
   is BLAT, but in the output file there is no information regarding the
   program version so it defaults to ‘<unknown version>’.

-  The query ID, description, and its sequence length. Notice here that
   these details are slightly different from the ones we saw in BLAST.
   The ID is ‘mystery_seq’ instead of 42991, there is no known
   description, but the query length is still 61. This is actually a
   difference introduced by the file formats themselves. BLAST sometimes
   creates its own query IDs and uses your original ID as the sequence
   description.

-  The target database is not known, as it is not stated in the BLAT
   output file.

-  And finally, the list of hits we have is completely different. Here,
   we see that our query sequence only hits the ‘chr19’ database entry,
   but in it we see 17 HSP regions. This should not be surprising
   however, given that we are using a different program, each with its
   own target database.

All the details you saw when invoking the ``print`` method can be
accessed individually using Python’s object attribute access notation
(a.k.a. the dot notation). There are also other format-specific
attributes that you can access using the same method.

.. cont-doctest

.. code:: pycon

   >>> print("%s %s" % (blast_qresult.program, blast_qresult.version))
   blastn 2.2.27+
   >>> print("%s %s" % (blat_qresult.program, blat_qresult.version))
   blat <unknown version>
   >>> blast_qresult.param_evalue_threshold  # blast-xml specific
   10.0

For a complete list of accessible attributes, you can check each
format-specific documentation. e.g. :py:mod:`Bio.SearchIO.BlastIO`
and :py:mod:`Bio.SearchIO.BlatIO`.

Having looked at using ``print`` on ``QueryResult`` objects, let’s drill
down deeper. What exactly is a ``QueryResult``? In terms of Python
objects, ``QueryResult`` is a hybrid between a list and a dictionary. In
other words, it is a container object with all the convenient features
of lists and dictionaries.

Like Python lists and dictionaries, ``QueryResult`` objects are
iterable. Each iteration returns a ``Hit`` object:

.. code:: pycon

   >>> for hit in blast_qresult:
   ...     hit
   ...
   Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)
   Hit(id='gi|301171311|ref|NR_035856.1|', query_id='42291', 1 hsps)
   Hit(id='gi|270133242|ref|NR_032573.1|', query_id='42291', 1 hsps)
   Hit(id='gi|301171322|ref|NR_035857.1|', query_id='42291', 2 hsps)
   Hit(id='gi|301171267|ref|NR_035851.1|', query_id='42291', 1 hsps)
   ...

To check how many items (hits) a ``QueryResult`` has, you can simply
invoke Python’s ``len`` method:

.. cont-doctest

.. code:: pycon

   >>> len(blast_qresult)
   100
   >>> len(blat_qresult)
   1

Like Python lists, you can retrieve items (hits) from a ``QueryResult``
using the slice notation:

.. cont-doctest

.. code:: pycon

   >>> blast_qresult[0]  # retrieves the top hit
   Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)
   >>> blast_qresult[-1]  # retrieves the last hit
   Hit(id='gi|397513516|ref|XM_003827011.1|', query_id='42291', 1 hsps)

To retrieve multiple hits, you can slice ``QueryResult`` objects using
the slice notation as well. In this case, the slice will return a new
``QueryResult`` object containing only the sliced hits:

.. cont-doctest

.. code:: pycon

   >>> blast_slice = blast_qresult[:3]  # slices the first three hits
   >>> print(blast_slice)
   Program: blastn (2.2.27+)
     Query: 42291 (61)
            mystery_seq
    Target: refseq_rna
      Hits: ----  -----  ----------------------------------------------------------
               #  # HSP  ID + description
            ----  -----  ----------------------------------------------------------
               0      1  gi|262205317|ref|NR_030195.1|  Homo sapiens microRNA 52...
               1      1  gi|301171311|ref|NR_035856.1|  Pan troglodytes microRNA...
               2      1  gi|270133242|ref|NR_032573.1|  Macaca mulatta microRNA ...

Like Python dictionaries, you can also retrieve hits using the hit’s ID.
This is particularly useful if you know a given hit ID exists within a
search query results:

.. cont-doctest

.. code:: pycon

   >>> blast_qresult["gi|262205317|ref|NR_030195.1|"]
   Hit(id='gi|262205317|ref|NR_030195.1|', query_id='42291', 1 hsps)

You can also get a full list of ``Hit`` objects using ``hits`` and a
full list of ``Hit`` IDs using ``hit_keys``:

.. code:: pycon

   >>> blast_qresult.hits
   [...]       # list of all hits
   >>> blast_qresult.hit_keys
   [...]       # list of all hit IDs

What if you just want to check whether a particular hit is present in
the query results? You can do a simple Python membership test using the
``in`` keyword:

.. cont-doctest

.. code:: pycon

   >>> "gi|262205317|ref|NR_030195.1|" in blast_qresult
   True
   >>> "gi|262205317|ref|NR_030194.1|" in blast_qresult
   False

Sometimes, knowing whether a hit is present is not enough; you also want
to know the rank of the hit. Here, the ``index`` method comes to the
rescue:

.. cont-doctest

.. code:: pycon

   >>> blast_qresult.index("gi|301171437|ref|NR_035870.1|")
   22

Remember that we’re using Python’s indexing style here, which is
zero-based. This means our hit above is ranked at no. 23, not 22.

Also, note that the hit rank you see here is based on the native hit
ordering present in the original search output file. Different search
tools may order these hits based on different criteria.

If the native hit ordering doesn’t suit your taste, you can use the
``sort`` method of the ``QueryResult`` object. It is very similar to
Python’s ``list.sort`` method, with the addition of an option to create
a new sorted ``QueryResult`` object or not.

Here is an example of using ``QueryResult.sort`` to sort the hits based
on each hit’s full sequence length. For this particular sort, we’ll set
the ``in_place`` flag to ``False`` so that sorting will return a new
``QueryResult`` object and leave our initial object unsorted. We’ll also
set the ``reverse`` flag to ``True`` so that we sort in descending
order.

.. cont-doctest

.. code:: pycon

   >>> for hit in blast_qresult[:5]:  # id and sequence length of the first five hits
   ...     print("%s %i" % (hit.id, hit.seq_len))
   ...
   gi|262205317|ref|NR_030195.1| 61
   gi|301171311|ref|NR_035856.1| 60
   gi|270133242|ref|NR_032573.1| 85
   gi|301171322|ref|NR_035857.1| 86
   gi|301171267|ref|NR_035851.1| 80

   >>> sort_key = lambda hit: hit.seq_len
   >>> sorted_qresult = blast_qresult.sort(key=sort_key, reverse=True, in_place=False)
   >>> for hit in sorted_qresult[:5]:
   ...     print("%s %i" % (hit.id, hit.seq_len))
   ...
   gi|397513516|ref|XM_003827011.1| 6002
   gi|390332045|ref|XM_776818.2| 4082
   gi|390332043|ref|XM_003723358.1| 4079
   gi|356517317|ref|XM_003527287.1| 3251
   gi|356543101|ref|XM_003539954.1| 2936

The advantage of having the ``in_place`` flag here is that we’re
preserving the native ordering, so we may use it again later. You should
note that this is not the default behavior of ``QueryResult.sort``,
however, which is why we needed to set the ``in_place`` flag to ``True``
explicitly.

At this point, you’ve known enough about ``QueryResult`` objects to make
it work for you. But before we go on to the next object in the
``Bio.SearchIO`` model, let’s take a look at two more sets of methods
that could make it even easier to work with ``QueryResult`` objects: the
``filter`` and ``map`` methods.

If you’re familiar with Python’s list comprehensions, generator
expressions or the built-in ``filter`` and ``map`` functions, you’ll
know how useful they are for working with list-like objects (if you’re
not, check them out!). You can use these built-in methods to manipulate
``QueryResult`` objects, but you’ll end up with regular Python lists and
lose the ability to do more interesting manipulations.

That’s why, ``QueryResult`` objects provide its own flavor of ``filter``
and ``map`` methods. Analogous to ``filter``, there are ``hit_filter``
and ``hsp_filter`` methods. As their name implies, these methods filter
its ``QueryResult`` object either on its ``Hit`` objects or ``HSP``
objects. Similarly, analogous to ``map``, ``QueryResult`` objects also
provide the ``hit_map`` and ``hsp_map`` methods. These methods apply a
given function to all hits or HSPs in a ``QueryResult`` object,
respectively.

Let’s see these methods in action, beginning with ``hit_filter``. This
method accepts a callback function that checks whether a given ``Hit``
object passes the condition you set or not. In other words, the function
must accept as its argument a single ``Hit`` object and returns ``True``
or ``False``.

Here is an example of using ``hit_filter`` to filter out ``Hit`` objects
that only have one HSP:

.. cont-doctest

.. code:: pycon

   >>> filter_func = lambda hit: len(hit.hsps) > 1  # the callback function
   >>> len(blast_qresult)  # no. of hits before filtering
   100
   >>> filtered_qresult = blast_qresult.hit_filter(filter_func)
   >>> len(filtered_qresult)  # no. of hits after filtering
   37
   >>> for hit in filtered_qresult[:5]:  # quick check for the hit lengths
   ...     print("%s %i" % (hit.id, len(hit.hsps)))
   ...
   gi|301171322|ref|NR_035857.1| 2
   gi|262205330|ref|NR_030198.1| 2
   gi|301171447|ref|NR_035871.1| 2
   gi|262205298|ref|NR_030190.1| 2
   gi|270132717|ref|NR_032716.1| 2

``hsp_filter`` works the same as ``hit_filter``, only instead of looking
at the ``Hit`` objects, it performs filtering on the ``HSP`` objects in
each hits.

As for the ``map`` methods, they too accept a callback function as their
arguments. However, instead of returning ``True`` or ``False``, the
callback function must return the modified ``Hit`` or ``HSP`` object
(depending on whether you’re using ``hit_map`` or ``hsp_map``).

Let’s see an example where we’re using ``hit_map`` to rename the hit
IDs:

.. cont-doctest

.. code:: pycon

   >>> def map_func(hit):
   ...     # renames "gi|301171322|ref|NR_035857.1|" to "NR_035857.1"
   ...     hit.id = hit.id.split("|")[3]
   ...     return hit
   ...
   >>> mapped_qresult = blast_qresult.hit_map(map_func)
   >>> for hit in mapped_qresult[:5]:
   ...     print(hit.id)
   ...
   NR_030195.1
   NR_035856.1
   NR_032573.1
   NR_035857.1
   NR_035851.1

Again, ``hsp_map`` works the same as ``hit_map``, but on ``HSP`` objects
instead of ``Hit`` objects.

.. _`sec:searchio-hit`:

Hit
~~~

``Hit`` objects represent all query results from a single database
entry. They are the second-level container in the ``Bio.SearchIO``
object hierarchy. You’ve seen that they are contained by ``QueryResult``
objects, but they themselves contain ``HSP`` objects.

Let’s see what they look like, beginning with our BLAST search:

.. doctest examples

.. code:: pycon

   >>> from Bio import SearchIO
   >>> blast_qresult = SearchIO.read("my_blast.xml", "blast-xml")
   >>> blast_hit = blast_qresult[3]  # fourth hit from the query result
   >>> print(blast_hit)
   Query: 42291
          mystery_seq
     Hit: gi|301171322|ref|NR_035857.1| (86)
          Pan troglodytes microRNA mir-520c (MIR520C), microRNA
    HSPs: ----  --------  ---------  ------  ---------------  ---------------------
             #   E-value  Bit score    Span      Query range              Hit range
          ----  --------  ---------  ------  ---------------  ---------------------
             0   8.9e-20     100.47      60           [1:61]                [13:73]
             1   3.3e-06      55.39      60           [0:60]                [13:73]

You see that we’ve got the essentials covered here:

-  The query ID and description is present. A hit is always tied to a
   query, so we want to keep track of the originating query as well.
   These values can be accessed from a hit using the ``query_id`` and
   ``query_description`` attributes.

-  We also have the unique hit ID, description, and full sequence
   lengths. They can be accessed using ``id``, ``description``, and
   ``seq_len``, respectively.

-  Finally, there’s a table containing quick information about the HSPs
   this hit contains. In each row, we’ve got the important HSP details
   listed: the HSP index, its e-value, its bit score, its span (the
   alignment length including gaps), its query coordinates, and its hit
   coordinates.

Now let’s contrast this with the BLAT search. Remember that in the BLAT
search we had one hit with 17 HSPs.

.. cont-doctest

.. code:: pycon

   >>> blat_qresult = SearchIO.read("my_blat.psl", "blat-psl")
   >>> blat_hit = blat_qresult[0]  # the only hit
   >>> print(blat_hit)
   Query: mystery_seq
          <unknown description>
     Hit: chr19 (59128983)
          <unknown description>
    HSPs: ----  --------  ---------  ------  ---------------  ---------------------
             #   E-value  Bit score    Span      Query range              Hit range
          ----  --------  ---------  ------  ---------------  ---------------------
             0         ?          ?       ?           [0:61]    [54204480:54204541]
             1         ?          ?       ?           [0:61]    [54233104:54264463]
             2         ?          ?       ?           [0:61]    [54254477:54260071]
             3         ?          ?       ?           [1:61]    [54210720:54210780]
             4         ?          ?       ?           [0:60]    [54198476:54198536]
             5         ?          ?       ?           [0:61]    [54265610:54265671]
             6         ?          ?       ?           [0:61]    [54238143:54240175]
             7         ?          ?       ?           [0:60]    [54189735:54189795]
             8         ?          ?       ?           [0:61]    [54185425:54185486]
             9         ?          ?       ?           [0:60]    [54197657:54197717]
            10         ?          ?       ?           [0:61]    [54255662:54255723]
            11         ?          ?       ?           [0:61]    [54201651:54201712]
            12         ?          ?       ?           [8:60]    [54206009:54206061]
            13         ?          ?       ?          [10:61]    [54178987:54179038]
            14         ?          ?       ?           [8:61]    [54212018:54212071]
            15         ?          ?       ?           [8:51]    [54234278:54234321]
            16         ?          ?       ?           [8:61]    [54238143:54238196]

Here, we’ve got a similar level of detail as with the BLAST hit we saw
earlier. There are some differences worth explaining, though:

-  The e-value and bit score column values. As BLAT HSPs do not have
   e-values and bit scores, the display defaults to ‘?’.

-  What about the span column? The span values is meant to display the
   complete alignment length, which consists of all residues and any
   gaps that may be present. The PSL format do not have this information
   readily available and ``Bio.SearchIO`` does not attempt to try guess
   what it is, so we get a ‘?’ similar to the e-value and bit score
   columns.

In terms of Python objects, ``Hit`` behaves almost the same as Python
lists, but contain ``HSP`` objects exclusively. If you’re familiar with
lists, you should encounter no difficulties working with the ``Hit``
object.

Just like Python lists, ``Hit`` objects are iterable, and each iteration
returns one ``HSP`` object it contains:

.. cont-doctest

.. code:: pycon

   >>> for hsp in blast_hit:
   ...     hsp
   ...
   HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='42291', 1 fragments)
   HSP(hit_id='gi|301171322|ref|NR_035857.1|', query_id='42291', 1 fragments)

You can invoke ``len`` on a ``Hit`` to see how many ``HSP`` objects it
has:

.. cont-doctest

.. code:: pycon

   >>> len(blast_hit)
   2
   >>> len(blat_hit)
   17

You can use the slice notation on ``Hit`` objects, whether to retrieve
single ``HSP`` or multiple ``HSP`` objects. Like ``QueryResult``, if you
slice for multiple ``HSP``, a new ``Hit`` object will be returned
containing only the sliced ``HSP`` objects:

.. cont-doctest

.. code:: pycon

   >>> blat_hit[0]  # retrieve single items
   HSP(hit_id='chr19', query_id='mystery_seq', 1 fragments)
   >>> sliced_hit = blat_hit[4:9]  # retrieve multiple items
   >>> len(sliced_hit)
   5
   >>> print(sliced_hit)
   Query: mystery_seq
          <unknown description>
     Hit: chr19 (59128983)
          <unknown description>
    HSPs: ----  --------  ---------  ------  ---------------  ---------------------
             #   E-value  Bit score    Span      Query range              Hit range
          ----  --------  ---------  ------  ---------------  ---------------------
             0         ?          ?       ?           [0:60]    [54198476:54198536]
             1         ?          ?       ?           [0:61]    [54265610:54265671]
             2         ?          ?       ?           [0:61]    [54238143:54240175]
             3         ?          ?       ?           [0:60]    [54189735:54189795]
             4         ?          ?       ?           [0:61]    [54185425:54185486]

You can also sort the ``HSP`` inside a ``Hit``, using the exact same
arguments like the sort method you saw in the ``QueryResult`` object.

Finally, there are also the ``filter`` and ``map`` methods you can use
on ``Hit`` objects. Unlike in the ``QueryResult`` object, ``Hit``
objects only have one variant of ``filter`` (``Hit.filter``) and one
variant of ``map`` (``Hit.map``). Both of ``Hit.filter`` and ``Hit.map``
work on the ``HSP`` objects a ``Hit`` has.

.. _`sec:searchio-hsp`:

HSP
~~~

``HSP`` (high-scoring pair) represents region(s) in the hit sequence
that contains significant alignment(s) to the query sequence. It
contains the actual match between your query sequence and a database
entry. As this match is determined by the sequence search tool’s
algorithms, the ``HSP`` object contains the bulk of the statistics
computed by the search tool. This also makes the distinction between
``HSP`` objects from different search tools more apparent compared to
the differences you’ve seen in ``QueryResult`` or ``Hit`` objects.

Let’s see some examples from our BLAST and BLAT searches. We’ll look at
the BLAST HSP first:

.. doctest examples

.. code:: pycon

   >>> from Bio import SearchIO
   >>> blast_qresult = SearchIO.read("my_blast.xml", "blast-xml")
   >>> blast_hsp = blast_qresult[0][0]  # first hit, first hsp
   >>> print(blast_hsp)
         Query: 42291 mystery_seq
           Hit: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520...
   Query range: [0:61] (1)
     Hit range: [0:61] (1)
   Quick stats: evalue 4.9e-23; bitscore 111.29
     Fragments: 1 (61 columns)
        Query - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          Hit - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

Just like ``QueryResult`` and ``Hit``, invoking ``print`` on an ``HSP``
shows its general details:

-  There are the query and hit IDs and descriptions. We need these to
   identify our ``HSP``.

-  We’ve also got the matching range of the query and hit sequences. The
   slice notation we’re using here is an indication that the range is
   displayed using Python’s indexing style (zero-based, half open). The
   number inside the parenthesis denotes the strand. In this case, both
   sequences have the plus strand.

-  Some quick statistics are available: the e-value and bitscore.

-  There is information about the HSP fragments. Ignore this for now; it
   will be explained later on.

-  And finally, we have the query and hit sequence alignment itself.

These details can be accessed on their own using the dot notation, just
like in ``QueryResult`` and ``Hit``:

.. cont-doctest

.. code:: pycon

   >>> blast_hsp.query_range
   (0, 61)

.. code:: pycon

   >>> blast_hsp.evalue
   4.91307e-23

They’re not the only attributes available, though. ``HSP`` objects come
with a default set of properties that makes it easy to probe their
various details. Here are some examples:

.. cont-doctest

.. code:: pycon

   >>> blast_hsp.hit_start  # start coordinate of the hit sequence
   0
   >>> blast_hsp.query_span  # how many residues in the query sequence
   61
   >>> blast_hsp.aln_span  # how long the alignment is
   61

Check out the ``HSP`` documentation under :py:mod:`Bio.SearchIO`
for a full list of these predefined properties.

Furthermore, each sequence search tool usually computes its own
statistics / details for its ``HSP`` objects. For example, an XML BLAST
search also outputs the number of gaps and identical residues. These
attributes can be accessed like so:

.. cont-doctest

.. code:: pycon

   >>> blast_hsp.gap_num  # number of gaps
   0
   >>> blast_hsp.ident_num  # number of identical residues
   61

These details are format-specific; they may not be present in other
formats. To see which details are available for a given sequence search
tool, you should check the format’s documentation in ``Bio.SearchIO``.
Alternatively, you may also use ``.__dict__.keys()`` for a quick list of
what’s available:

.. code:: pycon

   >>> blast_hsp.__dict__.keys()
   ['bitscore', 'evalue', 'ident_num', 'gap_num', 'bitscore_raw', 'pos_num', '_items']

Finally, you may have noticed that the ``query`` and ``hit`` attributes
of our HSP are not just regular strings:

.. cont-doctest

.. code:: pycon

   >>> blast_hsp.query
   SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG'), id='42291', name='aligned query sequence', description='mystery_seq', dbxrefs=[])
   >>> blast_hsp.hit
   SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG'), id='gi|262205317|ref|NR_030195.1|', name='aligned hit sequence', description='Homo sapiens microRNA 520b (MIR520B), microRNA', dbxrefs=[])

They are ``SeqRecord`` objects you saw earlier in
Section :ref:`chapter:seq_annot`! This means that
you can do all sorts of interesting things you can do with ``SeqRecord``
objects on ``HSP.query`` and/or ``HSP.hit``.

It should not surprise you now that the ``HSP`` object has an
``alignment`` property which is a ``MultipleSeqAlignment`` object:

.. cont-doctest

.. code:: pycon

   >>> print(blast_hsp.aln)
   Alignment with 2 rows and 61 columns
   CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG 42291
   CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAG...GGG gi|262205317|ref|NR_030195.1|

Having probed the BLAST HSP, let’s now take a look at HSPs from our BLAT
results for a different kind of HSP. As usual, we’ll begin by invoking
``print`` on it:

.. cont-doctest

.. code:: pycon

   >>> blat_qresult = SearchIO.read("my_blat.psl", "blat-psl")
   >>> blat_hsp = blat_qresult[0][0]  # first hit, first hsp
   >>> print(blat_hsp)
         Query: mystery_seq <unknown description>
           Hit: chr19 <unknown description>
   Query range: [0:61] (1)
     Hit range: [54204480:54204541] (1)
   Quick stats: evalue ?; bitscore ?
     Fragments: 1 (? columns)

Some of the outputs you may have already guessed. We have the query and
hit IDs and descriptions and the sequence coordinates. Values for evalue
and bitscore is ‘?’ as BLAT HSPs do not have these attributes. But The
biggest difference here is that you don’t see any sequence alignments
displayed. If you look closer, PSL formats themselves do not have any
hit or query sequences, so ``Bio.SearchIO`` won’t create any sequence or
alignment objects. What happens if you try to access ``HSP.query``,
``HSP.hit``, or ``HSP.aln``? You’ll get the default values for these
attributes, which is ``None``:

.. cont-doctest

.. code:: pycon

   >>> blat_hsp.hit is None
   True
   >>> blat_hsp.query is None
   True
   >>> blat_hsp.aln is None
   True

This does not affect other attributes, though. For example, you can
still access the length of the query or hit alignment. Despite not
displaying any attributes, the PSL format still have this information so
``Bio.SearchIO`` can extract them:

.. cont-doctest

.. code:: pycon

   >>> blat_hsp.query_span  # length of query match
   61
   >>> blat_hsp.hit_span  # length of hit match
   61

Other format-specific attributes are still present as well:

.. cont-doctest

.. code:: pycon

   >>> blat_hsp.score  # PSL score
   61
   >>> blat_hsp.mismatch_num  # the mismatch column
   0

So far so good? Things get more interesting when you look at another
‘variant’ of HSP present in our BLAT results. You might recall that in
BLAT searches, sometimes we get our results separated into ‘blocks’.
These blocks are essentially alignment fragments that may have some
intervening sequence between them.

Let’s take a look at a BLAT HSP that contains multiple blocks to see how
``Bio.SearchIO`` deals with this:

.. cont-doctest

.. code:: pycon

   >>> blat_hsp2 = blat_qresult[0][1]  # first hit, second hsp
   >>> print(blat_hsp2)
         Query: mystery_seq <unknown description>
           Hit: chr19 <unknown description>
   Query range: [0:61] (1)
     Hit range: [54233104:54264463] (1)
   Quick stats: evalue ?; bitscore ?
     Fragments: ---  --------------  ----------------------  ----------------------
                  #            Span             Query range               Hit range
                ---  --------------  ----------------------  ----------------------
                  0               ?                  [0:18]     [54233104:54233122]
                  1               ?                 [18:61]     [54264420:54264463]

What’s happening here? We still some essential details covered: the IDs
and descriptions, the coordinates, and the quick statistics are similar
to what you’ve seen before. But the fragments detail is all different.
Instead of showing ‘Fragments: 1’, we now have a table with two data
rows.

This is how ``Bio.SearchIO`` deals with HSPs having multiple fragments.
As mentioned before, an HSP alignment may be separated by intervening
sequences into fragments. The intervening sequences are not part of the
query-hit match, so they should not be considered part of query nor hit
sequence. However, they do affect how we deal with sequence coordinates,
so we can’t ignore them.

Take a look at the hit coordinate of the HSP above. In the
``Hit range:`` field, we see that the coordinate is
``[54233104:54264463]``. But looking at the table rows, we see that not
the entire region spanned by this coordinate matches our query.
Specifically, the intervening region spans from ``54233122`` to
``54264420``.

Why then, is the query coordinates seem to be contiguous, you ask? This
is perfectly fine. In this case it means that the query match is
contiguous (no intervening regions), while the hit match is not.

All these attributes are accessible from the HSP directly, by the way:

.. cont-doctest

.. code:: pycon

   >>> blat_hsp2.hit_range  # hit start and end coordinates of the entire HSP
   (54233104, 54264463)
   >>> blat_hsp2.hit_range_all  # hit start and end coordinates of each fragment
   [(54233104, 54233122), (54264420, 54264463)]
   >>> blat_hsp2.hit_span  # hit span of the entire HSP
   31359
   >>> blat_hsp2.hit_span_all  # hit span of each fragment
   [18, 43]
   >>> blat_hsp2.hit_inter_ranges  # start and end coordinates of intervening regions in the hit sequence
   [(54233122, 54264420)]
   >>> blat_hsp2.hit_inter_spans  # span of intervening regions in the hit sequence
   [31298]

Most of these attributes are not readily available from the PSL file we
have, but ``Bio.SearchIO`` calculates them for you on the fly when you
parse the PSL file. All it needs are the start and end coordinates of
each fragment.

What about the ``query``, ``hit``, and ``aln`` attributes? If the HSP
has multiple fragments, you won’t be able to use these attributes as
they only fetch single ``SeqRecord`` or ``MultipleSeqAlignment``
objects. However, you can use their ``*_all`` counterparts:
``query_all``, ``hit_all``, and ``aln_all``. These properties will
return a list containing ``SeqRecord`` or ``MultipleSeqAlignment``
objects from each of the HSP fragment. There are other attributes that
behave similarly, i.e. they only work for HSPs with one fragment. Check
out the ``HSP`` documentation under :py:mod:`Bio.SearchIO` for a full
list.

Finally, to check whether you have multiple fragments or not, you can
use the ``is_fragmented`` property like so:

.. cont-doctest

.. code:: pycon

   >>> blat_hsp2.is_fragmented  # BLAT HSP with 2 fragments
   True
   >>> blat_hsp.is_fragmented  # BLAT HSP from earlier, with one fragment
   False

Before we move on, you should also know that we can use the slice
notation on ``HSP`` objects, just like ``QueryResult`` or ``Hit``
objects. When you use this notation, you’ll get an ``HSPFragment``
object in return, the last component of the object model.

.. _`sec:searchio-hspfragment`:

HSPFragment
~~~~~~~~~~~

``HSPFragment`` represents a single, contiguous match between the query
and hit sequences. You could consider it the core of the object model
and search result, since it is the presence of these fragments that
determine whether your search have results or not.

In most cases, you don’t have to deal with ``HSPFragment`` objects
directly since not that many sequence search tools fragment their HSPs.
When you do have to deal with them, what you should remember is that
``HSPFragment`` objects were written with to be as compact as possible.
In most cases, they only contain attributes directly related to
sequences: strands, reading frames, molecule types, coordinates, the
sequences themselves, and their IDs and descriptions.

These attributes are readily shown when you invoke ``print`` on an
``HSPFragment``. Here’s an example, taken from our BLAST search:

.. doctest examples

.. code:: pycon

   >>> from Bio import SearchIO
   >>> blast_qresult = SearchIO.read("my_blast.xml", "blast-xml")
   >>> blast_frag = blast_qresult[0][0][0]  # first hit, first hsp, first fragment
   >>> print(blast_frag)
         Query: 42291 mystery_seq
           Hit: gi|262205317|ref|NR_030195.1| Homo sapiens microRNA 520b (MIR520...
   Query range: [0:61] (1)
     Hit range: [0:61] (1)
     Fragments: 1 (61 columns)
        Query - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG
                |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
          Hit - CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTTTAGAGGG

At this level, the BLAT fragment looks quite similar to the BLAST
fragment, save for the query and hit sequences which are not present:

.. cont-doctest

.. code:: pycon

   >>> blat_qresult = SearchIO.read("my_blat.psl", "blat-psl")
   >>> blat_frag = blat_qresult[0][0][0]  # first hit, first hsp, first fragment
   >>> print(blat_frag)
         Query: mystery_seq <unknown description>
           Hit: chr19 <unknown description>
   Query range: [0:61] (1)
     Hit range: [54204480:54204541] (1)
     Fragments: 1 (? columns)

In all cases, these attributes are accessible using our favorite dot
notation. Some examples:

.. cont-doctest

.. code:: pycon

   >>> blast_frag.query_start  # query start coordinate
   0
   >>> blast_frag.hit_strand  # hit sequence strand
   1
   >>> blast_frag.hit  # hit sequence, as a SeqRecord object
   SeqRecord(seq=Seq('CCCTCTACAGGGAAGCGCTTTCTGTTGTCTGAAAGAAAAGAAAGTGCTTCCTTT...GGG'), id='gi|262205317|ref|NR_030195.1|', name='aligned hit sequence', description='Homo sapiens microRNA 520b (MIR520B), microRNA', dbxrefs=[])

.. _`sec:searchio-standards`:

A note about standards and conventions
--------------------------------------

Before we move on to the main functions, there is something you ought to
know about the standards ``Bio.SearchIO`` uses. If you’ve worked with
multiple sequence search tools, you might have had to deal with the many
different ways each program deals with things like sequence coordinates.
It might not have been a pleasant experience as these search tools
usually have their own standards. For example, one tools might use
one-based coordinates, while the other uses zero-based coordinates. Or,
one program might reverse the start and end coordinates if the strand is
minus, while others don’t. In short, these often creates unnecessary
mess must be dealt with.

We realize this problem ourselves and we intend to address it in
``Bio.SearchIO``. After all, one of the goals of ``Bio.SearchIO`` is to
create a common, easy to use interface to deal with various search
output files. This means creating standards that extend beyond the
object model you just saw.

Now, you might complain, "Not another standard!". Well, eventually we
have to choose one convention or the other, so this is necessary. Plus,
we’re not creating something entirely new here; just adopting a standard
we think is best for a Python programmer (it is Biopython, after all).

There are three implicit standards that you can expect when working with
``Bio.SearchIO``:

-  The first one pertains to sequence coordinates. In ``Bio.SearchIO``,
   all sequence coordinates follows Python’s coordinate style:
   zero-based and half open. For example, if in a BLAST XML output file
   the start and end coordinates of an HSP are 10 and 28, they would
   become 9 and 28 in ``Bio.SearchIO``. The start coordinate becomes 9
   because Python indices start from zero, while the end coordinate
   remains 28 as Python slices omit the last item in an interval.

-  The second is on sequence coordinate orders. In ``Bio.SearchIO``,
   start coordinates are always less than or equal to end coordinates.
   This isn’t always the case with all sequence search tools, as some of
   them have larger start coordinates when the sequence strand is minus.

-  The last one is on strand and reading frame values. For strands,
   there are only four valid choices: ``1`` (plus strand), ``-1`` (minus
   strand), ``0`` (protein sequences), and ``None`` (no strand). For
   reading frames, the valid choices are integers from ``-3`` to ``3``
   and ``None``.

Note that these standards only exist in ``Bio.SearchIO`` objects. If you
write ``Bio.SearchIO`` objects into an output format, ``Bio.SearchIO``
will use the format’s standard for the output. It does not force its
standard over to your output file.

.. _`sec:searchio-input`:

Reading search output files
---------------------------

There are two functions you can use for reading search output files into
``Bio.SearchIO`` objects: ``read`` and ``parse``. They’re essentially
similar to ``read`` and ``parse`` functions in other submodules like
``Bio.SeqIO`` or ``Bio.AlignIO``. In both cases, you need to supply the
search output file name and the file format name, both as Python
strings. You can check the documentation for a list of format names
``Bio.SearchIO`` recognizes.

``Bio.SearchIO.read`` is used for reading search output files with only
one query and returns a ``QueryResult`` object. You’ve seen ``read``
used in our previous examples. What you haven’t seen is that ``read``
may also accept additional keyword arguments, depending on the file
format.

Here are some examples. In the first one, we use ``read`` just like
previously to read a BLAST tabular output file. In the second one, we
use a keyword argument to modify so it parses the BLAST tabular variant
with comments in it:

.. doctest ../Tests/Blast

.. code:: pycon

   >>> from Bio import SearchIO
   >>> qresult = SearchIO.read("tab_2226_tblastn_003.txt", "blast-tab")
   >>> qresult
   QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)
   >>> qresult2 = SearchIO.read("tab_2226_tblastn_007.txt", "blast-tab", comments=True)
   >>> qresult2
   QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)

These keyword arguments differs among file formats. Check the format
documentation to see if it has keyword arguments that modifies its
parser’s behavior.

As for the ``Bio.SearchIO.parse``, it is used for reading search output
files with any number of queries. The function returns a generator
object that yields a ``QueryResult`` object in each iteration. Like
``Bio.SearchIO.read``, it also accepts format-specific keyword
arguments:

.. doctest ../Tests/Blast

.. code:: pycon

   >>> from Bio import SearchIO
   >>> qresults = SearchIO.parse("tab_2226_tblastn_001.txt", "blast-tab")
   >>> for qresult in qresults:
   ...     print(qresult.id)
   ...
   gi|16080617|ref|NP_391444.1|
   gi|11464971:4-101
   >>> qresults2 = SearchIO.parse("tab_2226_tblastn_005.txt", "blast-tab", comments=True)
   >>> for qresult in qresults2:
   ...     print(qresult.id)
   ...
   random_s00
   gi|16080617|ref|NP_391444.1|
   gi|11464971:4-101

.. _`sec:searchio-index`:

Dealing with large search output files with indexing
----------------------------------------------------

Sometimes, you’re handed a search output file containing hundreds or
thousands of queries that you need to parse. You can of course use
``Bio.SearchIO.parse`` for this file, but that would be grossly
inefficient if you need to access only a few of the queries. This is
because ``parse`` will parse all queries it sees before it fetches your
query of interest.

In this case, the ideal choice would be to index the file using
``Bio.SearchIO.index`` or ``Bio.SearchIO.index_db``. If the names sound
familiar, it’s because you’ve seen them before in
Section :ref:`sec:SeqIO-index`. These functions also
behave similarly to their ``Bio.SeqIO`` counterparts, with the addition
of format-specific keyword arguments.

Here are some examples. You can use ``index`` with just the filename and
format name:

.. doctest ../Tests/Blast

.. code:: pycon

   >>> from Bio import SearchIO
   >>> idx = SearchIO.index("tab_2226_tblastn_001.txt", "blast-tab")
   >>> sorted(idx.keys())
   ['gi|11464971:4-101', 'gi|16080617|ref|NP_391444.1|']
   >>> idx["gi|16080617|ref|NP_391444.1|"]
   QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)
   >>> idx.close()

Or also with the format-specific keyword argument:

.. cont-doctest

.. code:: pycon

   >>> idx = SearchIO.index("tab_2226_tblastn_005.txt", "blast-tab", comments=True)
   >>> sorted(idx.keys())
   ['gi|11464971:4-101', 'gi|16080617|ref|NP_391444.1|', 'random_s00']
   >>> idx["gi|16080617|ref|NP_391444.1|"]
   QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)
   >>> idx.close()

Or with the ``key_function`` argument, as in ``Bio.SeqIO``:

.. cont-doctest

.. code:: pycon

   >>> key_function = lambda id: id.upper()  # capitalizes the keys
   >>> idx = SearchIO.index("tab_2226_tblastn_001.txt", "blast-tab", key_function=key_function)
   >>> sorted(idx.keys())
   ['GI|11464971:4-101', 'GI|16080617|REF|NP_391444.1|']
   >>> idx["GI|16080617|REF|NP_391444.1|"]
   QueryResult(id='gi|16080617|ref|NP_391444.1|', 3 hits)
   >>> idx.close()

``Bio.SearchIO.index_db`` works like as ``index``, only it writes the
query offsets into an SQLite database file.

.. _`sec:searchio-write`:

Writing and converting search output files
------------------------------------------

It is occasionally useful to be able to manipulate search results from
an output file and write it again to a new file. ``Bio.SearchIO``
provides a ``write`` function that lets you do exactly this. It takes as
its arguments an iterable returning ``QueryResult`` objects, the output
filename to write to, the format name to write to, and optionally some
format-specific keyword arguments. It returns a four-item tuple, which
denotes the number or ``QueryResult``, ``Hit``, ``HSP``, and
``HSPFragment`` objects that were written.

.. code:: pycon

   >>> from Bio import SearchIO
   >>> qresults = SearchIO.parse("mirna.xml", "blast-xml")  # read XML file
   >>> SearchIO.write(qresults, "results.tab", "blast-tab")  # write to tabular file
   (3, 239, 277, 277)

You should note different file formats require different attributes of
the ``QueryResult``, ``Hit``, ``HSP`` and ``HSPFragment`` objects. If
these attributes are not present, writing won’t work. In other words,
you can’t always write to the output format that you want. For example,
if you read a BLAST XML file, you wouldn’t be able to write the results
to a PSL file as PSL files require attributes not calculated by BLAST
(e.g. the number of repeat matches). You can always set these attributes
manually, if you really want to write to PSL, though.

Like ``read``, ``parse``, ``index``, and ``index_db``, ``write`` also
accepts format-specific keyword arguments. Check out the documentation
for a complete list of formats ``Bio.SearchIO`` can write to and their
arguments.

Finally, ``Bio.SearchIO`` also provides a ``convert`` function, which is
simply a shortcut for ``Bio.SearchIO.parse`` and ``Bio.SearchIO.write``.
Using the convert function, our example above would be:

.. code:: pycon

   >>> from Bio import SearchIO
   >>> SearchIO.convert("mirna.xml", "blast-xml", "results.tab", "blast-tab")
   (3, 239, 277, 277)

As ``convert`` uses ``write``, it is only limited to format conversions
that have all the required attributes. Here, the BLAST XML file provides
all the default values a BLAST tabular file requires, so it works just
fine. However, other format conversions are less likely to work since
you need to manually assign the required attributes first.
