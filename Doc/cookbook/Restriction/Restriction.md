Working with restriction enzymes
================================

## Table of contents

1. [The restriction enzymes classes](#1)
    1. [Importing the enzymes](#1.1)
    2. [Naming convention](#1.2)
    3. [Searching for restriction sites](#1.3)
    4. [Retrieving the sequences produced by a digestion](#1.4)
    5. [Analysing circular sequences](#1.5)
    6. [Comparing enzymes with each others](#1.6)
    7. [Other facilities provided by the enzyme classes](#1.7)
2. [The RestrictionBatch class: a class to deal with several enzymes](#2)
    1. [Creating a RestrictionBatch](#2.1)
    2. [Restricting a RestrictionBatch to a particular supplier](#2.2)
    3. [Adding enzymes to a RestrictionBatch](#2.3)
    4. [Removing enzymes from a RestrictionBatch](#2.4)
    5. [Manipulating RestrictionBatch](#2.5)
    6. [Analysing sequences with a RestrictionBatch](#2.6)
    7. [Other RestrictionBatch methods](#2.7)
3. [AllEnzymes and CommOnly: two preconfigured RestrictionBatches](#3)
4. [The Analysis class: even simpler restriction analysis](#4)
    1. [Setting up an Analysis](#4.1)
    2. [Full restriction analysis](#4.2)
    3. [Changing the title](#4.3)
    4. [Customising the output](#4.4)
    5. [Fancier restriction analysis](#4.5)
    6. [More complex analysis](#4.6)
5. [Advanced features: the FormattedSeq class](#5)
    1. [Creating a FormattedSeq](#5.1)
    2. [Unlike Bio.Seq, FormattedSeq retains information about their shape](#5.2)
    3. [Changing the shape of a FormattedSeq](#5.3)
    4. [Using / and // operators with FormattedSeq](#5.4)
6. [More advanced features](#6)
    1. [Updating the enzymes from Rebase](#6.1)
        1. [Fetching the recent enzyme files manually from Rebase](#6.1.1)
        2. [Fetching the recent enzyme files with rebase_update.py](#6.1.2)
        3. [Compiling a new dictionary with ranacompiler.py](#6.1.3)
    2. [Subclassing the class Analysis](#6.2)
7. [Limitation and caveat](#7)
    1. [All DNA are non methylated](#7.1)
    2. [No support for star activity](#7.2)
    3. [Safe to use with degenerated DNA](#7.3)
    4. [Non standard bases in DNA are not allowed](#7.4)
    5. [Sites found at the edge of linear DNA might not be accessible in a real digestion](#7.5)
    6. [Restriction reports cutting sites not enzyme recognition sites](#7.6)
8. [Annexe: modifying dir() to use with from Bio.Restriction import *](#8)

### <a name="1"></a>1. The restriction enzymes classes
The restriction enzyme package is situated in `Bio.Restriction`. This package
will allow you to work with restriction enzymes and realise restriction
analysis on your sequence. Restriction make use of the facilities offered
by **REBASE** and contains classes for more than 800 restriction enzymes.
This chapter will lead you through a quick overview of the facilities offered
by the `Restriction` package of Biopython. The chapter is constructed as an
interactive Python session and the best way to read it is with a Python shell
open alongside you.

#### <a name="1.1"></a> 1.1  Importing the enzymes
To import the enzymes, open a Python shell and type:

``` pycon
>>> from Bio import Restriction
>>> dir()
['Restriction', '__annotations__', '__builtins__', '__doc__', '__loader__', '__name__', '__package__', '__spec__']
>>> Restriction.EcoRI
EcoRI
>>> Restriction.EcoRI.site
'GAATTC'
>>>
```

You will certainly notice that the package is quite slow to load. This is normal
as each enzyme possess its own class and there is a lot of them. This will not
affect the speed of Python after the initial import.

I don't know for you but I find it quite cumbersome to have to prefix each
operation with `Restriction.`, so here is another way to import the package.

``` pycon
>>> from Bio.Restriction import *
>>> EcoRI
EcoRI
>>> EcoRI.site
'GAATTC'
>>>
```

However, this method has one big disadvantage:
It is almost impossible to use the command `dir()` anymore as there is so much
enzymes the results is hardly readable. A workaround is provided at the end of
this tutorial. I let you decide which method you prefer. But in this tutorial I
will use the second. If you prefer the first method you will need to prefix
each call to a restriction enzyme with `Restriction.` in the remaining of the
tutorial.

#### <a name="1.2"></a>1.2  Naming convention

To access an enzyme simply enter its name. You must respect the usual naming
convention with the upper case letters and Latin numbering (in upper case as
well):

``` pycon
>>> EcoRI
EcoRI
>>> ecori
Traceback (most recent call last):
  File "<pyshell#25>", line 1, in -toplevel-
    ecori
NameError: name 'ecori' is not defined
>>> EcoR1
Traceback (most recent call last):
  File "<pyshell#26>", line 1, in -toplevel-
    EcoR1
NameError: name 'EcoR1' is not defined
>>> KpnI
KpnI
>>>
```

`ecori` or `EcoR1` are not enzymes, `EcoRI` and `KpnI` are.

#### <a name="1.3"></a>1.3  Searching for restriction sites

So what can we do with these restriction enzymes? To see that we will need a
DNA sequence. Restriction enzymes support both `Bio.Seq.MutableSeq`and
`Bio.Seq.Seq` objects. Your sequence must comply with the IUPAC alphabet.
That means using A, C, G and T or U, plus N for any base, and various
other standard codes like S for C or G, and V for A, C or G.

``` pycon
>>> from Bio.Seq import Seq
>>> my_seq = Seq("AAAAAAAAAAAAAA")
```

Searching a sequence for the presence of restriction site for your preferred
enzyme is as simple as:

``` pycon
>>> EcoRI.search(my_seq)
[]
```

The results is a list. Here the list is empty since there is obviously no EcoRI
site in *my_seq*.  Let's try to get a sequence with an EcoRI site.

``` pycon
>>> ecoseq = my_seq + Seq(EcoRI.site) + my_seq
>>> ecoseq
Seq('AAAAAAAAAAAAAAGAATTCAAAAAAAAAAAAAA')
>>> EcoRI.search(ecoseq)
[16]
```

We therefore have a site at position 16 of the sequence *ecoseq*. The position
returned by the method search is the first base of the downstream segment
produced by a restriction (i.e. the first base after the position where the
enzyme will cut). The `Restriction` package follows biological convention (the
first base of a sequence is base 1). No need to make difficult conversions
between your recorded biological data and the results produced by the enzymes
in this package.

#### <a name="1.4"></a>1.4 Retrieving the sequences produced by a digestion

`Seq` objects as all Python sequences, have different conventions and the first
base of a sequence is base 0. Therefore to get the sequences produced by an
EcoRI digestion of *ecoseq*, one should do the following:

``` pycon
>>> ecoseq[:15], ecoseq[15:]
(Seq('AAAAAAAAAAAAAAG'), Seq('AATTCAAAAAAAAAAAAAA'))
```

I hear you thinking "this is a cumbersome and error prone method to get these
sequences".  To simplify your life, `Restriction` provides another method to get
these sequences without hassle: `catalyse`. This method will return a tuple
containing all the fragments produced by a complete digestion of the sequence.
Using it is as simple as before:

``` pycon
>>> EcoRI.catalyse(ecoseq)
(Seq('AAAAAAAAAAAAAAG'), Seq('AATTCAAAAAAAAAAAAAA'))
```

BTW, you can also use spell it the American way `catalyze`:

``` pycon
>>> EcoRI.catalyze(ecoseq)
(Seq('AAAAAAAAAAAAAAG'), Seq('AATTCAAAAAAAAAAAAAA'))
```

#### <a name="1.5"></a>1.5 Analysing circular sequences

Now, if you have entered the previous command in your shell you may have noticed
that both `search` and `catalyse` can take a second argument `linear` which
defaults to `True`. Using this will allow you to simulate circular sequences
such as plasmids. Setting `linear` to `False` informs the enzyme to make the
search over a circular sequence and to search for potential sites spanning over
the boundaries of the sequence.

``` pycon
>>> EcoRI.search(ecoseq, linear=False)
[16]
>>> EcoRI.catalyse(ecoseq, linear=False)
(Seq('AATTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAG'),)
>>> ecoseq  # for memory
Seq('AAAAAAAAAAAAAAGAATTCAAAAAAAAAAAAAA')
```

OK, this is quite a difference, we only get one fragment, which correspond to
the linearised sequence. The beginning sequence has been shifted to take this
fact into account. Moreover we can see another difference:

``` pycon
>>> new_seq = Seq("TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA")
>>> EcoRI.search(new_seq)
[]
>>> EcoRI.search(new_seq, linear=False)
[33]
```

As you can see using `linear=False`, make a site appearing in the sequence
*new_seq*. This site does not exist in a linear sequence as the EcoRI site is
split into two halves at the start and the end of the sequence. In a circular
sequence however, the site is effectively present when the beginning and end of
the sequence are joined.

#### <a name="1.6"></a>1.6 Comparing enzymes with each others

`Restriction` enzymes define 4 comparative operators `==`, `!=`, `>>` and `%`.
All these operator compares two enzymes together and either return `True` or
`False`.

`==` (test identity)
: It will return `True` if the two sides of the operator are the same. *Same"
is defined as: same name, same site, same overhang (i.e. the only thing which
is equal to `EcoRI` is `EcoRI`).

`!=` (test for different site or cutting)
: It will return `True` if the two sides of the operator are different. Two
enzymes are not different if the result produced by one enzyme will always be
the same as the result produced by the other (i.e. true isoschizomers will not
being the same enzymes, are not different since they are interchangeable).

`>>` (test for neoschizomer)
: `True` if the enzymes recognise the same site, but cut it in a different way
(i.e. the enzymes are neoschizomers).

`%` (test compatibility)
: Test the compatibility of the ending produced by the enzymes (will be `True`
if the fragments produced with one of the enzyme can directly be ligated to
fragments produced by the other).

Let's use `Acc65I` and its isoschizomers as example:

``` pycon
>>> Acc65I.isoschizomers()
[Asp718I, KpnI]
>>> Acc65I.elucidate()
'G^GTAC_C'
>>> Asp718I.elucidate()
'G^GTAC_C'
>>> KpnI.elucidate()
'G_GTAC^C'
>>> # Asp718I and Acc65I are true isoschizomers,
>>> # they recognise the same site and cut it the
>>> # same way.
>>> # KpnI is a neoschizomers of the 2 others.
>>> # Here are the results of the 4 operators
>>> # for each pair of enzymes:
>>>
>>> ############# x == y  (x is y)
>>> Acc65I == Acc65I  # same enzyme => True
True
>>> Acc65I == KpnI  # all other cases => False
False
>>> Acc65I == Asp718I
False
>>> Acc65I == EcoRI
False
>>> ############ x != y  (x and y are not true isoschizomers)
>>> Acc65I != Acc65I  # same enzyme => False
False
>>> Acc65I != Asp718I  # different enzymes, but cut same manner => False
False
>>> Acc65I != KpnI  # all other cases => True
True
>>> Acc65I != EcoRI
True
>>> ###########  x >> y (x is neoschizomer of y)
>>> Acc65I >> Acc65I  # same enzyme => False
False
>>> Acc65I >> Asp718I  # same site, same cut => False
False
>>> Acc65I >> EcoRI  # different site => False
False
>>> Acc65I >> KpnI  # same site, different cut => True
True
>>> ########### x % y   (fragments produced by x and fragments produced by y
>>> #            can be directly ligated to each other)
>>> Acc65I % Asp718I
True
>>> Acc65I % Acc65I
True
>>> Acc65I % KpnI  # KpnI -> '3 overhang, Acc65I-> 5' overhang => False
False
>>>
>>> SunI.elucidate()
'C^GTAC_G'
>>> SunI == Acc65I
False
>>> SunI != Acc65I
True
>>> SunI >> Acc65I
False
>>> SunI % Acc65I  # different site, same overhang (5' GTAC) => True
True
>>> SmaI % EcoRV  # 2 Blunt enzymes, all blunt enzymes are compatible => True
True
```

#### <a name="1.7"></a>1.7 Other facilities provided by the enzyme classes

The `Restriction` class provides quite a number of others methods. We will not
go through all of them, but only have a quick look to the most useful ones.

Not all enzymes possess the same properties when it comes to the way they digest
a DNA. If you want to know more about the way a particular enzyme cut you can
use the three following methods. They are fairly straightforward to understand
and refer to the ends that the enzyme produces: blunt, 5' overhanging (also
called 3' recessed) sticky end and 3' overhanging (or 5' recessed) sticky end.

``` pycon
>>> EcoRI.is_blunt()
False
>>> EcoRI.is_5overhang()
True
>>> EcoRI.is_3overhang()
False
```

A more detailed view of the restriction site can be produced using the
 `elucidate()` method. The `^` refers to the position of the cut in the sense
strand of the sequence, `_` to the cut on the antisense or complementary strand.
`^_` means blunt.

``` pycon
>>> EcoRI.elucidate()
'G^AATT_C'
>>> KpnI.elucidate()
'G_GTAC^C'
>>> EcoRV.elucidate()
'GAT^_ATC'
```

The method `frequency()` will give you the statistical frequency of the enzyme
site.

``` pycon
>>> EcoRI.frequency()
4096
>>> XhoII.elucidate()
'R^GATC_Y'
>>> XhoII.frequency()
1024
```

To get the length of a the recognition sequence of an enzyme use the built-in
function `len()`:

``` pycon
>>> len(EcoRI)
6
>>> BstXI.elucidate()
'CCAN_NNNN^NTGG'
>>> len(BstXI)
12
>>> FokI.site
'GGATG'
>>> FokI.elucidate()  # FokI cut well outside its recognition site
'GGATGNNNNNNNNN^NNNN_N'
>>> len(FokI)  # its length is the length of the recognition site
5
```

Also interesting are the methods dealing with isoschizomers. For memory, two
enzymes are *isoschizomers* if they share a same recognition site.
A further division is made between isoschizomers (same name, recognise the same
sequence and cut the same way) and *neoschizomers* which cut at different
positions. *Equischizomer* is an arbitrary choice to design
"isoschizomers_that_are_not_neoschizomers" as this last one was a bit long.
Another set of method `one_enzyme.is_*schizomers(one_other_enzyme)`, allow to
test 2 enzymes against each other.

``` pycon
>>> Acc65I.isoschizomers()
[Asp718I, KpnI]
>>> Acc65I.neoschizomers()
[KpnI]
>>> Acc65I.equischizomers()
[Asp718I]
>>> KpnI.elucidate()
'G_GTAC^C'
>>> Acc65I.elucidate()
'G^GTAC_C'
>>> KpnI.is_neoschizomer(Acc65I)
True
>>> KpnI.is_neoschizomer(KpnI)
False
>>> KpnI.is_isoschizomer(Acc65I)
True
>>> KpnI.is_isoschizomer(KpnI)
True
>>> KpnI.is_equischizomer(Acc65I)
False
>>> KpnI.is_equischizomer(KpnI)
True
```

`suppliers()` will get you the list of all the suppliers of the enzyme.
`all_suppliers()` will give you all the suppliers in the database.

### <a name="2"></a>2. The RestrictionBatch class: a class to deal with several enzymes

If you want to make a restriction map of a sequence, using individual enzymes
can become tedious and will endures a big overhead due to the repetitive
conversion of the sequence to a `FormattedSeq` (see [Chapter 5](#5)).
`Restriction` provides a class to make easier the use of large number of enzymes
in one go: `RestrictionBatch`.
`RestrictionBatch` will help you to manipulate lots of enzymes with a single
command. Moreover all the enzymes in the restriction batch will share the same
converted sequence, reducing the overhead.

#### <a name="2.1"></a><span class="mozTocH4"></span>2.1 Creating a RestrictionBatch

You can initiate a restriction batch by passing it a list of enzymes or enzyme
names as argument.

``` pycon
>>> rb = RestrictionBatch([EcoRI])
>>> rb
RestrictionBatch(['EcoRI'])
>>> rb2 = RestrictionBatch(["EcoRI"])
>>> rb2
RestrictionBatch(['EcoRI'])
>>> rb == rb2
True
```

Adding a new enzyme to a restriction batch is easy:

``` pycon
>>> rb.add(KpnI)
>>> rb
RestrictionBatch(['EcoRI', 'KpnI'])
>>> rb += EcoRV
>>> rb
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI'])])
```

Another way to create a RestrictionBatch is by simply adding restriction enzymes
together, this is particularly useful for small batches:

``` pycon
>>> rb3 = EcoRI + KpnI + EcoRV
>>> rb3
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI'])
```

#### <a name="2.2"></a>2.2 Restricting a RestrictionBatch to a particular supplier

The Restriction package is based upon the **REBASE** database. This database
gives a list of suppliers for each enzyme. It would be a shame not to make use
of this facility. You can produce a `RestrictionBatch` containing only enzymes
from one or a few supplier(s). Here is how to do it:

``` pycon
>>> rb_supp = RestrictionBatch(
...     first=[],
...     suppliers=[
...         "C",
...         "B",
...         "E",
...         "I",
...         "K",
...         "J",
...         "M",
...         "O",
...         "N",
...         "Q",
...         "S",
...         "R",
...         "V",
...         "Y",
...         "X",
...     ],
... )
>>> # This will create a RestrictionBatch with the all enzymes which possess a s
upplier.
>>> len(rb_supp)  # May 2020
621
```

The argument `suppliers` take a list of one or several single letter codes
corresponding to the supplier(s). The codes are the same as defined in REBASE.
As it would be a pain to have to remember each supplier code, `RestrictionBatch`
provides a method which show the pair code <=> supplier:

``` pycon
>>> RestrictionBatch.show_codes()  # as of May 2016 REBASE release.
C = Minotech Biotechnology
B = Life Technologies
E = Agilent Technologies
I = SibEnzyme Ltd.
K = Takara Bio Inc.
J = Nippon Gene Co., Ltd.
M = Roche Applied Science
O = Toyobo Biochemicals
N = New England Biolabs
Q = Molecular Biology Resources - CHIMERx
S = Sigma Chemical Corporation
R = Promega Corporation
V = Vivantis Technologies
Y = SinaClon BioScience Co.
X = EURx Ltd.
>>> # You can now choose a code and built your RestrictionBatch
```

This way of producing a `RestrictionBatch` can drastically reduce the amount of
useless output from a restriction analysis, limiting the search to enzymes that
you can get hold of and limiting the risks of nervous breakdown. Nothing is more
frustrating than to get the perfect enzyme for a sub-cloning only to find it's
not commercially available.

#### <a name="2.3"></a>2.3 Adding enzymes to a RestrictionBatch

Adding an enzyme to a batch if the enzyme is already present will not raise an
exception, but will have no effects. Sometimes you want to get an enzyme from a
`RestrictionBatch` or add it to the batch if it is not present.
You will use the `get` method setting the second argument `add` to `True`.

``` pycon
>>> rb3
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI'])
>>> rb3.add(EcoRI)
>>> rb3
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI'])
>>> rb3.get(EcoRI)
EcoRI
>>> rb3.get(SmaI)

Traceback (most recent call last):
  File "<pyshell#4>", line 1, in -toplevel-
    rb3.get(SmaI)
  File "/usr/lib/Python2.3/site-packages/Bio/Restriction/Restriction.py", line 1800, in get
    raise ValueError, 'enzyme %s is not in RestrictionBatch'%e.__name__
ValueError: enzyme SmaI is not in RestrictionBatch
>>> rb3.get(SmaI, add=True)
SmaI
>>> rb3
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI', 'SmaI'])
```

#### <a name="2.4"></a>2.4 Removing enzymes from a RestrictionBatch

Removing enzymes from a batch is done using the `remove()` method. If the enzyme
is not present in the batch this will raise a `KeyError`. If the value you want
to remove is not an enzyme this will raise a `ValueError`.

``` pycon
>>> rb3.remove(EcoRI)
>>> rb3
RestrictionBatch(['EcoRV', 'KpnI', 'SmaI'])
>>> rb3.remove(EcoRI)

Traceback (most recent call last):
  File "<pyshell#14>", line 1, in -toplevel-
    rb3.remove('EcoRI')
  File "/usr/lib/Python2.3/site-packages/Bio/Restriction/Restriction.py", line 1839, in remove
    return Set.remove(self, self.format(other))
  File "/usr/lib/Python2.3/sets.py", line 534, in remove
    del self._data[element]
KeyError: EcoRI
>>> rb3 += EcoRI
>>> rb3
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI', 'SmaI'])
>>> rb3.remove("EcoRI")
>>> rb3
RestrictionBatch(['EcoRV', 'KpnI', 'SmaI'])
>>> rb3.remove("spam")

Traceback (most recent call last):
  File "<pyshell#18>", line 1, in -toplevel-
    rb3.remove('spam')
  File "/usr/lib/Python2.3/site-packages/Bio/Restriction/Restriction.py", line 1839, in remove
    return Set.remove(self, self.format(other))
  File "/usr/lib/Python2.3/site-packages/Bio/Restriction/Restriction.py", line 1871, in format
    raise ValueError, '%s is not a RestrictionType'%y.__class__
ValueError: <type 'str'> is not a RestrictionType
```

#### <a name="2.5"></a>2.5 Manipulating RestrictionBatch

You can not, however, add batches together, as they are Python `sets`. You must
use the pipe operator `|` instead. You can find the intersection between 2
batches using `&` (see the Python documentation about `sets` for more
information.

``` pycon
>>> rb3 = EcoRI + KpnI + EcoRV
>>> rb3
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI'])
>>> rb4 = SmaI + PstI
>>> rb4
RestrictionBatch(['PstI', 'SmaI'])
>>> rb3 + rb4

Traceback (most recent call last):
  File "<pyshell#23>", line 1, in -toplevel-
    rb3 + rb4
  File "/usr/lib/Python2.3/site-packages/Bio/Restriction/Restriction.py", line 1829, in __add__
    new.add(other)
  File "/usr/lib/Python2.3/site-packages/Bio/Restriction/Restriction.py", line 1848, in add
    return Set.add(self, self.format(other))
  File "/usr/lib/Python2.3/site-packages/Bio/Restriction/Restriction.py", line 1871, in format
    raise ValueError, '%s is not a RestrictionType'%y.__class__
ValueError: <class 'Bio.Restriction.Restriction.RestrictionBatch'> is not a RestrictionType
>>> rb3 | rb4
RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI', 'PstI', 'SmaI'])
>>> rb3 & rb4
RestrictionBatch([])
>>> rb4 += EcoRI
>>> rb4
RestrictionBatch(['EcoRI', 'PstI', 'SmaI'])
>>> rb3 & rb4
RestrictionBatch(['EcoRI'])
```

#### <a name="2.6"></a>2.6 Analysing sequences with a RestrictionBatch

To analyse a sequence for potential site, you can use the `search` method of the
batch, the same way you did for restriction enzymes. The results is no longer a
list however, but a dictionary. The keys of the dictionary are the names of the
enzymes and the value a list of position site. `RestrictionBatch` does not
implement a `catalyse` method, as it would not have a real meaning when used
with large batch.

``` pycon
>>> new_seq = Seq("TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA")
>>> rb.search(new_seq)
{'KpnI': [], 'EcoRV': [], 'EcoRI': []}
>>> rb.search(new_seq, linear=False)
{'KpnI': [], 'EcoRV': [], 'EcoRI': [33]}
```

#### <a name="2.7"></a>2.7 Other RestrictionBatch methods

Amongst the other methods provided by `RestrictionBatch`, `elements()` which
return a list of all the element names alphabetically sorted, is certainly the
most useful.

``` pycon
>>> rb = EcoRI + KpnI + EcoRV
>>> rb.elements()
['EcoRI', 'EcoRV', 'KpnI']
```

If you don't care about the alphabetical order use the method `as_string()`, to
get the same thing a bit faster. The list is not sorted. The order is random as
Python sets are dictionary.

``` pycon
>>> rb = EcoRI + KpnI + EcoRV
>>> rb.as_string()
['EcoRI', 'KpnI', 'EcoRV']
```

Other `RestrictionBatch` methods are generally used for particular purposes and
will not be discussed here. See the
[source](https://github.com/biopython/biopython/tree/master/Bio/Restriction)
if you are interested.

### <a name="3"></a>3. AllEnzymes and CommOnly: two preconfigured RestrictionBatches

While it is sometime practical to produce a `RestrictionBatch` of your own you
will certainly more frequently use the two batches provided with the
`Restriction` packages: `AllEnzymes` and `CommOnly`. These two batches contain
respectively all the enzymes in the database and only the enzymes which have a
commercial supplier. They are rather big, but that's what make them useful. With
these batch you can produce a full description of a sequence with a single
command. You can use these two batch as any other batch.

``` pycon
>>> len(AllEnzymes)
778
>>> len(CommOnly)
622
>>> AllEnzymes.search(new_seq)
```

There is not a lot to say about them apart the fact that they are present. They
are really normal batches, and you can use them as any other batch.

### <a name="4"></a>4. The Analysis class: even simpler restriction analysis

`RestrictionBatch` can give you a dictionary with the sites for all the enzymes
in a batch. However, it is sometime nice to get something a bit easier to read
than a Python dictionary. Complex restriction analysis are not easy with
`RestrictionBatch`. Some refinements in the way to search a sequence for
restriction sites will help. `Analysis` provides a series of command to
customise the results obtained from a pair restriction batch/sequence and some
facilities to make the output slightly more human readable.

#### <a name="4.1"></a>4.1 Setting up an Analysis

To build a restriction analysis you will need a `RestrictionBatch` and a
sequence and to tell it if the sequence is linear or circular. The first
argument `Analysis` takes is the restriction batch, the second is the sequence.
If the third argument is not provided, `Analysis` will assume the sequence is
linear.

``` pycon
>>> new_seq = Seq("TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA")
>>> rb = RestrictionBatch([EcoRI, KpnI, EcoRV])
>>> Ana = Analysis(rb, new_seq, linear=False)
>>> Ana
Analysis(RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI']),Seq('TTCAAAAAAAAAAAAAAAAAA
AAAAAAAAAAGAA'),False)
```

#### <a name="4.2"></a>4.2 Full restriction analysis

Once you have created your new `Analysis`, you can use it to get a restriction
analysis of your sequence. The way to make a full restriction analysis of the
sequence is:

``` pycon
>>> Ana.full()
{'KpnI': [], 'EcoRV': [], 'EcoRI': [33]}
```

This is much the same as the output of a `RestrictionBatch.search` method. You
will get a more easy to read output with `print_that` used without argument:

``` pycon
>>> # let's create a something a bit more complex to analyse.
>>>
>>> rb = RestrictionBatch([], ["C"])  # we will explain the meaning of the
>>> # double list argument later.
>>>
>>> multi_site = Seq.Seq(
...     "AAA"
...     + EcoRI.site
...     + "G"
...     + KpnI.site
...     + EcoRV.site
...     + "CT"
...     + SmaI.site
...     + "GT"
...     + FokI.site
...     + "GAAAGGGC"
...     + EcoRI.site
...     + "ACGT"
... )
>>> Analong = Analysis(rb, multi_site)
>>> Analong.full()
{BglI: [], BstEII: [], AsuII: [], HinfI: [], SfiI: [], PspPI: [], BsiSI: [27], S
alI: [], SlaI: [], NcoI: [], NotI: [], PstI: [], StyI: [], BseBI: [], PvuII: [],
HindIII: [], BglII: [], ApaLI: [], TaqI: [], BssAI: [], AluI: [], SstI: [], Bse
CI: [], Sau3AI: [], HpaI: [], SnaBI: [], NheI: [], BclI: [], KpnI: [16], NruI: [
], MspCI: [], BshFI: [], CspAI: [], RsaI: [14], EcoRV: [20], SphI: [], BamHI: []
, MboI: [], SgrBI: [], SspI: [], ScaI: [], XbaI: [], SseBI: [], NaeI: [], EcoRI:
[5, 47], SmaI: [28], BseAI: []}
>>>
>>> # The results are here but it is difficult to read. let's try print_that
>>>
>>> Analong.print_that()

BsiSI      :  27.
RsaI       :  14.
EcoRI      :  5, 47.
EcoRV      :  20.
KpnI       :  16.
SmaI       :  28.

   Enzymes which do not cut the sequence.

AluI      BshFI     MboI      Sau3AI    TaqI      BseBI     HinfI     PspPI
ApaLI     AsuII     BamHI     BclI      BglII     BseAI     BseCI     BssAI
CspAI     HindIII   HpaI      MspCI     NaeI      NcoI      NheI      NruI
PstI      PvuII     SalI      ScaI      SgrBI     SlaI      SnaBI     SphI
SseBI     SspI      SstI      StyI      XbaI      BstEII    NotI      BglI
SfiI
```

Much clearer, isn't ? The output is optimised for a shell 80 columns wide. If
the output seems odd, check that the width of your shell is at least 80 columns.

#### <a name="4.3"></a>4.3 Changing the title

You can provide a title to the analysis and modify the sentence 'Enzymes which
do not cut the sequence', by setting the two optional arguments of `print_that`,
`title` and `s1`. No formatting will be done on these strings so if you have to
include the newline (`\n`) as you see fit:

``` pycon
>>> Analong.print_that(None, title="sequence = multi_site\n\n")

sequence = multi_site

BsiSI      :  27.
RsaI       :  14.
EcoRI      :  5, 47.
EcoRV      :  20.
KpnI       :  16.
SmaI       :  28.

   Enzymes which do not cut the sequence.

AluI      BshFI     MboI      Sau3AI    TaqI      BseBI     HinfI     PspPI
ApaLI     AsuII     BamHI     BclI      BglII     BseAI     BseCI     BssAI
CspAI     HindIII   HpaI      MspCI     NaeI      NcoI      NheI      NruI
PstI      PvuII     SalI      ScaI      SgrBI     SlaI      SnaBI     SphI
SseBI     SspI      SstI      StyI      XbaI      BstEII    NotI      BglI
SfiI

>>> Analong.print_that(None, title="sequence = multi_site\n\n", s1="\n no site:\n\n")

sequence = multi_site

BsiSI      :  27.
RsaI       :  14.
EcoRI      :  5, 47.
EcoRV      :  20.
KpnI       :  16.
SmaI       :  28.

 no site:

AluI      BshFI     MboI      Sau3AI    TaqI      BseBI     HinfI     PspPI
ApaLI     AsuII     BamHI     BclI      BglII     BseAI     BseCI     BssAI
CspAI     HindIII   HpaI      MspCI     NaeI      NcoI      NheI      NruI
PstI      PvuII     SalI      ScaI      SgrBI     SlaI      SnaBI     SphI
SseBI     SspI      SstI      StyI      XbaI      BstEII    NotI      BglI
SfiI
```

#### <a name="4.4"></a>4.4 Customising the output

You can modify some aspects of the output interactively. There is three main
type of output, two listing types (alphabetically sorted and sorted by number
of site) and map-like type. To change the output, use the method `print_as()` of
`Analysis`. The change of output is permanent for the instance of `Analysis`
(that is until the next time you use `print_as()`). The argument of `print_as()`
are strings: `'map'`, `'number'` or `'alpha'`. As you have seen previously the
default behaviour is an alphabetical list (`'alpha'`).

``` pycon
>>> Analong.print_as("map")
>>> Analong.print_that()

    5 EcoRI
    |
    |        14 RsaI
    |        |
    |        | 16 KpnI
    |        | |
    |        | |   20 EcoRV
    |        | |   |
    |        | |   |      27 BsiSI
    |        | |   |      |
    |        | |   |      |28 SmaI
    |        | |   |      ||
    |        | |   |      ||                  47 EcoRI
    |        | |   |      ||                  |
AAAGAATTCGGGTACCGATATCCTCCCGGGGTGGATGGAAAGGGCGAATTCACGT
|||||||||||||||||||||||||||||||||||||||||||||||||||||||
TTTCTTAAGCCCATGGCTATAGGAGGGCCCCACCTACCTTTCCCGCTTAAGTGCA
1                                                    55


   Enzymes which do not cut the sequence.

AluI      BshFI     MboI      Sau3AI    TaqI      BseBI     HinfI     PspPI
ApaLI     AsuII     BamHI     BclI      BglII     BseAI     BseCI     BssAI
CspAI     HindIII   HpaI      MspCI     NaeI      NcoI      NheI      NruI
PstI      PvuII     SalI      ScaI      SgrBI     SlaI      SnaBI     SphI
SseBI     SspI      SstI      StyI      XbaI      BstEII    NotI      BglI
SfiI

>>> Analong.print_as("number")
>>> Analong.print_that()



enzymes which cut 1 times :

BsiSI      :  27.
RsaI       :  14.
EcoRV      :  20.
KpnI       :  16.
SmaI       :  28.


enzymes which cut 2 times :

EcoRI      :  5, 47.

   Enzymes which do not cut the sequence.

AluI      BshFI     MboI      Sau3AI    TaqI      BseBI     HinfI     PspPI
ApaLI     AsuII     BamHI     BclI      BglII     BseAI     BseCI     BssAI
CspAI     HindIII   HpaI      MspCI     NaeI      NcoI      NheI      NruI
PstI      PvuII     SalI      ScaI      SgrBI     SlaI      SnaBI     SphI
SseBI     SspI      SstI      StyI      XbaI      BstEII    NotI      BglI
SfiI

>>>
```

To come back to the previous behaviour:

``` pycon
>>> Analong.print_as("alpha")
>>> Analong.print_that()

BsiSI      :  27.
RsaI       :  14.
EcoRI      :  5, 47.
EcoRV      :  20.
etc ...
```

#### <a name="4.5"></a>4.5 Fancier restriction analysis

I will not go into the detail for each single method, here are all the functions
that are available. Most are perfectly self explanatory and the others are
fairly well documented (use `help('Analysis.command_name')`). The methods are:

``` python
full(self, linear=True)
blunt(self, dct=None)
overhang5(self, dct=None)
overhang3(self, dct=None)
defined(self, dct=None)
with_sites(self, dct=None)
without_site(self, dct=None)
with_N_sites(self, N, dct=None)
with_number_list(self, list, dct=None)
with_name(self, names, dct=None)
with_site_size(self, site_size, dct=None)
only_between(self, start, end, dct=None)
between(self, start, end, dct=None)
show_only_between(self, start, end, dct=None)
only_outside(self, start, end, dct=None)
outside(self, start, end, dct=None)
do_not_cut(self, start, end, dct=None)
```

Using these methods is simple:

``` pycon
>>> new_seq = Seq("TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA")
>>> rb = RestrictionBatch([EcoRI, KpnI, EcoRV])
>>> Ana = Analysis(rb, new_seq, linear=False)
>>> Ana
Analysis(RestrictionBatch(['EcoRI', 'EcoRV', 'KpnI']),Seq('TTCAAAAAAAAAAAAAAAAAA
AAAAAAAAAAGAA'),False)
>>> Ana.blunt()  # output only the result for enzymes which cut blunt
{'EcoRV': []}
>>> Ana.full()  # all the enzymes in the RestrictionBatch
{'KpnI': [], 'EcoRV': [], 'EcoRI': [33]}
>>> Ana.with_sites()  # output only the result for enzymes which have a site
{'EcoRI': [33]}
>>> Ana.without_site()  # output only the enzymes which have no site
{'KpnI': [], 'EcoRV': []}
>>> Ana.only_between(1, 20)  # the enzymes which cut between position 1 and 20
{}
>>> Ana.only_between(20, 34)  # etc...
{'EcoRI': [33]}
>>> Ana.only_outside(20, 34)
{}
>>> Ana.with_name([EcoRI])
{'EcoRI': [33]}
>>>
```

To get a nice output, you still use `print_that` but this time with the command
you want executed as argument.

``` pycon
>>> Ana.print_that(Ana.blunt())

   Enzymes which do not cut the sequence.

EcoRV

>>> pt = Ana.print_that
>>> pt(Ana.with_sites())

EcoRI      :  33.

>>> pt(Ana.without_site())

   Enzymes which do not cut the sequence.

EcoRV     KpnI

>>> # etc ...
```

#### <a name="4.6"></a>4.6 More complex analysis

All of these methods (except `full()` which, well ... do a full restriction
analysis) can be supplied with an additional dictionary. If no dictionary is
supplied a full restriction analysis is used as starting point. Otherwise the
dictionary provided by the argument `dct` is used. The dictionary must be
formatted as the result of `RestrictionBatch.search`. Therefore of the form
`{'enzyme_name': [position1, position2],...}`, where *position1* and
*position2* are integers. All methods list previously output such dictionaries
and can be used as starting point.

Using this method you can build really complex query by chaining several method
one after the other. For example if you want all the enzymes which are 5'
overhang and cut the sequence only once, you have two ways to go:

The hard way consist to build a restriction batch containing only 5' overhang
enzymes and use this batch to create a new `Analysis` instance and then use the
method `with_N_sites()` as follow:

``` pycon
>>> rbov5 = RestrictionBatch([x for x in rb if x.is_5overhang()])
>>> Anaov5 = Analysis(rbov5, new_seq, linear=False)
>>> Anaov5.with_N_sites(1)
{'EcoRI' : [33]}
```

The easy solution is to chain several `Analysis` methods. This is possible since
each method return a dictionary as results and is able to take a dictionary as
input:

``` pycon
>>> Ana.with_N_sites(1, Ana.overhang5())
{'EcoRI': [33]}
```

The dictionary is always the last argument whatever the command you use.

The way to prefer certainly depends of the conditions you will use your
`Analysis` instance. If you are likely to frequently reuse the same batch with
different sequences, using a dedicated `RestrictionBatch` might be faster as the
batch is likely to be smaller. Chaining methods is generally quicker when
working with an interactive shell. In a script, the extended syntax may be
easier to understand in a few months.

### <a name="5"></a>5. Advanced features: the FormattedSeq class

Restriction enzymes require a much more strict formatting of the DNA sequences
than `Bio.Seq` object provides. For example, the restriction enzymes expect to
find an ungapped (no space) upper-case sequence, while `Bio.Seq` object allow
sequences to be in lower-case separated by spaces. Therefore when a restriction
enzyme analyse a `Bio.Seq` object (be it a `Seq` or a `MutableSeq`), the object
undergoes a conversion. The class `FormattedSeq` ensure the smooth conversion
from a `Bio.Seq` object to something which can be safely be used by the enzyme.

While this conversion is done automatically by the enzymes if you provide them
with a `Seq` or a `MutableSeq`, there is time where it will be more efficient to
realise the conversion before hand. Each time a `Seq` object is passed to an
enzyme for analysis you pay a overhead due to the conversion. When analysing the
same sequence over and over, it will be faster to convert the sequence, store
the conversion and then use only the converted sequence.

#### <a name="5.1"></a>5.1 Creating a FormattedSeq

Creating a `FormattedSeq` from a `Bio.Seq` object is simple. The first argument
of `FormattedSeq` is the sequence you wish to convert. You can specify a shape
with the second argument `linear`, if you don't the `FormattedSeq` will be
linear:

``` pycon
>>> from Bio.Restriction import *
>>> from Bio.Seq import Seq
>>> seq = Seq("TTCAAAAAAAAAAGAATTCAAAAGAA")
>>> linear_fseq = FormattedSeq(seq, linear=True)
>>> default_fseq = FormattedSeq(seq)
>>> circular_fseq = FormattedSeq(seq, linear=False)
>>> linear_fseq
FormattedSeq(Seq('TTCAAAAAAAAAAGAATTCAAAAGAA'), linear=True)
>>> linear_fseq.is_linear()
True
>>> default_fseq.is_linear()
True
>>> circular_fseq.is_linear()
False
>>> circular_fseq
FormattedSeq(Seq('TTCAAAAAAAAAAGAATTCAAAAGAA'), linear=False)
```

#### <a name="5.2"></a>5.2 Unlike Bio.Seq, FormattedSeq retains information about their shape

`FormattedSeq` retains information about the shape of the sequence. Therefore
unlike with `Seq` and `MutableSeq` you don't need to specify the shape of the
sequence when using `search()` or `catalyse()`:

``` pycon
>>> EcoRI.search(linear_fseq)
[15]
>>> EcoRI.search(circular_fseq)  # no need to specify the shape
[15, 25]
```

In fact, the shape of a FormattedSeq is not altered by the second argument of
the commands `search()` and `catalyse()`:

``` pycon
>>> # In fact the shape is blocked.
>>> # The 3 following commands give the same results
>>> # which correspond to a circular sequence
>>> EcoRI.search(circular_fseq)
[15, 25]
>>> EcoRI.search(circular_fseq, linear=True)
[15, 25]
>>> EcoRI.search(circular_fseq, linear=False)
[15, 25]
```

#### <a name="5.3"></a>5.3 Changing the shape of a FormattedSeq

You can however change the shape of the `FormattedSeq`. The command to use are:

``` python
FormattedSeq.to_circular()  # new FormattedSeq, shape will be circular.
FormattedSeq.to_linear()  # new FormattedSeq, shape will be linear
FormattedSeq.circularise()  # change the shape of FormattedShape to circular
FormattedSeq.linearise()  # change the shape of FormattedShape to linear
```

``` pycon
>>> circular_fseq
FormatedSeq(Seq('TTCAAAAAAAAAAGAATTCAAAAGAA'), linear=False)
>>> circular_fseq.is_linear()
False
>>> circular_fseq == linear_fseq
False
>>> newseq = circular_fseq.to_linear()
>>> circular_fseq
FormatedSeq(Seq('TTCAAAAAAAAAAGAATTCAAAAGAA'), linear=False)
>>> newseq
FormatedSeq(Seq('TTCAAAAAAAAAAGAATTCAAAAGAA'), linear=True)
>>> circular_fseq.linearise()
>>> circular_fseq
FormatedSeq(Seq('TTCAAAAAAAAAAGAATTCAAAAGAA'), linear=True)
>>> circular_fseq.is_linear()
True
>>> circular_fseq == linear_fseq
True
>>> EcoRI.search(circular_fseq)  # which is now linear
[15]
```

#### <a name="5.4"></a>5.4 Using / and // operators with FormattedSeq

Not having to specify the shape of the sequence to analyse gives you the
opportunity to use the shorthand '/' and '//' with restriction enzymes:

``` pycon
>>> EcoRI / linear_fseq  # <=> EcoRI.search(linear_fseq)
[15]
>>> linear_fseq / EcoRI  # <=> EcoRI.search(linear_fseq)
[15]
>>> EcoRI // linear_fseq  # <=> linear_fseq//EcoRI <=> EcoRI.catalyse(linear_fseq)
(Seq('TTCAAAAAAAAAAG'), Seq('AATTCAAAAGAA'))
```

Another way to avoid the overhead due to a repetitive conversion from a `Seq`
object to a `FormattedSeq` is to use a [`RestrictionBatch`](#2).

To conclude, the performance gain achieved when using a `FormattedSeq` instead
of a `Seq` is not huge. The analysis of a 10 kb sequence by all the enzymes in
`AllEnzymes` (`for x in AllEnzymes: x.search(seq)`, 867 enzymes) is 7 % faster
when using a `FormattedSeq` than a `Seq`. Using a `RestrictionBatch`
(`AllEnzymes.search(seq)`) is about as fast as using a `FormattedSeq` the first
time the search is run. This however is dramatically reduced in subsequent runs
with the same sequence (`RestrictionBatch` keeps in memory the result of their
last run while the sequence is not changed).

### <a name="6"></a>6.  More advanced features

This chapter addresses some more advanced features of the packages, most users
can safely ignore it.

#### <a name="6.1"></a>6.1 Updating the enzymes from REBASE

Most people will certainly not need to update the enzymes. The restriction
enzyme package will be updated in with each new release of Biopython. But if you
wish to get an update in between Biopython-releases here is how to do it.

First, you have to download the two scripts `rebase_update.py` and
`ranacompile.py`:
Go to https://github.com/biopython/biopython/tree/master/Scripts/Restriction,
click on the respective file and press the '**Raw**' button on the top right of
the code window. Then, with right-click, save the file. Both scripts must be in
the same directory.

##### <a name="6.1.1"></a>6.1.1 Fetching the recent enzyme files manually from REBASE

Each month, [REBASE](http://rebase.neb.com/rebase/rebase.html) release a new
compilation of data about restriction enzymes. While the enzymes do not change
so frequently, you may wish to update the restriction enzymes classes. The first
thing to do is to get the last rebase file. You can find the release of REBASE
at http://rebase.neb.com/rebase/rebase.files.html. The file you are interested
in are in the EMBOSS format. You can download the files directly from the REBASE
ftp server using your browser. The file are situated at
ftp://ftp.neb.com/pub/rebase.
You will have to download 3 files: `emboss_e.###`, `emboss_r.###`, and
`emboss_s.###`.
The `###` is a 3 digit number corresponding to the year and month of the
release. The first digit is the year, the two last are the month: so July 2015
will be: 507; October 2016: 610, etc... Download the three file corresponding to
the current month and place them in the same folder as your `rebase_update.py`
and `ranacompiler.py` scripts.

##### <a name="6.1.2"></a>6.1.2 Fetching the recent enzyme files with rebase_update.py

Another way to do the same thing is to use the `rebase_update.py` script. It
will connect directly to the rebase ftp server and download the last batch of
emboss files. From a DOS or Unix shell do the following:

``` bash
$ cd path_to_the_update_script
$ rebase_update.py -p http://www.somewhere.com:8000

Please wait, trying to connect to Rebase

copying ftp://ftp.neb.com/pub/rebase/emboss_e.407
to /cvsroot/bioPython/Bio/Restriction/Scripts/emboss_e.407
copying ftp://ftp.neb.com/pub/rebase/emboss_s.407
to /cvsroot/bioPython/Bio/Restriction/Scripts/emboss_s.407
copying ftp://ftp.neb.com/pub/rebase/emboss_r.407
to /cvsroot/bioPython/Bio/Restriction/Scripts/emboss_r.407
```

Some explanation are needed: `-p` is the switch to indicate to the script you
are using a proxy. If you use a ftp proxy enter its address and the connection
port after the '`:`'.


##### <a name="6.1.3"></a>6.1.3 Compiling a new dictionary with ranacompiler.py

Once you have got the recent emboss files you can compile a new module
containing the data necessary to create restriction enzyme.

Note: if the emboss files are not present in the current directory or if they
are not up to date, `ranacompiler.py` will invoke the script
[`rebase_update.py`](#6.1.2), which needs to be installed in the same folder.
You will need to use the same options as before (ie `-m` and `-p`). See the
previous paragraph on [`rebase_update.py`](#6.1.2) for more details.

For simplicity let's assume we have put the emboss files in the same folder as
the files which contains the script `ranacompiler.py`. You may have the change
the mode of the file to make it executable:

``` bash
$ cd path_to_the_ranacompiler_script
$ chmod '+x' ranacompiler.py
```

Now execute the script:

``` bash
$ Python ranacompiler.py  # or ./ranacompiler.py
```

You get normally the following message:

``` bash
$ ./ranacompiler.py

 Using the files : emboss_e.407, emboss_r.407, emboss_s.407

WARNING : HaeIV cut twice with different overhang length each time.
        Unable to deal with this behaviour.
        This enzyme will not be included in the database. Sorry.
        Checking :
Anyway, HaeIV is not commercially available.

WARNING : HpyUM037X has two different sites.


The new database contains 867 enzymes.

Writing the dictionary containing the new Restriction classes...
OK.

Writing the dictionary containing the suppliers data...
OK.

Writing the dictionary containing the Restriction types....
OK.

 ******************************************************************************

                Compilation of the new dictionary : OK.
                Installation : No.

 You will find the newly created 'Restriction_Dictionary.py' file
 in the  :

        /path/where/you/run/ranacompiler.py

 Make a copy of 'Restriction_Dictionary.py' and place it with
 the other Restriction libraries.

 note :
 This folder should be :

        path_to_python/site-packages/Bio/Restriction

 ******************************************************************************
```

The first line indicate which emboss files have been used for the present
compilation. You can safely ignore the warnings as long as the
`compilation of the new dictionary : OK.` is present in the last part of the
output. They are here for debugging purpose. The number of enzymes in the new
module is indicated as well as a list of the dictionary which have been
compiled. The last part indicate that the module has been successfully created
but not installed. To finish the update you must copy the file
`Restriction_Dictionary.py` into the folder
`/your_python_path/site-packages/Bio/Restriction/` as indicated by the script.
Looking into the present folder, you will see to new files: the newly created
dictionary `Restriction_Dictionary.py` and `Restriction_Dictionary.old`. This
last file containing the old dictionary to which you can revert in case anything
the new file is corrupted (this should not happen since the script is happy
enough the new dictionary is good, but if there is a problem it is always nice
to know you can revert to the previous setting without having to reinstall the
whole thing.

If you wish, the script may install the folder for you as well, but you will
have to run it as root if your normal user has no write access to your Python
installation (and it shouldn't). Use the command `ranacompiler.py -i` or
`ranacompiler.py --install` for this.

If anything goes wrong (you have no write access to the destination folder for
example) the script will let you know it did not perform the installation. It
will however still save the new module in the current directory.

As you can see the script is not very bright and will redo the compilation each
time it is invoked, no matter if a previous version of the module is already
present.

#### <a name="6.2"></a>6.2 Subclassing the class Analysis

As seen previously, you can modify some aspects of the `Analysis` output
interactively. However if you want to write your own `Analysis` class, you may
wish to provide others output facilities than is given in this package.
Depending on what you want to do you may get away with simply changing the
`make_format` method of your derived class or you will need to provide new
methods. Rather than get into a long explanation, here is the implementation of
a rather useless `Analysis` class:

``` pycon
>>> class UselessAnalysis(Analysis):
...     def __init__(self, rb=RestrictionBatch(), seq=Seq(""), lin=True):
...         """UselessAnalysis -> A class that waste your time"""
...         #
...         #    Unless you want to do something more fancy all
...         #    you need to do here is instantiate Analysis.
...         #    Don't forget the self in __init__
...         #
...         Analysis.__init__(self, rb, seq, lin)
...     def make_format(self, cut=[], t="", nc=[], s=""):
...         """not funny"""
...         #
...         #    Generally, you don't need to do anything else here
...         #    This will tell to your new class to default to the
...         #    _make_joke format.
...         #
...         return self._make_joke(cut, t, nc, s)
...     def print_as(self, what="joke"):
...         """Somebody might want to change the behaviour of this class."""
...         #
...         #    add your new option to print_as
...         #
...         if what == "joke":
...             self.make_format = self._make_joke
...             return
...         else:
...             #
...             #   The other options will be treated as before
...             #
...             return Analysis.print_as(self, what)
...     def _make_joke(self, cut=[], title="", nc=[], s1=""):
...         """UA._make_joke(cut, t, nc, s) -> new analysis output"""
...         #
...         #    starting your new method with '_make_'
...         #    will give a hint to what it is suppose to do
...         #
...         #    We will not process the non-cutting enzymes
...         #    Their names are in nc
...         #    s1 is the string printed before them
...         #
...         if not title:
...             title = "\nYou have guessed right the following enzymes:\n\n"
...         for name, sites in cut:
...             #
...             #    cut contains:
...             #    - the name of the enzymes which cut the sequence (name)
...             #    - a list of the site positions (sites)
...             guess = raw_input("next enzyme is %s, Guess how many sites ?\n>>> " % name)
...             try:
...                 guess = int(guess)
...             except:
...                 guess = None
...             if guess == len(sites):
...                 print("You did guess right. Good. Next.")
...             result = "%i site" % guess
...             if guess > 1:
...                 result += "s"
...             #
...             #    now we format the line. See the PrintFormat module
...             #    for some examples
...             #   PrintFormat.__section_list and _make_map are good start.
...             #
...             title = "".join(
...                 (title, str(name).ljust(self.NameWidth), " :  ", result, ".\n")
...             )
...         print("\nNo more enzyme.")
...
..          return title
...         #
...         #    I you want to print the non cutting enzymes use
...         #    the following return instead of the previous one:
...         #
...         # return  title + t + self._make_nocut_only(nc,s1)
...
>>> # You initiate and use it as before
>>> rb = RestrictionBatch([], ["A"])
>>> multi_site = Seq(
...     "AAA"
...     + EcoRI.site
...     + "G"
...     + KpnI.site
...     + EcoRV.site
...     + "CT"
...     + SmaI.site
...     + "GT"
...     + FokI.site
...     + "GAAAGGGC"
...     + EcoRI.site
...     + "ACGT"
... )
>>>
>>> b = UselessAnalysis(rb, multi_site)
>>> b.print_that()  # Well, I let you discover if you haven't already guessed
```

Using this example, as a template you should now be able to subclass `Analysis`
as you wish. You will found more implementation details in the module
`Bio.Restriction.PrintFormat` which contains the class providing all the
`_make_*` methods.

### <a name="7"></a>7. Limitation and caveat

Particularly, the class `Analysis` is a quick and dirty implementation based on
the facilities furnished by the package. Please check your results and report
any fault.

On a more general basis, `Restriction` have some other limitations:

#### <a name="7.1"></a>7.1 All DNA are non methylated

No facility to work with methylated DNA has been implemented yet. As far as the
enzyme classes are concerned all DNA is non methylated DNA. Implementation of
methylation sensibility will possibly occur in the future. But for now, if your
sequence is methylated, you will have to check if the site is methylated using
other means.

#### <a name="7.2"></a>7.2 No support for star activity

As before no support has been yet implemented to find site mis-recognised by
enzymes under high salt concentration conditions, the so-called star activity.
This will be implemented as soon as I can get a good source of information for
that.

#### <a name="7.3"></a>7.3 Safe to use with degenerated DNA

It is safe to use degenerated DNA as input for the query. You will not be
flooded with meaningless results. But this come at a price: GAA***N***TC will
not be recognised as a potential EcoRI site for example, in fact it will not be
recognised at all. Degenerated sequences will not be analysed. If your sequence
is not fully sequenced, you will certainly miss restriction sites:

``` pycon
>>> a = Seq("nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnGAATTCrrrrrrrrrrr")
>>> EcoRI.search(a)
[36]
>>> b = Seq("nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnGAAnTCrrrrrrrrrrr")
>>> EcoRI.search(b)
[]
```

#### <a name="7.4"></a>7.4 Non standard bases in DNA are not allowed

While you can use degenerated DNA, using non standard base alphabet will make
the enzymes choke, even if `Bio.Seq.Seq` accepts them. However, space-like
characters (' ', '\n', '\t', ...) and digit will be removed but will not stop
the enzyme analysing the sequence. You can use them but the fragments produced
by `catalyse` will have lost any formatting. `catalyse` tries to keep the
original case of the sequence (i.e lower case sequences will generate lower case
fragments, upper case sequences upper case fragments), but mixed case will
return upper case fragments:

``` pycon
>>> c = Seq("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxGAANTCrrrrrrrrrrr")
>>> EcoRI.search(c)

Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "/usr/lib/python3.6/site-packages/Bio/Restriction/Restriction.py", line 553, in search
    cls.dna = FormattedSeq(dna, linear)
  File "/usr/lib/python3.6/site-packages/Bio/Restriction/Restriction.py", line 171, in __init__
    self.data = _check_bases(stringy)
  File "/usr/lib/python3.6/site-packages/Bio/Restriction/Restriction.py", line 122, in _check_bases
    raise TypeError("Invalid character found in %s" % repr(seq_string))
TypeError: Invalid character found in 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXGAANTCRRRRRRRRRRR'
>>> d = Seq(
...     "1 nnnnn nnnnn nnnnn nnnnn nnnnn \n" "26 nnnnn nnnnG AATTC rrrrr rrrrr \n" "51 r"
... )
>>> d
Seq('1 nnnnn nnnnn nnnnn nnnnn nnnnn \n26 nnnnn nnnnG AATTC rrrrr rrrrr \n51 r')
>>> EcoRI.search(d)
[36]
>>> EcoRI.catalyse(d)
(Seq('AATTCRRRRRRRRRRR'), Seq('NNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNG'))
>>> e = Seq("nnnnGAATTCrr")
>>> f = Seq("NNNNGAATTCRR")
>>> g = Seq("nnnngaattcrr")
>>> EcoRI.catalyse(e)
(Seq('NNNNG'), Seq('AATTCRR'))
>>> EcoRI.catalyse(f)
(Seq('NNNNG'), Seq('AATTCRR'))
>>> EcoRI.catalyse(g)
(Seq('nnnng'), Seq('aattcrr'))
```

Not allowing other letters than IUPAC might seems drastic but this is really to
limit errors. It is not totally fool proof but it does help.

#### <a name="7.5"></a>7.5 Sites found at the edge of linear DNA might not be accessible in a real digestion

While sites clearly outsides a sequence will not be reported, nothing has been
done to try to determine if a restriction site at the end of a linear sequence
is valid:

``` pycon
>>> d = Seq("GAATTCAAAAAAAAAAAAAAAAAAAAAAAAAAGGATG")
>>> FokI.site  # site present
'GGATG'
>>> FokI.elucidate()  # but cut outside the sequence
'GGATGNNNNNNNNN^NNNN_N'
>>> FokI.search(d)  # therefore no site found
[]
>>> EcoRI.search(d)
[2]
```

`EcoRI` finds a site at position 2 even if it is highly unlikely that EcoRI
accepts to cut this site in a tube. It is generally considered that at about 5
nucleotides must separate the site from the edge of the sequence to be
reasonably sure the enzyme will work correctly. This "security margin" is
variable from one enzyme to the other. In doubt consult the documentation for
the enzyme.

#### <a name="7.6"></a>7.6 Restriction reports cutting sites not enzyme recognition sites

Some enzymes will cut twice each time they encounter a restriction site. The
enzymes in this package report both cut not the site. Other software may only
reports restriction sites. Therefore the output given for some enzymes might
seems to be the double when compared with the results of these software. It is
not a bug.

``` pycon
>>> AloI.cut_twice()
True
>>> AloI.fst5  # first cut
-7
>>> AloI.scd5  # second cut
25
>>> AloI.site
'GAACNNNNNNTCC'
>>> b = Seq("AAAAAAAAAAA" + AloI.site + "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
>>> b
Seq('AAAAAAAAAAAGAACNNNNNNTCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
>>> AloI.search(b)  # one site, two cuts -> two positions
[5, 37]
```

### <a name="8"></a>8.  Annex: modifying dir() to use with from Bio.Restriction import *

Having all the enzymes imported directly in the shell is useful when working in
an interactive shell (even if it is not recommended by the purists). Here is a
little hack to get some sanity back when using dir() in those conditions:

``` pycon
>>> # we will change the builtin dir() function to get ride of the enzyme names.
>>> import sys
>>> def dir(object=None):
...     """dir([object]) -> list of string.
...     over-ride the built-in function to get some clarity."""
...     if object:
...         # we only want to modify dir(),
...         # so here we return the result of the builtin function.
...         return __builtins__.dir(object)
...     else:
...         # now the part we want to modify.
...         # All the enzymes are in a RestrictionBatch (we will talk about
...         # that later, for the moment simply believe me).
...         # So if we remove from the results of dir() everything which is
...         # in AllEnzymes we will get a much shorter list when we do dir()
...         #
...         # the current level is __main__ ie dir() is equivalent to
...         # ask what's in __main__ at the moment.
...         # we can't access __main__ directly.
...         # so we will use sys.modules['__main__'] to reach it.
...         # the following list comprehension remove from the result of
...         # dir() everything which is also present in AllEnzymes.
...         #
...         return [
...             x for x in __builtins__.dir(sys.modules["__main__"]) if not x in AllEnzymes
...         ]
...
>>> # now let's see if it works.
>>> dir()
['AllEnzymes', 'Analysis', 'CommOnly', 'NonComm', 'PrintFormat', 'RanaConfig',
 'Restriction', 'RestrictionBatch', 'Restriction_Dictionary', '__builtins__',
 '__doc__', '__name__', 'dir', 'sys']
>>> # ok that's much better.
>>> # The enzymes are still there
>>> EcoRI.site
'GAATTC'
```
