.. _`chapter:phenotype`:

Bio.phenotype: analyze phenotypic data
======================================

This chapter gives an overview of the functionalities of the
``Bio.phenotype`` package included in Biopython. The scope of this
package is the analysis of phenotypic data, which means parsing and
analyzing growth measurements of cell cultures. In its current state the
package is focused on the analysis of high-throughput phenotypic
experiments produced by the `Phenotype Microarray
technology <https://en.wikipedia.org/wiki/Phenotype_microarray>`__, but
future developments may include other platforms and formats.

.. _`sec:phenotypemicroarrays`:

Phenotype Microarrays
---------------------

The `Phenotype
Microarray <https://en.wikipedia.org/wiki/Phenotype_microarray>`__ is a
technology that measures the metabolism of bacterial and eukaryotic
cells on roughly 2000 chemicals, divided in twenty 96-well plates. The
technology measures the reduction of a tetrazolium dye by NADH, whose
production by the cell is used as a proxy for cell metabolism; color
development due to the reduction of this dye is typically measured once
every 15 minutes. When cells are grown in a media that sustains cell
metabolism, the recorded phenotypic data resembles a sigmoid growth
curve, from which a series of growth parameters can be retrieved.

Parsing Phenotype Microarray data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``Bio.phenotype`` package can parse two different formats of
Phenotype Microarray data: the
`CSV <https://en.wikipedia.org/wiki/Comma-separated_values>`__ (comma
separated values) files produced by the machine’s proprietary software
and `JSON <https://en.wikipedia.org/wiki/JSON>`__ files produced by
analysis software, like
`opm <https://www.dsmz.de/research/microorganisms/projects/analysis-of-omnilog-phenotype-microarray-data.html>`__
or `DuctApe <https://combogenomics.github.io/DuctApe/>`__. The parser
will return one or a generator of PlateRecord objects, depending on
whether the read or parse method is being used. You can test the parse
function by using the
`Plates.csv <https://github.com/biopython/biopython/blob/master/Doc/examples/Plates.csv>`__
file provided with the Biopython source code.

.. doctest examples lib:numpy

.. code:: pycon

   >>> from Bio import phenotype
   >>> for record in phenotype.parse("Plates.csv", "pm-csv"):
   ...     print("%s %i" % (record.id, len(record)))
   ...
   PM01 96
   PM01 96
   PM09 96
   PM09 96

The parser returns a series of PlateRecord objects, each one containing
a series of WellRecord objects (holding each well’s experimental data)
arranged in 8 rows and 12 columns; each row is indicated by a uppercase
character from A to H, while columns are indicated by a two digit
number, from 01 to 12. There are several ways to access WellRecord
objects from a PlateRecord objects:

Well identifier
   If you know the well identifier (row + column identifiers) you can
   access the desired well directly.

   .. cont-doctest

   .. code:: pycon

      >>> record["A02"]
      WellRecord('(0.0, 12.0), (0.25, 18.0), (0.5, 27.0), (0.75, 35.0), (1.0, 37.0), ..., (71.75, 143.0)')

Well plate coordinates
   The same well can be retrieved by using the row and columns numbers
   (0-based index).

   .. doctest examples lib:numpy

   .. code:: pycon

      >>> from Bio import phenotype
      >>> record = list(phenotype.parse("Plates.csv", "pm-csv"))[-1]
      >>> print(record[0, 1].id)
      A02

Row or column coordinates
   A series of WellRecord objects contiguous to each other in the plate
   can be retrieved in bulk by using the python list slicing syntax on
   PlateRecord objects; rows and columns are numbered with a 0-based
   index.

   .. cont-doctest

   .. code:: pycon

      >>> print(record[0])
      Plate ID: PM09
      Well: 12
      Rows: 1
      Columns: 12
      PlateRecord('WellRecord['A01'], WellRecord['A02'], WellRecord['A03'], ..., WellRecord['A12']')
      >>> print(record[:, 0])
      Plate ID: PM09
      Well: 8
      Rows: 8
      Columns: 1
      PlateRecord('WellRecord['A01'], WellRecord['B01'], WellRecord['C01'], ..., WellRecord['H01']')
      >>> print(record[:3, :3])
      Plate ID: PM09
      Well: 9
      Rows: 3
      Columns: 3
      PlateRecord('WellRecord['A01'], WellRecord['A02'], WellRecord['A03'], ..., WellRecord['C03']')

Manipulating Phenotype Microarray data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Accessing raw data
^^^^^^^^^^^^^^^^^^

The raw data extracted from the PM files is comprised of a series of
tuples for each well, containing the time (in hours) and the
colorimetric measure (in arbitrary units). Usually the instrument
collects data every fifteen minutes, but that can vary between
experiments. The raw data can be accessed by iterating on a WellRecord
object; in the example below only the first ten time points are shown.

.. doctest examples lib:numpy

.. code:: pycon

   >>> from Bio import phenotype
   >>> record = list(phenotype.parse("Plates.csv", "pm-csv"))[-1]
   >>> well = record["A02"]

.. code:: pycon

   >>> for time, signal in well:
   ...     print(time, signal)
   ...
   (0.0, 12.0)
   (0.25, 18.0)
   (0.5, 27.0)
   (0.75, 35.0)
   (1.0, 37.0)
   (1.25, 41.0)
   (1.5, 44.0)
   (1.75, 44.0)
   (2.0, 44.0)
   (2.25, 44.0)
   [...]

This method, while providing a way to access the raw data, doesn’t allow
a direct comparison between different WellRecord objects, which may have
measurements at different time points.

Accessing interpolated data
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To make it easier to compare different experiments and in general to
allow a more intuitive handling of the phenotypic data, the module
allows to define a custom slicing of the time points that are present in
the WellRecord object. Colorimetric data for time points that have not
been directly measured are derived through a linear interpolation of the
available data, otherwise a NaN is returned. This method only works in
the time interval where actual data is available. Time intervals can be
defined with the same syntax as list indexing; the default time interval
is therefore one hour.

.. cont-doctest

.. code:: pycon

   >>> well[:10]
   [12.0, 37.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0]
   >>> well[9.55]
   44.0

Different time intervals can be used, for instance five minutes:

.. cont-doctest

.. code:: pycon

   >>> for value in well[63 : 64 : 5 / 60]:
   ...     print(f"{value:0.2f}")
   ...
   110.00
   111.00
   112.00
   113.00
   113.33
   113.67
   114.00
   114.33
   114.67
   115.00
   115.00
   115.00
   >>> for value in well[63.33:73.33]:
   ...     print(f"{value:0.2f}")
   ...
   113.32
   117.00
   120.32
   128.00
   129.64
   132.96
   136.96
   140.00
   142.00
   nan

Control well subtraction
^^^^^^^^^^^^^^^^^^^^^^^^

Many Phenotype Microarray plates contain a control well (usually A01),
that is a well where the media shouldn’t support any growth; the low
signal produced by this well can be subtracted from the other wells. The
PlateRecord objects have a dedicated function for that, which returns
another PlateRecord object with the corrected data.

.. cont-doctest

.. code:: pycon

   >>> corrected = record.subtract_control(control="A01")
   >>> record["A01"][63]
   336.0
   >>> corrected["A01"][63]
   0.0

Parameters extraction
^^^^^^^^^^^^^^^^^^^^^

Those wells where metabolic activity is observed show a sigmoid behavior
for the colorimetric data. To allow an easier way to compare different
experiments a sigmoid curve can be fitted onto the data, so that a
series of summary parameters can be extracted and used for comparisons.
The parameters that can be extracted from the curve are:

-  Minimum (**min**) and maximum (**max**) signal;

-  Average height (**average_height**);

-  Area under the curve (**area**);

-  Curve plateau point (**plateau**);

-  Curve slope during exponential metabolic activity (**slope**);

-  Curve lag time (**lag**).

All the parameters (except **min**, **max** and **average_height**)
require the `scipy library <https://www.scipy.org/>`__ to be installed.

The fit function uses three sigmoid functions:

Gompertz
   :math:`Ae^{-e^{(\frac{\mu_{m}e}{A}(\lambda - t) + 1)}} + y0`

Logistic
   :math:`\frac{A}{1+e^{(\frac{4\mu_{m}}{A}(\lambda - t) + 2)}} + y_{0}`

Richards
   :math:`A(1 + ve^{1 + v} + e^{\frac{\mu_{m}}{A}(1 + v)(1 + \frac{1}{v})(\lambda - t)})^{-\frac{1}{v}} + y0`

Where:

-  corresponds to the **plateau**

-  corresponds to the **slope**

-  corresponds to the **lag**

These functions have been derived from `this
publication <https://www.ncbi.nlm.nih.gov/pubmed/16348228>`__. The fit
method by default tries first to fit the gompertz function: if it fails
it will then try to fit the logistic and then the richards function. The
user can also specify one of the three functions to be applied.

.. code:: pycon

   >>> from Bio import phenotype
   >>> record = list(phenotype.parse("Plates.csv", "pm-csv"))[-1]
   >>> well = record["A02"]
   >>> well.fit()
   >>> print("Function fitted: %s" % well.model)
   Function fitted: gompertz
   >>> for param in ["area", "average_height", "lag", "max", "min", "plateau", "slope"]:
   ...     print("%s\t%.2f" % (param, getattr(well, param)))
   ...
   area    4414.38
   average_height  61.58
   lag     48.60
   max     143.00
   min     12.00
   plateau 120.02
   slope   4.99

Writing Phenotype Microarray data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PlateRecord objects can be written to file in the form of
`JSON <https://en.wikipedia.org/wiki/JSON>`__ files, a format compatible
with other software packages such as
`opm <https://www.dsmz.de/research/microorganisms/projects/analysis-of-omnilog-phenotype-microarray-data.html>`__
or `DuctApe <https://combogenomics.github.io/DuctApe/>`__.

.. code:: pycon

   >>> phenotype.write(record, "out.json", "pm-json")
   1
