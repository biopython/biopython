# Copyright 2014-2016 by Marco Galardini.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Classes to work with Phenotype Microarray data.

More information on the single plates can be found here: http://www.biolog.com/

Classes:
 - PlateRecord - Object that contain time course data on each well of the
   plate, as well as metadata (if any).
 - WellRecord - Object that contains the time course data of a single well
 - JsonWriter - Writer of PlateRecord objects in JSON format.

Functions:
 - JsonIterator -  Incremental PM JSON parser, this is an iterator that returns
   PlateRecord objects.
 - CsvIterator - Incremental PM CSV parser, this is an iterator that returns
   PlateRecord objects.
 - _toOPM - Used internally by JsonWriter, converts PlateRecord objects in
   dictionaries ready to be serialized in JSON format.

"""

import warnings
import json
import csv
import numpy as np

from Bio import BiopythonParserWarning

# Private csv headers - hardcoded because this are supposedly never changed
_datafile = "Data File"
_plate = "Plate Type"
_strainType = "Strain Type"
_sample = "Sample Number"
_strainName = "Strain Name"
_strainNumber = "Strain Number"
_other = "Other"
_hour = "Hour"
_file = "File"
_position = "Position"
_setupTime = "Setup Time"

_platesPrefix = "PM"
_platesPrefixMammalian = "PM-M"
#

# Json identifiers - hardcoded as they are set by the creators of opm
_csvData = "csv_data"
_measurements = "measurements"
#


class PlateRecord:
    """PlateRecord object for storing Phenotype Microarray plates data.

    A PlateRecord stores all the wells of a particular phenotype
    Microarray plate, along with metadata (if any). The single wells can be
    accessed calling their id as an index or iterating on the PlateRecord:

    >>> from Bio import phenotype
    >>> plate = phenotype.read("phenotype/Plate.json", "pm-json")
    >>> well = plate['A05']
    >>> for well in plate:
    ...    print(well.id)
    ...
    A01
    ...

    The plate rows and columns can be queried with an indexing system similar
    to NumPy and other matrices:

    >>> print(plate[1])
    Plate ID: PM01
    Well: 12
    Rows: 1
    Columns: 12
    PlateRecord('WellRecord['B01'], WellRecord['B02'], WellRecord['B03'], ..., WellRecord['B12']')

    >>> print(plate[:,1])
    Plate ID: PM01
    Well: 8
    Rows: 8
    Columns: 1
    PlateRecord('WellRecord['A02'], WellRecord['B02'], WellRecord['C02'], ..., WellRecord['H02']')

    Single WellRecord objects can be accessed using this indexing system:

    >>> print(plate[1,2])
    Plate ID: PM01
    Well ID: B03
    Time points: 384
    Minum signal 0.00 at time 11.00
    Maximum signal 76.25 at time 18.00
    WellRecord('(0.0, 11.0), (0.25, 11.0), (0.5, 11.0), (0.75, 11.0), (1.0, 11.0), ..., (95.75, 11.0)')

    The presence of a particular well can be inspected with the "in" keyword:
    >>> 'A01' in plate
    True

    All the wells belonging to a "row" (identified by the first character of
    the well id) in the plate can be obtained:

    >>> for well in plate.get_row('H'):
    ...     print(well.id)
    ...
    H01
    H02
    H03
    ...

    All the wells belonging to a "column" (identified by the number of the well)
    in the plate can be obtained:

    >>> for well in plate.get_column(12):
    ...     print(well.id)
    ...
    A12
    B12
    C12
    ...

    Two PlateRecord objects can be compared: if all their wells are equal the
    two plates are considered equal:

    >>> plate2 = phenotype.read("phenotype/Plate.json", "pm-json")
    >>> plate == plate2
    True

    Two PlateRecord object can be summed up or subtracted from each other: the
    the signals of each well will be summed up or subtracted. The id of the
    left operand will be kept:

    >>> plate3 = plate + plate2
    >>> print(plate3.id)
    PM01

    Many Phenotype Microarray plate have a "negative control" well, which can
    be subtracted to all wells:

    >>> subplate = plate.subtract_control()

    """

    def __init__(self, plateid, wells=None):
        """Initialize the class."""
        self.id = plateid

        if wells is None:
            wells = []

        # Similar behaviour as GenBank
        # Contains all the attributes
        self.qualifiers = {}

        # Well_id --> WellRecord objects
        self._wells = {}
        try:
            for w in wells:
                self._is_well(w)
                self[w.id] = w
        except TypeError:
            raise TypeError(
                "You must provide an iterator-like object containing the single wells"
            )

        self._update()

    def _update(self):
        """Update the rows and columns string identifiers (PRIVATE)."""
        self._rows = sorted({x[0] for x in self._wells})
        self._columns = sorted({x[1:] for x in self._wells})

    def _is_well(self, obj):
        """Check if the given object is a WellRecord object (PRIVATE).

        Used both for the class constructor and the __setitem__ method
        """
        # Value should be of WellRecord type
        if not isinstance(obj, WellRecord):
            raise ValueError(
                f"A WellRecord type object is needed as value (got {type(obj)})"
            )

    def __getitem__(self, index):
        """Access part of the plate.

        Depending on the indices, you can get a WellRecord object
        (representing a single well of the plate),
        or another plate
        (representing some part or all of the original plate).

        plate[wid] gives a WellRecord (if wid is a WellRecord id)
        plate[r,c] gives a WellRecord
        plate[r] gives a row as a PlateRecord
        plate[r,:] gives a row as a PlateRecord
        plate[:,c] gives a column as a PlateRecord

        plate[:] and plate[:,:] give a copy of the plate

        Anything else gives a subset of the original plate, e.g.
        plate[0:2] or plate[0:2,:] uses only row 0 and 1
        plate[:,1:3] uses only columns 1 and 2
        plate[0:2,1:3] uses only rows 0 & 1 and only cols 1 & 2

        >>> from Bio import phenotype
        >>> plate = phenotype.read("phenotype/Plate.json", "pm-json")

        You can access a well of the plate, using its id.

        >>> w = plate['A01']

        You can access a row of the plate as a PlateRecord using an integer
        index:

        >>> first_row = plate[0]
        >>> print(first_row)
        Plate ID: PM01
        Well: 12
        Rows: 1
        Columns: 12
        PlateRecord('WellRecord['A01'], WellRecord['A02'], WellRecord['A03'], ..., WellRecord['A12']')
        >>> last_row = plate[-1]
        >>> print(last_row)
        Plate ID: PM01
        Well: 12
        Rows: 1
        Columns: 12
        PlateRecord('WellRecord['H01'], WellRecord['H02'], WellRecord['H03'], ..., WellRecord['H12']')

        You can also access use python's slice notation to sub-plates
        containing only some of the plate rows:

        >>> sub_plate = plate[2:5]
        >>> print(sub_plate)
        Plate ID: PM01
        Well: 36
        Rows: 3
        Columns: 12
        PlateRecord('WellRecord['C01'], WellRecord['C02'], WellRecord['C03'], ..., WellRecord['E12']')

        This includes support for a step, i.e. plate[start:end:step], which
        can be used to select every second row:

        >>> sub_plate = plate[::2]

        You can also use two indices to specify both rows and columns.
        Using simple integers gives you the single wells. e.g.

        >>> w = plate[3, 4]
        >>> print(w.id)
        D05

        To get a single column use this syntax:

        >>> sub_plate = plate[:, 4]
        >>> print(sub_plate)
        Plate ID: PM01
        Well: 8
        Rows: 8
        Columns: 1
        PlateRecord('WellRecord['A05'], WellRecord['B05'], WellRecord['C05'], ..., WellRecord['H05']')

        Or, to get part of a column,

        >>> sub_plate = plate[1:3, 4]
        >>> print(sub_plate)
        Plate ID: PM01
        Well: 2
        Rows: 2
        Columns: 1
        PlateRecord(WellRecord['B05'], WellRecord['C05'])

        However, in general you get a sub-plate,

        >>> print(plate[1:5, 3:6])
        Plate ID: PM01
        Well: 12
        Rows: 4
        Columns: 3
        PlateRecord('WellRecord['B04'], WellRecord['B05'], WellRecord['B06'], ..., WellRecord['E06']')

        This should all seem familiar to anyone who has used the NumPy
        array or matrix objects.
        """
        # Well identifier access
        if isinstance(index, str):
            try:
                return self._wells[index]
            except KeyError:
                raise KeyError(f"Well {index} not found!")

        # Integer index
        elif isinstance(index, int):
            try:
                row = self._rows[index]
            except IndexError:
                raise IndexError("Row %d not found!" % index)
            return PlateRecord(
                self.id, filter(lambda x: x.id.startswith(row), self._wells.values())
            )

        # Slice
        elif isinstance(index, slice):
            rows = self._rows[index]
            return PlateRecord(
                self.id, filter(lambda x: x.id[0] in rows, self._wells.values())
            )

        # Other access
        elif len(index) != 2:
            raise TypeError("Invalid index type.")

        row_index, col_index = index
        if isinstance(row_index, int) and isinstance(col_index, int):
            # Return a single WellRecord
            try:
                row = self._rows[row_index]
            except IndexError:
                raise IndexError("Row %d not found!" % row_index)
            try:
                col = self._columns[col_index]
            except IndexError:
                raise IndexError("Column %d not found!" % col_index)

            return self._wells[row + col]

        elif isinstance(row_index, int):
            try:
                row = self._rows[row_index]
            except IndexError:
                raise IndexError("Row %d not found!" % row_index)
            cols = self._columns[col_index]

            return PlateRecord(
                self.id,
                filter(
                    lambda x: x.id.startswith(row) and x.id[1:] in cols,
                    self._wells.values(),
                ),
            )

        elif isinstance(col_index, int):
            try:
                col = self._columns[col_index]
            except IndexError:
                raise IndexError("Columns %d not found!" % col_index)
            rows = self._rows[row_index]

            return PlateRecord(
                self.id,
                filter(
                    lambda x: x.id.endswith(col) and x.id[0] in rows,
                    self._wells.values(),
                ),
            )

        else:
            rows = self._rows[row_index]
            cols = self._columns[col_index]

            return PlateRecord(
                self.id,
                filter(
                    lambda x: x.id[0] in rows and x.id[1:] in cols, self._wells.values()
                ),
            )

    def __setitem__(self, key, value):
        if not isinstance(key, str):
            raise ValueError("Well identifier should be string-like")
        self._is_well(value)
        # Provided key and well ID should be the same
        if value.id != key:
            raise ValueError(
                "WellRecord ID and provided key are different (got '%s' and '%s')"
                % (type(value.id), type(key))
            )
        self._wells[key] = value

        self._update()

    def __delitem__(self, key):
        if not isinstance(key, str):
            raise ValueError("Well identifier should be string-like")
        del self._wells[key]

        self._update()

    def __iter__(self):
        for well in sorted(self._wells):
            yield self._wells[well]

    def __contains__(self, wellid):
        if wellid in self._wells:
            return True
        return False

    def __len__(self):
        """Return the number of wells in this plate."""
        return len(self._wells)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self._wells == other._wells
        else:
            return False

    def __add__(self, plate):
        """Add another PlateRecord object.

        The wells in both plates must be the same

        A new PlateRecord object is returned, having the same id as the
        left operand.
        """
        if not isinstance(plate, PlateRecord):
            raise TypeError("Expecting a PlateRecord object")

        if {x.id for x in self} != {x.id for x in plate}:
            raise ValueError("The two plates have different wells")

        wells = []

        for w in self:
            wells.append(w + plate[w.id])

        newp = PlateRecord(self.id, wells=wells)

        return newp

    def __sub__(self, plate):
        """Subtract another PlateRecord object.

        The wells in both plates must be the same

        A new PlateRecord object is returned, having the same id as the
        left operand.
        """
        if not isinstance(plate, PlateRecord):
            raise TypeError("Expecting a PlateRecord object")

        if {x.id for x in self} != {x.id for x in plate}:
            raise ValueError("The two plates have different wells")

        wells = []

        for w in self:
            wells.append(w - plate[w.id])

        newp = PlateRecord(self.id, wells=wells)

        return newp

    def get_row(self, row):
        """Get all the wells of a given row.

        A row is identified with a letter (e.g. 'A')
        """
        # Key is casted to str implicitly
        try:
            row = str(row)
        except Exception:
            # Is it even possible to get an exception here?
            raise ValueError("Row identifier should be string-like")
        if len(row) > 1:
            raise ValueError("Row identifier must be of maximum one letter")

        for w in sorted(filter(lambda x: x.startswith(row), self._wells)):
            yield self._wells[w]

    def get_column(self, column):
        """Get all the wells of a given column.

        A column is identified with a number (e.g. '6')
        """
        # Column is casted to int implicitly
        try:
            column = int(column)
        except Exception:
            raise ValueError("Column identifier should be a number")

        # A 96-well plate has well numbers in two digits
        for w in sorted(filter(lambda x: x.endswith("%02d" % column), self._wells)):
            yield self._wells[w]

    def subtract_control(self, control="A01", wells=None):
        """Subtract a 'control' well from the other plates wells.

        By default the control is subtracted to all wells, unless
        a list of well ID is provided

        The control well should belong to the plate
        A new PlateRecord object is returned
        """
        if control not in self:
            raise ValueError("Control well not present in plate")
        wcontrol = self[control]

        if wells is None:
            wells = self._wells.keys()

        missing = {w for w in wells if w not in self}
        if missing:
            raise ValueError("Some wells to be subtracted are not present")

        nwells = []

        for w in self:
            if w.id in wells:
                nwells.append(w - wcontrol)
            else:
                nwells.append(w)

        newp = PlateRecord(self.id, wells=nwells)

        return newp

    def __repr__(self):
        """Return a (truncated) representation of the plate for debugging."""
        if len(self._wells) > 4:
            # Show the last well and the first three
            return "%s('%s, ..., %s')" % (
                self.__class__.__name__,
                ", ".join(
                    [
                        "%s['%s']" % (self[x].__class__.__name__, self[x].id)
                        for x in sorted(self._wells.keys())[:3]
                    ]
                ),
                "%s['%s']"
                % (
                    self[sorted(self._wells.keys())[-1]].__class__.__name__,
                    self[sorted(self._wells.keys())[-1]].id,
                ),
            )
        else:
            return "%s(%s)" % (
                self.__class__.__name__,
                ", ".join(
                    [
                        "%s['%s']" % (self[x].__class__.__name__, self[x].id)
                        for x in sorted(self._wells.keys())
                    ]
                ),
            )

    def __str__(self):
        """Return a human readable summary of the record (string).

        The python built in function str works by calling the object's __str__
        method.  e.g.

        >>> from Bio import phenotype
        >>> record = next(phenotype.parse("phenotype/Plates.csv", "pm-csv"))
        >>> print(record)
        Plate ID: PM01
        Well: 96
        Rows: 8
        Columns: 12
        PlateRecord('WellRecord['A01'], WellRecord['A02'], WellRecord['A03'], ..., WellRecord['H12']')

        Note that long well lists are shown truncated.
        """
        lines = []
        if self.id:
            lines.append(f"Plate ID: {self.id}")
        lines.append("Well: %i" % len(self))
        # Here we assume that all well ID start with a char
        lines.append("Rows: %d" % len({x.id[0] for x in self}))
        # Here we assume that well number is a two-digit number
        lines.append("Columns: %d" % len({x.id[1:3] for x in self}))
        lines.append(repr(self))
        return "\n".join(lines)


class WellRecord:
    """WellRecord stores all time course signals of a phenotype Microarray well.

    The single time points and signals can be accessed iterating on the
    WellRecord or using lists indexes or slices:

    >>> from Bio import phenotype
    >>> plate = phenotype.read("phenotype/Plate.json", "pm-json")
    >>> well = plate['A05']
    >>> for time, signal in well:
    ...    print("Time: %f, Signal: %f" % (time, signal)) # doctest:+ELLIPSIS
    ...
    Time: 0.000000, Signal: 14.000000
    Time: 0.250000, Signal: 13.000000
    Time: 0.500000, Signal: 15.000000
    Time: 0.750000, Signal: 15.000000
    ...
    >>> well[1]
    16.0
    >>> well[1:5]
    [16.0, 20.0, 18.0, 15.0]
    >>> well[1:5:0.5]
    [16.0, 19.0, 20.0, 18.0, 18.0, 18.0, 15.0, 18.0]

    If a time point was not present in the input file but it's between the
    minimum and maximum time point, the interpolated signal is returned,
    otherwise a nan value:

    >>> well[1.3]
    19.0
    >>> well[1250]
    nan

    Two WellRecord objects can be compared: if their input time/signal pairs
    are exactly the same, the two records are considered equal:

    >>> well2 = plate['H12']
    >>> well == well2
    False

    Two WellRecord objects can be summed up or subtracted from each other: a new
    WellRecord object is returned, having the left operand id.

    >>> well1 = plate['A05']
    >>> well2 = well + well1
    >>> print(well2.id)
    A05

    If SciPy is installed, a sigmoid function can be fitted to the PM curve,
    in order to extract some parameters; three sigmoid functions are available:
    * gompertz
    * logistic
    * richards
    The functions are described in Zwietering et al., 1990 (PMID: 16348228)

    For example::

        well.fit()
        print(well.slope, well.model)
        (61.853516785566917, 'logistic')

    If not sigmoid function is specified, the first one that is successfully
    fitted is used. The user can also specify a specific function.

    To specify gompertz::

        well.fit('gompertz')
        print(well.slope, well.model)
        (127.94630059171354, 'gompertz')

    If no function can be fitted, the parameters are left as None, except for
    the max, min, average_height and area.
    """

    def __init__(self, wellid, plate=None, signals=None):
        """Initialize the class."""
        if plate is None:
            self.plate = PlateRecord(None)
        else:
            self.plate = plate

        self.id = wellid

        # Curve parameters (to be calculated with the "fit" function)
        # Parameters that don't need scipy
        self.max = None
        self.min = None
        self.average_height = None

        # Parameters that need scipy
        self.area = None
        self.plateau = None
        self.slope = None
        self.lag = None
        self.v = None
        self.y0 = None
        self.model = None

        # Original signals (private)
        if signals is None:
            self._signals = {}
        else:
            self._signals = signals

    def _interpolate(self, time):
        """Linear interpolation of the signals at certain time points (PRIVATE)."""
        times = sorted(self._signals.keys())

        return np.interp(
            time, times, [self._signals[x] for x in times], left=np.nan, right=np.nan
        )

    def __setitem__(self, time, signal):
        """Assign a signal at a certain time point."""
        try:
            time = float(time)
        except ValueError:
            raise ValueError("Time point should be a number")
        try:
            signal = float(signal)
        except ValueError:
            raise ValueError("Signal should be a number")

        self._signals[time] = signal

    def __getitem__(self, time):
        """Return a subset of signals or a single signal."""
        if isinstance(time, slice):
            # Fix the missing values in the slice
            if time.start is None:
                start = 0
            else:
                start = time.start

            if time.stop is None:
                stop = max(self.get_times())
            else:
                stop = time.stop

            time = np.arange(start, stop, time.step)
            return list(self._interpolate(time))

        elif isinstance(time, (float, int)):
            return self._interpolate(time)

        raise ValueError("Invalid index")

    def __iter__(self):
        for time in sorted(self._signals.keys()):
            yield time, self._signals[time]

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            if list(self._signals.keys()) != list(other._signals.keys()):
                return False
            # Account for the presence of NaNs
            for k in self._signals:
                if np.isnan(self[k]) and np.isnan(other[k]):
                    continue
                elif self[k] != other[k]:
                    return False
            return True
        else:
            return False

    def __add__(self, well):
        """Add another WellRecord object.

        A new WellRecord object is returned, having the same id as the
        left operand
        """
        if not isinstance(well, WellRecord):
            raise TypeError("Expecting a WellRecord object")

        signals = {}

        times = set(self._signals.keys()).union(set(well._signals.keys()))
        for t in sorted(times):
            signals[t] = self[t] + well[t]

        neww = WellRecord(self.id, signals=signals)

        return neww

    def __sub__(self, well):
        """Subtract another WellRecord object.

        A new WellRecord object is returned, having the same id as the
        left operand
        """
        if not isinstance(well, WellRecord):
            raise TypeError("Expecting a WellRecord object")

        signals = {}

        times = set(self._signals.keys()).union(set(well._signals.keys()))
        for t in sorted(times):
            signals[t] = self[t] - well[t]

        neww = WellRecord(self.id, signals=signals)

        return neww

    def __len__(self):
        """Return the number of time points sampled."""
        return len(self._signals)

    def __repr__(self):
        """Return a (truncated) representation of the signals for debugging."""
        if len(self) > 7:
            # Shows the last time point and the first five
            return "%s('%s, ..., %s')" % (
                self.__class__.__name__,
                ", ".join([str(x) for x in self.get_raw()[:5]]),
                str(self.get_raw()[-1]),
            )
        else:
            return "%s(%s)" % (
                self.__class__.__name__,
                ", ".join([str(x) for x in self.get_raw()]),
            )

    def __str__(self):
        """Return a human readable summary of the record (string).

        The python built-in function str works by calling the object's __str__
        method.  e.g.

        >>> from Bio import phenotype
        >>> plate = phenotype.read("phenotype/Plate.json", "pm-json")
        >>> record = plate['A05']
        >>> print(record)
        Plate ID: PM01
        Well ID: A05
        Time points: 384
        Minum signal 0.25 at time 13.00
        Maximum signal 19.50 at time 23.00
        WellRecord('(0.0, 14.0), (0.25, 13.0), (0.5, 15.0), (0.75, 15.0), (1.0, 16.0), ..., (95.75, 16.0)')

        Note that long time spans are shown truncated.
        """
        lines = []
        if self.plate and self.plate.id:
            lines.append(f"Plate ID: {self.plate.id}")
        if self.id:
            lines.append(f"Well ID: {self.id}")
        lines.append("Time points: %i" % len(self))
        lines.append("Minum signal %.2f at time %.2f" % min(self, key=lambda x: x[1]))
        lines.append("Maximum signal %.2f at time %.2f" % max(self, key=lambda x: x[1]))
        lines.append(repr(self))
        return "\n".join(lines)

    def get_raw(self):
        """Get a list of time/signal pairs."""
        return [(t, self._signals[t]) for t in sorted(self._signals.keys())]

    def get_times(self):
        """Get a list of the recorded time points."""
        return sorted(self._signals.keys())

    def get_signals(self):
        """Get a list of the recorded signals (ordered by collection time)."""
        return [self._signals[t] for t in sorted(self._signals.keys())]

    def fit(self, function=("gompertz", "logistic", "richards")):
        """Fit a sigmoid function to this well and extract curve parameters.

        If function is None or an empty tuple/list, then no fitting is done.
        Only the object's ``.min``, ``.max`` and ``.average_height`` are
        calculated.

        By default the following fitting functions will be used in order:
         - gompertz
         - logistic
         - richards

        The first function that is successfully fitted to the signals will
        be used to extract the curve parameters and update ``.area`` and
        ``.model``. If no function can be fitted an exception is raised.

        The function argument should be a tuple or list of any of these three
        function names as strings.

        There is no return value.
        """
        avail_func = ("gompertz", "logistic", "richards")

        # Parameters not dependent on curve fitting
        self.max = max(self, key=lambda x: x[1])[1]
        self.min = min(self, key=lambda x: x[1])[1]

        self.average_height = np.array(self.get_signals()).mean()

        if not function:
            self.area = None
            self.model = None
            return
        for sigmoid_func in function:
            if sigmoid_func not in avail_func:
                raise ValueError(f"Fitting function {sigmoid_func!r} not supported")

        # Parameters that depend on scipy curve_fit
        from .pm_fitting import fit, get_area
        from .pm_fitting import logistic, gompertz, richards

        function_map = {
            "logistic": logistic,
            "gompertz": gompertz,
            "richards": richards,
        }

        self.area = get_area(self.get_signals(), self.get_times())

        self.model = None
        for sigmoid_func in function:
            func = function_map[sigmoid_func]
            try:
                (self.plateau, self.slope, self.lag, self.v, self.y0), pcov = fit(
                    func, self.get_times(), self.get_signals()
                )

                self.model = sigmoid_func
                return
            except RuntimeError:
                continue
        raise RuntimeError("Could not fit any sigmoid function")


def JsonIterator(handle):
    """Iterate over PM json records as PlateRecord objects.

    Arguments:
     - handle - input file

    """
    try:
        data = json.load(handle)
    except ValueError:
        raise ValueError("Could not parse JSON file")

    # We can have one single plate or several
    # we need to discriminate
    if hasattr(data, "keys"):
        data = [data]

    for pobj in data:
        try:
            plateID = pobj[_csvData][_plate]
        except TypeError:
            raise TypeError("Malformed JSON input")
        except KeyError:
            raise KeyError("Could not retrieve plate id")

        # Parse also non-standard plate IDs
        if not plateID.startswith(_platesPrefix) and not plateID.startswith(
            _platesPrefixMammalian
        ):
            warnings.warn(
                f"Non-standard plate ID found ({plateID})", BiopythonParserWarning
            )
        else:
            # Simplify the plates IDs, removing letters, as opm does
            if plateID.startswith(_platesPrefixMammalian):
                pID = plateID[len(_platesPrefixMammalian) :]
            else:
                pID = plateID[len(_platesPrefix) :]
            while len(pID) > 0:
                try:
                    int(pID)
                    break
                except ValueError:
                    pID = pID[:-1]

            # No luck
            if len(pID) == 0:
                warnings.warn(
                    f"Non-standard plate ID found ({plateID})", BiopythonParserWarning
                )
            elif int(pID) < 0:
                warnings.warn(
                    f"Non-standard plate ID found ({plateID}), using {_platesPrefix}{abs(int(pID))}"
                )
                plateID = _platesPrefix + str(abs(int(pID)))
            else:
                if plateID.startswith(_platesPrefixMammalian):
                    plateID = _platesPrefixMammalian + "%02d" % int(pID)
                else:
                    plateID = _platesPrefix + "%02d" % int(pID)

        try:
            times = pobj[_measurements][_hour]
        except KeyError:
            raise KeyError("Could not retrieve the time points")

        plate = PlateRecord(plateID)

        for k in pobj[_measurements]:
            # Skip the time points
            if k == _hour:
                continue

            plate[k] = WellRecord(
                k,
                plate=plate,
                signals={
                    times[i]: pobj[_measurements][k][i] for i in range(len(times))
                },
            )

        # Remove the measurements and assign the other qualifiers
        del pobj["measurements"]
        plate.qualifiers = pobj

        yield plate


def CsvIterator(handle):
    """Iterate over PM csv records as PlateRecord objects.

    Arguments:
     - handle - input file

    """
    plate = None
    data = False
    qualifiers = {}
    idx = {}
    wells = {}

    tblreader = csv.reader(handle, delimiter=",", quotechar='"')
    for line in tblreader:
        if len(line) < 2:
            continue

        elif _datafile in line[0].strip():
            # Do we have a previous plate?
            if plate is not None:
                qualifiers[_csvData][_datafile] = line[1].strip()
                plate = PlateRecord(plate.id)
                for k, v in wells.items():
                    plate[k] = WellRecord(k, plate, v)
                plate.qualifiers = qualifiers
                yield plate
            plate = PlateRecord(None)
            data = False
            qualifiers[_csvData] = {}
            idx = {}
            wells = {}

        elif _plate in line[0].strip():
            plateID = line[1].strip()

            qualifiers[_csvData][_plate] = plateID

            # Parse also non-standard plate IDs
            if not plateID.startswith(_platesPrefix) and not plateID.startswith(
                _platesPrefixMammalian
            ):
                warnings.warn(
                    f"Non-standard plate ID found ({plateID})", BiopythonParserWarning
                )
            else:
                # Simplify the plates IDs, removing letters, as opm does
                if plateID.startswith(_platesPrefixMammalian):
                    pID = plateID[len(_platesPrefixMammalian) :]
                else:
                    pID = plateID[len(_platesPrefix) :]
                while len(pID) > 0:
                    try:
                        int(pID)
                        break
                    except ValueError:
                        pID = pID[:-1]

                # No luck
                if len(pID) == 0:
                    warnings.warn(
                        f"Non-standard plate ID found ({plateID})",
                        BiopythonParserWarning,
                    )
                elif int(pID) < 0:
                    warnings.warn(
                        f"Non-standard plate ID found ({plateID}), using {_platesPrefix}{abs(int(pID))}"
                    )
                    plateID = _platesPrefix + str(abs(int(pID)))
                else:
                    if plateID.startswith(_platesPrefixMammalian):
                        plateID = _platesPrefixMammalian + "%02d" % int(pID)
                    else:
                        plateID = _platesPrefix + "%02d" % int(pID)

            plate.id = plateID

        elif _strainType in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_strainType] = line[1].strip()

        elif _sample in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_sample] = line[1].strip()

        elif _strainNumber in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_strainNumber] = line[1].strip()

        elif _strainName in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_strainName] = line[1].strip()

        elif _other in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_other] = line[1].strip()

        elif _file in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_file] = line[1].strip()

        elif _position in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_position] = line[1].strip()

        elif _setupTime in line[0].strip():
            if plate is None:
                continue
            qualifiers[_csvData][_setupTime] = line[1].strip()

        elif _hour in line[0].strip():
            if plate is None:
                continue
            data = True
            for i in range(1, len(line)):
                x = line[i]
                if x == "":
                    continue
                wells[x.strip()] = {}
                idx[i] = x.strip()

        elif data:
            if plate is None:
                continue

            # Workaround for bad-formatted files
            try:
                float(line[0])
            except ValueError:
                continue

            time = float(line[0])

            for i in range(1, len(line)):
                x = line[i]

                try:
                    signal = float(x)
                except ValueError:
                    continue

                well = idx[i]
                wells[well][time] = signal

    if plate is not None and plate.id is not None:
        plate = PlateRecord(plate.id)
        for k, v in wells.items():
            plate[k] = WellRecord(k, plate, v)
        plate.qualifiers = qualifiers
        yield plate


def _toOPM(plate):
    """Transform a PlateRecord object into a dictionary (PRIVATE)."""
    d = dict(plate.qualifiers.items())

    d[_csvData] = {}
    d[_csvData][_plate] = plate.id
    d[_measurements] = {}
    d[_measurements][_hour] = []
    times = set()
    for wid, w in plate._wells.items():
        d[_measurements][wid] = []
        for hour in w._signals:
            times.add(hour)

    for hour in sorted(times):
        d[_measurements][_hour].append(hour)
        for wid, w in plate._wells.items():
            if hour in w._signals:
                d[_measurements][wid].append(w[hour])
            # This shouldn't happen
            else:
                d[_measurements][wid].append(float("nan"))

    return d


class JsonWriter:
    """Class to write PM Json format files."""

    def __init__(self, plates):
        """Initialize the class."""
        self.plates = plates

    def write(self, handle):
        """Write this instance's plates to a file handle."""
        out = []
        for plate in self.plates:
            try:
                out.append(_toOPM(plate))
            except ValueError:
                raise ValueError("Could not export plate(s) in JSON format")

        handle.write(json.dumps(out) + "\n")

        return len(out)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
