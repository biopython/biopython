# Copyright 2004 by Harry Zuzan. All rights reserved.
# Copyright 2016 by Adam Kurkiewicz. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Reading information from Affymetrix CEL files version 3 and 4."""

from __future__ import print_function

import struct

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError
    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile")


class ParserError(ValueError):
    """Affymetrix parser error."""

    def __init__(self, *args):
        """Initialise class."""
        super(ParserError, self).__init__(*args)


_modeError = ParserError("You're trying to open an Affymetrix v4"
                         " CEL file. You have to use a read binary mode,"
                         " like this: open(filename, 'rb')")

# for debugging
# import pprint
# pp = pprint.PrettyPrinter(indent=4)


class Record(object):
    """Stores the information in a cel file.

    Example usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v3_example.CEL", "r") as handle:
    ...     c = CelFile.read(handle)
    ...
    >>> print(c.ncols, c.nrows)
    5 5
    >>> print(c.intensities)
    [[   234.    170.  22177.    164.  22104.]
     [   188.    188.  21871.    168.  21883.]
     [   188.    193.  21455.    198.  21300.]
     [   188.    182.  21438.    188.  20945.]
     [   193.  20370.    174.  20605.    168.]]
    >>> print(c.stdevs)
    [[   24.     34.5  2669.     19.7  3661.2]
     [   29.8    29.8  2795.9    67.9  2792.4]
     [   29.8    88.7  2976.5    62.   2914.5]
     [   29.8    76.2  2759.5    49.2  2762. ]
     [   38.8  2611.8    26.6  2810.7    24.1]]
    >>> print(c.npix)
    [[25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]
     [25 25 25 25 25]]

    """

    def __init__(self):
        """Initialize class."""
        self.version = None
        self.GridCornerUL = None
        self.GridCornerUR = None
        self.GridCornerLR = None
        self.GridCornerLL = None
        self.DatHeader = None
        self.Algorithm = None
        self.AlgorithmParameters = None
        self.NumberCells = None
        self.intensities = None
        self.stdevs = None
        self.npix = None
        self.nrows = None
        self.ncols = None
        self.nmask = None
        self.mask = None
        self.noutliers = None
        self.outliers = None
        self.modified = None


def read(handle):
    """Read Affymetrix CEL file and return Record object.

    CEL files version 3 and 4 are supported, and the parser attempts version detection.

    Example Usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v4_example.CEL", "rb") as handle:
    ...     c = CelFile.read(handle)
    ...
    >>> c.version == 4
    True

    """
    # If we fail to read the magic number, then it will remain None, and thus
    # we will invoke read_v3 (if mode is not strict), or raise IOError if mode
    # is strict.
    magicNumber = None
    # We check if the handle is a file-like object. If it isn't, and the mode
    # is strict, we raise an error. If it isn't and the mode isn't strict, we
    # continue (perhaps somebody has got a CEL-file-like iterable, which used
    # to work with previous versions of biopython and we don't want to maintain
    # backwards compatibility).
    try:
        mode = handle.mode
        # By definition an Affymetrix v4 CEL file has 64 as the first 4 bytes.
        # Note that we use little-endian irrespective of the platform, again by
        # definition.
        position = handle.tell()
        magicNumber = struct.unpack("<i", handle.read(4))[0]
    except (AttributeError, TypeError):
        pass
    except UnicodeDecodeError:
        raise _modeError
    finally:
        try:
            # reset the offset, to avoid breaking either v3 or v4.
            handle.seek(position)
        except AttributeError:
            pass
    if magicNumber != 64:
        return read_v3(handle)

    else:
        return read_v4(handle)


# read Affymetrix files version 4.
def read_v4(f):
    """Read verion 4 Affymetrix CEL file, returns corresponding Record object.

    Most importantly record.intensities correspond to intensities from the CEL
    file.

    record.mask and record.outliers are not set.

    Example Usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v4_example.CEL", "rb") as handle:
    ...     c = CelFile.read_v4(handle)
    ...
    >>> c.version == 4
    True
    >>> print("%i by %i array" % c.intensities.shape)
    5 by 5 array

    """
    # We follow the documentation here:
    # http://www.affymetrix.com/estore/support/developer/powertools/changelog/gcos-agcc/cel.html.affx
    record = Record()
    preHeaders = ["magic", "version", "columns", "rows", "cellNo", "headerLen"]
    preHeadersMap = {}
    headersMap = {}

    # load pre-headers
    try:
        for name in preHeaders:
            preHeadersMap[name] = struct.unpack("<i", f.read(4))[0]
    except UnicodeDecodeError:
        raise _modeError

    char = f.read(preHeadersMap["headerLen"])
    header = char.decode("ascii", "ignore")
    for header in header.split("\n"):
        if "=" in header:
            header = header.split("=")
            headersMap[header[0]] = "=".join(header[1:])

    # for debugging
    # pp.pprint("preHeadersMap")
    # pp.pprint(preHeadersMap)
    # pp.pprint("headersMap")
    # pp.pprint(headersMap)

    record.version = preHeadersMap["version"]
    if record.version != 4:
        raise ParserError("You are trying to parse CEL file version 4. This"
                          " file violates the structure expected from CEL file"
                          " version 4")
    record.GridCornerUL = headersMap["GridCornerUL"]
    record.GridCornerUR = headersMap["GridCornerUR"]
    record.GridCornerLR = headersMap["GridCornerLR"]
    record.GridCornerLL = headersMap["GridCornerLL"]
    record.DatHeader = headersMap["DatHeader"]
    record.Algorithm = headersMap["Algorithm"]
    record.AlgorithmParameters = headersMap["AlgorithmParameters"]
    record.NumberCells = preHeadersMap["cellNo"]
    # record.intensities are set below
    # record.stdevs are set below
    # record.npix are set below
    record.nrows = int(headersMap["Rows"])
    record.ncols = int(headersMap["Cols"])

    # These cannot be reliably set in v4, because of discrepancies between real
    # data and the documented format.
    record.nmask = None
    record.mask = None
    record.noutliers = None
    record.outliers = None
    record.modified = None

    # Real data never seems to have anything but zeros here, but we don't want
    # to take chances. Raising an error is better than returning unreliable
    # data.
    def raiseBadHeader(field, expected):
        actual = int(headersMap[field])
        message = "The header {field} is expected to be 0, not {value}".format(value=actual, field=field)
        if actual != expected:
            raise ParserError(message)

    raiseBadHeader("Axis-invertX", 0)

    raiseBadHeader("AxisInvertY", 0)

    raiseBadHeader("OffsetX", 0)

    raiseBadHeader("OffsetY", 0)

    # This is unfortunately undocumented, but it turns out that real data has
    # the record.AlgorithmParameters repeated in the data section, until an
    # EOF, i.e. b"\x04".
    char = b"\x00"
    safetyValve = 10**4
    for i in range(safetyValve):
        char = f.read(1)
        # For debugging
        # print([i for i in char], end="")
        if char == b"\x04":
            break
        if i == safetyValve:
            raise ParserError("Parse Error. The parser expects a short, "
                              "undocumented binary blob terminating with "
                              "ASCII EOF, x04")

    # After that there are precisely 15 bytes padded. Again, undocumented.
    padding = f.read(15)

    # That's how we pull out the values (triplets of the form float, float,
    # signed short).
    structa = struct.Struct("< f f h")

    # There are 10 bytes in our struct.
    structSize = 10

    # We initialize the most important: intensities, stdevs and npixs.
    record.intensities = numpy.empty(record.NumberCells, dtype=float)
    record.stdevs = numpy.empty(record.NumberCells, dtype=float)
    record.npix = numpy.empty(record.NumberCells, dtype=int)

    b = f.read(structSize * record.NumberCells)
    for i in range(record.NumberCells):
        binaryFragment = b[i * structSize: (i + 1) * structSize]
        intensity, stdevs, npix = structa.unpack(binaryFragment)
        record.intensities[i] = intensity
        record.stdevs[i] = stdevs
        record.npix[i] = npix

    # reshape without copying.
    def reshape(array):
        view = array.view()
        view.shape = (record.nrows, record.ncols)
        return view

    record.intensities = reshape(record.intensities)
    record.stdevs = reshape(record.stdevs)
    record.npix = reshape(record.npix)

    return record


def read_v3(handle):
    """Read version 3 Affymetrix CEL file, and return corresponding Record object.

    Example Usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v3_example.CEL", "r") as handle:
    ...     c = CelFile.read_v3(handle)
    ...
    >>> c.version == 3
    True

    """
    # Needs error handling.
    # Needs to know the chip design.
    record = Record()
    section = ""
    for line in handle:
        if not line.strip():
            continue
        # Set current section
        if line[:5] == "[CEL]":
            section = "CEL"
        elif line[:8] == "[HEADER]":
            section = "HEADER"
        elif line[:11] == "[INTENSITY]":
            section = "INTENSITY"
            record.intensities = numpy.zeros((record.nrows, record.ncols))
            record.stdevs = numpy.zeros((record.nrows, record.ncols))
            record.npix = numpy.zeros((record.nrows, record.ncols), int)
        elif line[:7] == "[MASKS]":
            section = "MASKS"
            record.mask = numpy.zeros((record.nrows, record.ncols))
        elif line[:10] == "[OUTLIERS]":
            section = "OUTLIERS"
            record.outliers = numpy.zeros((record.nrows, record.ncols))
        elif line[:10] == "[MODIFIED]":
            section = "MODIFIED"
            record.modified = numpy.zeros((record.nrows, record.ncols))
        elif line[0] == "[":
            # This would be an unknown section
            section = ""
        elif section == "CEL":
            keyword, value = line.split("=", 1)
            if keyword == "Version":
                record.version = int(value)
        elif section == "HEADER":
            # Set record.ncols and record.nrows, remaining data goes into
            # record.header dict
            keyword, value = line.split("=", 1)
            if keyword == "Cols":
                record.ncols = int(value)
            elif keyword == "Rows":
                record.nrows = int(value)
            elif keyword == "GridCornerUL":
                x, y = value.split()
                record.GridCornerUL = (int(x), int(y))
            elif keyword == "GridCornerUR":
                x, y = value.split()
                record.GridCornerUR = (int(x), int(y))
            elif keyword == "GridCornerLR":
                x, y = value.split()
                record.GridCornerLR = (int(x), int(y))
            elif keyword == "GridCornerLL":
                x, y = value.split()
                record.GridCornerLL = (int(x), int(y))
            elif keyword == "DatHeader":
                record.DatHeader = value.strip("\n\r")
            elif keyword == "Algorithm":
                record.Algorithm = value.strip("\n\r")
            elif keyword == "AlgorithmParameters":
                record.AlgorithmParameters = value.strip("\n\r")
        elif section == "INTENSITY":
            if "NumberCells" in line:
                record.NumberCells = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.intensities[x, y] = float(words[2])
                record.stdevs[x, y] = float(words[3])
                record.npix[x, y] = int(words[4])
        elif section == "MASKS":
            if "NumberCells" in line:
                record.nmask = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.mask[x, y] = int(1)
        elif section == "OUTLIERS":
            if "NumberCells" in line:
                record.noutliers = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.outliers[x, y] = int(1)
        elif section == "MODIFIED":
            if "NumberCells" in line:
                record.nmodified = int(line.split("=", 1)[1])
            elif "CellHeader" in line:
                pass
            else:
                words = line.split()
                y = int(words[0])
                x = int(words[1])
                record.modified[x, y] = float(words[2])
        else:
            continue
    return record


if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()
