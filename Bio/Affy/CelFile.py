# Copyright 2004 by Harry Zuzan. All rights reserved.
# Copyright 2016 by Adam Kurkiewicz. All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Reading information from Affymetrix CEL files version 3 and 4."""


import struct

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.Affy.CelFile"
    ) from None


class ParserError(ValueError):
    """Affymetrix parser error."""

    def __init__(self, *args):
        """Initialise class."""
        super().__init__(*args)


class Record:
    """Stores the information in a cel file.

    Example usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v3_example.CEL") as handle:
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


def read(handle, version=None):
    """Read Affymetrix CEL file and return Record object.

    CEL files format versions 3 and 4 are supported.
    Please specify the CEL file format as 3 or 4 if known for the version
    argument. If the version number is not specified, the parser will attempt
    to detect the version from the file contents.

    The Record object returned by this function stores the intensities from
    the CEL file in record.intensities.
    Currently, record.mask and record.outliers are not set in when parsing
    version 4 CEL files.

    Example Usage:

    >>> from Bio.Affy import CelFile
    >>> with open("Affy/affy_v3_example.CEL") as handle:
    ...     record = CelFile.read(handle)
    ...
    >>> record.version == 3
    True
    >>> print("%i by %i array" % record.intensities.shape)
    5 by 5 array

    >>> with open("Affy/affy_v4_example.CEL", "rb") as handle:
    ...     record = CelFile.read(handle, version=4)
    ...
    >>> record.version == 4
    True
    >>> print("%i by %i array" % record.intensities.shape)
    5 by 5 array

    """
    try:
        data = handle.read(0)
    except AttributeError:
        raise ValueError("handle should be a file handle") from None
    data = handle.read(4)
    if not data:
        raise ValueError("Empty file.")
    if data == b"[CEL":
        raise ValueError("CEL file in version 3 format should be opened in text mode")
    if data == "[CEL":
        # Version 3 format. Continue to read the header here before passing
        # control to _read_v3 to avoid having to seek to the beginning of
        # the file.
        data += next(handle)
        if data.strip() != "[CEL]":
            raise ValueError("Failed to parse Affy Version 3 CEL file.")
        line = next(handle)
        keyword, value = line.split("=", 1)
        if keyword != "Version":
            raise ValueError("Failed to parse Affy Version 3 CEL file.")
        version = int(value)
        if version != 3:
            raise ValueError("Incorrect version number in Affy Version 3 CEL file.")
        return _read_v3(handle)
    try:
        magicNumber = struct.unpack("<i", data)
    except TypeError:
        raise ValueError(
            "CEL file in version 4 format should be opened in binary mode"
        ) from None
    except struct.error:
        raise ValueError(
            "Failed to read magic number from Affy Version 4 CEL file"
        ) from None
    if magicNumber != (64,):
        raise ValueError("Incorrect magic number in Affy Version 4 CEL file")
    return _read_v4(handle)


# read Affymetrix files version 4.
def read_v4(f):
    """Read version 4 Affymetrix CEL file, and return a Record object (DEPRECATED)."""
    raise Exception(
        "The read_v4 function in Bio.Affy.CelFile is deprecated. "
        "Instead, please use the read function in Bio.Affy.CelFile "
        "specifying version=4."
    )


def read_v3(handle):
    """Read version 3 Affymetrix CEL file, and return a Record object (DEPRECATED)."""
    raise Exception(
        "The read_v3 function in Bio.Affy.CelFile is deprecated. "
        "Instead, please use the read function in Bio.Affy.CelFile "
        "specifying version=3."
    )


def _read_v4(f):
    # We follow the documentation here:
    # http://www.affymetrix.com/estore/support/developer/powertools/changelog/gcos-agcc/cel.html.affx
    record = Record()
    preHeaders = ["version", "columns", "rows", "cellNo", "headerLen"]
    preHeadersMap = {}
    headersMap = {}

    # Load pre-headers. The magic number was already parsed in the read
    # function calling _read_v4.
    preHeadersMap["magic"] = 64
    try:
        for name in preHeaders:
            preHeadersMap[name] = struct.unpack("<i", f.read(4))[0]
    except struct.error:
        raise ParserError("Failed to parse CEL version 4 file") from None

    char = f.read(preHeadersMap["headerLen"])
    header = char.decode("ascii", "ignore")
    for header in header.split("\n"):
        if "=" in header:
            header = header.split("=")
            headersMap[header[0]] = "=".join(header[1:])

    record.version = preHeadersMap["version"]
    if record.version != 4:
        raise ParserError("Incorrect version number in CEL version 4 file")

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
        message = f"The header {field} is expected to be 0, not {actual}"
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
    safetyValve = 10 ** 4
    for i in range(safetyValve):
        char = f.read(1)
        # For debugging
        # print([i for i in char], end="")
        if char == b"\x04":
            break
        if i == safetyValve:
            raise ParserError(
                "Parse Error. The parser expects a short, "
                "undocumented binary blob terminating with "
                "ASCII EOF, x04"
            )

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
        binaryFragment = b[i * structSize : (i + 1) * structSize]
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


def _read_v3(handle):
    # Needs error handling.
    # Needs to know the chip design.
    record = Record()
    # The version number was already obtained when the read function calling
    # _read_v3 parsed the CEL section.
    record.version = 3
    section = ""
    for line in handle:
        line = line.rstrip("\r\n")
        if not line:
            continue
        # Set current section
        if line.startswith("[HEADER]"):
            section = "HEADER"
        elif line.startswith("[INTENSITY]"):
            section = "INTENSITY"
            record.intensities = numpy.zeros((record.nrows, record.ncols))
            record.stdevs = numpy.zeros((record.nrows, record.ncols))
            record.npix = numpy.zeros((record.nrows, record.ncols), int)
        elif line.startswith("[MASKS]"):
            section = "MASKS"
            record.mask = numpy.zeros((record.nrows, record.ncols), bool)
        elif line.startswith("[OUTLIERS]"):
            section = "OUTLIERS"
            record.outliers = numpy.zeros((record.nrows, record.ncols), bool)
        elif line.startswith("[MODIFIED]"):
            section = "MODIFIED"
            record.modified = numpy.zeros((record.nrows, record.ncols))
        elif line.startswith("["):
            raise ParserError("Unknown section found in version 3 CEL file")
        else:  # read the data in a section
            if section == "HEADER":
                # Set record.ncols and record.nrows, remaining data goes into
                # record.header dict
                key, value = line.split("=", 1)
                if key == "Cols":
                    record.ncols = int(value)
                elif key == "Rows":
                    record.nrows = int(value)
                elif key == "GridCornerUL":
                    x, y = value.split()
                    record.GridCornerUL = (int(x), int(y))
                elif key == "GridCornerUR":
                    x, y = value.split()
                    record.GridCornerUR = (int(x), int(y))
                elif key == "GridCornerLR":
                    x, y = value.split()
                    record.GridCornerLR = (int(x), int(y))
                elif key == "GridCornerLL":
                    x, y = value.split()
                    record.GridCornerLL = (int(x), int(y))
                elif key == "DatHeader":
                    # not sure if all parameters here are interpreted correctly
                    record.DatHeader = {}
                    index = line.find(":")
                    _, filename = line[:index].split()
                    record.DatHeader["filename"] = filename
                    index += 1
                    field = line[index : index + 9]
                    assert field[:4] == "CLS="
                    assert field[8] == " "
                    record.DatHeader["CLS"] = int(field[4:8])
                    index += 9
                    field = line[index : index + 9]
                    assert field[:4] == "RWS="
                    assert field[8] == " "
                    record.DatHeader["RWS"] = int(field[4:8])
                    index += 9
                    field = line[index : index + 7]
                    assert field[:4] == "XIN="
                    assert field[6] == " "
                    record.DatHeader["XIN"] = int(field[4:6])
                    index += 7
                    field = line[index : index + 7]
                    assert field[:4] == "YIN="
                    assert field[6] == " "
                    record.DatHeader["YIN"] = int(field[4:6])
                    index += 7
                    field = line[index : index + 6]
                    assert field[:3] == "VE="
                    assert field[5] == " "
                    record.DatHeader["VE"] = int(field[3:5])
                    index += 6
                    field = line[index : index + 7]
                    assert field[6] == " "
                    temperature = field[:6].strip()
                    if temperature:
                        record.DatHeader["temperature"] = int(temperature)
                    else:
                        record.DatHeader["temperature"] = None
                    index += 7
                    field = line[index : index + 4]
                    assert field.endswith(" ")
                    record.DatHeader["laser-power"] = float(field)
                    index += 4
                    field = line[index : index + 18]
                    assert field[8] == " "
                    record.DatHeader["scan-date"] = field[:8]
                    assert field[17] == " "
                    record.DatHeader["scan-time"] = field[9:17]
                    index += 18
                    field = line[index:]
                    subfields = field.split(" \x14 ")
                    assert len(subfields) == 12
                    subfield = subfields[0]
                    try:
                        scanner_id, scanner_type = subfield.split()
                    except ValueError:
                        scanner_id = subfield.strip()
                    record.DatHeader["scanner-id"] = scanner_id
                    record.DatHeader["scanner-type"] = scanner_type
                    record.DatHeader["array-type"] = subfields[2]
                    record.DatHeader["image-orientation"] = int(subfields[11])
                elif key == "Algorithm":
                    record.Algorithm = value
                elif key == "AlgorithmParameters":
                    parameters = value.split(";")
                    values = {}
                    for parameter in parameters:
                        key, value = parameter.split(":", 1)
                        if key in (
                            "Percentile",
                            "CellMargin",
                            "FullFeatureWidth",
                            "FullFeatureHeight",
                            "PoolWidthExtenstion",
                            "PoolHeightExtension",
                        ):
                            values[key] = int(value)
                        elif key in ("OutlierHigh", "OutlierLow", "StdMult"):
                            values[key] = float(value)
                        elif key in (
                            "FixedCellSize",
                            "IgnoreOutliersInShiftRows",
                            "FeatureExtraction",
                            "UseSubgrids",
                            "RandomizePixels",
                        ):
                            if value == "TRUE":
                                value = True
                            elif value == "FALSE":
                                value = False
                            else:
                                raise ValueError("Unexpected boolean value")
                            values[key] = value
                        elif key in ("AlgVersion", "ErrorBasis"):
                            values[key] = value
                        else:
                            raise ValueError("Unexpected tag in AlgorithmParameters")
                    record.AlgorithmParameters = values
            elif section == "INTENSITY":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.NumberCells = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y", "MEAN", "STDV", "NPIXELS"]:
                        raise ParserError(
                            "Unexpected CellHeader in INTENSITY "
                            "section CEL version 3 file"
                        )
                else:
                    words = line.split()
                    y = int(words[0])
                    x = int(words[1])
                    record.intensities[x, y] = float(words[2])
                    record.stdevs[x, y] = float(words[3])
                    record.npix[x, y] = int(words[4])
            elif section == "MASKS":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.nmask = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y"]:
                        raise ParserError(
                            "Unexpected CellHeader in MASKS "
                            "section in CEL version 3 file"
                        )
                else:
                    words = line.split()
                    y = int(words[0])
                    x = int(words[1])
                    record.mask[x, y] = True
            elif section == "OUTLIERS":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.noutliers = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y"]:
                        raise ParserError(
                            "Unexpected CellHeader in OUTLIERS "
                            "section in CEL version 3 file"
                        )
                else:
                    words = line.split()
                    y = int(words[0])
                    x = int(words[1])
                    record.outliers[x, y] = True
            elif section == "MODIFIED":
                if line.startswith("NumberCells="):
                    key, value = line.split("=", 1)
                    record.nmodified = int(value)
                elif line.startswith("CellHeader="):
                    key, value = line.split("=", 1)
                    if value.split() != ["X", "Y", "ORIGMEAN"]:
                        raise ParserError(
                            "Unexpected CellHeader in MODIFIED "
                            "section in CEL version 3 file"
                        )
                else:
                    words = line.split()
                    y = int(words[0])
                    x = int(words[1])
                    record.modified[x, y] = float(words[2])
    return record


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
