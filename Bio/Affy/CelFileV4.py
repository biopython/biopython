from collections import namedtuple
from numpy import zeros
from struct import calcsize, unpack
import sys

# TODO(hammer): define formats in a variable (e.g. int = '<i') to facilitate reuse
# TODO(hammer): use struct.calcsize to determine number of bytes to read

class Record(object):
  """
  Stores the information in a cel file
  """
  def __init__(self):
    self.HEADER = Header()
    self.INTENSITY = Intensity()
    self.MASKS = None
    self.OUTLIERS = None

  def __repr__(self):
    state = ['HEADER=%s' % self.HEADER,
             'INTENSITY=%s' % self.INTENSITY,
             'MASKS=%s' % self.MASKS,
             'OUTLIERS=%s' % self.OUTLIERS,
            ]
    return 'Record(%s)' % ','.join(state)

class Header(object):
  """
  Stores the information from the header of the CEL file.
  """
  def __init__(self):
    self.ncols = None
    self.nrows = None
    self.ncells = None
    self.header_entries = None
    self.algorithm_name = None
    self.algorithm_parameters = None
    self.cell_margin = None
    self.noutlier_cells = None
    self.nmasked_cells = None
    self.nsubgrids = None

  def __repr__(self):
    state = ['ncols=%d' % self.ncols, 
             'nrows=%d' % self.nrows,
             'ncells=%d' % self.ncells,
             'header_entries=%s' % self.header_entries,
             'algorithm_name=%s' % self.algorithm_name,
             'algorithm_parameters=%s' % self.algorithm_parameters,
             'cell_margin=%d' % self.cell_margin,
             'noutlier_cells=%d' % self.noutlier_cells,
             'nmasked_cells=%d' % self.nmasked_cells,
             'nsugrids=%d' % self.nsubgrids,
            ]
    return 'Header(%s)' % ','.join(state)

class Intensity(object):
  """
  Stores the information from the intensities section of the CEL file.
  """
  def __init__(self):
    self.mean = None
    self.stdev = None
    self.npixels = None

  def __repr__(self):
    state = ['mean=%s' % self.mean,
             'stdev=%s' % self.stdev,
             'npixels=%s' % self.npixels,
            ]
    return 'Intensity(%s)' % ','.join(state)


# Built using documentation from http://www.stat.lsa.umich.edu/~kshedden/Courses/Stat545/Notes/AffxFileFormats/
def read(handle):
  # Magic number. Always set to 64.
  if (unpack('<i', handle.read(4))[0] != 64):
    sys.exit("Invalid magic number. Expected 64.")

  # Version number. Always set to 4.
  if (unpack('<i', handle.read(4))[0] != 4):
    sys.exit("Invalid version number. Expected 4.")

  record = Record()

  # Number of columns.
  ncols = unpack('<i', handle.read(4))[0]
  record.HEADER.ncols = ncols

  # Number of rows.
  nrows = unpack('<i', handle.read(4))[0]
  record.HEADER.nrows = nrows

  # Number of cells (rows*cols).
  ncells = unpack('<i', handle.read(4))[0]
  if (ncells != (nrows * ncols)):
    sys.exit("Number of cells does not match computed number of cells, %d." % nrows * ncells)
  record.HEADER.ncells = ncells

  # Header length
  header_length = unpack('<i', handle.read(4))[0]

  # Header as defined in the HEADER section of the version 3 CEL files. The string contains TAG=VALUE separated by a space where the TAG names are defined in the version 3 HEADER section.
  # TODO(hammer): parse the DAT header
  header = unpack('<%dc' % header_length, handle.read(header_length))
  header_entries = [line.split('=', 1) for line in ''.join(header).splitlines()]
  record.HEADER.header_entries = header_entries

  # Algorithm name length.
  algorithm_name_length = unpack('<i', handle.read(4))[0]

  # The algorithm name used to create the CEL file.
  algorithm_name = unpack('<%dc' % algorithm_name_length, handle.read(algorithm_name_length))
  algorithm_name = ''.join(algorithm_name)
  record.HEADER.algorithm_name = algorithm_name

  # Algorithm parameters length.
  algorithm_parameters_length = unpack('<i', handle.read(4))[0]

  # The parameters used by the algorithm. The format is TAG:VALUE pairs separated by semi-colons or TAG=VALUE pairs separated by spaces.
  # TODO(hammer): detect delimiter format and parse
  algorithm_parameters = unpack('<%dc' % algorithm_parameters_length, handle.read(algorithm_parameters_length))
  record.HEADER.algorithm_parameters = ''.join(algorithm_parameters)

  # Cell margin used for computing the cells intensity value.
  cell_margin = unpack('<i', handle.read(4))[0]
  record.HEADER.cell_margin = cell_margin

  # Number of outlier cells.
  noutlier_cells = unpack('<I', handle.read(4))[0]
  record.HEADER.noutlier_cells = noutlier_cells

  # Number of masked cells.
  nmasked_cells = unpack('<I', handle.read(4))[0]
  record.HEADER.nmasked_cells = nmasked_cells

  # Number of sub-grids.
  nsubgrids = unpack('<i', handle.read(4))[0]
  record.HEADER.nsubgrids = nsubgrids

  # Cell entries - this consists of an intensity value, standard deviation value and pixel count for each cell in the array.
  # The values are stored by row then column starting with the X=0, Y=0 cell. As an example, the first five entries are for cells defined by XY coordinates: (0,0), (1,0), (2,0), (3,0), (4,0).< /p>
  record.INTENSITY.mean = zeros((record.HEADER.nrows, record.HEADER.ncols))
  record.INTENSITY.stdev = zeros((record.HEADER.nrows, record.HEADER.ncols))
  record.INTENSITY.npixels = zeros((record.HEADER.nrows, record.HEADER.ncols), int)
  for x in range(nrows):
    for y in range(ncols):
      record.INTENSITY.mean[x, y], record.INTENSITY.stdev[x, y], record.INTENSITY.npixels[x, y] = unpack('<ffh', handle.read(10))

  # Masked entries - this consists of the XY coordinates of those cells masked by the user.
  record.MASKS = []
  for i in range(nmasked_cells):
    record.MASKS.append(unpack('<hh', handle.read(4)))

  # Outlier entries - this consists of the XY coordinates of those cells called outliers by the software.
  record.OUTLIERS = []
  for i in range(noutlier_cells):
    record.OUTLIERS.append(unpack('<hh', handle.read(4)))

  # TODO(hammer): do something with subgrid entries
  # Sub-grid entries - This is the sub-grid definition. There are as many sub-grids in the file as defined by the number of sub-grids above. Each sub-grid is defined as:
  # - row number (integer)
  # - column number (integer)
  # - upper left x coordinate in pixels (float)
  # - upper left y coordinate in pixels (float)
  # - upper right x coordinate in pixels (float)
  # - upper right x coordinate in pixels (float) (sic)
  # - lower left x coordinate in pixels (float)
  # - lower left y coordinate in pixels (float)
  # - lower right x coordinate in pixels (float)
  # - lower right x coordinate in pixels (float) (sic)
  # - left cell position (integer)
  # - top cell position (integer)
  # - right cell position (integer)
  # - bottom cell position (integer)
  subgrid_entries = []
  SubgridEntry = namedtuple('SubgridEntry', 'rnum cnum ulx uly urx ury llx lly lrx lry left top right bottom')
  for i in range(nsubgrids):
    subgrid_entries.append(SubgridEntry._make(unpack('<iiffffffffiiii', handle.read(56))))

  return record


if __name__ == "__main__":
  parsed_record = read(open(sys.argv[1], "rb"))
  print parsed_record
