from Numeric import *

from Bio.SVDSuperimposer import *

# start with two coordinate sets (Nx3 arrays - Float0)

x=array([[51.65, -1.90, 50.07],
	 [50.40, -1.23, 50.65],
	 [50.68, -0.04, 51.54],
	 [50.22, -0.02, 52.85]], Float0)

y=array([[51.30, -2.99, 46.54],
	 [51.09, -1.88, 47.58],
	 [52.36, -1.20, 48.03],
	 [52.71, -1.18, 49.38]], Float0)

sup=SVDSuperimposer()

# set the coords
# y will be rotated and translated on x
sup.set(x, y)

# do the lsq fit
sup.run()

# get the rmsd
rms=sup.get_rms()

# get rotation (right multiplying!) and the translation
rot, tran=sup.get_rotran()

# rotate y on x manually
y_on_x1=matrixmultiply(y, rot)+tran

# same thing
y_on_x2=sup.get_transformed()

# output results
print y_on_x1
print
print y_on_x2
print
print "%.2f" % rms


	
