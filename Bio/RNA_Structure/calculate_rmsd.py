#!/usr/bin/env python

"""
Calculate RMSD between two XYZ files

by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Andersen Bratholm <larsbratholm@gmail.com>
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE

"""

import numpy as np
import re


def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P = rotate(P, Q)
    return(rmsd(P, Q))


def rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm
    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return(P)


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.

    Using the Kabsch algorithm with two sets of paired point P and Q,
    centered around the center-of-mass.
    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters:
    P -- (N, number of points)x(D, dimension) matrix
    Q -- (N, number of points)x(D, dimension) matrix

    Returns:
    U -- Rotation matrix

    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U = np.dot(V, W)

    return(U)

# based on doi:10.1016/1049-9660(91)90036-O
def quaternion_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    """
    P_rot = quaternion_rotate(P, Q)
    return(rmsd(P_rot, Q))

# get optimal rotation (translation will be zero when the centroids of each molecule are the same)
def quaternion_transform(r):
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3,:3]
    return(rot)

# matrix involved in quaternion rotation
def makeW(r1,r2,r3,r4=0):
    W = np.asarray([
             [r4, r3, -r2, r1],
             [-r3, r4, r1, r2],
             [r2, -r1, r4, r3],
             [-r1, -r2, -r3, r4] ])
    return(W)

# matrix involved in quaternion rotation
def makeQ(r1,r2,r3,r4=0):
    Q = np.asarray([
             [r4, -r3, r2, r1],
             [r3, r4, -r1, r2],
             [-r2, r1, r4, r3],
             [-r1, -r2, -r3, r4] ])
    return(Q)

# calculate the rotation
def quaternion_rotate(X, Y):
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T,W[k]) for k in range(N)])
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    C1 = -np.sum(Qt_dot_W,axis=0)
    C2 = 0.5*N
    C3 = np.sum(W_minus_Q,axis=0)
    A = np.dot(C3.T,C3)*C2-C1
    eigen = np.linalg.eigh(A)
    r = eigen[1][:,eigen[0].argmax()]
    rot = quaternion_transform(r)
    return(np.dot(X,rot))

def centroid(X):
    """
    Calculate the centroid from a vectorset X
    """
    C = sum(X)/len(X)
    return(C)


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return(np.sqrt(rmsd/N))


def write_coordinates(atoms, V):
    """
    Print coordinates V
    """
    N, D = V.shape

    print(str(N))
    print

    for i in xrange(N):
        line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}".format(atoms[i], V[i, 0], V[i, 1], V[i, 2])
        print(line)

def get_coordinates(filename, fmt, ignore_hydrogens):
    """
    Get coordinates from filename.

    """
    if fmt == "xyz":
        return get_coordinates_xyz(filename, ignore_hydrogens)
    elif fmt == "pdb":
        return get_coordinates_pdb(filename, ignore_hydrogens)
    exit("Could not recognize file format: {:s}".format(fmt))


def get_coordinates_pdb(filename, ignore_hydrogens):
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.

    """
    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the above
    # column indices as a fallback.
    x_column = None
    V = []
    # Same with atoms and atom naming. The most robust way to do this is probably
    # to assume that the atomtype is given in column 3.
    atoms = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
                tokens = line.split()
                # Try to get the atomtype
                try:
                    atom = tokens[2][0]
                    if ignore_hydrogens and atom == "H":
                        continue
                    elif atom in ["H", "C", "N", "O", "S", "P"]:
                        atoms.append(atom)
                except:
                        exit("Error parsing atomtype for the following line: \n%s" % line)

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i+1] and "." in tokens[i+2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("Error parsing coordinates for the following line: \n%s" % line)
                # Try to read the coordinates
                try:
                    V.append(np.asarray(tokens[x_column:x_column+3],dtype=float))
                except:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        V.append(np.asarray([x,y,z],dtype=float))
                    except:
                        exit("Error parsing input for the following line: \n%s" % line)


    V = np.asarray(V)
    return(atoms, V)


def get_coordinates_xyz(filename, ignore_hydrogens):
    """
    Get coordinates from a filename.xyz and return a vectorset with all the
    coordinates.

    This function has been written to parse XYZ files, but can easily be
    written to parse others.

    """
    f = open(filename, 'r')
    V = []
    atoms = []
    n_atoms = 0
    lines_read = 0

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(f.next())
    except ValueError:
        exit("Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    f.next()

    # Use the number of atoms to not read beyond the end of a file
    for line in f:

        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
        numbers = [float(number) for number in numbers]

        # ignore hydrogens
        if ignore_hydrogens and atom.lower() == "h":
            continue

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) == 3:
            V.append(np.array(numbers))
            atoms.append(atom)
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

        lines_read += 1

    f.close()
    V = np.array(V)
    return(atoms, V)


if __name__ == "__main__":

    import argparse
    import sys

    description = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format. The order of the atoms *must* be the same for both structures.

citation:

 - Kabsch algorithm:
   Kabsch W., 1976, A solution for the best rotation to relate two sets of
   vectors, Acta Crystallographica, A32:922-923, doi:10.1107/S0567739476001873

 - Quaternion algorithm:
   Michael W. Walker and Lejun Shao and Richard A. Volz, 1991, Estimating 3-D
   location parameters using dual number quaternions, CVGIP: Image Understanding,
   54:358-367, doi: 10.1016/1049-9660(91)90036-o

 - Implementation:
   Calculate RMSD for two XYZ structures, GitHub,
   http://github.com/charnley/rmsd

"""

    epilog = """
The script will return three RMSD values:
Normal: The RMSD calculated the straight-forward way.
Kabsch: RMSD after coordinates are translated and rotated using Kabsch.
Quater: RMSD after coordinates are translated and rotated using quaternions.
"""

    parser = argparse.ArgumentParser(
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('structure_a', metavar='structure_a.xyz', type=str)
    parser.add_argument('structure_b', metavar='structure_b.xyz', type=str)
    parser.add_argument('-o', '--output', action='store_true', help='print out structure A, centered and rotated unto structure B\'s coordinates in XYZ format')
    parser.add_argument('-n', '--no-hydrogen', action='store_true', help='ignore hydrogens when calculating RMSD')
    parser.add_argument('-f', '--format', action='store', help='Format of input files. Supports xyz or pdb.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # As default, load the extension as format
    if args.format == None:
        args.format = args.structure_a.split('.')[-1]

    atomsP, P = get_coordinates(args.structure_a, args.format, args.no_hydrogen)
    atomsQ, Q = get_coordinates(args.structure_b, args.format, args.no_hydrogen)

    # Calculate 'dumb' RMSD
    normal_rmsd = rmsd(P, Q)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    if args.output:
        V = rotate(P, Q)
        V += Qc
        write_coordinates(atomsP, V)
        quit()

    print("Normal RMSD:", normal_rmsd)
    print("Kabsch RMSD:", kabsch_rmsd(P, Q))
    print("Quater RMSD:", quaternion_rmsd(P, Q))
