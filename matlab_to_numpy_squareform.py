#!/usr/bin/python
"""Convert Matrix of .mat to squareform .npy format."""
import sys
from numpy import *
from scipy.spatial.distance import squareform
import scipy.io as sio


fname = sys.argv[1]
mname = sys.argv[2]
M = sio.loadmat(fname)
Q = squareform(M[mname], checks=False)
R = argsort(Q)

save(fname+".npy", Q)
save(fname+"-sorted.npy", R)

