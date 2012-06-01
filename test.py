#!/usr/bin/python
"""Argsort a numpy matrix."""
import sys
from py_symmetric_matrix import *
import numpy as np

def main(matrix_filename):
  fp = open(matrix_filename)
  M=np.load(fp)
  np.save("argsorted.matrix.npy", M.argsort())
  
if __name__ == "__main__":
  main(sys.argv[1])


def invert_sym_idx(idx, n):
  # do stuff
  return x, y
