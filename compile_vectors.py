#!/usr/bin/python
"""Compile PCC2, Spearman2, dCOR, and MIC into max vectors.

SAMPLE USE:
python compile_vectors.py pcc2=pcc2.npy spearman2=spearman2.npy dcor2=
"""
import numpy as np
import sys
import os
from itertools import *


def main(pcc2, spearman2, dcor2, mic, outdir=""):
  # for each combination of these four measures (with at least 2)
  assert os.path.exists(outdir)
  assert os.path.isdir(outdir)
  
  vectors = {
    'pcc2': np.load(pcc2),
    'spearman2': np.load(spearman2),
    'dcor2': np.load(dcor2),
    'mic': np.load(mic),
  }
  n = np.size(vectors['pcc2'])
  # all vectors are same size
  assert all([n == np.size(q) for q in vectors.values()])
  
  for keys in powerset(vectors.keys(), 2):
    name = "+".join(keys)
    print "Max of %s..." % name
    v = multimax([vectors[k] for k in keys])
    path = os.path.join(outdir, name+".np")
    print "Saving %s..." % path
    np.save(path, v)

def powerset(iterable, min_r=1):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(min_r,len(s)+1))

def multimax(vs):
  """multimax(x,y,z...) ---> element-wise max of x,y,z..."""
  v = vs[0]
  for i in range(1,len(vs)):
    v = np.maximum(v, v[i])
  return v


if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))  
