#!/usr/bin/python
"""EXAMPLE USE:

$ python convert_mine_to_matrix.py ~/Dropbox/biostat/eqtl_data/GSE25935/gse25935.1.out GSE25935 write_varlist= rel_type="PCC2"
"""
import sys
import re
import numpy as np
from py_symmetric_matrix import *

def clean(s):
  # leave spaces; handle them specifically
  return re.sub('[^a-zA-Z0-9 ]', '', s.upper())

def main(mine_filename, study_name, **kwds):
  fp = open(mine_filename)
  save_mine_to_matrix(fp, name=study_name, **kwds)
  fp.close()
  

def save_mine_to_matrix(fp, name="MINE", write_varlist=True, write_matrix=True, rel_type="MIC"):
  """Load MINE rankings from MINE.jar results file. This is a bit hacky.
  HACK: save variable name list and matrix to file

  Args:
    fp: [*str] of csv MINE.jar results ordered by MIC
  Returns:
    [(str,str)] of ranked pairs
  """
  # First read pass: Get list of all variables. Save it.
  print "Reading variables..."
  if write_varlist:
    varset = set()
    for line in fp:
      row = line.split(',')
      # skip header lines
      if row[0] == "X var": continue
      varset.add(clean(row[0]))
      varset.add(clean(row[1]))
      
    print "Writing %d variables." % len(varset)
    varlist = list(sorted(varset))
    if write_varlist:
      fp_out = open("%s.varlist.txt" % name, 'w')
      for varname in varlist:
        fp_out.write("%s\n" % varname)
      fp_out.close()
  else:
    # Read from file
    varlist = []
    fp_in = open("%s.varlist.txt" % (name), 'r')
    for line in fp_in:
      varlist.append(line[:-1])

  # Second read pass: load into named pair matrix. Save it.
  print "Reading matrix... getting relation type %s" % rel_type
  M = NamedSymmetricMatrix(var_list=varlist, n=len(varlist), store_diagonal=False)
  fp.seek(0)
  for line in fp:
    row = line.split(',')
    # skip header lines
    if row[0] == "X var": continue
    if rel_type == "PCC2":
      v = float(row[7]) ** 2
    else:
      v = float(row[2])
    x, y = clean(row[0]), clean(row[1])
    M.set(x,y,v)

  print "Writing matrix..."
  if write_matrix:
    np.save("%s.%s.matrix.npy" % (name, rel_type), M._m)


if __name__ == "__main__":
  main(mine_filename=sys.argv[1], study_name=sys.argv[2], **dict([s.split('=') for s in sys.argv[3:]]))

  
