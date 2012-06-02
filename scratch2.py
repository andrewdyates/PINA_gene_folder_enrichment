#!/usr/bin/python
"""Copied from test.py
Interactive inspection of ranks.

$ python matrix_filename pina_filename varlist_filename presorted_filename


$ python test.py gse25935.mic.matrix.npy PINA_Homo_sapiens-20110628_minitab.txt gse25935.varlist.txt gse25935.mic.matrix.argsorted.npy

$ python test.py GSE25935.PCC2.matrix.npy PINA_Homo_sapiens-20110628_minitab.txt gse25935.varlist.txt GSE25935.PCC2.matrix.sorted.npy
"""
import sys
from py_symmetric_matrix import *
import numpy as np
import re

RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")

GENE_RENAMES = {
  'MO25 ALPHA': 'MO25',
  'STRAD BETA': 'STRAB',
}


def compare_ranks(matrix_filename, pina_filename, varlist_filename, presorted_filename=None):

  # Load computed matrices and ranks
  M = np.load(matrix_filename)
  # sort by increasing order
  if not presorted_filename:
    Q = M.argsort()
  else:
    Q = np.load(presorted_filename)

  # Load variable list associate with matrix
  varlist = [s.strip() for s in open(varlist_filename)]
  

  # Load PINA known binary relationships
  pina_list, pina_gene_set = load_pina_minitab(pina_filename)
  # for membership testing
  pina_hash = set(["%s,%s" % (a,b) for a,b in pina_list])

  # Get intersection of gene names
  shared_set = pina_gene_set & set(varlist)
  print "Num genes in PINA: %d, num genes in matrix %d, overlapping: %d" % \
    (len(pina_gene_set), len(varlist), len(shared_set))
  
  top_names = []
  top_scores = []
  in_pina = []

  for i in xrange(1,2000000):
    x, y = inv_sym_idx(Q[-i], len(varlist))
    pair = sorted((varlist[x], varlist[y]))
    v = M[Q[-i]]
    top_names.append(pair)
    top_scores.append(v)
    if "%s,%s" % (pair[0], pair[1]) in pina_hash:
      in_pina.append(True)
      #print "%s, %s: %.3f" % (pair[0], pair[1], v)
    else:
      in_pina.append(False)
    
    
  print "size PINA set", len(pina_list)
  print "PINA gene list", len(pina_gene_set)
  print "size of varlist", len(varlist)

  print "out of %d, %d in list." % (len(in_pina), sum([1 for i in filter(None, in_pina)]))

  

def clean(s):
  # leave spaces; handle them specifically
  return re.sub('[^a-zA-Z0-9 ]', '', s.upper())



def load_pina_minitab(filename, name_col_a=2, name_col_b=3):
  """Load a PINA minitab file as a list of genes and dependency matrix.

  Args:
    fp: [*str] of open file pointer to PINA file
    name_col_1: int of column index (from zero) of gene name (a)
    name_col_b: int of column index (from zero) of gene name (b)
  """
  fp = open(filename)
  fp.next()   # skip first header line
  # List of pairs and unique genes.
  pairs = []
  genes = set()
  for line_num, line in enumerate(fp):
    row = line[:-1].split('\t')
    s_a, s_b = (row[name_col_a], row[name_col_b])
    m_a, m_b = RX_GENE_NAME.match(s_a), RX_GENE_NAME.match(s_b)
    if not (m_a and m_b):
      print "Gene names [%s, %s] do not match pattern." % (s_a, s_b)
      continue
    else:
      a, b = clean(m_a.group(1)), clean(m_b.group(1))
    # skip genes missing symbols
    if not (a and b):
      continue
    # manual remap of gene names
    if a in GENE_RENAMES:
      a = GENE_RENAMES[a]
    if b in GENE_RENAMES:
      b = GENE_RENAMES[b]
    if a == b:
      continue
    genes.add(a)
    genes.add(b)
    pairs.append(sorted((a,b)))

  return pairs, genes


  
if __name__ == "__main__":
  compare_ranks(*sys.argv[1:])
