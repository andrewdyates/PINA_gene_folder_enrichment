#!/usr/bin/python
"""Compute top K enrichments and print results from numpy similarity matrices.

SAMPLE USE

python script.py enriched_file=PINA_Homo_sapiens-20110628_minitab.txt tabdata_file=GSE25935.GPL4133.eQTL.tab limit=1000 pcc2=gse25935-matlab-correlations-all.npy,gse25935-matlab-correlations-all-argsorted.npy
"""
from __init__ import *
import sys


RECOGNIZED_MATRICES = set(['pcc2', 'mic', 'spearman'])

def print_enrichments(enriched_file=None, tabdata_file=None, limit=None, np_matrices=None):
  """Compute top K enriched pairs from numpy matrices and known pair list.

  Args:
    enriched_file: str of path to PINA enriched filename
    tabdata_file: str of path to tab-delimited data 
      VARIABLE NAMES MUST BE IN SAME ORDER AS IN NP_MATRICES!
    limit: int of 
    kwds: {str: (str,str)} of matrix_type => (values_filepath, orders_filepath)
      matrix_type: dependency measure in RECOGNIZED_MATRICES
      values_filepath: path to numpy squareform matrix of pairwise dependencies
        (.npy array, VARIABLE NAMES IN SAME ORDER as rows in `tabdata_file`)
      orders_filepath: path to argsorted matrix in values_filepath
        (.npy array, VARIABLE NAMES IN SAME ORDER as rows in `tabdata_file`)
  """
  if type(limit) == str:
    limit = int(limit)

  # At least one similarity matrix in the set of recognized_matrices is required
  assert enriched_file
  assert tabdata_file
  assert np_matrices
  assert len(np_matrices) >= 1
  assert len(set(np_matrices.keys()) - RECOGNIZED_MATRICES) == 0
  # assert that np_matrices values all ordered pairs
  assert all([len(q) == 2 for q in np_matrices.values()])

  # Load enriched set of pairs.
  enriched = EnrichedSet(enriched_file)

  # Load variable list in order from tab data row variable names.
  depends = DependencySet(enriched, tabdata_file)
  # Add each dependency matrix.
  for k, v in np_matrices.items():
    depends.add(k, v[0], v[1])

  # For each dependency matrix, compare rankings
  results = []
  for k in np_matrices:
    # r is [{str:var}] of enriched dependencies in decreasing dependencies order
    #  where each element is a dict including var names, enrichment, etc.
    r = depends.compare(k, limit=limit)
    results.append(r)

  # TODO XXX HANDLE RESULTS
  print results

  # print "merged pcc", len(r_merged)
  # print "merged mic", len(r_merged_mic)
  # print "all", len(r_all)
  # print
  # print r_merged
  # print
  # print r_merged_mic
  # print
  # print r_all
  # print 
  # all_pairs = set(["%s,%s" % (d['name1'], d['name2']) for d in r_all])
  # merged_pairs = set(["%s,%s" % (d['name1'], d['name2']) for d in r_merged])
  # merged_pairs_mic = set(["%s,%s" % (d['name1'], d['name2']) for d in r_merged_mic])
  # print 'all_pairs', all_pairs
  # print 'merged_pairs', merged_pairs
  # print 'merged_pairs_mic', merged_pairs_mic
  
  

if __name__ == "__main__":
  kwds = dict([ s.split('=') for s in sys.argv[1:]])
  # hack for pairs of dependencies
  np_matrices = {}
  for k, v in kwds.items():
    if k in RECOGNIZED_MATRICES:
      np_matrices[k] = v.split(',')
  for k in RECOGNIZED_MATRICES:
    if k in kwds:
      del kwds[k]
  print_enrichments(np_matrices=np_matrices, **kwds)
  
