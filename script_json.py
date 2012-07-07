#!/usr/bin/python
"""Compute top K enrichments and print results from numpy similarity matrices.
Use JSON configuration file; output JSON.

SAMPLE USE

python script.py enriched_file=PINA_Homo_sapiens-20110628_minitab.txt tabdata_file=GSE25935.GPL4133.eQTL.tab limit=1000 pcc2=gse25935-matlab-correlations-all.npy,gse25935-matlab-correlations-all-argsorted.npy


python script.py enriched_file=PINA_Homo_sapiens-20110628_minitab.txt tabdata_file=GSE25935.GPL4133.eQTL.tab limit=1000 mic=gse25935.mic.matrix.npy,gse25935.mic.matrix.argsorted.npy,gse25935.varlist.txt

python script.py enriched_file=PINA_Homo_sapiens-20110628_minitab.txt tabdata_file=GSE25935.GPL4133.eQTL.tab limit=1000 pcc2=GSE25935.PCC2.matrix.npy,GSE25935.PCC2.matrix.sorted.npy,gse25935.varlist.txt mic=gse25935.mic.matrix.npy,gse25935.mic.matrix.argsorted.npy,gse25935.varlist.txt

python script.py enriched_file=PINA_Homo_sapiens-20110628_minitab.txt limit=1000 tabdata_file= spearman=/Users/qq/Desktop/yang_ppi_enrichment/yang_gse2034_spearman_squareform.npy,/Users/qq/Desktop/yang_ppi_enrichment/yang_gse2034_spearman_squareform_sorted.npy,/Users/qq/Desktop/yang_ppi_enrichment/GeneCoexpWang_genelist.txt
"""
from __init__ import *
import sys
import json
import os


def print_enrichments(enriched_file=None, tabdata_file=None, limit=None, np_matrices_json=None):
  """Compute top K enriched pairs from numpy matrices and known pair list.

  Args:
    enriched_file: str of path to PINA enriched filename
    tabdata_file: str of path to tab-delimited data 
      VARIABLE NAMES MUST BE IN SAME ORDER AS IN NP_MATRICES!
    limit: int of number of pairs to consider
    np_matrices_json: str of filepath to dependency matrix definitions
  """
  if type(limit) == str:
    limit = int(limit)
  assert limit > 0
  assert enriched_file and os.path.exists(enriched_file), enriched_file
  assert tabdata_file and os.path.exists(tabdata_file), tabdata_file
  assert np_matrices_json and os.path.exists(np_matrices_json), np_matrices_json

  # Load enriched set of pairs.
  enriched = EnrichedSet(enriched_file)

  # Load variable list in order from tab data row variable names.
  depends = DependencySet(enriched, varlist_filename=tabdata_file)

  # Add each dependency matrix.
  D = json.load(open(np_matrices_json))

  #R = [{"measure": "name", "size": size, [other stuff], "list": []},]
  All_R = []

  for dep_name, d in D['dependencies'].items():
    depends.add(dep_name, d['values_file'], d['argsorted_file'])

  # For each dependency matrix, compare rankings
  for dep_name in D['dependencies']:

    # r is [{str:var}] of enriched dependencies in decreasing dependencies order
    #  where each element is a dict including var names, enrichment, etc.
    r = depends.compare(dep_name, limit=limit)
    # print results
    R = {}
    R['measure'] = dep_name
    R['size'] = len(r)
    R['limit'] = limit
    R['meta'] = depends.meta[dep_name]
    R['list'] = r
    All_R.append(R)

  # OUTPUT AT END OF LOOP
  print json.dumps(All_R)
  

if __name__ == "__main__":
  kwds = dict([ s.split('=') for s in sys.argv[1:]])
  print_enrichments(**kwds)
  
