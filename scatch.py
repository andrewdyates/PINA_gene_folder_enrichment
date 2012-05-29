#!/usr/bin/python
"""Implementation of Dr. Yang Xiang's PPI Folder Enrichment pseudocode.

Writen to resemble pseudocode as much as possible in Python.

Download Enrichment list:
  http://cbg.garvan.unsw.edu.au/pina/download/Homo%20sapiens-20110628.txt
"""
from __future__ import division
import sys
from py_symmetric_matrix import *

def PPIFolderEnrichment(M_PPI, CorrelationRank):
  """Return gene enrichment score given dependency matrix and a ranking.

  Args:
    M_PPI: SymmetricMatrix boolean matrix of given gene/gene interactions.
    CorrelationRank: [{'a': str, 'b': str}] list of decreasing rank relations.
  Returns:
    float of folder_enrich score
  """
  total_confirmed = 0
  for i in range(len(CorrelationRank)):
    # are both genes in the pair in the list of genes considered?
    if CorrelationRank[i]['a'] in GList_PPI and CorrelationRank[i]['b'] in GList_PPI:
      total_confirmed += 1
      
  if total_confirmed == 0:
    raise Exception, "PPIFolderEnrichment fails due to insufficient information."
  
  folder_enrich = -1
  confirmed = 0
  
  for i in range(len(CorrelationRank)):
    if CorrelationRank[i]['a'] in GList_PPI and CorrelationRank[i]['b'] in GList_PPI:
      if M_PPI.get(CorrelationRank[i]['a'], CorrelationRank[i]['b']):
        confirmed += 1
        tmp = (confirmed / total_confirmed) / (i+1 / len(CorrelationRank))
        if tmp > folder_enrich:
  	folder_enrich = tmp
  
  return folder_enrich


def main(relations_filepath, dependency_rank_filepath):
  """Read files of known relations anddependency rankings; print FolderEnrich.

  Args:
    relations_filepath: str of filepath
    dependency_rank_filepath: str of filepath
  """
  # load M_PPI
  fp = open(relations_filepath)
  # TODO: Load in file
  M_PPI = []
  fp.close()

  # load CorrelationRank
  fp = open(relations_filepath)
  # TODO: Load in file
  CorrelationRank = []
  fp.close()
  
  score = PPIFolderEnrichment(M_PPI, CorrelationRank)
  print "PPI Folder Enrichment Score: %f" % score

  
if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2])
