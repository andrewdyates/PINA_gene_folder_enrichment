#!/usr/bin/python
"""Implementation of Dr. Yang Xiang's PPI Folder Enrichment pseudocode.
WARNING: HACKY

Writen to resemble pseudocode as much as possible in Python.

Download Enrichment list:
  http://cbg.garvan.unsw.edu.au/pina/download/Homo%20sapiens-20110628.txt

Currently, compares on Gene Symbol (name).
to help limit comparison problems, strip whitespace and punctuation and convert to all caps
NOTE: Should use first two names. I need a good gene name tool!

NOTE: some rows have the exact? same gene listed as dependent with itself. ignore these.


def main(relations_fname, varlist_fname, matrix_filename):
USE:
  python scratch.py PINA_Homo_sapiens-20110628_minitab.txt gse25935.varlist.txt argsorted.matrix.npy
"""
from __future__ import division
import sys
from py_symmetric_matrix import *
import re
import numpy as np
import gzip
import zipfile


RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")

GENE_RENAMES = {
  'MO25 ALPHA': 'MO25',
  'STRAD BETA': 'STRAB',
}

def clean(s):
  # leave spaces; handle them specifically
  return re.sub('[^a-zA-Z0-9 ]', '', s.upper())

def open_with_zip(filepath):
  """Open a file with .gzip, .zip, or no compression.

  Returns:
    [*str] read filepointer of uncompressed file stream
  """
  ext = filepath.lower().rpartition('.')[-1]
  if ext == "gzip":
    return gzip.open(filepath, "rb")
  elif ext == "zip":
    z = zipfile.ZipFile(filepath, "r")
    # Return first file in zip.
    return z.open(z.namelist()[0], "r")
  else:
    return open(filepath, "rb")



def load_pina_minitab(fp, name_col_a=2, name_col_b=3):
  """Load a PINA minitab file as a list of genes and dependency matrix.

  Args:
    fp: [*str] of open file pointer to PINA file
    name_col_1: int of column index (from zero) of gene name (a)
    name_col_b: int of column index (from zero) of gene name (b)
  Returns:
    ([str], NamedSymmetricMatrix) of variable list and dependency matrix
    i.e., (GList_PPI, M_PPI) in the parlance of Dr. Yang's pseudocode.
  """
  fp.next()   # skip first header line
  # List of pairs and unique genes.
  pairs = []
  GList_PPI = set()
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
    pairs.append((a,b))
    GList_PPI.add(a)
    GList_PPI.add(b)

  # Generate Dependency Matrix from list of pairs
  M_PPI = NamedSymmetricMatrix(var_list=GList_PPI,
    n=len(GList_PPI), store_diagonal=False, dtype=np.bool)
  for a, b in pairs:
    M_PPI.set(a, b, 1)

  return (GList_PPI, M_PPI)


def load_varlist(varlist_fname):
  # load variable list
  fp = open(varlist_fname)
  varlist = [s.strip() for s in fp]
  fp.close()
  return varlist


def load_mine_results(matrix_filename, varlist, order='increasing'):
  """Load MINE rankings from sorted numpy array.

  Args:
    matrix_filename: str of filename of sorted npy matrix
    varlist: [str] of variable names corresponding to matrix 
    order: str in ('increasing', 'decreasing') of order of 
  Returns:
    [(str,str)] of ranked pairs in descending order
  """
  M = np.load(matrix_filename)
  CorrelationRank = []
  n = len(varlist)
  if order == 'increasing':
    for i in range(1,len(M)+1):
      x,y = inv_sym_idx(M[-i], n)
      try:
        CorrelationRank.append((varlist[x], varlist[y]))
      except IndexError:
        print "Unknown error. i=%d, x=%d, y=%d, n=%d." % (i,x,y,n)
      
  return CorrelationRank


def yield_top_pairs(M, varlist, order='increasing'):
  """YIELD variable name pairs in rank order from matrix

  Args:
    matrix_filename: str of filename of sorted npy matrix
    varlist: [str] of variable names corresponding to matrix 
    order: str in ('increasing', 'decreasing') of order of 
  Returns:
    [(str,str)] of ranked pairs in descending order
  """
  n = len(varlist)
  if order == 'increasing':
    for i in range(1,len(M)+1):
      x,y = inv_sym_idx(M[-i], n)
      yield (varlist[x], varlist[y])
  


def get_total(GList_PPI, CorrelationRank):
  total_confirmed = 0
  total = 0
  for a, b in CorrelationRank:
    total += 1
    # are both genes in the pair in the list of genes considered?
    if a in GList_PPI and b in GList_PPI:
      total_confirmed += 1
    if total % 100000 == 0:
      print "counting up to %d..." % total
      
  if total_confirmed == 0:
    raise Exception, "PPIFolderEnrichment fails due to insufficient information."
  return total_confirmed, total

      
def PPIFolderEnrichment(GList_PPI, M_PPI, CorrelationRank, total_confirmed, corr_len):
  """Return gene enrichment score given dependency matrix and a ranking.

  Args:
    GList_PPI: set([str]) of all variable names considered
    M_PPI: NamedSymmetricMatrix boolean matrix of given gene/gene interactions.
    CorrelationRank: [*(str, str)] iter of decreasing rank relations.
  Returns:
    float of folder_enrich score
  """
  
  folder_enrich = -1
  confirmed = 0
  enrich_seq = []
  total = 0
  
  for i, pair in enumerate(CorrelationRank):
    a, b = pair
    total += 1
    if a in GList_PPI and b in GList_PPI:
      if M_PPI.get(a, b):
        confirmed += 1
        tmp = (confirmed / total_confirmed) / (i+1 / corr_len)
        enrich_seq.append(tmp)
        if tmp > folder_enrich:
          folder_enrich = tmp
    if total % 100000 == 0:
      print "counting up to %d..." % total
  
  return folder_enrich, enrich_seq




def main(relations_fname, varlist_fname, matrix_filename, total=109675455):
  """Read files of known relations anddependency rankings; print FolderEnrich.

  Args:
    relations_filepath: str of filepath
    dependency_rank_filepath: str of filepath
  """

  varlist = load_varlist(varlist_fname)
  
  fp = open(relations_fname)
  GList_PPI, M_PPI = load_pina_minitab(fp)
  fp.close()
  
  #CorrelationRank = yield_top_pairs
  print "counting total..."
  M = np.load(matrix_filename)
  CorrIter = yield_top_pairs(M, varlist)
  total_confirmed,q = get_total(GList_PPI, CorrIter)
  print "total confirmed: %d; total: %d (%d)" % (total_confirmed, total, (total-total_confirmed))
  
  # write CorrelationRank to file!
  CorrIter = yield_top_pairs(M, varlist)
  score, enrich_seq = PPIFolderEnrichment(GList_PPI, M_PPI, CorrIter, total_confirmed, total)

  print "PPI Folder Enrichment Score: %f" % score
  fp = open('enrich_seq.txt', 'w')
  for x in enrich_seq:
    fp.write("%.8f\n" % x)
  fp.close()

  
if __name__ == "__main__":
  main(sys.argv[1], sys.argv[2], sys.argv[3])
