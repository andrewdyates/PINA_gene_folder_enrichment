#!/usr/bin/python
from __future__ import division
"""Interactive inspection of ranks.
Include both MIC and PCC^2 scores.
WARNING: SUPER HACKY!!!!
WARNING: SUPER HACKY!!!!
WARNING: SUPER HACKY!!!!
WARNING: SUPER HACKY!!!!

n_total all entries in matrix 109675455
n_total_entries all non-zero entries 109655762
n_total_confirmed confirmed as in PPI Enrichment 82393400

USE
$ python counts.py mic_matrix_filename pcc2_matrix_filename mic_presorted_filename pcc2_presorted_filename pina_filename varlist_filename gse25935.varlist.txt

EXAMPLE
$ python counts.py gse25935.mic.matrix.npy GSE25935.PCC2.matrix.npy gse25935.mic.matrix.argsorted.npy GSE25935.PCC2.matrix.sorted.npy PINA_Homo_sapiens-20110628_minitab.txt gse25935.varlist.txt

HACK
$ python counts.py PINA_Homo_sapiens-20110628_minitab.txt gse25935.varlist.txt gse25935.mic.matrix.npy
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

# XXX SUPER HACKY!
TOTAL_CONFIRMED = 82393400
SIZE_CORRELATION_RANK = 109655762



def compute_total_confirmed(pina_filename, varlist_filename, mic_matrix_filename):
  # exploit the fact that MIC is never 0

  # Load PINA known binary relationships
  pina_list, pina_gene_set = load_pina_minitab(pina_filename)
  # Load MIC value matrix
  M = np.load(mic_matrix_filename)
  # Load variable list associate with matrix
  varlist = [s.strip() for s in open(varlist_filename)]

  # loop through each entry in M. 
  # If the var combination is in pina_gene_set, then check to see if the MIC value exists.
  n_total = 0            # all entries in matrix
  n_total_entries = 0    # all non-zero entries
  n_total_confirmed = 0  # confirmed as in PPI Enrichment
  n = len(varlist)
  for i in xrange(n):
    for j in xrange(i+1,n):
      n_total += 1
      if n_total % 100000 == 0:
        print "%d..." % n_total
      idx = sym_idx(i,j,n)
      # check for zero entries
      if not M[idx] or M[idx] <= .000001:
        continue
      n_total_entries += 1
      x_name, y_name = varlist[i], varlist[j]
      if not (x_name in pina_gene_set and y_name in pina_gene_set):
        continue
      n_total_confirmed += 1

  print "n_total all entries in matrix", n_total            
  print "n_total_entries all non-zero entries", n_total_entries    
  print "n_total_confirmed confirmed as in PPI Enrichment", n_total_confirmed  
  






  
def compare_ranks(mic_matrix_filename, pcc2_matrix_filename, mic_presorted_filename, pcc2_presorted_filename, pina_filename, varlist_filename, data_matrix_file=None):
  """
  1) load both the MIC and PCC^2 matrix.
  2) compute ranks for both
  3) print out both rank lists with both MIC and PCC^2 and MIC-PCC^2 scores
  """

  # Load computed matrices and ranks
  M_mic = np.load(mic_matrix_filename)
  Q_mic = np.load(mic_presorted_filename)
  M_pcc2 = np.load(pcc2_matrix_filename)
  Q_pcc2 = np.load(pcc2_presorted_filename)

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

  


  N = 5000000
  # ====================
  # RANK MIC
  # ====================

  top_mic_enriched = [] #(name, rank, mic, pcc)
  n_counted = 0 # total pairs considered in both matrix and enriched set
  n_confirmed = 0
  folder_enrich = -1
  
  for i in xrange(1,N+1):
    # get gene pair
    idx = Q_mic[-i]
    x, y = inv_sym_idx(idx, len(varlist))
    pair = sorted((varlist[x], varlist[y]))
    # if a gene name is not in the enriched list, skip
    if not (pair[0] in pina_gene_set and pair[1] in pina_gene_set):
      continue
    n_counted += 1

    # if pair is in enriched set, then compute its values and return
    if "%s,%s" % (pair[0], pair[1]) in pina_hash:
      n_confirmed += 1
      mic = M_mic[idx]
      pcc2 = M_pcc2[idx]

      # update enrichment
      tmp = (n_confirmed / TOTAL_CONFIRMED) / (i / SIZE_CORRELATION_RANK)
      if tmp > folder_enrich:
        folder_enrich = tmp

      # save stats
      top_mic_enriched.append(\
        (pair, mic, pcc2, n_confirmed, n_counted, tmp, folder_enrich))

        
  print "MIC Ranked"
  print "Matrix entries considered=%d" % N
  print "counted=%d" % n_counted
  print "confirmed=%d" % n_confirmed
  print 
  print "Gene1, Gene2, MIC, PCC^2, rank_confirmed, rank_total, curr_enrich, folder_enrich"
  for pair, mic, pcc2, n_confirmed, n_counted, tmp, folder_enrich in top_mic_enriched:
    print "(%s, %s), %.3f, %.3f, %d, %d, %f, %f" % \
      (pair[0], pair[1], mic, pcc2, n_confirmed, n_counted, tmp, folder_enrich)





  # WARNING: COPY PASTA!!!!
  # ====================
  # RANK PCC^2
  # ====================
  top_mic_enriched = [] #(name, rank, mic, pcc)
  n_counted = 0 # total pairs considered in both matrix and enriched set
  n_confirmed = 0
  folder_enrich = -1
    
  for i in xrange(1,N+1):
    # get gene pair
    idx = Q_pcc2[-i]
    x, y = inv_sym_idx(idx, len(varlist))
    pair = sorted((varlist[x], varlist[y]))
    # if a gene name is not in the enriched list, skip
    if not (pair[0] in pina_gene_set and pair[1] in pina_gene_set):
      continue
    n_counted += 1

    if "%s,%s" % (pair[0], pair[1]) in pina_hash:
      n_confirmed += 1
      mic = M_mic[idx]
      pcc2 = M_pcc2[idx]

      # update enrichment
      tmp = (n_confirmed / TOTAL_CONFIRMED) / (i / SIZE_CORRELATION_RANK)
      if tmp > folder_enrich:
        folder_enrich = tmp

      # save stats
      top_mic_enriched.append(\
        (pair, mic, pcc2, n_confirmed, n_counted, tmp, folder_enrich))
    else:
      print pair
      sys.exit(1)

  print
  print
  print "PCC^2 Ranked"
  print "Matrix entries considered=%d" % N
  print "counted=%d" % n_counted
  print "confirmed=%d" % n_confirmed
  print 
  print "Gene1, Gene2, MIC, PCC^2, rank_confirmed, rank_total, curr_enrich, folder_enrich"
  for pair, mic, pcc2, n_confirmed, n_counted, tmp, folder_enrich in top_mic_enriched:
    print "(%s, %s), %.3f, %.3f, %d, %d, %f, %f" % \
      (pair[0], pair[1], mic, pcc2, n_confirmed, n_counted, tmp, folder_enrich)
  



    
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
  # XXX SUPER HACKY
  # compare_ranks(*sys.argv[1:])
  compute_total_confirmed(*sys.argv[1:])
