#!/usr/bin/python
from __future__ import division
"""

"""
import numpy as np
from py_symmetric_matrix import *
import re


# pattern to select gene name from minitab file
RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")
# manual renaming
GENE_RENAMES = {
  'MO25 ALPHA': 'MO25',
  'STRAD BETA': 'STRAB',
}
def clean(s):
  # leave spaces; handle them specifically
  s = re.sub('[^a-zA-Z0-9 ]', '', s.upper())
  if s in GENE_RENAMES:
    s = GENE_RENAMES[s]
  return s


class EnrichedSet(object):
  """Confirmed protein relationships.

  Attributes:
    pairs: [(str,str)] list of sorted ordered pairs of cleaned gene names
    genes: set(str): set of all gene names found in enrichment set
    mismatched: set(str): set of strings unmatched to gene names
    pairs_hash: set(str): set of "%s,%s" of `self.pairs` Use to confirm pair.
  """

  def __init__(self, filename, name_col_a=2, name_col_b=3):
    """Load a PINA minitab file as a list of genes and dependency matrix.
  
    Args:
      filename: str of filepath to PINA minitab
      name_col_1: int of column index (from zero) of gene name (a)
      name_col_b: int of column index (from zero) of gene name (b)
    """
    self.pairs = []
    self.genes = set()
    self.mismatched = set()
      
    fp = open(filename)
    fp.next()   # skip first header line
    # List of pairs and unique genes.
  
    for line_num, line in enumerate(fp):
      row = line[:-1].split('\t')
      s_a, s_b = (row[name_col_a], row[name_col_b])
      m_a, m_b = RX_GENE_NAME.match(s_a), RX_GENE_NAME.match(s_b)
      
      # skip missing symbols, mismatches, and same-pairs
      if not (m_a and m_b):
        if not m_a: self.mismatched.add(s_a)
        if not m_b: self.mismatched.add(s_b)
        continue
  
      # get cleaned variable names
      a, b = clean(m_a.group(1)), clean(m_b.group(1))
      if a == b or not a or not b:
        continue
      
      self.genes.add(a)
      self.genes.add(b)
      self.pairs.append(tuple(sorted((a,b))))
  
    # make all-pairs hash
    self.pairs_hash = set(["%s,%s" % (pair[0], pair[1]) for pair in self.pairs])


class DependencySet(object):
  """Dependecy matrices with rank order and variable lists.
  NOTE: all associated variable names must be sorted

  Attributes:
    varlist: [str] of sorted variable names corresponding to matrices.
    n: int of number of variables (size of matrix in rows)
    dependencies: {str=>(np.array(float), np.array(float)} of 
      relation=>(Value_Matrix, Rank_Matrix)
      matrices are symmetric representations; use sym_idx
  """
  PLOT_TOP_K = 200
  # XXX SUPER HACKY!
  # FOR GSE25935
  N_TOTAL_CONFIRMED = 41902435 # 9155 choose 2, 9155 intersecting genes
  N_CORRELATION_RANKS = 109655762

  def __init__(self, varlist_filename):
    self.varlist = [s.strip() for s in open(varlist_filename)]
    self.dependencies = {}
    self.n = len(self.varlist)
  
  def add(self, name, matrix_filename, rank_filename):
    """Add a dependency."""
    M = np.load(matrix_filename)
    Q = np.load(rank_filename)
    self.dependencies[name] = (M, Q)

  def get_overlap(self, enriched):
    """Return intersection of genes in data and genes in enrichment.

    Returns:
      set([str]) of shared genes
    """
    shared_set = enriched.genes & set(self.varlist)
    return shared_set

  def compare(self, name, enriched, limit=None):
    """Compare a dependency ranking with an enriched set of pairs.

    Args:
      name: str in self.dependencies
      enriched: EnrichedSet 
      limit: int of pairs to check. None=check all pairs
      plot: bool if to save top K plots of confirmed pairs.
    Returns:
      [{str:var}] of enriched dependencies in decreasing dependencies order
    """
    top_enriched = []  # list of dictionaries
    n_counted = 0
    n_confirmed = 0
    folder_enrich = -1
    
    # M=values, Q=ranks
    M, Q = self.dependencies[name]

    if limit is None:
      limit = len(M)
      
    for i in xrange(1,limit+1):
      idx = Q[-i]
      x, y = inv_sym_idx(idx, self.n)
      pair = sorted((self.varlist[x], self.varlist[y]))
      
      # if both gene name are not in the enriched list, skip
      if not (pair[0] in enriched.genes and pair[1] in enriched.genes):
	continue
      n_counted += 1

      assert pair[0] < pair[1]
      if "%s,%s" % (pair[0], pair[1]) in enriched.pairs_hash:
	n_confirmed += 1
	# build dict of confirmation record
	r = {
	  'name1': pair[0],
	  'name2': pair[1],
	  'counted': n_counted,
	  'confirmed': n_confirmed,
	}
	# add all dependency scores
	for name, R in self.dependencies.items():
	  r[name] = R[0][idx]
	# add nonlinearity if possible
	if 'mic' in self.dependencies and 'pcc2' in self.dependencies:
	  r['mic-pcc2'] = r['mic'][idx] - r['pcc2'][idx]

	# update enrichment
	tmp = (n_confirmed / self.N_TOTAL_CONFIRMED) / \
	  (i / self.N_CORRELATION_RANKS)
	if tmp > folder_enrich:
	  folder_enrich = tmp
	r['curr_folder'] = tmp
	r['max_folder'] = folder_enrich
	
	# append record to list of matches
	top_enriched.append(r)
        
    return top_enriched

    
