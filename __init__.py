#!/usr/bin/python
from __future__ import division
"""Perform gene enrichment.

TODO: 
  - convert all files to "array objects"

USEFUL COMMANDS:

Truncate symmetric matrix to squareform array.
from scipy.spatial.distance import squareform
>>> Q = squareform(M, checks=False)
http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html

Load matlab data
M = sio.loadmat('gse25935-matlab-correlations-all.mat') #WARNING! REMOVE COLS FROM OUTLIER PATIENT IN CORRESPONDING DATA!
NOTE: this may load M as a dictionary if .mat contains multiple files.
"""
from scipy.spatial.distance import squareform
import numpy as np
from py_symmetric_matrix import *
import re


# pattern to select gene name from minitab file
RX_GENE_NAME = re.compile("uniprotkb:([^)]*)\(gene name\)")
# manual renaming
GENE_RENAMES = {
  'MO25 ALPHA': 'MO25A',
  'MO25 BETA': 'MO25B',
  'STRAD BETA': 'STRADB',
  'STRAD ALPHA': 'STRADA',
}
def clean(s):
  """Standardize gene names to be all caps, alphanumeric."""
  s = re.sub('[^a-zA-Z0-9 ]', '', s.upper())
  if s in GENE_RENAMES:
    s = GENE_RENAMES[s]
  # Loudly complain about spaces so that they can be handled specially
  assert ' ' not in s, s
  return s

def matlab_to_squareform(M):
  """Truncate symmetric matrix (as numpy.ndarray) to numpy.squareform."""
  return squareform(M, checks=False)

def tab_to_varlist(tab_fname):
  """Transform a .tab gene expression value matrix into a list of variables.

  Args:
    tab_fname: str of path to .tab MINE.jar ready matrix file
  Returns:
    [str] of cleaned variable names IN FILE ORDER
  """ 
  fp = open(tab_fname)
  varlist = []
  for line in fp:
    if line[0] == "#": continue
    var,c,cc = line.partition('\t')
    varlist.append(clean(var))
  return varlist


class EnrichedSet(object):
  """Confirmed protein relationships.

  NOTE: all variable names are "cleaned"!

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

  def exists(self, x, y):
    """Return if x, y variable pair is in enriched set.

    Args:
      x: str of variable name
      y: str of variable name
    Returns: 
      bool if (x,y) variable pair is in enriched set.
    """
    assert x in self.genes
    assert y in self.genes
    assert not x == y
    if x > y:
      x, y = y, x
    return "%s,%s" % (x, y) in self.pairs_hash


class DependencySet(object):
  """Dependecy matrices with rank order and variable lists.
  
  REQUIRED: row and column order of numpy dependecy matrices must
    be in the same order as the row order of .tab variable list.

  Attributes:
    varlist: [str] of sorted variable names corresponding to matrices.
    n: int of number of variables (size of matrix in rows)
    dependencies: {str=>(np.array(float), np.array(float)} of 
      relation=>(Value_Matrix, Rank_Matrix)
      matrices are symmetric representations; use sym_idx
    enriched: corresponding `EnrichedSet` object
  """
  PLOT_TOP_K = 200

  def __init__(self, enriched, varlist_filename=None, varlist=None):
    """Initialize dependency set with either filename or varlist.
    Filename can either be a .tab matrix or [str] list but not both.

    Args:
      enriched: obj of `EnrichedSet`
      varlist_filename: str of path to ordered list of variables
      varlist: [str] of variable names
    """
    assert bool(varlist_filename) != bool(varlist)
    if varlist_filename:
      self.varlist = tab_to_varlist(varlist_filename)
    else:
      self.varlist = varlist
    self.enriched = enriched
    
    self.dependencies = {}
    self.n = len(self.varlist)
  
  def add(self, name, similarity_fname, rank_fname):
    """Add a dependency matrix.

    Args:
      name: str of name of dependency type. 
      similarity_fname: str of path to similarity matrix as numpy array
        array is length (n choose 2), use sym_idx functions to index
      rank_fname: str of path to ranks of similarity matrix as numpy array
        array is length (n choose 2), use sym_idx functions to index
    """
    M = np.load(similarity_fname)
    Q = np.load(rank_fname)
    self.dependencies[name] = (M, Q)

  def get_overlap(self):
    """Return intersection of genes in data and genes in enrichment.

    Returns:
      set([str]) of shared genes
    """
    shared_set = self.enriched.genes & set(self.varlist)
    return shared_set

  def compare(self, name, limit=None):
    """Compare a dependency ranking with an enriched set of pairs.

    Args:
      name: str in self.dependencies
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

    Named_M = NamedSymmetricMatrix(var_list=self.varlist, matrix=M)

    # Do not count values ranked exactly zero for MIC
    if name == 'MIC':
      n_correlation_rank = len(filter(None, M))
    else:
      n_correlation_rank = len(M)
      
    # All enriched pairs in overlapping set are confirmed
    n_total_confirmed = 0
    overlapping_genes_set = self.get_overlap()
    for x,y in self.enriched.pairs:
      if x in overlapping_genes_set and y in overlapping_genes_set:
        # Only confirm MIC pairs with values greater than 0
        if name == 'MIC' and Named_M.get(x,y) != 0:
          n_total_confirmed += 1
        else:
          n_total_confirmed += 1

    print n_total_confirmed
    print len(overlapping_genes_set)
    print len(self.enriched.pairs)

    # For the top `limit` pairs, compute enrichment
    for i in xrange(1,len(M)+1):
      # Break if n_counted is equal to limit. If no limit, count all pairs
      if limit is not None and n_counted >= limit:
        break
      
      idx = Q[-i]
      x, y = inv_sym_idx(idx, self.n)
      pair = sorted((self.varlist[x], self.varlist[y]))
      
      # if both gene name are not in the enriched list, skip
      if not (pair[0] in self.enriched.genes and pair[1] in self.enriched.genes):
	continue
      n_counted += 1

      # Count pairs in enrichment set.
      assert pair[0] < pair[1]
      if "%s,%s" % (pair[0], pair[1]) in self.enriched.pairs_hash:
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
	  r['mic-pcc2'] = r['mic'] - r['pcc2']

	# update enrichment
	tmp = (n_confirmed / n_total_confirmed) / (i / n_correlation_rank)
	if tmp > folder_enrich:
	  folder_enrich = tmp
	r['curr_folder'] = tmp
	r['max_folder'] = folder_enrich
	
	# append record to list of matches
	top_enriched.append(r)
        
    return top_enriched


def to_floats(s):
  try:
    x = float(s)
  except ValueError:
    x = None
  return x

class TabData(object):
  """GSE dataset .tab file representation with an index to rows.

  Attributes:
    col_titles = [str] of GSM sample name column titles in col order
    rows = {str=>[float|None] of variable name to values in col order
      var name cleaned using `clean`
  """

  def __init__(self, filename):
    """Load .tab matrix file into self.

    Args:
      filename: str of filepath to .tab file.
      clean: bool if to 'normalize' gene name
    """
    fp = open(filename)

    self.rows = {}
    
    for line in fp:
      line = line.strip('\n')
      if not line: continue

      row = line.split('\t')
      # handle headers
      if line[0] == '#':
        # get column headers
        if "#GENE_SYMBOL" == line[0:len('#GENE_SYMBOL')]:
          self.col_titles = row[1:]
        continue
      # add line
      self.rows[clean(row[0])] = map(to_floats, row[1:])

  def output(self, row_id):
    row_id = clean(row_id)
    try:
      row = self.rows[row_id]
    except KeyError:
      return None
    else:
      return "\t".join([row_id] + row_to_strs(row))
    
    
def row_to_strs(row):
  line = []
  for x in row:
    if x is None:
      line.append("")
    else:
      line.append("%.6f"%x)
  return line
