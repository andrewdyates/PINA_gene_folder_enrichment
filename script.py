#!/usr/bin/python
from __init__ import *
import sys

def compare_all_w_avg():
  enriched = EnrichedSet('PINA_Homo_sapiens-20110628_minitab.txt')

  depends_all = DependencySet('GSE25935.GPL4133.eQTL.tab')
  depends_all.add('pcc2', 'gse25935-matlab-correlations-all.npy', 'gse25935-matlab-correlations-all-argsorted.npy')

  depends_merged = DependencySet('gse25935.varlist.txt')
  depends_merged.add('pcc2', 'GSE25935.PCC2.matrix.npy', 'GSE25935.PCC2.matrix.sorted.npy')
  depends_merged.add('mic', 'gse25935.mic.matrix.npy', 'gse25935.mic.matrix.argsorted.npy')

  r_merged = depends_merged.compare('pcc2', enriched, limit=1000000)
  r_merged_mic = depends_merged.compare('mic', enriched, limit=1000000)
  r_all = depends_all.compare('pcc2', enriched, limit=1000000)

  print "merged pcc", len(r_merged)
  print "merged mic", len(r_merged_mic)
  print "all", len(r_all)
  print
  print r_merged
  print
  print r_merged_mic
  print
  print r_all
  print 
  all_pairs = set(["%s,%s" % (d['name1'], d['name2']) for d in r_all])
  merged_pairs = set(["%s,%s" % (d['name1'], d['name2']) for d in r_merged])
  merged_pairs_mic = set(["%s,%s" % (d['name1'], d['name2']) for d in r_merged_mic])
  print 'all_pairs', all_pairs
  print 'merged_pairs', merged_pairs
  print 'merged_pairs_mic', merged_pairs_mic
  
  
  

if __name__ == "__main__":
  compare_all_w_avg()
  
