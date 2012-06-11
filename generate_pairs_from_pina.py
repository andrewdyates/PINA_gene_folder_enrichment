#!/usr/bin/python
"""Given a PINA and .tab matrix, output known pairs as adjacent vectors.

USE:

$ python generate_pairs_from_pina.py data_file=GSE25935.GPL4133.eQTL.nooutliers.tab enriched_file=PINA_Homo_sapiens-20110628_minitab.txt
"""
import sys
from __init__ import *

def main(data_file, enriched_file):
  
  eset = EnrichedSet(filename=enriched_file)
  data = TabData(data_file)

  for x, y in eset.pairs:
    x_line = data.output(x)
    y_line = data.output(y)
    if x_line and y_line:
      print x_line
      print y_line

  
if __name__ == "__main__":
  main(**dict([s.split('=') for s in sys.argv[1:]]))
