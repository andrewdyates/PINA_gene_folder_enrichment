import sys
fp = open(sys.argv[1])
for line in fp:
  row = line.strip('\n').split('\t')
  # remove
  # GSM637026 (46_rep1, col #84) w/o variable name
  # GSM637056 (46_rep2, col #114) 
  print '\t'.join(row[0:84] + row[85:114] + row[115:])
  
