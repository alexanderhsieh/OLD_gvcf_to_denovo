'''
Script to take as input sample_map and sample_id and output the file path corresponding to sample_id
Usage: python parse_sample_map.py <sample_map> <sample_id>
'''

import sys


####################################################################################################
## handle arguments
####################################################################################################


## parse sample_map 
pathd = {} # {id : path}
with open(sys.argv[1], 'r') as smapf:
  for line in smapf:
    tmp = line.strip().split('\t')
    id = tmp[0]
    path = tmp[1]
    pathd[id] = path


## return path
print(pathd[sys.argv[2]])