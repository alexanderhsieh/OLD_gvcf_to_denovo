'''
Script to take as input sample_map and sample_id and output the file path corresponding to sample_id
'''

import sys
from optparse import OptionParser


####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-m', '--smap', dest='sample_map',help='sample map (picard)')
parser.add_option('-s', '--sid', dest='sample_id',help='sample id')
(options, args) = parser.parse_args()

## check all arguments present
if (options.sample_id == None or options.sample_map == None):
	print '\n' + '## ERROR: missing arguments' + '\n'
	parser.print_help()
	print '\n'
	sys.exit()


sample_id = options.sample_id
sample_map = options.sample_map

## parse sample_map
pathd = {} # {id : path}
with open(sample_map, 'r') as smapf:
  for line in smapf:
    tmp = line.strip().split('\t')
    id = tmp[0]
    path = tmp[1]
    pathd[id] = path


## return path
print(pathd[sample_id])