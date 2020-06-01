'''
Script to take as input sample_map and sample_id and output the file path corresponding to sample_id
produces 3 output files:
- tmp.pb_path.txt
- tmp.fa_path.txt
- tmp.mo_path.txt
'''

import sys
from optparse import OptionParser


####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-m', '--smap', dest='sample_map',help='sample map (picard)')
parser.add_option('-p', '--ped', dest='ped', help='pedigree file')
parser.add_option('-s', '--sid', dest='sample_id',help='sample id')
(options, args) = parser.parse_args()

## check all arguments present
if (options.sample_map == None or options.ped == None or options.sample_id == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()


sample_id = options.sample_id
sample_map = options.sample_map
ped = options.ped

## parse ped file
pedd = {} # { id : {'fa': father_id, 'mo': mother_id} }
with open(ped, 'r') as pedf:
  for line in pedf:
    tmp = line.strip().split('\t')
    fid = tmp[0]
    sid = tmp[1]
    faid = tmp[2]
    moid = tmp[3]
    sex = tmp[4] # 1=male, 2=female, other=unknown
    aff = tmp[5] # 1=unaffect, 2=affect, -9=missing, 0=missing
    pedd[sid] = {'fa': faid, 'mo': moid}

## parse sample_map
pathd = {} # {id : path}
with open(sample_map, 'r') as smapf:
  for line in smapf:
    tmp = line.strip().split('\t')
    id = tmp[0]
    path = tmp[1]
    pathd[id] = path


## return path
pb_path = '.'
pb_path = pathd[sample_id]
## if parent, return default values
if pedd[sample_id]['fa'] == '0' or pedd[sample_id]['mo'] == '0':
	err_msg = '## ERROR! PARENT SAMPLE, NO PARENTAL GVCF PATH'
	fa_path, mo_path = '.', '.'
else:
	fa_path = pathd[pedd[sample_id]['fa']]
	mo_path = pathd[pedd[sample_id]['mo']]

print(pb_path)
print(fa_path)
print(mo_path)


pbout = open('tmp.pb_path.txt', 'w')
pbout.write(pb_path)
pbout.close()

faout = open('tmp.fa_path.txt', 'w')
faout.write(fa_path)
faout.close()

moout = open('tmp.mo_path.txt', 'w')
moout.write(mo_path)
moout.close()




