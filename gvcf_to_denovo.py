#!/usr/bin/python3
## Purpose: call de novos from a single sample gVCF
'''
Usage: gvcf_to_denovo.py -s <sample id> \
                         -m <sample map> \
                         -p <ped file> \
                         -x <proband min altdp> \
                         -y <parent max altdp> \
                         =z <parent min dp> \
                         -o <output filename> 

Options:
  -h, --help            show this help message and exit
  -s sample_id, --sid=sample_id
                        sample ID
  -m sample_map, --smap=sample_map
                        sample map
  -p ped, --ped=ped
                    	  pedigree file in plink format
  -x pb_min_alt, --min_alt=pb_min_alt
                        proband minimum alt allele read depth
  -y par_max_alt, --max_alt=par_max_alt
                        parent maximum alt allele read depth
  -z par_min_dp, --min_dp=par_min_dp
                        parent minimum read depth
  -o output_file, --output=output_file
                    	output tab-separated variants file with additional columns
## CAVEATS:
# -assumes parent gvcfs are tabix indexed and .tbi files are present in same directory as gvcf
# -SNVs only, no indels or SVs
# -splits multiallelic sites into individual lines
# -ignores missing genotypes ('./.') and sites without AD or DP information

# Output format:
# chr, pos, ref, alt, refdp, altdp, dp, adfref, adfalt, adrref, adralt, <proband VCF information repeated>,
#                <father FORMAT field>, <father GT information>,
#                <mother FORMAT field>, <mother GT information>
'''
import sys
from optparse import OptionParser
import subprocess
import os
import gzip

####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-s', '--sid', dest='sample_id',help='sample id')
parser.add_option('-m', '--smap', dest='sample_map',help='sample map (picard)')
parser.add_option('-p', '--ped', dest='ped', help='pedigree file (plink)')
parser.add_option('-x', '--min_alt', dest='pb_min_alt',help='proband minimum alternate allele read depth')
parser.add_option('-y', '--max_alt', dest='par_max_alt',help='parent maximum alternate allele read depth')
parser.add_option('-z', '--min_dp', dest='par_min_dp',help='parent minimum read depth')
parser.add_option('-o', '--output', dest='output_file',help='output tab-separated variants file')
(options, args) = parser.parse_args()

## check all arguments present
if (options.sample_id == None or options.sample_map == None or options.ped == None or options.pb_min_alt == None or options.par_max_alt == None or options.par_min_dp == None):
	print '\n' + '## ERROR: missing arguments' + '\n'
	parser.print_help()
	print '\n'
	sys.exit()


sample_id = options.sample_id
sample_map = options.sample_map
ped = options.ped
pb_min_alt = options.pb_min_alt
par_max_alt = options.par_max_alt
par_min_dp = options.par_min_dp
output_file = options.output_file

####################################################################################################
## Function that, given parent gvcf and variant information, checks if variant is present
## if multi-line return, iterate over tabix result and compare chr:pos:ref:alt
## return parent dp, parent altdp, parent FORMAT, parent gt
####################################################################################################
def parse_parent(gvcf, region, chr, pos, ref, alt):

  tmp_head = subprocess.Popen(('tabix', '-H', gvcf), stdout=subprocess.PIPE) # get header lines
  tmp_cols = subprocess.check_output(('grep', '^#CHROM'), stdin=tmp_head.stdout).strip().split('\t') # parse column names
  tmp_head.wait()

  tmp = subprocess.check_output(['tabix',gvcf, region]).strip().split('\n')

  # intialize values
  # handles the case of empty tabix return
  dp = 0
  altdp = 0
  fmt = 'NA'
  gt = 'NA'
  inf = 'NA'
  
  if len(tmp) == 1 and (not tmp[0] == ''): # single variant return, non-empty tabix result
    tmpv = tmp[0].split('\t')
    tmp_d = dict(zip(tmp_cols, tmpv))
    tmp_gtd = dict(zip(tmp_d['FORMAT'].split(':'), tmpv[-1].split(':'))) # dictionary of FORMAT : GT value mapping

    fmt = tmp_d['FORMAT']
    gt = tmpv[-1]
    inf = tmp_d['INFO']
    ## check coordinates and parse AD, DP information 
    if ('END=' in tmp_d['INFO']): # if non-variant block, altdp = 0
      try:
        dp = int(tmp_gtd['DP'])
      except:
        dp = 0
    elif ('AS_RAW' in tmp_d['INFO']): # if variant block
      if ('GT' in tmp_d['FORMAT']) and ('AD' in tmp_d['FORMAT']) and ('DP' in tmp_d['FORMAT']): # check that necessary info is there
        
        if tmp_d['#CHROM'] == chr and tmp_d['POS'] == pos: # check that coordinates are correct
          if alt in tmp_d['ALT'].strip(',<NON_REF>').split(','): # check that proband variant allele present among ALT in parent
            altidx = tmp_d['ALT'].strip(',<NON_REF>').split(',').index(alt) + 1 # Get index of alternate allele
            if tmp_d['ALT'].strip(',<NON_REF>').split(',')[altidx-1] == alt: # note: alt_idx is for parsing GT field.  For parsing ALT col, use alt_idx-1
              if not './.' in tmp_gtd['GT']:
                altdp = int(tmp_gtd['AD'].split(',')[altidx-1])
                dp = int(tmp_gtd['DP'])
          else: # if variant allele not present, parent has effective altdp = 0
            try: # handle missing DP cases e.g. tabix 11003-mo.g.vcf.gz chr8:144530986-144530986
              dp = int(tmp_gtd['DP'])
            except:
              dp = 0

  elif len(tmp) > 1 and (not tmp[0] == ''): # multiple variant return, non-empty tabix result
    for v in tmp:
      tmpv = v.split('\t')
      tmp_d = dict(zip(tmp_cols, tmpv))
      tmp_gtd = dict(zip(tmp_d['FORMAT'].split(':'), tmpv[-1].split(':'))) # dictionary of FORMAT : GT value mapping
      
      fmt = tmp_d['FORMAT']
      gt = tmpv[-1]
      inf = tmp_d['INFO']
      ## check coordinates and parse AD, DP information 
      if ('END=' in tmp_d['INFO']): # if non-variant block, altdp = 0
        altdp = 0
        try:
          dp = int(tmp_gtd['DP'])
        except:
          dp = 0
      elif ('AS_RAW' in tmp_d['INFO']): # if variant block
        if ('GT' in tmp_d['FORMAT']) and ('AD' in tmp_d['FORMAT']) and ('DP' in tmp_d['FORMAT']): # check that necessary info is there

          if tmp_d['#CHROM'] == chr and tmp_d['POS'] == pos: # check that coordinates are correct
            if alt in tmp_d['ALT'].strip(',<NON_REF>').split(','): # check that proband variant allele present among ALT in parent
              altidx = tmp_d['ALT'].strip(',<NON_REF>').split(',').index(alt) + 1 # Get index of alternate allele
              if tmp_d['ALT'].strip(',<NON_REF>').split(',')[altidx-1] == alt: # note: alt_idx is for parsing GT field.  For parsing ALT col, use alt_idx-1
                if not './.' in tmp_gtd['GT']:
                  altdp = int(tmp_gtd['AD'].split(',')[altidx-1])
                  dp = int(tmp_gtd['DP'])
            else: # if variant allele not present, parent has effective altdp = 0
              altdp = 0
              try: # handle missing DP cases e.g. tabix 11003-mo.g.vcf.gz chr8:144530986-144530986
                dp = int(tmp_gtd['DP'])
              except:
                dp = 0

  outd = {'altdp': altdp, 'dp': dp, 'fmt': fmt, 'gt':gt, 'info': inf}

  return(outd)

####################################################################################################
## read sample map and create dictionary
####################################################################################################
print('')
print('## READING SAMPLE MAP')
print('')
pathd = {} # {id : path}
with open(sample_map, 'r') as smapf:
  for line in smapf:
    tmp = line.strip().split('\t')
    id = tmp[0]
    path = tmp[1]
    pathd[id] = path

####################################################################################################
## read pedigree file and create dictionaary
####################################################################################################
print('')
print('## READING PEDIGREE FILE')
print('')
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

####################################################################################################
## iterate over proband gVCF and 
####################################################################################################
bufsize=1
outf = open(output_file, 'w', buffering=bufsize)


##
## IF NOT PROBAND OR SIBLING, WRITE ERROR MESSAGE AND EXIT
##
if pedd[sample_id]['fa'] == '0' or pedd[sample_id]['mo'] == '0':
  err_msg = '## ERROR! PARENT SAMPLE, UNABLE TO CALL DE NOVOS'
  
  print(err_msg)
  
  outf.write(err_msg)
  outf.close()
  
  sys.exit()



## get gvcf paths
sample_gvcf = pathd[sample_id]
fa_gvcf = pathd[pedd[sample_id]['fa']]
mo_gvcf = pathd[pedd[sample_id]['mo']]



cmd2 = 'zcat < %s | grep -v "#"| wc -l'%(sample_gvcf)
tot = int(subprocess.check_output(cmd2, shell=True).strip().split(' ')[0])
print('## TOTAL VARIANT LINES: %s'%(str(tot)))


print('')
print('## ITERATING OVER VARIANT LINES')
print('')

head = ['id', 'chr', 'pos', 'ref', 'alt', 'dp', 'adfref', 'adfalt', 'adrref', 'adralt', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'S_GT', 'FA_FORMAT', 'FA_GT', 'MO_FORMAT', 'MO_GT']
#print '\t'.join(head)
outf.write('\t'.join(head) + '\n')

i = 0


## iterate over proband gVCF
with gzip.open(sample_gvcf, 'r') as f:  
#with open(tmp_var_only, 'r') as f:  
  for line in f:


    tmp = line.strip().split('\t')

    ## handle vcf header information
    if line.startswith('##'):
      continue
    
    ## handle vcf header containing column information
    if line.startswith('#CHROM'):
      idx = {col:index for index, col in enumerate(tmp)}

    ## handle variant lines
    else:
      i += 1
      if i%1000 == 0:
        print('## %d/%d lines processed ... '%(i, tot))

      ## initialize values to avoid iteration bugs
      chr, pos, ref, alt = '', '', '',''
      region = ''
      info = ''
      fmt, gt = [], []
      gtd = {}


      ## ignore non-variant blocks 
      if not 'END=' in line.strip():
        
        ## get proband variant information
        chr, pos, ref = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']]

        ## parse alt allele
        alt = tmp[idx['ALT']].strip(',<NON_REF>')

        
        ## how to handle multiallelic sites? e.g. chr1    1646352 .       A       C,G,<NON_REF>
        ## iterate over all alternate alleles present in ALT
        for a in alt.split(','):


          if not a == '*': # ignore point deletions for now; messy when matching alleles with parents
            if len(ref) == 1 and len(a) == 1: # ignore indels for now;
              
              
              # Get index of current alternate allele
              pb_altidx = alt.split(',').index(a) + 1 

              # get region for tabixing parents
              region = chr + ':' + pos + '-' + pos

              # save INFO field
              info = tmp[idx['INFO']]

              # create dictionary of FORMAT:GT mapping
              fmt = tmp[idx['FORMAT']].split(':')
              gt = tmp[-1].split(':') ## ASSUMES THAT SAMPLE GENOTYPE INFORMATION IS IN THE LAST COLUMN; didn't use ID since column ID differs from sample id....

              gtd = dict(zip(fmt, gt)) # e.g. {'GT': '0/1', 'AD': '5,7,0', 'GQ': '99', 'PL': '157,0,104,172,125,297', 'SB': '5,0,7,0', 'DP': '12'}
              
              print(region)
              #print(gtd)


              if not gtd['GT'] == './.': ## ignore sites with missing genotypes
                if ('AD'in gtd) and ('DP' in gtd): # ignore sites with no AD or DP information

                  ## parse strand-specific allelic depth information
                  adf = gtd['F1R2']
                  adr = gtd['F2R1']

                  adfref = adf.split(',')[0]
                  adrref = adr.split(',')[0]

                  adfalt = adf.split(',')[pb_altidx]
                  adralt = adr.split(',')[pb_altidx]

                  pb_refdp = int(gtd['AD'].split(',')[0])
                  pb_altdp = int(gtd['AD'].split(',')[pb_altidx])
                  pb_dp = int(gtd['DP'])


                  fa_d = parse_parent(fa_gvcf, region, chr, pos, ref, a)

                  fa_altdp = fa_d['altdp']
                  fa_dp = fa_d['dp']
                  fa_fmt = fa_d['fmt']
                  fa_gt = fa_d['gt']


                  mo_d = parse_parent(mo_gvcf, region, chr, pos, ref, a)

                  mo_altdp = mo_d['altdp']
                  mo_dp = mo_d['dp']
                  mo_fmt = mo_d['fmt']
                  mo_gt = mo_d['gt']

                  

                  ## APPLY DE NOVO CALLING CRITERIA
                  if not (int(pb_refdp) == 0): # ignore hom alt sites
                    if int(pb_altdp) >= int(pb_min_alt):
                      if int(fa_altdp) <= int(par_max_alt) and int(mo_altdp) <= int(par_max_alt):
                        if int(fa_dp) >= int(par_min_dp) and int(mo_dp) >= int(par_min_dp):
                          out = map(str, [sample_id, chr.strip('chr'), pos, ref, a, pb_refdp, pb_altdp, pb_dp, adfref, adfalt, adrref, adralt])

                          #print '\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + tmpfa_d['FORMAT'] + '\t' + tmpfa[-1] + '\t' + tmpmo_d['FORMAT'] + '\t' + tmpmo[-1]
                          #outf.write('\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + tmpfa_d['FORMAT'] + '\t' + tmpfa[-1] + '\t' + tmpmo_d['FORMAT'] + '\t' + tmpmo[-1] + '\n')
                          outstring = '\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + fa_fmt + '\t' + fa_gt + '\t' + mo_fmt + '\t' + mo_gt
                          print outstring
                          outf.write(outstring + '\n')





outf.close()





