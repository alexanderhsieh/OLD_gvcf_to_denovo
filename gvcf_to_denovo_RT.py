#!/usr/bin/python3
## Purpose: call de novos from a single sample gVCF
'''
Usage: gvcf_to_denovo.py -s <sample id> \
                         -p <proband gvcf> \
                         -f <father gvcf> \
                         -m <mother gvcf> \
                         -r <relations in pedigree format>
                         -x <proband min altdp> \
                         -y <parent max altdp> \
                         =z <parent min dp> \
                         -o <output filename> 

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
import io


####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-s', '--sid', dest='sample_id',help='sample id')
parser.add_option('-p', '--pb', dest='sample_gvcf', help='sample gvcf')
parser.add_option('-f', '--fa', dest='fa_gvcf',help='father gvcf path')
parser.add_option('-m', '--mo', dest='mo_gvcf',help='mother gvcf path')
parser.add_option('-r', '--ped', dest='ped', help='ped file')
parser.add_option('-x', '--min_vaf', dest='pb_min_vaf',help='proband minimum variant allele frequency')
parser.add_option('-y', '--max_alt', dest='par_max_alt',help='parent maximum alternate allele read depth')
parser.add_option('-z', '--min_dp', dest='par_min_dp',help='parent minimum read depth')
parser.add_option('-o', '--output', dest='output_file',help='output tab-separated variants file')
(options, args) = parser.parse_args()

## check all arguments present
if (options.sample_id == None or options.sample_gvcf == None or options.fa_gvcf == None or options.mo_gvcf == None or options.ped == None or options.pb_min_vaf == None or options.par_max_alt == None or options.par_min_dp == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()


sample_id = options.sample_id
sample_gvcf = options.sample_gvcf
fa_gvcf = options.fa_gvcf
mo_gvcf = options.mo_gvcf
ped = options.ped
pb_min_vaf = options.pb_min_vaf
par_max_alt = options.par_max_alt
par_min_dp = options.par_min_dp
output_file = options.output_file

####################################################################################################
## Function that, given parent gvcf and variant information, checks if variant is present
## if multi-line return, iterate over tabix result and compare chr:pos:ref:alt
## return parent dp, parent altdp, parent FORMAT, parent gt
####################################################################################################
def parse_parent(gvcf, region, chr, pos, ref, alt):
  tabix_cmd = 'tabix -H %s'%(gvcf)
  #tmp_head = subprocess.run(tabix_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8', text=True) # get header lines
  tmp_head = subprocess.run(tabix_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8') # get header lines
  tmp_cols = tmp_head.stdout.strip().split('\n')[-1].split('\t')

  # GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
  # tabix -R region.merged.bed ~{splitfile} | bgzip -c > SR_1.txt.gz
  # gs://fc-8f6b542d-a30e-4164-92c6-40eb8887a470/3f0e5d0a-4da6-48d2-91ce-e675b23336cb/ExomeGermlineSingleSample/59f64ccd-2014-4d13-b3a9-edf6b6754b42/call-BamToGvcf/VariantCalling/b8dea78a-0037-473f-b3e2-153168339273/call-MergeVCFs/11003-fa.g.vcf.gz
  gcs_cmd = 'GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`'
  tmp_gcs = subprocess.run(gcs_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8')
  tmp = subprocess.run('tabix %s %s'%(gvcf, region), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout


  #tmp = subprocess.run('tabix %s %s'%(gvcf, region), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').stdout
  


  #print('##')
  #print('## %s'%(tabix_cmd))
  #print('## %s'%('tabix %s %s'%(gvcf, region)))
  #print('## %s'%(tmp))

  # intialize values
  # handles the case of empty tabix return
  dp = 0
  altdp = 0
  fmt = 'NA'
  gt = 'NA'
  inf = 'NA'
  
  if len([tmp]) == 1 and (not [tmp][0] == ''): # single variant return, non-empty tabix result
    tmpv = tmp.strip().split('\t')
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
                #print('# parent variant found')
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
                  #print('# parent variant found')
            else: # if variant allele not present, parent has effective altdp = 0
              altdp = 0
              try: # handle missing DP cases e.g. tabix 11003-mo.g.vcf.gz chr8:144530986-144530986
                dp = int(tmp_gtd['DP'])
              except:
                dp = 0

  outd = {'altdp': altdp, 'dp': dp, 'fmt': fmt, 'gt':gt, 'info': inf}

  return(outd)


####################################################################################################
## iterate over proband gVCF and 
####################################################################################################
print('## GVCF PATHS:')
print('## PROBAND: %s'%(sample_gvcf))
print('## FATHER: %s'%(fa_gvcf))
print('## MOTHER: %s'%(mo_gvcf))

#bufsize=1
#outf = open(output_file, 'w', buffering=bufsize)
outf = open(output_file, 'w')



cmd2 = 'cat %s | grep -v "#"| wc -l'%(sample_gvcf)
#tot = int(subprocess.check_output(cmd2, shell=True, encoding='utf8').strip().split(' ')[0])
tot = int(subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf8').communicate()[0].strip().split(' ')[0])
print('## TOTAL VARIANT LINES: %s'%(str(tot)))


print('')
print('## ITERATING OVER VARIANT LINES')
print('')

head = ['id', 'chr', 'pos', 'ref', 'alt', 'refdp', 'altdp', 'dp', 'adfref', 'adfalt', 'adrref', 'adralt', 'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'S_GT', 'FA_FORMAT', 'FA_GT', 'MO_FORMAT', 'MO_GT']
#print '\t'.join(head)
outf.write('\t'.join(head) + '\n')

i = 0
dnct = 0

## list of potentially problematic parental genotypes
gt_blacklist = ['0/1', '0|1', '0/2', '0|2', '1/1', '1|1', '1/2', '1|2', '2/2', '2|2', '1/0', '1|0', '1/3', '1|3', '2/3', '2|3', '3/3', '3|3']

## iterate over proband gVCF
#with gzip.open(sample_gvcf, 'rb') as f:
#  with io.TextIOWrapper(f, encoding='utf-8') as decodef:  
#    for line in decodef:

with open(sample_gvcf, 'r') as f:

  for line in f:
    #print(line)

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
              
              #print(region)
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

                  if pb_dp > 0:
                    pb_vaf = float(pb_altdp)/float(pb_dp)
                  else:
                    pb_vaf = 0.0

                  ## parse father gvcf using tabix
                  fa_d = parse_parent(fa_gvcf, region, chr, pos, ref, a)

                  fa_altdp = fa_d['altdp']
                  fa_dp = fa_d['dp']
                  fa_fmt = fa_d['fmt']
                  fa_gt = fa_d['gt']

                  ## parse mother gvcf using tabix
                  mo_d = parse_parent(mo_gvcf, region, chr, pos, ref, a)

                  mo_altdp = mo_d['altdp']
                  mo_dp = mo_d['dp']
                  mo_fmt = mo_d['fmt']
                  mo_gt = mo_d['gt']

                  

                  ## APPLY DE NOVO CALLING CRITERIA
                  if not (int(pb_refdp) == 0): # ignore hom alt sites
                    if float(pb_vaf) >= float(pb_min_vaf):
                      if int(fa_altdp) <= int(par_max_alt) and int(mo_altdp) <= int(par_max_alt):
                        if int(fa_dp) >= int(par_min_dp) and int(mo_dp) >= int(par_min_dp):
                          if not (fa_gt in gt_blacklist or mo_gt in gt_blacklist):
                            out = map(str, [sample_id, chr.strip('chr'), pos, ref, a, pb_refdp, pb_altdp, pb_dp, adfref, adfalt, adrref, adralt])

                            #print '\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + tmpfa_d['FORMAT'] + '\t' + tmpfa[-1] + '\t' + tmpmo_d['FORMAT'] + '\t' + tmpmo[-1]
                            #outf.write('\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + tmpfa_d['FORMAT'] + '\t' + tmpfa[-1] + '\t' + tmpmo_d['FORMAT'] + '\t' + tmpmo[-1] + '\n')
                            outstring = '\t'.join(out) + '\t' + '\t'.join(tmp) + '\t' + fa_fmt + '\t' + fa_gt + '\t' + mo_fmt + '\t' + mo_gt
                            #print(outstring)
                            outf.write(outstring + '\n')
                            # deal with empty output?
                            outf.flush()
                            os.fsync(outf)

                            dnct += 1

                            print('## %d de novo variants found ...'%(dnct))





outf.close()


