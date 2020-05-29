## Copyright Broad Institute, 2020
## This workflow calls de novo SNVs using sample gVCF + paternal/maternal gVCFs
## Requires (1) sample id, (2) sample map (picard), (3) pedigree file (plink), (4) options for de novo calling criteria
## 
##  NOTE: uses docker: "gatksv/sv-base-mini:cbb1fc" to borrow samtools, tabix functionality
##
## TESTED: 
## Versions of other tools on this image at the time of testing:
##
## LICENSING : This script is released under the WDL source code license (BSD-3) (see LICENSE in https://github.com/broadinstitute/wdl). 
## Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. 
## Please see the docker for detailed licensing information pertaining to the included programs.
##


###########################################################################
#WORKFLOW DEFINITION
###########################################################################
workflow gvcf_to_denovo {
  
  File localize_script
  File dn_script
  String sample_id 
  File sample_map
  File ped
  Int pb_min_alt
  Int par_max_alt
  Int par_min_dp
  String output_suffix


  parameter_meta{
    localize_script: "parse_sample_map.py"
    dn_script: "gvcf_to_denovo.py"
    sample_id: "sample ID for which to call de novo SNVs"
    sample_map: "sample map containing id:gvcf_path mapping; generated via Picard"
    ped: "pedigree file containing relatedness information; plink format"
    pb_min_alt: "proband; minimum number of reads supporting the variant allele"
    par_max_alt: "parent; maximum number of reads supporting the variant allele"
    par_min_dp: "parent; minimum read depth at the variant position"
    output_suffix: "output de novo SNVs filename suffix"
  }
  meta{
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }

  call localize_path{
    input:
    script = localize_script,
    sample_map = sample_map,
    sample_id = sample_id

  }

  call call_denovos {
    input:
    script = dn_script,
    sample_id = sample_id,

    sample_gvcf = localize_path.local_gvcf,

    sample_map = sample_map,
    ped = ped,
    pb_min_alt = pb_min_alt,
    par_max_alt = par_max_alt,
    par_min_dp = par_min_dp,
    output_suffix = output_suffix
  }


  #Outputs a .txt file containing de novo SNVs
  output {
    File denovos = call_denovos.outfile
      
  }

}


###########################################################################
#Task Definitions
###########################################################################
# Takes 
task localize_path {
  File script
  File sample_map
  String sample_id

  command{
    
    TMP_PATH=`python -m ${sample_map} -s ${sample_id}`
    
    echo "## BUCKET PATH: "$TMP_PATH

    gsutil cp $TMP_PATH ./tmp.g.vcf.gz

    echo `ls -hl tmp.g.vcf.gz`

  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"

  }

  output {
    File local_gvcf = "tmp.g.vcf.gz"
  }
}

#Calls denovos from proband gvcf + parent paths
# NOTE: currently runs gsutil cp to localize proband gvcf
task call_denovos {
  File script
  String sample_id
  File sample_gvcf
  File sample_map
  File ped
  Int pb_min_alt
  Int par_max_alt
  Int par_min_dp
  String output_suffix
  String output_file = "${sample_id}${output_suffix}"

  command {

    python -u ${script} -s ${sample_id} -g ${sample_gvcf} -m ${sample_map} -p ${ped} -x ${pb_min_alt} -y ${par_max_alt} -z ${par_min_dp} -o ${output_file}

  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"

  }

  output {
    File outfile = "${output_file}"
  }
}