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
  
  File script
  String sample_id 
  Array[File] sample_map
  File ped
  Int pb_min_alt
  Int par_max_alt
  Int par_min_dp
  String output_suffix


  parameter_meta{
    script: "gvcf_to_denovo.py"
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

  call localize_sample_map {
    input:
    in_map = sample_map
  }



  #Outputs a .txt file containing de novo SNVs
  output {
    File localout = localize_sample_map.outfile
      
  }

}


###########################################################################
#Task Definitions
###########################################################################

#Generates new sample map using localized file paths
task localize_sample_map {
  Array[File] in_map

  command{

    echo ${sep=" " in_map}

    echo ${sep=" " in_map} > tmp.stdout.txt

  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"

  }

  output {
    File outfile = "tmp.stdout.txt"
  }
}

#Scores variants in dnSNVs file
task call_denovos {
  File script
  String sample_id
  Array[File] sample_map
  File ped
  Int pb_min_alt
  Int par_max_alt
  Int par_min_dp
  String output_suffix
  String output_file = "${sample_id}${output_suffix}"

  command {

    python -u ${script} -s ${sample_id} -m ${sample_map} -p ${ped} -x ${pb_min_alt} -y ${par_max_alt} -z ${par_min_dp} -o ${output_file}
  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"

  }

  output {
    File outfile = "${output_file}"
  }
}