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
  
  File dn_script

  Array[File] trio_gvcf_array
  Array[File] trio_gvcf_index_array
  Array[String] trio_readgroup_ids

  String sample_id 

  Float pb_min_vaf
  Int par_max_alt
  Int par_min_dp
  String output_suffix

  File ref_fasta
  File ref_fasta_index

  ## TO UPDATE
  parameter_meta{
    dn_script: "gvcf_to_denovo_ALT.py"
    sample_id: "sample ID for which to call de novo SNVs"
    sample_map: "sample map containing id:gvcf_path mapping; generated via Picard"
    ped: "pedigree file containing relatedness information; plink format"
    pb_min_vaf: "proband; minimum variant allele frequency to be called"
    par_max_alt: "parent; maximum number of reads supporting the variant allele"
    par_min_dp: "parent; minimum read depth at the variant position"
    output_suffix: "output de novo SNVs filename suffix"
  }
  meta{
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }


  call merge_trio_gvcf{
    input:
    sample_id = sample_id,
    trio_gvcf_array = trio_gvcf_array,
    trio_gvcf_index_array = trio_gvcf_index_array,
    trio_readgroup_ids = trio_readgroup_ids,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index
  }

  call call_denovos {
    input:
    script = dn_script,
    sample_id = merge_trio_gvcf.out_pb_id,
    father_id = merge_trio_gvcf.out_fa_id,
    mother_id = merge_trio_gvcf.out_mo_id,

    gvcf = merge_trio_gvcf.out_gvcf,

    pb_min_vaf = pb_min_vaf,
    par_max_alt = par_max_alt,
    par_min_dp = par_min_dp,

    output_suffix = output_suffix

  }


  
  #Outputs a .txt file containing de novo SNVs
  output {
    File denovos = call_denovos.out
      
  }

}


###########################################################################
#Task Definitions
###########################################################################

## merges proband, father, mother gvcfs into trio gvcf
task merge_trio_gvcf {
  
  String sample_id

  Array[File] trio_gvcf_array
  Array[File] trio_gvcf_index_array
  Array[String] trio_readgroup_ids

  File ref_fasta
  File ref_fasta_index

  String outfname = "${sample_id}.TRIO.g.vcf.gz"

  command {

    bcftools merge -g ${ref_fasta} -l ${write_lines(trio_gvcf_array)} -o ${outfname} -O z -m all

    tabix -p vcf ${outfname}

    cat ${write_lines(trio_readgroup_ids)} > id_list.txt
    sed -n '1p' id_list.txt = pb_id.txt
    sed -n '2p' id_list.txt = fa_id.txt
    sed -n '3p' id_list.txt = mo_id.txt
    

  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
    memory: "8G"
    preemptible: 3
    maxRetries: 3
  }

  output {
    File out_gvcf = "${outfname}"
    File out_idx = "${outfname}.tbi"

    File out_pb_id = "pb_id.txt"
    File out_fa_id = "fa_id.txt"
    File out_mo_id = "mo_id.txt"
  }

}


#Calls denovos from proband gvcf + parent paths
task call_denovos {
  File script
  
  File sample_id
  File father_id
  File mother_id

  File gvcf

  Float pb_min_vaf
  Int par_max_alt
  Int par_min_dp

  String output_suffix
  
  String output_file = "${sample_id}${output_suffix}"

  command {

    S_ID=`cat ${sample_id}`
    FA_ID=`cat ${father_id}`
    MO_ID=`cat ${mother_id}`

    python ${script} -s "$S_ID" -f "$FA_ID" -m "$MO_ID" -g ${gvcf} -x ${pb_min_vaf} -y ${par_max_alt} -z ${par_min_dp} -o ${output_file}

  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
    preemptible: 3
    maxRetries: 3
  }

  output {
    File out = "${output_file}"
  }
}



