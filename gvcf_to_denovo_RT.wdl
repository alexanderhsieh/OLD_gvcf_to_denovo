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
  File sample_gvcf
  File sample_gvcf_index
  File sample_map
  File ped
  Float pb_min_vaf
  Int par_max_alt
  Int par_min_dp
  String output_suffix


  parameter_meta{
    localize_script: "parse_sample_map.py"
    dn_script: "gvcf_to_denovo_v4.py"
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

  call parse_parent_paths {
    input:
    script = localize_script,
    sample_map = sample_map,
    ped = ped,
    sample_id = sample_id

  }

  # Step *: split gvcf by chromosome
  call split_gvcf {
    input:
    gvcf = sample_gvcf,
    index = sample_gvcf_index
  }

  # for each chr vcf, call de novos
  scatter (idx in range(length(split_gvcf.out))) {
    
    call call_denovos {
      input:
      script = dn_script,
      sample_id = sample_id,

      sample_gvcf = split_gvcf.out[idx],
      father_gvcf_path = parse_parent_paths.father_path,
      mother_gvcf_path = parse_parent_paths.mother_path,

      sample_map = sample_map,
      ped = ped,
      pb_min_vaf = pb_min_vaf,
      par_max_alt = par_max_alt,
      par_min_dp = par_min_dp,

      shard = "${idx}"

    }


  }

  # Step 3: gather shards into final output 
  call gather_shards {
    input:
    shards = call_denovos.outfile,
    headers = call_denovos.header,
    prefix = sample_id,
    suffix = output_suffix
  }


  
  #Outputs a .txt file containing de novo SNVs
  output {

    File denovos = gather_shards.out
      
  }

}


###########################################################################
#Task Definitions
###########################################################################
# Takes sample map as input, generates 3 temporary text files containing the 
# google bucket path, then copies the gvcfs and tabix indices to the current
# workflow instance. 
# NOTE: As this is for calling de novos, requires trio to be present
# if from the ped no parents are listed, print status message and exit
task parse_parent_paths {
  File script
  File sample_map
  File ped
  String sample_id

  command{
    
    ## PARSE SAMPLE MAP GOOGLE BUCKET PATHS
    python ${script} -m ${sample_map} -p ${ped} -s ${sample_id}
    
    ## LOCALIZE PROBAND
    PB_PATH=`cat tmp.pb_path.txt`
    
    echo "## PROBAND BUCKET PATH: "$PB_PATH

    ## LOCALIZE FATHER AND MOTHER
    FA_PATH=`cat tmp.fa_path.txt`
    MO_PATH=`cat tmp.mo_path.txt`

    echo "## FATHER BUCKET PATH: "$FA_PATH
    echo "## MOTHER BUCKET PATH: "$MO_PATH

  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1" 
    preemptible: 3
    maxRetries: 3
  }

  output {
    File sample_path = "tmp.pb_path.txt"
    File father_path = "tmp.fa_path.txt"
    File mother_path = "tmp.mo_path.txt"

  }
}

## splits vcf by chromosome
task split_gvcf {

  File gvcf # input gvcf
  File index # input gvcf index
  String outprefix = basename(gvcf, '.g.vcf.gz')

  Int disk_size = 100 # start with 100G

  command {

    tabix -H ${gvcf} > "header.txt"

    # split vcf by chromosome - use tabix -l to get all contig names from tabix index
    for i in $(tabix -l ${gvcf})
    do 
      (cat header.txt; tabix ${gvcf} $i)  > "${outprefix}.$i.vcf"
    done

    ## get full directory paths
    readlink -f *.vcf > file_full_paths.txt

  }

  output {
    Array[File] out = glob("*.vcf") 
    File filepaths = "file_full_paths.txt"

  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: disk_size
  }

}

#Calls denovos from proband gvcf + parent paths
# NOTE: currently runs gsutil cp to localize proband gvcf
task call_denovos {
  File script
  String sample_id

  File sample_gvcf

  File father_gvcf_path
  File mother_gvcf_path

  File sample_map
  File ped
  Float pb_min_vaf
  Int par_max_alt
  Int par_min_dp

  String shard
  
  String output_file = "${sample_id}.${shard}.denovo.txt"

  command {

    FA_PATH=`cat ${father_gvcf_path}`
    MO_PATH=`cat ${mother_gvcf_path}`

    python ${script} -s ${sample_id} -p ${sample_gvcf} -f $FA_PATH -m $MO_PATH -r ${ped} -x ${pb_min_vaf} -y ${par_max_alt} -z ${par_min_dp} -o ${output_file}

    head -n 1 ${output_file} > "header.txt"
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
    preemptible: 3
    maxRetries: 3
  }

  output {
    File outfile = "${output_file}"
    File header = "header.txt"

  }
}

#Gathers shards of raw de novo call files into a single call set
task gather_shards {

  Array[File] shards 
  Array[File] headers
  String prefix
  String suffix

  File head = select_first(headers)

  command {

    set -eou pipefail

    while read file; do
      cat $file | grep -v "^id" >> "tmp.cat.denovo.raw.txt"
    done < ${write_lines(shards)};

    (cat ${head}; cat "tmp.cat.denovo.raw.txt") > "${prefix}${suffix}"

  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
    preemptible: 3
    maxRetries: 3
  }

  output {
    File out = "${prefix}${suffix}"
  }
}

