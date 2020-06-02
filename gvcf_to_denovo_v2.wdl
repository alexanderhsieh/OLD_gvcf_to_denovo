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
    dn_script: "gvcf_to_denovo_v2.py"
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
    ped = ped,
    sample_id = sample_id

  }

  # Step *: split gvcf by chromosome
  call split_gvcf {
    input:
    gvcf = localize_path.local_pb_gvcf,
    index = localize_path.local_pb_gvcf_index
  }

  # for each chr gvcf, call de novos
  
  scatter (idx in range(length(split_gvcf.out))) {
    
    call call_denovos {
      input:
      script = dn_script,
      sample_id = sample_id,

      sample_gvcf = split_gvcf.out[idx],
      sample_gvcf_index = split_gvcf.out_index[idx],
      father_gvcf = localize_path.local_fa_gvcf,
      father_gvcf_index = localize_path.local_fa_gvcf_index,
      mother_gvcf = localize_path.local_mo_gvcf,
      mother_gvcf_index = localize_path.local_mo_gvcf_index,

      sample_map = sample_map,
      ped = ped,
      pb_min_alt = pb_min_alt,
      par_max_alt = par_max_alt,
      par_min_dp = par_min_dp,

      shard = "${idx}"

    }


  }

  # Step 3: gather shards into final output 
  call gather_shards {
    input:
    shards = call_denovos.outfile,
    prefix = sample_id,
    suffix = output_suffix
  }


  


  #Outputs a .txt file containing de novo SNVs
  output {

    File denovos = gather_shards.out

    File split_gvcf_header = split_gvcf.header
    Array[File] call_denovos_header = call_denovos.header
      
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
task localize_path {
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
    
    gsutil cp $PB_PATH ./tmp.pb.g.vcf.gz
    gsutil cp $PB_PATH".tbi" ./tmp.pb.g.vcf.gz.tbi


    ## LOCALIZE FATHER AND MOTHER
    FA_PATH=`cat tmp.fa_path.txt`
    MO_PATH=`cat tmp.mo_path.txt`

    ## if no father listed in ped
    if [ "$FA_PATH" = "." ]
    then
      touch ./tmp.fa.g.vcf.gz
      touch ./tmp.fa.g.vcf.gz.tbi
      echo "## ERROR: MISSING FATHER GVCF PATH"; exit
    
    ## if no mother listed in ped
    elif [ "$MO_PATH" = "." ]
    then
      touch ./tmp.mo.g.vcf.gz
      touch ./tmp.mo.g.vcf.gz.tbi
      echo "## ERROR: MISSING MOTHER GVCF PATH"; exit
    
    ## if both parents found
    else
      echo "## FATHER BUCKET PATH: "$FA_PATH
      echo "## MOTHER BUCKET PATH: "$MO_PATH

      gsutil cp $FA_PATH ./tmp.fa.g.vcf.gz
      gsutil cp $FA_PATH".tbi" ./tmp.fa.g.vcf.gz.tbi
      
      gsutil cp $MO_PATH ./tmp.mo.g.vcf.gz
      gsutil cp $MO_PATH".tbi" ./tmp.mo.g.vcf.gz.tbi
    fi


  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"

  }

  output {
    File local_pb_gvcf = "tmp.pb.g.vcf.gz"
    File local_pb_gvcf_index = "tmp.pb.g.vcf.gz.tbi"
    File local_fa_gvcf = "tmp.fa.g.vcf.gz"
    File local_fa_gvcf_index = "tmp.fa.g.vcf.gz.tbi"
    File local_mo_gvcf = "tmp.mo.g.vcf.gz"
    File local_mo_gvcf_index = "tmp.mo.g.vcf.gz.tbi"
  }
}

## splits vcf by chromosome
task split_gvcf {

  File gvcf # input gvcf
  File index # input gvcf index
  String outprefix = basename(gvcf, '.g.vcf.gz')

  command {
    # pull header lines
    zcat < ${gvcf} | grep "^#" ${gvcf} > header.txt

    # split vcf by chromosome - use tabix -l to get all contig names from tabix index
    for i in $(tabix -l ${gvcf})
    do 
      (cat header.txt; tabix ${gvcf} $i) | bgzip -c > "${outprefix}.$i.g.vcf.gz"
    done

    ## get full directory paths
    readlink -f *.g.vcf.gz > file_full_paths.txt

  }

  output {
    Array[File] out = glob("*.g.vcf.gz") 
    Array[File] out_index = glob("*.g.vcf.gz.tbi")
    File filepaths = "file_full_paths.txt"

    File header = "header.txt"
  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
  }

}

#Calls denovos from proband gvcf + parent paths
# NOTE: currently runs gsutil cp to localize proband gvcf
task call_denovos {
  File script
  String sample_id

  File sample_gvcf
  File sample_gvcf_index

  File father_gvcf
  File father_gvcf_index
  File mother_gvcf
  File mother_gvcf_index


  File sample_map
  File ped
  Int pb_min_alt
  Int par_max_alt
  Int par_min_dp

  String shard
  
  String output_file = "${sample_id}.${shard}.denovo.txt"

  command {

    zcat ${sample_gvcf} | head -n 100 > header.${shard}.txt

    python ${script} -s ${sample_id} -p ${sample_gvcf} -f ${father_gvcf} -m ${mother_gvcf} -r ${ped} -x ${pb_min_alt} -y ${par_max_alt} -z ${par_min_dp} -o ${output_file}
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"

  }

  output {
    File outfile = "${output_file}"
    
    File header = "header.${shard}.txt"
  }
}

#Gathers shards of raw de novo call files into a single call set
task gather_shards {

  Array[File] shards 
  String prefix
  String suffix

  command {

    while read file; do
      cat $file >> "tmp.cat.denovo.raw.txt"
    done < ${write_lines(shards)};

    grep "^id" "tmp.cat.denovo.raw.txt" > "header.txt" 

    (cat header.txt; grep -v "^id" "tmp.cat.denovo.raw.txt") > "${prefix}${suffix}"

  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
  }

  output {
    File out = "${prefix}${suffix}"
  }
}

