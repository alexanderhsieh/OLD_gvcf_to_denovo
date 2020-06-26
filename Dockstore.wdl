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
  Float pb_min_vaf
  Int par_max_alt
  Int par_min_dp
  String output_suffix

  File ref_fasta
  File ref_fasta_index

  ## TO UPDATE
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

  call localize_path{
    input:
    script = localize_script,
    sample_map = sample_map,
    ped = ped,
    sample_id = sample_id

  }

  call merge_trio_gvcf{
    input:
    sample_id = sample_id,
    pb_gvcf = localize_path.local_pb_gvcf,
    pb_idx = localize_path.local_pb_gvcf_index,
    fa_gvcf = localize_path.local_fa_gvcf,
    fa_idx = localize_path.local_fa_gvcf_index,
    mo_gvcf = localize_path.local_mo_gvcf,
    mo_idx = localize_path.local_mo_gvcf_index,
    ref_fasta = ref_fasta,
    ref_fasta_index = ref_fasta_index
  }

  # Step *: split gvcf by chromosome
  call split_gvcf {
    input:
    gvcf = localize_path.local_pb_gvcf,
    index = localize_path.local_pb_gvcf_index,
    header = localize_path.header
  }

  # for each chr vcf, call de novos
  scatter (idx in range(length(split_gvcf.out))) {
    
    call call_denovos {
      input:
      script = dn_script,
      sample_id = localize_path.new_sample_id,
      father_id = localize_path.new_father_id,
      mother_id = localize_path.new_mother_id,

      gvcf = split_gvcf.out[idx],
  
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

    ## PARSE HEADER LINE
    zgrep "^#" ./tmp.pb.g.vcf.gz > header.txt

    ## LOCALIZE FATHER AND MOTHER
    FA_PATH=`cat tmp.fa_path.txt`
    MO_PATH=`cat tmp.mo_path.txt`

    ## if no father or mother listed in ped
    if [[ "$FA_PATH" == "." ]] | [[ "$MO_PATH" == "." ]]
    then
      touch ./tmp.fa.g.vcf.gz
      touch ./tmp.fa.g.vcf.gz.tbi
      echo "## ERROR: MISSING FATHER GVCF PATH"

      touch ./tmp.mo.g.vcf.gz
      touch ./tmp.mo.g.vcf.gz.tbi
      echo "## ERROR: MISSING MOTHER GVCF PATH"
    ## if both parents found
    else
      echo "## FATHER BUCKET PATH: "$FA_PATH
      echo "## MOTHER BUCKET PATH: "$MO_PATH

      gsutil cp $FA_PATH ./tmp.fa.g.vcf.gz
      gsutil cp $FA_PATH".tbi" ./tmp.fa.g.vcf.gz.tbi
      
      gsutil cp $MO_PATH ./tmp.mo.g.vcf.gz
      gsutil cp $MO_PATH".tbi" ./tmp.mo.g.vcf.gz.tbi
    fi

    ## parse sample name from vcf path to be passed downstream
    basename $PB_PATH '.g.vcf.gz' > pb_id.txt
    basename $FA_PATH '.g.vcf.gz' > fa_id.txt
    basename $MO_PATH '.g.vcf.gz' > mo_id.txt

  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1" 
    preemptible: 3
    maxRetries: 3
  }

  output {
    File local_pb_gvcf = "tmp.pb.g.vcf.gz"
    File local_pb_gvcf_index = "tmp.pb.g.vcf.gz.tbi"
    File local_fa_gvcf = "tmp.fa.g.vcf.gz"
    File local_fa_gvcf_index = "tmp.fa.g.vcf.gz.tbi"
    File local_mo_gvcf = "tmp.mo.g.vcf.gz"
    File local_mo_gvcf_index = "tmp.mo.g.vcf.gz.tbi"

    File header = "header.txt"

    File new_sample_id = "pb_id.txt"
    File new_father_id = "fa_id.txt"
    File new_mother_id = "mo_id.txt"

  }
}

## merges proband, father, mother gvcfs into trio gvcf
task merge_trio_gvcf {
  
  String sample_id

  File pb_gvcf
  File pb_idx
  File fa_gvcf
  File fa_idx
  File mo_gvcf
  File mo_idx

  File ref_fasta
  File ref_fasta_index

  String outfname = "${sample_id}.TRIO.g.vcf.gz"

  command {
    bcftools merge -g ${ref_fasta} ${pb_gvcf} ${fa_gvcf} ${mo_gvcf} -o ${outfname} -O z

    tabix -p vcf ${outfname}


  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
  }

  output {
    File out_gvcf = "${outfname}"
    File out_idx = "${outfname}.tbi"
  }

}

## splits vcf by chromosome
task split_gvcf {

  File gvcf # input gvcf
  File index # input gvcf index
  String outprefix = basename(gvcf, '.g.vcf.gz')
  File header # header from localize_paths step

  Int disk_size = 100 # start with 100G

  command {

    # split vcf by chromosome - use tabix -l to get all contig names from tabix index
    for i in $(tabix -l ${gvcf})
    do 
      (cat ${header}; tabix ${gvcf} $i)  > "${outprefix}.$i.vcf"
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
task call_denovos {
  File script
  
  File sample_id
  File father_id
  File mother_id

  File gvcf

  Float pb_min_vaf
  Int par_max_alt
  Int par_min_dp

  String shard
  
  String output_file = "${sample_id}.${shard}.denovo.txt"

  command {

    python ${script} -s ${read_string(sample_id)} -f ${read_string(father_id)} -m ${read_string(mother_id)} -g ${gvcf} -x ${pb_min_vaf} -y ${par_max_alt} -z ${par_min_dp} -o ${output_file}

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

