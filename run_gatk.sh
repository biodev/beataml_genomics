#!/bin/bash

cur_host=`hostname`
JOB_RUNNER="-jobRunner Condor"
RSCRIPT="Rscript"
GW_PATH="beataml_genomics"
JAVA_CMD="java"

collect_vals (){
    
    input_files=""
    
    for i in ${@:2}
    do
        input_files=${input_files}" --$1 "${i}
    done
    
    echo $input_files
}

if ( `echo "$cur_host" | grep -q "exa"` )
then
    CLUST_BASE="/home/exacloud/lustre1/NGSdev/resources/gatk_v3.3"
    JAR_BASE="/home/exacloud/lustre1/NGSdev/"
    JOB_RUNNER="-jobRunner Slurm"
    ENSEMBL_TOOLS=${CLUST_BASE}/vep/ensembl-tools
    SUBREAD_VERSION="subread-1.5.0-p2"
    PINDEL="/home/exacloud/lustre1/NGSdev/resources/gatk_v3.3/programs/pindel/"

    #add tabix and bgzip to path for vep
    export PATH=/usr/lib/jvm/java-1.7.0/bin/:/home/exacloud/lustre1/NGSdev/resources/gatk_v3.3/programs/htslib/:/home/exacloud/lustre1/NGSdev/resources/gatk_v3.3/programs/samtools-1.1/:$PATH

    #to set up the custom R library paths
    export LD_LIBRARY_PATH=/home/exacloud/lustre1/NGSdev/resources/gatk_v3.3/programs/SlurmRlibs/system_libs/:/home/exacloud/lustre1/NGSdev/resources/gatk_v3.3/programs/slurm-libs/lib/:$LD_LIBRARY_PATH
    export R_LIBS_USER=/home/exacloud/lustre1/NGSdev/resources/gatk_v3.3/programs/SlurmRlibs/x86_64-redhat-linux-gnu-library/3.4/

    PERLFAIL="yes"

else
    echo "ERROR: Unknown host, please modify script"
    exit 1
fi

echo "Running on $cur_host"
echo "Using JobRunner: $JOB_RUNNER"

export PERL5LIB=${CLUST_BASE}/vep/:$PERL5LIB

export PATH=`ls -d ${CLUST_BASE}/programs/${SUBREAD_VERSION}*/bin/`:$PATH

if [ $# -lt 1 ]
then
    echo "Need to supply the name of the command to run and its corresponding arguments if any (e.g. ./run_gatk.sh mutect_tumor_only out_dir test1.bam test2.bam)"
	exit 1  
fi

cur_scala=""

if [ -f "$1.scala" ]
then
    cur_scala="$1.scala"
elif [ -f "${GW_PATH}/$1.scala" ]
then
    cur_scala="${GW_PATH}/$1.scala"
else
    echo "ERROR: Cannot find scala file $1.scala"
    exit 1
    
fi

echo "Found script: $cur_scala"

run_cmd="${JAVA_CMD} -Djava.io.tmpdir=sum_geno_tmp \
         -jar ${JAR_BASE}/Queue.jar \
         -S ${cur_scala} \
	 -retry 3 \
         -l DEBUG --disableJobReport \
         ${JOB_RUNNER} \
         --graphviz ${cur_scala}.dot \
         --base_dir ${CLUST_BASE}/"

#base arguments that are common to the specified sets

if [ $1 == "subjunc_rnaseq" ]
then
    
    args_bad="no"
    should_trim=""
    
    if [ $# -lt 5 ]
    then
        args_bad="yes"
        
    elif (( $# == 5 )) && [ "$5" == "trim" ]
    then
        args_bad="yes"
    
    elif (( $# > 5 )) && [ "$5" == "trim" ]
    then
        input_files=`collect_vals inp_fastqs ${@:6}`
        should_trim="--should_trim"
    else
        input_files=`collect_vals inp_fastqs ${@:5}`
    fi
    
    if [ "$args_bad" == "yes" ]
    then
        echo "$1: ERROR: Need to supply the output directory, read orientation, strand specificity, optionally [trim] and desired fastq files"
        echo "$1: EXAMPLE1: ./run_gatk.sh subjunc_rnaseq test_out fr 2 /Users/bottomly/Desktop/github/gatk_tests/test_rnaseq/151210_D00735_0076_AC8DKVANXX/RNA150311JT/*.fastq.gz"
        echo "$1: EXAMPLE2: ./run_gatk.sh subjunc_rnaseq test_out fr 2 trim /Users/bottomly/Desktop/github/gatk_tests/test_rnaseq/151210_D00735_0076_AC8DKVANXX/RNA150311JT/*.fastq.gz"
        exit 1  
    fi
    
    
    run_cmd="$run_cmd \
             --out_base_dir $2 \
             --num_threads 4 \
             --mem_limit 10 \
             --subjunc_indexes ${CLUST_BASE}/pipeline/${SUBREAD_VERSION}/human_g1k_v37 \
             --annotation_file ${CLUST_BASE}/pipeline/Homo_sapiens.GRCh37.75.gtf \
             --subjunc_readsorientation $3 \
             --fc_orientation $4 \
             ${should_trim} \
             ${input_files}"

elif [ $1 == "filter_varscan" ]
then

    if [ $# -lt 3 ]
    then
        echo "$1: ERROR: Need to provide the output directory and the VCF files to process"
        echo "$1: EXAMPLE1: ./run_gatk.sh filter_varscan snp_out test_muts/*snp.vcf"
        echo "$1: EXAMPLE2: ./run_gatk.sh filter_varscan indel_out test_muts/*.flag.fixed.la.vcf"
        exit 1
    fi
    
    if  ! perl -V | grep 'version 14 subversion 2' &> /dev/null;
    then
        echo "$1 ERROR: Need to use perl 5.14.2.  Easiest way is utilize perlbrew: perlbrew use 5.14.2"
        
        if [ $PERLFAIL == "yes" ]
        then
            exit 1
        fi
    fi

    input_files=`collect_vals input_vcfs ${@:3}`
    
    run_cmd="$run_cmd \
            --out_base ${2} \
            --mem_limit 10 \
            --ens_tools_path ${ENSEMBL_TOOLS} \ 
            ${input_files}"

elif [ $1 == "tumor_only_varscan" ]
then

    if [ $# -lt 5 ]
    then
        echo "$1: ERROR: Need to provide the output directory, type of run [dna, rna] , a file of mpileup paths processed to have blank entries for a normal and a file of corresponding tumor-only processed BAM paths"
        echo "$1: EXAMPLE1: ./run_gatk.sh tumor_only_varscan test_out dna mpup.txt bam.txt"
        exit 1
    fi
    
    if  ! perl -V | grep 'version 14 subversion 2' &> /dev/null;
    then
        echo "$1 ERROR: Need to use perl 5.14.2.  Easiest way is utilize perlbrew: perlbrew use 5.14.2"
        
        if [ $PERLFAIL == "yes" ]
        then
            exit 1
        fi
    fi
    
    run_cmd="$run_cmd \
            --out_base ${2} \
            --mem_limit 10 \
            --ens_tools_path ${ENSEMBL_TOOLS} \
            --filter_type ${3} \
            --mpup_file ${4} \
            --bam_file ${5}"


elif [ $1 == "tophat_fusion_search" ]
then
    
    if [ $# -lt 4 ]
    then
        echo "$1: ERROR: Need to provide the output directory, the number of samples per tophat run and the 'processed' FASTQ or 'tophat_SampleX/fusion.out' files to run"
        echo "$1: EXAMPLE1: ./run_gatk.sh tophat_fusion_search test_out 5 processed_fastq/*.fastq.gz"
        echo "$1: EXAMPLE2: ./run_gatk.sh tophat_fusion_search test_out 5 tophat_alignments/tophat_Sample*/fusions.out"
        exit 1
    fi
    
    #make sure bowtie, tophat, blast are added to the path
    export PATH=${CLUST_BASE}/programs/bowtie-1.1.1/:${CLUST_BASE}/programs/ncbi-blast-2.2.30+/bin/:${CLUST_BASE}/programs/samtools-1.1/:${CLUST_BASE}/programs/tophat-2.0.14-build/bin/:$PATH
    
    input_files=`collect_vals inp_fastqs ${@:4}`
    
    run_cmd="$run_cmd \
            --out_base_dir ${2} \
            --tf_num ${3} \
            --mem_limit 30 \
            --num_threads 5 \
            --rscript_path ${RSCRIPT} \
            ${input_files}"      
    
elif [ $1 == "pindel" ]
then
    
    if [ $# -lt 6 ]
    then
        echo "$1: ERROR: Need to provide an insert size file, chromosome, pairs file, output directory and paths to BAM files"
        echo "$1: EXAMPLE1: ./run_gatk.sh pindel insert_size_summary.txt 13 pairs.txt test_out /home/exacloud/lustre1/NGSdev/dnaseq_out/BeatAML_Release1/final/somatic/inputs/*.bam"
        exit 1
    fi
    
    if  ! perl -V | grep 'version 14 subversion 2' &> /dev/null;
    then
        echo "$1 ERROR: Need to use perl 5.14.2.  Easiest way is utilize perlbrew: perlbrew use 5.14.2"
        
        if [ $PERLFAIL == "yes" ]
        then
            exit 1
        fi
    fi
    
    input_files=`collect_vals bam_files ${@:6}`
    
    run_cmd="$run_cmd \
            --mem_limit 10 \
            --num_threads 5 \
            --ins_file ${2} \
            --pindel_chr ${3} \
            --tn_file ${4} \
            --out_dir ${5} \
            --rscript_path ${RSCRIPT} \
            --pindel_path ${PINDEL} \
            --ens_tools_path ${ENSEMBL_TOOLS} \
            ${input_files}"

elif [ $1 == "gatk_lane_level" ]
then
    
    if [ "$#" -lt 3 ]
    then
        echo "$1: ERROR: Need to supply the 'CaptureGroup' directory, and the relative path to the genomic intervals (in bed format) to be used by GATK"
        echo "$1: EXAMPLE1: ./run_gatk.sh gatk_lane_level /home/exacloud/lustre1/NGSdev/dnaseq_out/SampleGroup20/CaptureGroup1 bundle_2_8/nextera_v1_2_for_gatk.bed"
        echo "$1: EXAMPLE2: ./run_gatk.sh gatk_lane_level /home/exacloud/lustre1/NGSdev/dnaseq_out/SampleGroup20/CaptureGroup1 bundle_2_8/Nimblegen_SeqCap_EZ_v3.bed"
        exit 1
    fi


    alignments_basedir=$2
    interval_bed=$3
    
    run_cmd="$run_cmd \
        --input_align_basedir ${2} \
        --interval_coords ${3} \
        --scatter_count 5 \
        --mem_limit 20 \
        --num_threads 4"

elif [ $1 == "summarize_align_lane" ]
then
    
    if [ $# -lt 4 ]
    then
      echo "$1: ERROR: Need to supply the 'FlowCell' directory, the overall output directory and the SampleGroup label"
      echo "$1: EXAMPLE1: ./run_gatk.sh summarize_align_lane /home/exacloud/lustre1/NGSdev/seqcap/SampleGroup20/FlowCell2 /home/exacloud/lustre1/NGSdev/dnaseq_out BeatAML10"
      exit 1
    fi
    
    run_cmd="$run_cmd \
        --num_threads 4 \
        --flowcell_path ${2} \
        --out_base_dir ${3} \
        --id_pref ${4}"


elif [ $1 == "gatk_picard_align_stats" ]
then
    
    if [ $# -lt 3 ]
    then
      echo "$1: ERROR: Need to supply the output directory and paths to the BAM files to process"
      echo "$1: EXAMPLE1: ./run_gatk.sh gatk_picard_align_stats lane_level_test /home/exacloud/lustre1/NGSdev/dnaseq_out/SampleGroup21/FlowCell2/alignments/*.bam"
      echo "$1: EXAMPLE2: ./run_gatk.sh gatk_picard_align_stats sample_level_test /home/exacloud/lustre1/NGSdev/dnaseq_out/SampleGroup21/SomaticGenotyping1/paired_genotypes/*.dedup.realign.bam"
      exit 1
    fi
    
    input_files=`collect_vals align_files ${@:3}`
    
    run_cmd="$run_cmd --mem_limit 30 \
            --rscript_path ${RSCRIPT} \
            --out_dir ${2} \
            ${input_files}"

elif [ $1 == "tumor_only_genotyping" ]
then
    
    if (( $# >= 5 )) && [ "$3" == "skip" ]
    then
        skip_str="--skip_realignment skip"
    elif (( $# >= 5 )) && [ "$3" == "noskip" ]
    then
        skip_str="--skip_realignment noskip"
    elif (( $# >= 5 )) && [ "$3" == "nogeno" ]
    then
        skip_str="--skip_realignment nogeno"
    else
        echo "$1: ERROR: Need to supply the output directory, whether to skip realignment/genotyping [skip,noskip,nogeno], the basename of GATK BED file to use and the path to one or more BAM files"
        echo "$1: EXAMPLE1: ./run_gatk.sh tumor_only_genotyping TumorOnly noskip nextera_v1_2_for_gatk.bed ../test_bams/Sample_13-005*.dedup.bam"
        exit 1
    fi
    
    input_files=`collect_vals inp_bams ${@:5}`
    
    run_cmd="$run_cmd \
             --out_dir $2 \
             --mem_limit 10 \
             --interval_file $4 \
             $skip_str \
             ${input_files}"


elif [ $1 == "somatic_genotyping2" ]
then
    
    if [ $# -lt 6 ]
    then
        
        echo "$1: ERROR: Need to supply the input directory, output directory, One of [realign, dont_realign], the basename of GATK BED file to use and the path to a file containing the tumor/normal assignments"
        echo "$1: EXAMPLE: ./run_gatk.sh somatic_genotyping2 SomaticGenotyping1/paired_alignments test_out realign nextera_v1_2_for_gatk.bed tn_file.txt"
        exit 1
        
    fi
    
    realign_str=""
    bam_suffix=".dedup.bam"
    
    if [ $4 == "dont_realign" ]
    then
        realign_str="--dont_realign"
        bam_suffix=".realign.bam"
    elif [ $4 != "realign" ]
    then
        echo "Only a value of [realign, dont_realign] is valid"
    fi
    
    
    run_cmd="$run_cmd \
    --mem_limit 30 \
    --scatter_count 1 \
    --num_threads 7 \
    --genotyper mutect \
    --genotyper varscan \
    --alignment_path $2 \
    $realign_str \
    --bam_suffix $bam_suffix \
    --out_base $3 \
    --interval_file $5 \
    --tn_file $6"


elif [ $1 == "summarize_genotyping2" ]
then
    if [ $# -lt 4 ]
    then
        echo "$1 ERROR: Need to supply the output directory, annotation type for the VCFs [mafs, both, csq_summary, summary_only] and the mutect and/or varscan VCF files"
        echo "$1 EXAMPLE1 (release): ./run_gatk.sh summarize_genotyping2 test_out mafs single_genotypes/*filter.vcf"
        echo "$1 EXAMPLE2 (process previous runs for vizome): ./run_gatk.sh summarize_genotyping2 test_out csq_summary single_genotypes/*filter.vep.vcf"
        echo "$1 EXAMPLE3 (don't include annotation): ./run_gatk.sh summarize_genotyping2 test_out summary_only single_genotypes/*filter.vcf"
        exit 1
    fi
    
     #not doing this directly due to some strangeness with perlbrew and exacloud...
    
    if  ! perl -V | grep 'version 14 subversion 2' &> /dev/null;
    then
        echo "$1 ERROR: Need to use perl 5.14.2.  Easiest way is utilize perlbrew: perlbrew use 5.14.2"
        
        if [ $PERLFAIL == "yes" ]
        then
            exit 1
        fi
    fi
    
    input_files=`collect_vals input_vcfs ${@:4}`
    
    cur_script=""

    if [ -f "format_vcf.R" ]
    then
        cur_script="format_vcf.R"
    elif [ -f "${GW_PATH}/format_vcf.R" ]
    then
        cur_script="${GW_PATH}/format_vcf.R"
    else
        echo "ERROR: Cannot find annotation script file 'format_vcf.R'"
        exit 1
        
    fi
    
    run_cmd="$run_cmd \
             --out_base ${2} \
             --annotation_type ${3} \
             --ens_tools_path ${ENSEMBL_TOOLS} \
             --rscript_path ${RSCRIPT} \
             --summ_script ${cur_script} \
             ${input_files}"
    

else
  echo "ERROR: Unknown command: $1"
  exit 1
fi

#additional arguments that are specific to each program


fin_cmd="$run_cmd \
     -run"

echo $fin_cmd  
eval $fin_cmd
 
