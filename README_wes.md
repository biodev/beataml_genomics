# BeatAML WES protocol

## NOTICE: This protocol is provided for informational purposes only.  Although, originally based off of GATK and others' best practices at the time of development (~4-5 years ago), this may not be the case now.  Any re-genotyping efforts should make use of current tools and best practices.

**Those interested in the processed BeatAML data should retrieve the harmonized versions from the GDC (coming soon)**

## Environment setup

Files used for GATK preprocessing and genotyping were originally downloaded from ftp://ftp.broadinstitute.org/bundle/2.8/b37/ on June 2nd 2015.

Below is the expected setup for the programs and metadata files:

* resource_dir
  * bundle_2_8
    * 1000G_phase1.indels.b37.vcf
    * dbsnp_138.b37.vcf
    * human_g1k_v37.fasta (Including .dict and BWA input files)
    * Mills_and_1000G_gold_standard.indels.b37.vcf
    * Seperate files specific to the capture library:
      * nextera_v1_2_for_gatk.bed--A BED file containing the Nextera v1.2 capture intervals 
      * nextera_v1_2_for_gatk_UNION_500.bed A version of the Nextera BED file containing the collapsed intervals with a padding of 500 bases
        * A BGZIP'd version and TABIX index should be available as well.
  * mutect_bundle
    * b37_cosmic_v54_120711.vcf
    * dbsnp_132_b37.leftAligned.vcf
  * programs
    * htslib/ (v1.6)
    * samtools/ (v1.1)
    * pindel/ (v0.2.5b9)
    * bwa-0.7.10/
    * Trimmomatic-0.36/
    * picard-tools-1.124/
    * GenomeAnalysisTK.jar (v3.3-3-g3e87c07)
    * mutect-src/ (v1.1.7)
    * VarScan.v2.4.1p1.jar (After application of changes in varscan.patch)
    * other_varscan/v2_4_1/bam-readcount-bin/bin/bam-readcount (v0.7.4)
    * bcftools/ (v1.6)
    * vep/ (see below)

The Variant Effect Predictor (VEP) is a more complex setup, roughly this script expects a top-level vep directory with the corresponding ensembl-tools/ (release/83) and vcf2maf/ (v1.6.6).  The relevant installed cache was `homo_sapiens_vep_83_GRCh37.tar.gz` and utilizes the ExAC plugin and corresponding VCF file per `vcf2maf`.  The `vep/vep/current.fasta` file referenced in the scala scripts is a symlink to the corresponding VEP reference file. 

This protocol utilizes GATK Queue to manage and parallelize the workflow.  The main `Queue.jar` file (v3.3-3-g3e87c07) can reside in `resource_dir/programs` above or elsewhere.

### 1. Stage Fastq data

The following is assumed about the Fastq data:

1. The reads are 100 bases in length and paired-end from an Illumina sequencer
2. The Nextera v1.2 whole exome sequencing protocol was used.

Create an appropriate folder indicating the samples that are to be run for instance SampleGroup10/FlowCell5.  In this directory the fastq data should be placed with the top level directory set as the flowcell identifier (e.g. SampleGroup10/FlowCell5/160607_D00735_0122_AC9561ANXX).

The below workflow further assumes that the fastq files are in a subdirectory prefixed by 'DNA' and have the read 1 files named like `^DNA.*_R1_.*fastq.gz$`.  Samples are assumed to be named like '\d{2}-\d{5}'.  There should be one file per lane and read number.  Lanes should be indicated by _L00[1-8]_.

### 2. Summarize the FastQC results provided by the core.  

From within the SampleGroup folder, an HTML summary of the Fastqc files provided by the sequencing core can be produced.

```
bash$ R
> source("/path/to/beataml_genomics/baml_helpers.R")
> make.fastqc.report("FlowCell5/160607_D00735_0122_AC9561ANXX/FastQC/", report.name="SG10_FC5_FastQC.html")
```

This generates an HTML document called `SG10_FC5_FastQC.html` which contains a table of all the FastQC plots and an overall summary.

### 3. FASTQ Summarization and Alignment

First the lane-level FASTQ files are first trimmed with the resulting files placed in the `processed_fastq` folder created as a subdirectory in the input FlowCell directory. These reads are then aligned using BWA with the resulting SAM files sorted and converted to corresponding BAM files. 

```
(time sh run_gatk.sh summarize_align_lane SampleGroupX/FlowCellY dnaseq_out BeatAMLX) >> SGX_FLY_summarize.log 2>&1
```

From this point forward the paths are relative to the `dnaseq_out` folder.  So for instance, the resulting BAM files are output in `SampleGroupX/FlowCellY/alignments/` in the `dnaseq_out/` folder. 


### 4. Collect Lane-level Picard Alignment Metrics

Summary tables of the results of Picard's `CollectAlignmentSummaryMetrics` and `CollectInsertSizeMetrics` are produced in the specified output directory.

```
(time sh run_gatk.sh gatk_picard_align_stats FlowCellY/statistics FlowCellY/alignments/*.bam) >> SGX_FLY_stats.log 2>&1
```

### 5. Set up samples into CaptureGroups folders

The samples are divided into `CaptureGroup` folders assuming they have the expected number of lanes.  CaptureGroups are numbered sequentially and samples are skipped if already assigned to a CaptureGroup.  The `CaptureGroupX/alignments` folder contains symbolic links to the corresponding BAM files in the `SampleGroupX/FlowCellY/alignments/` folder(s).

The `sample_group_X_pairs.csv` file below describes the patient samples and should be formated as: `patient_id, AML, Normal` where AML and Normal refer to the tumor and normal samples respectively and patient_id is an arbitrary grouping identifier.

```
bash$ R
> source("/path/to/beataml_genomics/baml_helpers.R")
> pair.file <- "sample_group_X_pairs.csv"
> setup.capturegroups(directory=".", pair.file, expected.lanes=c(AML=5, Normal=3))
```

### 6. GATK Preprocessing

Carries out the GATK pre-processing for a given CaptureGroup.  This involves lane-level duplicate marking, indel realignment and base quality score recalibration.  The resulting lane-level BAMs are then combined into sample-level BAM files output to the `CaptureGroupX/sample_alignments` folder.

```
bash$ (time sh run_gatk.sh gatk_lane_level CaptureGroupX bundle_2_8/nextera_v1_2_for_gatk.bed ) >> SGX_CG1_preprocessing.log 2>&1
```

### 7. Categorize the Pre-processed Data for Genotyping

Using the pairing `.csv` file from step 5, divide the sample-level BAM files into Tumor/Normal pairs for somatic genotyping.  Similar to the formation of the CaptureGroups, the SomaticGenotyping folders are numbered sequentially skipping previous samples.  Each SomaticGenotypingX folder contains a `paired_aligments` subdirectory with symbolic links to the BAM files in the correspoding CaptureGroup. The `setup.sgs` function also divides samples into `TumorOnly` and `NormalOnly` if necessary.

```
bash$ R
> setup.sgs(directory=".", pair.file="sample_group_X_pairs.csv")
```

### 8a. Initial Somatic Genotyping

For each tumor/normal pair, another round of indel realignment is performed followed by somatic genotyping using MuTect and VarScan2.  In addition to the original and intermediate VCF files, the final filtered VCF files are placed in the `SomaticGenotypingX/paired_genotypes` folder.

```
bash$ (time sh run_gatk.sh somatic_genotyping2 SomaticGenotyping1/paired_alignments SomaticGenotyping1/ realign nextera_v1_2_for_gatk.bed sample_group_X_pairs.csv ) >> SGX_SomaticGenotyping1.log 2>&1
```

### 8b. Initial Tumor Only Genotyping

Similar to the SomaticGenotyping, each tumor-only sample underwent indel realignment and was genotyped using MuTect and VarScan2 with the final and all intermediate VCFs placed in the TumorOnlyX/single_genotypes folder.

```
bash$ (time sh run_gatk.sh tumor_only_genotyping TumorOnly1 noskip nextera_v1_2_for_gatk.bed TumorOnly1/single_alignments/*.dedup.bam) >> SGX_tumor_only.log 2>&1
```

### 9. Generating MAFs

For the MAF, we use the `vcf2maf` program to create MAFs for release.  The example below is for tumor-only runs but applies
similarly to somatic runs.  It will create a subfolder 'mafs' within the 'TumorOnly1' folder (in this case) containing the MAF files.

NOTE: It will also create `*.filt.vep.vcf` files in the original input directory.

```
bash$ (time sh run_gatk.sh summarize_genotyping2 TumorOnly1 mafs TumorOnly1/single_genotypes/*.filter.vcf) 2>&1 | tee tumor_only_annotation.log
```

### 10 BeatAML manuscript variant parameters

To use the VarScan2 parameterization as described in the flagship BeatAML manuscript, the following steps can be taken.

### 10a. Refiltering of Somatic Genotyping for BeatAML Manuscript

Given the 'unfiltered' versions of the VarScan2 SNV and indel calls, apply alternative filtering, annotate them and output them to the `vcf` and `maf` subdirectories of the supplied
output directory `varscan_out` in the below example.

```
bash$ ls SomaticGenotyping1/paired_genotypes/*snp.vcf > som_vcfs.txt
bash$ ls SomaticGenotyping1/paired_genotypes/*.flag.fixed.la.vcf >> som_vcfs.txt
bash$ (time sh run_gatk.sh filter_varscan varscan_out "$(< som_vcfs.txt)" ) > filtered_varscan.log 2>&1
```

### 10b. Rerunning Tumor-only Varscan for BeatAML Manuscript

Re-call the tumor-only variants as opposed to simply re-filtering the somatic ones above.  In this case provide the previously generated mpileup and BAM files are provided to the script and annotated final VCFs and MAFs along with their intermediates are output to the `vcf` and `maf` subdirectories of the supplied output directory.

```
bash$ ls TumorOnly1/single_genotypes/*No_Skin.mpileup > to_mpups.txt
bash$ ls TumorOnly1/single_alignments/*.realign.bam > to_bams.txt
bash$ (time sh run_gatk.sh tumor_only_varscan to_var_out dna to_mpups.txt to_bams.txt) |& tee tumor_only_var.log
```

### 10c. Parse the annotated VCFs to create per-sample tabular input files.

An alternative tabular form of the variants can be generated (similar to that at <vizome.org>).  Here, use VCFs from the MuTect runs and the analysis versions of the VarScan2 runs.  They should be `*.filter.vep.vcf` files in both cases.  Note this will create seperate folders for each program.

```
bash$ (time sh run_gatk.sh summarize_genotyping2 release_test csq_summary inputs/*.filter.vep.vcf) 2>&1 | tee reannotation.log
```

Once this is done, the per-sample tabular files can then be `rbind'ed` together to form a single `data.frame` and further filtered or transformed for analysis.

### 11.  Pindel

In order to run pindel several files are needed:

1. A pairing CSV file as described above in section 5.
2. Summarized output from the Picard insert size metrics as output in the 'insert_size_summary.txt' file from `gatk_picard_align_stats` above.
3. Paths to all the BAM files referenced by the pairing file.  Simplest is probably to put the paths in a file, say 'current_bams.txt', and supply the file to the `run_gatk.sh` script as is shown below.

Pindel is run and the resulting VCF files are annotated using `vcf2maf` producing annotated VCF files and MAF files in the output directory.

```
bash$ (time sh run_gatk.sh pindel insert_size_summary.txt ALL sample_pairs.csv output_dir "$( < current_bams.txt)" ) 2>&1 | tee pindel.log 
```






