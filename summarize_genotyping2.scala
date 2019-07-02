import java.io.{File, PrintWriter}

import org.broadinstitute.gatk.queue.extensions.gatk.{LeftAlignAndTrimVariants, SelectVariants}
import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.util.Logging

import scala.io.Source


case class SampleFilename (path: File){
  val split_str = path.getName.split("[_\\.]")

  val case_samp = split_str.slice(1,3).mkString("_")
  val control_samp = split_str.slice(3,5).mkString("_")
}

class summarize_genotyping extends QScript with Logging {
  qscript =>

  @Argument(doc = "VCFs to annotate, can be specifed multiple times", required = true)
  var input_vcfs: List[File] = Nil


  @Argument(doc = "Base directory", required = true)
  var base_dir: String = _

  @Argument(doc = "Number of threads", required = false)
  var num_threads: Int = 5

  @Argument(doc = "Memory Limit (GB)", required = false)
  var mem_limit: Int = 20

  @Argument(doc = "Directory to outfile files", required = false)
  var out_base: String = ""

  @Argument(doc = "Ensembl tools path")
  var ens_tools_path: String = _

  @Argument(doc = "Path to Rscript", required = false)
  var rscript_path: String = "Rscript"

  @Argument(doc="summarization R script")
  var summ_script: File = "beataml_genomics/format_vcf.R"

  @Argument(doc="type of annotation run", required = false)
  var annotation_type: String = "release"


  class Vcf2MafRunner extends CommandLineFunction {

    @Input(doc = "Input VCF file") var vcf_input: File = _
    @Output(doc = "Output VEP file") var vep_output: File = _

    @Argument(doc = "Output MAF file", required = true)
    var maf_output: String = _

    @Argument(doc = "tumor sample name", required = true)
    var tumor_id: String = _

    @Argument(doc = "normal sample name", required = false)
    var norm_id: String = _

    //basic usage: perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --tumor-id WD1309 --normal-id NB1308

    this.memoryLimit = qscript.mem_limit

    def commandLine = required("perl") + required(ens_tools_path.getParent() + "/vcf2maf/vcf2maf.pl") +
      required("--ref-fasta", qscript.base_dir + "vep/vep/current.fasta") +
      required("--vep-path", ens_tools_path + "/scripts/variant_effect_predictor/") +
      required("--vep-data", qscript.base_dir + "vep/vep") +
      required("--input-vcf", vcf_input) +
      required("--output-maf", maf_output) +
      required("--tumor-id", tumor_id) +
      optional("--normal-id", norm_id)

  }

  class VarscanVcf2MafRunner extends Vcf2MafRunner {

    //for varscan: perl vcf2maf.pl --input-vcf data/test_varscan.vcf --output-maf data/test_varscan.maf --tumor-id WD1309 --normal-id NB1308 --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

    override def commandLine = super.commandLine +
      required("--vcf-tumor-id", "TUMOR") +
      required("--vcf-normal-id", "NORMAL")

  }


  class Rscript extends CommandLineFunction {

    @Input(doc = "")
    var inp_files: List[File] = _

    @Argument(doc="")
    var addl_args: List[String] = _

    @Output(doc = "")
    var out_files: List[File] = _

    var source_code: String = ""

    override def commandLine = required("echo", source_code) + required("|", escape = false) +
      required(qscript.rscript_path) + required("--vanilla") + required("-", escape = false) + repeat(inp_files) + repeat(addl_args)

  }


  class MutectToTable extends Rscript {

    source_code =
      """
        |
        |  annot.cols <- c("consequence", "impact", "symbol", "gene", "feature_type", "feature", "biotype", "exon", "intron", "hgvsc", "hgvsp", "cdna_position", "cds_position", "protein_position", "amino_acids",
        |                  "codons", "existing_variation", "allele_num", "distance", "strand", "pick", "variant_class", "symbol_source", "hgnc_id", "canonical", "refseq", "sift", "polyphen", "exac_af", "filter",
        |                  "variant_classification", "hgvsp_short")
        |
        |  args <- commandArgs(trailingOnly = T)
        |
        |  source(args[1])
        |
        |  x <- args[2]
        |
        |  split.samps <- strsplit(basename(x), "[_\\.]")[[1]]
        |
        |  include.vec <- setNames(c("samples", "normal"), c(paste(split.samps[2:3], collapse="_"), paste(split.samps[4:5], collapse="_")))
        |
        |  both.sample <- regmatches(x, gregexpr("\\d{2}-\\d{5}", x))[[1]]
        |
        |  is.to <- 0
        |
        |  if (length(both.sample) == 1){
        |    is.to <- 1
        |  }
        |
        |  if (args[4] == "yes"){
        |
        |  temp <- suppressWarnings(vcf.to.dt(x, samples.to.include=include.vec, geno.to.include=ug.geno.func, output.type="full", parse.csq=T))
        |
        |  names(temp) <- tolower(names(temp))
        |
        |  annot.temp <- temp[,annot.cols, with=F]
        |
        |  for(i in colnames(annot.temp)){
        |    if (any(class(annot.temp[[i]]) %in% c("factor", "logical"))){
        |      annot.temp[[i]] <- as.character(annot.temp[[i]])
        |    }
        |  }
        |
        |  res.temp <- temp[,.(seqnames=as.character(chromosome), pos_start=start_position, pos_end=end_position, ref=ref, alt=alt, samples=both.sample[1], allele_count=samples_genotype, total_reads=samples_total_reads, allele_reads=samples_allele_reads,
        |                      normal_total_reads, normal_allele_reads)]
        |
        |  mut.dt <- cbind(res.temp, annot.temp)
        |
        | }else{
        |
        |   temp <- suppressWarnings(vcf.to.dt(x, samples.to.include=include.vec, geno.to.include=ug.geno.func, output.type="full", parse.csq=F))
        |   names(temp) <- tolower(names(temp))
        |
        |   mut.dt <- temp[,.(seqnames=as.character(chromosome), pos_start=start_position, pos_end=end_position, ref=ref, alt=alt, samples=both.sample[1], allele_count=samples_genotype, total_reads=samples_total_reads, allele_reads=samples_allele_reads,
        |                      normal_total_reads, normal_allele_reads)]
        |
        | }
        |
        |  mut.dt[,`:=`(genotyper="mutect", t_vaf=allele_reads/total_reads, n_vaf=ifelse(normal_total_reads == 0, 0, normal_allele_reads/normal_total_reads), tumor_only=is.to)]
        |
        |  out.file <- file.path(args[3], sub("\\.vcf", ".RData", basename(x)))
        |
        |  save(mut.dt, file=out.file)
        |
        |
      """.stripMargin



  }

  class VarscanToTable extends Rscript {

    source_code =
      """
        |
        | annot.cols <- c("consequence", "impact", "symbol", "gene", "feature_type", "feature", "biotype", "exon", "intron", "hgvsc", "hgvsp", "cdna_position", "cds_position", "protein_position", "amino_acids",
        |                  "codons", "existing_variation", "allele_num", "distance", "strand", "pick", "variant_class", "symbol_source", "hgnc_id", "canonical", "refseq", "sift", "polyphen", "exac_af", "filter",
        |                  "variant_classification", "hgvsp_short")
        |
        |
        |  args <- commandArgs(trailingOnly = T)
        |
        |  source(args[1])
        |
        |  x <- args[2]
        |
        |  include.vec <- setNames(c("samples", "normal"), c("TUMOR", "NORMAL"))
        |
        |  both.sample <- regmatches(x, gregexpr("\\d{2}-\\d{5}", x))[[1]]
        |
        |  if (args[4] == "yes"){
        |
        |  temp <- suppressWarnings(vcf.to.dt(x, samples.to.include=include.vec, geno.to.include=varscan.gen.func, output.type="full", parse.csq=T))
        |
        |  names(temp) <- tolower(names(temp))
        |
        |  annot.temp <- temp[,annot.cols, with=F]
        |
        |  for(i in colnames(annot.temp)){
        |    if (any(class(annot.temp[[i]]) %in% c("factor", "logical"))){
        |      annot.temp[[i]] <- as.character(annot.temp[[i]])
        |    }
        |  }
        |
        |  res.temp <- temp[,.(seqnames=as.character(chromosome), pos_start=start_position, pos_end=end_position, ref=ref, alt=alt, samples=both.sample[1], allele_count=samples_genotype, total_reads=samples_total_reads, allele_reads=samples_allele_reads,
        |                      normal_total_reads, normal_allele_reads)]
        |
        |  mut.dt <- cbind(res.temp, annot.temp)
        |
        | }else{
        |
        |   temp <- suppressWarnings(vcf.to.dt(x, samples.to.include=include.vec, geno.to.include=varscan.gen.func, output.type="full", parse.csq=F))
        |   names(temp) <- tolower(names(temp))
        |
        |   mut.dt <- temp[,.(seqnames=as.character(chromosome), pos_start=start_position, pos_end=end_position, ref=ref, alt=alt, samples=both.sample[1], allele_count=samples_genotype, total_reads=samples_total_reads, allele_reads=samples_allele_reads,
        |                      normal_total_reads, normal_allele_reads)]
        | }
        |
        |  is.to <- 0
        |
        |  if (length(both.sample) == 1){
        |    is.to <- 1
        |    mut.dt[,`:=`(normal_total_reads=NULL, normal_allele_reads=NULL )]
        |    mut.dt[,`:=`(normal_total_reads=0L, normal_allele_reads=0L )]
        |  }
        |
        |  mut.dt[,`:=`(genotyper="varscan", t_vaf=allele_reads/total_reads, n_vaf=ifelse(normal_total_reads == 0, 0, normal_allele_reads/normal_total_reads), tumor_only=is.to)]
        |
        |  out.file <- file.path(args[3], sub("\\.vcf", ".RData", basename(x)))
        |
        |  save(mut.dt, file=out.file)
        |
        |
      """.stripMargin

  }

  def script() {


    val out_base_file = new File(qscript.out_base)

    if (out_base_file.exists() == false) {
      out_base_file.mkdirs()
    }

    val maf_out = new File(out_base_file.getPath() + "/mafs")

    val mutect_out = new File(out_base_file.getPath() + "/mutect_processed_tables")
    val varscan_out = new File(out_base_file.getPath() + "/varscan_processed_tables")

    var should_annotate = "yes"

    for (i <- input_vcfs) {

      //determine the sample names

      val sample_names = SampleFilename(i)

      var process_vep = new Rscript
      var vep_out_folder = mutect_out

      var mut_vep_out = new File(i.getPath)

      if (Set("mafs", "both") contains qscript.annotation_type){

        if (maf_out.exists() == false){
          maf_out.mkdirs()
        }

        var mut_maf = new Vcf2MafRunner

        //vcf2maf runs vep its own way now...

        if ((i.getName contains "varscan")) {
          mut_maf = new VarscanVcf2MafRunner
          process_vep = new VarscanToTable
          vep_out_folder = varscan_out

        } else if (i.getName contains "mutect" == false) {

          throw new QException("ERROR: Unknown type of somatic file")

        }else{
          process_vep = new MutectToTable
        }

        val mut_maf_out = maf_out + "/" + i.getName.replace(".vcf", ".maf")

        mut_vep_out = new File(i.getPath.replace(".vcf", ".vep.vcf"))

        mut_maf.jobOutputFile = mut_maf_out + ".out"
        mut_maf.jobErrorFile = mut_maf_out + ".err"

        mut_maf.vcf_input = i
        mut_maf.maf_output = new File(mut_maf_out)
        mut_maf.vep_output = mut_vep_out

        mut_maf.norm_id = sample_names.control_samp
        mut_maf.tumor_id = sample_names.case_samp
        mut_maf.wallTime = 240
        mut_maf.memoryLimit = qscript.mem_limit
        mut_maf.nCoresRequest = qscript.num_threads
        add(mut_maf)
      }else if (Set("csq_summary", "summary_only") contains qscript.annotation_type){

        if ((i.getName contains "varscan")) {

          if (varscan_out.exists() == false){
            varscan_out.mkdirs()
          }

          process_vep = new VarscanToTable
          vep_out_folder = varscan_out

        } else if (i.getName contains "mutect" == false) {

          throw new QException("ERROR: Unknown type of somatic file")

        }else{
          process_vep = new MutectToTable

          if (mutect_out.exists() == false){
            mutect_out.mkdirs()
          }
        }

      }else{
        throw new QException("expect annotation_type to be one of 'mafs', 'both', 'csq_summary', 'summary_only'")
      }

      if (Set("csq_summary", "summary_only", "both") contains qscript.annotation_type){
        val process_vep_out = vep_out_folder + "/" + mut_vep_out.getName.replace(".vcf", ".RData")

        process_vep.jobOutputFile = process_vep_out + ".out"
        process_vep.jobErrorFile = process_vep_out + ".err"

        //inputs: script, input file and output directory

        if (qscript.annotation_type == "summary_only"){
          should_annotate = "no"
        }

        process_vep.inp_files = List(qscript.summ_script, mut_vep_out, vep_out_folder)
        process_vep.addl_args = List(should_annotate)
        process_vep.out_files = List(process_vep_out)

        process_vep.wallTime = 60
        process_vep.memoryLimit = 10

        add(process_vep)
      }



    }

  }

}

