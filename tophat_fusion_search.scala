/**
 * Created by bottomly on 1/21/15.
  * $Version$
 */

import java.io.{File, PrintWriter}

import scala.util.matching.Regex
import scala.io._
import scala.collection.mutable._
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.util.Logging

case class SampleFile(line: File) {
  val data = line.getName.split("_")
  val sample = data(1)
  val file = line
  val isR1 = line.getName.matches(".*_R1_.*")

}

class tophat_fusion_search extends QScript with Logging {
  qscript =>

  @Argument(doc = "Memory Limit (GB)")
  var mem_limit: Int = _

  @Argument(doc = "Number of Threads")
  var num_threads: Int = _

  @Argument(doc="Resource Base directory")
  var base_dir: String = _

  @Argument(doc="fastq files")
  var inp_fastqs: List[File] = _

  @Argument(doc = "Output directory where TopHat results are to be placed")
  var out_base_dir : File = _

  @Argument(doc="Number samples per tophat fusion runs")
  var tf_num: Int = _

  @Argument(doc="Path to appropriate Rscript")
  var rscript_path: File = _

  // TopHat Class
  class TopHat extends CommandLineFunction {
    /*
      == command line ==

      tophat -o tophat_testData -p 8 --fusion-search --keep-fasta-order
      --bowtie1 --no-coverage-search -r 100 --mate-std-dev 80
      --max-intron-length 100000 --fusion-min-dist 100000
      --fusion-anchor-length 13 --fusion-ignore-chromosomes M
      /path/to/test_ref
      /path/to/reads_1.fq
      /path/to/reads_2.fq
    */

    //@Input(doc = "")
    @Argument(doc = "number of threads")
    var num_threads: Int = 8

    @Argument(doc = "mate inner distance")
    var mate_inner_distance: Int = 100

    @Argument(doc="output directory")
    var out_dir: File = _

    @Output(doc = "name of fusion file")
    var out_fusions: File = _

    @Input(doc = "path to bowtie indexes")
    var bowtie_base: File = _

    @Input(doc = "path to first fastq file")
    var fastq1: File = _

    @Input(doc = "path to second fastq file")
    var fastq2: File = _

    //@Output(doc="List with paths samples directories created by tophat")
    //var tophat_outputs : List[File] = _

    def commandLine = required("tophat") +
      required("-o", out_dir) +
      required("-p", num_threads) +
      required("--fusion-search") +
      required("--bowtie1") +
      required("--no-coverage-search") +
      required("-r", mate_inner_distance) +
      required("--mate-std-dev", 80) +
      required("--max-intron-length", 100000) +
      required("--fusion-min-dist", 100000) +
      required("--fusion-anchor-length", 13) +
      required("--fusion-ignore-chromosomes", "chrM") +
      required(bowtie_base) +
      required(fastq1) +
      required(fastq2)

  }

  class TopHatFusionPost extends CommandLineFunction {
    /*
      == command-line ==
      tophat-fusion-post -p 8
      --num-fusion-reads 1
      --num-fusion-pairs 2 --num-fusion-both 5
      --skip-blast /path/to/test_ref
    */

    //@Input(doc="List of Tophat output directories")
    //var tophat_dirs: ListBuffer[File] = _

    @Input(doc="")
    var bowtie_base: File = _

    @Input(doc="")
    var group_dir: File = _

    @Output(doc="")
    var res_file: File = _

    @Argument(doc="number of threads")
    var num_threads: Int = 8

    @Argument(doc="number of reads that should span a fusion")
    var num_fusion_reads: Int = 1

    @Argument(doc="number of mate pairs should span a fusion")
    var num_fusion_pairs: Int = 2

    @Argument(doc="sum of supporting reads and pairs is at least this number for a fusion to be reported")
    var num_fusion_both: Int = 5

    @Argument(doc="flag to whether run or not run blast")
    var blast_skip : Boolean = _

    def commandLine = required("cd", group_dir+";", escape=false) + required("tophat-fusion-post") + required("-p", num_threads) +
    required("--num-fusion-reads", num_fusion_reads) + required("--num-fusion-pairs", num_fusion_pairs) +
      required("--num-fusion-both", num_fusion_both) + conditional(blast_skip, "--skip-blast") +
      required(bowtie_base)
  }


  class PrepareTophatFusion extends CommandLineFunction{

    @Output(doc="group directory")
    var holder: File = _

    @Input(doc="")
    var group_files: List[File] = _

    @Output(doc="")
    var group_dir: File = _

    def commandLine = required("mkdir -p", group_dir+";", escape=false) +
                      required("ln -s ", escape=false) + repeat(group_files.map(_.getParentFile)) + required(group_dir + "/;", escape=false) +
                      required("ln -s " + base_dir.getPath() + "/pipeline/tophat/ensGene.txt " + group_dir + "/ensGene.txt;", escape=false) +
                      required("ln -s " + base_dir.getPath() + "/pipeline/tophat/refGene.txt " + group_dir + "/refGene.txt;", escape=false) +
                      required("ln -s " + base_dir.getPath() + "/pipeline/tophat/blast " + group_dir + "/blast;", escape=false) +
                      required("echo 'temp' > ", holder, escape=false)

  }

  class Rscript extends CommandLineFunction{

    @Input(doc="")
    var source_file: File = _

    @Output(doc="")
    var out_file: File = _

    @Input(doc="")
    var inp_dirs: List[File] = _

    def commandLine = required(qscript.rscript_path) + required("--vanilla") + required(source_file) + "format.thfp.results" + required(out_file) + repeat(inp_dirs)

  }

  def script() {

    if (out_base_dir.exists() == false){
      out_base_dir.mkdirs()
    }

    val log_dir = new File(out_base_dir.getPath + "/logs")

    if (log_dir.exists() == false){
      log_dir.mkdirs()
    }

    val baml_helpers_path = new File("baml_helpers.R")

    if (baml_helpers_path.exists() == false){
      throw new QException("baml_helpers.R needs to be in the current directory")
    }

    val fusion_files = scala.collection.mutable.ListBuffer.empty[File]

    if (inp_fastqs(0).getName() == "fusions.out"){

      for (cur_samp <- inp_fastqs) {

        fusion_files.append(cur_samp)

      }

    }else{

      val th_out = new File(out_base_dir.getPath + "/alignments")

      if (th_out.exists() == false){
        th_out.mkdirs()
      }

      val sample_file_list = scala.collection.mutable.ListBuffer.empty[SampleFile]

      for (cur_samp <- inp_fastqs) {

        sample_file_list.append(SampleFile(cur_samp))

      }

      val samp_map = sample_file_list.groupBy(_.sample)

      //first carry out tophat on all the samples in parallel

      for((sample, vals) <- samp_map){

        val read_type = vals.toList.groupBy(_.isR1) //(if (_.isR1) "R1" else "R2")

        var fastq_map = Map(("R1", read_type(true)(0).file), ("R2", read_type(false)(0).file))

        val tophat = new TopHat

        val tophat_outdir = new File(th_out.getPath + "/tophat_Sample_" + sample + "/")
        val tophat_out = new File(tophat_outdir.getPath + "/fusions.out")

        tophat.fastq1 = fastq_map("R1")
        tophat.fastq2 = fastq_map("R2")
        tophat.bowtie_base = new File(qscript.base_dir + "pipeline/bowtie1_hg19_tophatIndexes/hg19")
        tophat.out_dir = tophat_outdir
        tophat.out_fusions = tophat_out
        tophat.jobOutputFile = out_base_dir + "/logs/tophat_Sample_" + sample + ".out"
        tophat.jobErrorFile = out_base_dir + "/logs/tophat_Sample_" + sample + ".err"
        tophat.memoryLimit = qscript.mem_limit
        tophat.nCoresRequest = num_threads
        tophat.wallTime = 1400
        add(tophat)

        fusion_files.append(tophat_out)

      }

    }

    //next, group the samples by the specified value and run tophat fusion post, after setting up the environment

    val samp_group = fusion_files.grouped(qscript.tf_num)

    //val samp_files = (real_files zip real_samples) cross List("SNP", "INDEL")

    val result_list = scala.collection.mutable.ListBuffer.empty[File]

    for((cur_group, group_ind) <- samp_group.zipWithIndex){

      val group_dir = new File(out_base_dir.getPath + "/tmp" + group_ind)

      //val thfp_holder = new File(group_dir.getPath() + "/holder.txt")

      val prepare_tfp = new PrepareTophatFusion
      prepare_tfp.group_files = cur_group.toList
      prepare_tfp.group_dir = group_dir
      prepare_tfp.jobErrorFile = out_base_dir + "/logs/prepthfp."+group_ind+".err"
      prepare_tfp.jobOutputFile = out_base_dir + "/logs/prepthfp."+group_ind+".out"
      add(prepare_tfp)

      val thfp = new TopHatFusionPost

      val thfp_out = new File(group_dir.getPath + "/tophatfusion_out/result.txt")

      thfp.bowtie_base = new File(qscript.base_dir + "pipeline/bowtie1_hg19_tophatIndexes/hg19")
      thfp.group_dir = group_dir
      thfp.num_threads = qscript.num_threads
      thfp.blast_skip = false
      thfp.jobOutputFile =  out_base_dir + "/logs/tophatFusionPost."+group_ind+".out"
      thfp.jobErrorFile = out_base_dir + "/logs/tophatFusionPost."+group_ind+".err"
      thfp.memoryLimit = qscript.mem_limit + 30
      thfp.nCoresRequest = num_threads
      thfp.res_file = thfp_out
      thfp.wallTime = 1400
      add(thfp)

      result_list.append(thfp_out)

    }

    val summary_run = new Rscript

    val summary_run_out = new File(out_base_dir.getPath + "/tophat_fusion_results_with_blast.txt")

    summary_run.jobOutputFile = out_base_dir.getPath + "/logs/summary.out"
    summary_run.jobErrorFile = out_base_dir.getPath + "/logs/summary.err"

    summary_run.source_file = baml_helpers_path
    summary_run.inp_dirs = result_list.toList
    summary_run.out_file = summary_run_out
    summary_run.memoryLimit = 5
    summary_run.wallTime = 300

    add(summary_run)

  }

}
