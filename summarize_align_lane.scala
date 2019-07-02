//adapted from http://gatkforums.broadinstitute.org/discussion/3060/which-data-processing-steps-should-i-do-per-lane-vs-per-sample/

//Note, the Condor Engine, can only output to either STDOUT or STDERR.
//$Version$

import java.io.{File, PrintWriter}
import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.utils.interval.IntervalSetRule
import org.broadinstitute.gatk.queue.util.Logging
import scala.util.matching.Regex
import scala.collection.mutable._
import scala.collection.mutable.Set


case class SampleFilename (path: File){
  val split_str = path.getName.split("[_\\.]")

  val sample_name = split_str(1)
  val lane = split_str(3)
}

class GATK_pipeline extends QScript with Logging {
  qscript =>

  @Argument(doc = "Number of Threads", required = false)
  var num_threads: Int = 7

  @Argument(doc = "Path to flowcell to process")
  var flowcell_path: String = _

  @Argument(doc = "Base directory", required = false)
  var base_dir: String = "/share/data/resources/gatk_v3.3/"

  @Argument(doc = "Base output directory")
  var out_base_dir: String = _

  @Argument(doc="BAM ID prefix")
  var id_pref: String = "BeatAML"

  var bam_suffix: String = ".sam"

  class TrimReads extends CommandLineFunction {

    @Input(doc="")
    var R1: File = _

    @Input(doc="")
    var R2: File = _

    @Output(doc="")
    var R1_out: File = _

    @Output(doc="")
    var R2_out: File = _

    //five.trim=3, three.trim=5
    //required("-e", "library(seqtools)") +
    override def commandLine = required("java", "-jar") + required(base_dir+"programs/Trimmomatic-0.36/trimmomatic-0.36.jar") +
      required("PE") + required("-phred33") + required("-threads", qscript.num_threads) +
      required(R1) + required(R2) +
      required(R1_out) + required(new File(R1_out.getPath.replace(".trim.fastq.gz", ".unp.fastq.gz"))) +
      required(R2_out) + required(new File(R2_out.getPath.replace(".trim.fastq.gz", ".unp.fastq.gz"))) +
      required("HEADCROP:3") + required("CROP:92")
  }

  class BWAlign extends CommandLineFunction {

    //mem -M -t 7 -R '@RG\tID:BeatAML_FlowCell1_L007_13-00218\tSM:13-00218\tPL:ILLUMINA\tLB:LIB_13-00218' /share/data/resources/gatk_v3.3//bundle_2_8/human_g1k_v37.fasta processed_fastq_lane/Sample_DNA140626JT_9_13-00218_R1_L007.trim.fastq.gz processed_fastq_lane/Sample_DNA140626JT_9_13-00218_R2_L007.trim.fastq.gz"
    //Output=alignments/Sample_L007_13-00218.sam

    @Input(doc="")
    var input_rs: List[File] = Nil

    @Argument(doc="")
    var sample_name: String = _

    @Argument(doc="")
    var flowcell_name: String = _

    @Argument(doc="")
    var lane: String = _

    @Output(doc="")
    var out_sam: File = _


    override def commandLine = required(qscript.base_dir + "programs/bwa-0.7.10/bwa") +
      required("mem") +
      required("-R", "'@RG\tID:"+id_pref+"_"+flowcell_name+"_L"+lane+"_"+sample_name+
        "\tSM:"+sample_name+"\tPL:ILLUMINA\tLB:LIB_"+sample_name+"'", escape=false) +
      required(qscript.base_dir + "bundle_2_8/human_g1k_v37.fasta") +
      repeat(input_rs) +
      required(">", escape=false) +
      required(out_sam)

  }

  def script() {

    // extract flowcell and sample group
    var tmp_inst = new File(flowcell_path)
    var flowcell = tmp_inst.getName()
    var sampleGroup = tmp_inst.getParent().getName()
    var file_separator = System.getProperty("file.separator")
    var out_flowcell_base = qscript.out_base_dir + file_separator + sampleGroup + file_separator + flowcell

    def addBWA(flowcell_path: String, sample_name: String, fastq_list: List[File], cur_lane: String, id_prefix: String): Unit = {

      val bwa = new BWAlign
      val flowcell_name = new File (flowcell_path).getCanonicalFile().getName()
      val bwa_out = new File(out_flowcell_base + "/alignments/Sample_" + cur_lane + "_" + sample_name + "_" + id_prefix + "_" + flowcell_name + ".sam")

      bwa.jobOutputFile = bwa_out.getPath() + ".out"
      bwa.jobErrorFile = bwa_out.getPath() + ".err"
      bwa.memoryLimit = 6
      bwa.nCoresRequest = qscript.num_threads

      bwa.sample_name = sample_name
      bwa.flowcell_name = flowcell_name
      bwa.input_rs = fastq_list
      bwa.lane = cur_lane
      bwa.out_sam = bwa_out
      bwa.wallTime = 470
      add(bwa)

      val ssam = new SortSam

      val ssam_out = new File(bwa_out.replace(".sam", ".bam"))

      ssam.jobOutputFile = ssam_out.getPath() + ".out"
      ssam.jobErrorFile = ssam_out.getPath() + ".err"
      ssam.input = Seq(bwa_out)
      ssam.output = ssam_out
      ssam.maxRecordsInRam=5000000
      ssam.memoryLimit= 10
      ssam.wallTime = 420
      add(ssam)

    }


    //iterate over the sequencing core flowcell directories (usually one) which always contains a string starting with a
    //digit.

    val flowcell_files = new File(flowcell_path).listFiles
    for (use_dir <- flowcell_files if use_dir.getName.matches("^\\d.*$")) {

      val dna_files = new File(use_dir.getAbsolutePath).listFiles

      for(use_dna_dir <- dna_files if use_dna_dir.getName.matches("^DNA.*$")){

        val sample_files = new File(use_dna_dir.getAbsolutePath).listFiles

        //for each sample

        for (cur_samp <- sample_files if cur_samp.getName.matches("^DNA.*_R1_.*fastq.gz$")) {

          val samp_file = SampleFilename(cur_samp)

          val is_pe = sample_files.exists(_.getName == cur_samp.getName.replace("_R1_", "_R2_"))

          if (is_pe == false){
            throw new QException("Only supports paired end data for now...")
          }

          val sum_reads = new TrimReads
          val sum_reads_out = qscript.flowcell_path + "/processed_fastq/" + cur_samp.getName().replace(".fastq.gz", ".trim.fastq.gz")

          val fastq_list = List(new File(sum_reads_out), new File(sum_reads_out.getPath.replace("_R1_", "_R2_")))

          sum_reads.jobOutputFile = sum_reads_out + ".out"
          sum_reads.jobErrorFile = sum_reads_out + ".err"
          sum_reads.R1 = cur_samp
          sum_reads.R2 = new File(cur_samp.getPath.replace("_R1_", "_R2_"))
          sum_reads.R1_out = fastq_list(0)
          sum_reads.R2_out = fastq_list(1)
          sum_reads.memoryLimit = 4
          sum_reads.nCoresRequest = qscript.num_threads
          sum_reads.wallTime = 600
          add(sum_reads)

          //carry out the aligments at the lane level
          addBWA(qscript.flowcell_path, samp_file.sample_name, fastq_list, samp_file.lane, qscript.id_pref)


          }

        }

      }
    }
}
