//adapted from http://gatkforums.broadinstitute.org/discussion/3060/which-data-processing-steps-should-i-do-per-lane-vs-per-sample/

//Note, the Condor Engine, can only output to either STDOUT or STDERR.
//$Version$

import java.io.{File, PrintWriter}

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import org.broadinstitute.gatk.queue.util.Logging

import scala.collection.mutable.Set


class GetAlignStats_pipeline extends QScript with Logging {

  qscript =>

  @Argument(doc="Memory Limit (GB)", required=false)
  var mem_limit: Int = 30

  @Argument(doc="BAM files to summarize")
  var align_files: List[File] = _

  @Argument(doc="Output directory")
  var out_dir: File =  _

  @Argument(doc = "Base directory", required = false)
  var base_dir: String = "/share/data/resources/gatk_v3.3/"

  @Argument(doc="Path to Rscript", required=false)
  var rscript_path: String = "Rscript"

  class PicardCommandLineFunction extends JavaCommandLineFunction{

    jarFile = new File(base_dir + "/programs/picard-tools-1.124/picard.jar")

    @Input(doc="BAM file", required=true)
    var input: File = _

    @Input(doc="Reference sequence")
    var reference: File = _

    @Output(doc="output file")
    var output: File = _

    var prog_type = ""

    override def commandLine = super.commandLine + required(prog_type) +
      required("INPUT=", input) +
      required("OUTPUT=", output) +
      required("REFERENCE_SEQUENCE=", reference) +
      required("METRIC_ACCUMULATION_LEVEL=", "null") +
      required("METRIC_ACCUMULATION_LEVEL=", "SAMPLE")

  }

  class CollectAlignmentSummaryMetrics extends PicardCommandLineFunction{

    prog_type = "CollectAlignmentSummaryMetrics"

  }

  class CollectInsertSizeMetrics extends PicardCommandLineFunction{

    prog_type = "CollectInsertSizeMetrics"

    override def commandLine = super.commandLine + required("HISTOGRAM_FILE=", output+".histogram.pdf")
  }

  class RscriptSummary extends CommandLineFunction{

    @Input(doc="")
    var source_file: File = _

    @Input(doc="")
    var output_file: File = _

    @Input(doc="")
    var listing: List[File] = _

    @Argument(doc="")
    var skip_lines: Int = _

    def commandLine = required(qscript.rscript_path) + required("--vanilla") + required(source_file) + required(skip_lines) +
      required(output_file) + repeat(listing)

  }

  var summary_script = List(
    "args <- commandArgs(TRUE)",
    "",
    "res.dta <- do.call(rbind, lapply(args[-c(1:2)], function(x){",
      " temp.dta <- read.delim(x, comment.char=\"#\", sep=\"\\t\", header=T, stringsAsFactors=F, nrows=as.numeric(args[1]))",
      " lane.presence <- regmatches(x, regexec(\"_L(00\\\\d+)_\", x))[[1]][2]",
      " sg.presence <- regmatches(x, regexec(\"SampleGroup(\\\\d+)\", x))[[1]][2]",
      " fc.presence <- regmatches(x, regexec(\"FlowCell(\\\\d+)\",x))[[1]][2]",
      " temp.dta$Lane <- lane.presence",
      " temp.dta$SampleGroup <- sg.presence",
      " temp.dta$FlowCell <- fc.presence",
      " temp.dta",
    "}))",
    "",
    "write.table(res.dta, file=args[2], sep=\"\\t\", col.names=T, row.names=F, quote=F)"
  )


  def script(){

    if (!out_dir.exists()) out_dir.mkdir();

    val summarization_script = "temp_script.R"

    val scr_pw = new PrintWriter(summarization_script)

    for (i <- summary_script){
      scr_pw.write(i + "\n")
    }

    scr_pw.close()

    val align_listing = scala.collection.mutable.ListBuffer.empty[File]
    val ins_listing = scala.collection.mutable.ListBuffer.empty[File]

    for (file <- align_files){

      // collect alignment stats
      val aln_stats = new CollectAlignmentSummaryMetrics

      val aln_stats_out = new File(out_dir.getPath() + "/" + file.getName().replace(".bam", ".alignment_summary_metrics"))

      aln_stats.jobOutputFile = aln_stats_out.getPath() + ".out"
      aln_stats.jobErrorFile =aln_stats_out.getPath() + ".err"
      aln_stats.input = file
      aln_stats.output = aln_stats_out
      aln_stats.reference = new File(qscript.base_dir + "bundle_2_8/human_g1k_v37.fasta")
      aln_stats.memoryLimit=qscript.mem_limit
      aln_stats.wallTime=420
      add(aln_stats)

      align_listing.append(aln_stats_out)

      val ins_stats = new CollectInsertSizeMetrics

      val ins_stats_out = new File(out_dir.getPath() + "/" + file.getName().replace(".bam", ".insert_size_metrics"))

      ins_stats.jobOutputFile = ins_stats_out.getPath() + ".out"
      ins_stats.jobErrorFile =ins_stats_out.getPath() + ".err"
      ins_stats.input = file
      ins_stats.output = ins_stats_out
      ins_stats.reference = new File(qscript.base_dir + "bundle_2_8/human_g1k_v37.fasta")
      ins_stats.memoryLimit=qscript.mem_limit
      ins_stats.wallTime=420
      add(ins_stats)

      ins_listing.append(ins_stats_out)
    }

    if (align_listing.length > 0){

      val al_summary = new RscriptSummary
      val al_summary_out = new File(out_dir + "/alignments_summary.txt")

      al_summary.jobOutputFile = al_summary_out.getPath + ".out"
      al_summary.jobErrorFile = al_summary_out.getPath + ".err"
      al_summary.source_file = new File(summarization_script)
      al_summary.listing = align_listing.toList
      al_summary.output_file = al_summary_out
      al_summary.wallTime = 300
      al_summary.memoryLimit = 10
      al_summary.skip_lines = -1
      add(al_summary)
    }

    if (ins_listing.length > 0){

      val ins_summary = new RscriptSummary
      val ins_summary_out = new File(out_dir + "/insert_size_summary.txt")

      ins_summary.jobOutputFile = ins_summary_out.getPath + ".out"
      ins_summary.jobErrorFile = ins_summary_out.getPath + ".err"
      ins_summary.source_file = new File(summarization_script)
      ins_summary.listing = ins_listing.toList
      ins_summary.output_file = ins_summary_out
      ins_summary.wallTime = 300
      ins_summary.memoryLimit = 10
      ins_summary.skip_lines = 1
      add(ins_summary)

    }

  }

}
