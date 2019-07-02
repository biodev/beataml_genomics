/*

  Functions adapted from gatk_workflows@beataml_rnaseq2:rnaseq.scala

 */

import java.io.{File, PrintWriter}

import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction

import scala.io._
import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.util.Logging
import scala.collection.mutable

case class SampleFile(line: File) {
  val data = line.getName.split("_")
  val sample = data(1)
  val file = line
  val isR1 = line.getName.matches(".*_R1_.*")

  val flowcell_re = ".*([Ff]low[cC]ell\\d+).*".r
  val flowcell = line.getAbsolutePath match {
    case flowcell_re(m) => m;
    case _ => "FlowCellY"
  }

  val sg_re = ".*(SampleGroup\\d+).*".r
  val samplegroup = line.getAbsolutePath match {
    case sg_re(m) => m;
    case _ => "SampleGroupX"
  }

}

class rnaseq_pipeline extends QScript with Logging {
  qscript =>

  @Argument(doc = "Memory Limit (GB)")
  var mem_limit: Int = _

  @Argument(doc = "Number of Threads")
  var num_threads: Int = _

  @Argument(doc = "Fastq paths")
  var inp_fastqs: List[File] = _

  @Argument(doc = "Resource Base directory", required = false)
  var base_dir: String = "/mnt/lustre1/NGSdev/resources/"

  @Argument(doc = "Path to Subjunc indexes")
  var subjunc_indexes: File = _

  @Argument(doc = "Output directory where Subjunc and TopHat results are to be placed")
  var out_base_dir: String = _

  @Argument(doc = "Subjunc reads orientation")
  var subjunc_readsorientation: String = _

  @Argument(doc = "annotation file for featureCounts", required = true)
  var annotation_file: File = _

  @Argument(doc = "strand specificity option of featureCounts")
  var fc_orientation: Int = _

  @Argument(doc = "Should the inputs be trimmed (default is false)", required = false)
  var should_trim: Boolean = false

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


  class Subjunc extends CommandLineFunction {

    /*
   *  subjunc -i path/to/reference/ \
   *  -u -r fastq1 -R fastq2 --gzFASTQinput \
   *  -o outputSAMFilename -I 5 -n 10 -T 7 \
   *  -H -d 50 -D 600 -S fr -m 3 -p 1 --allJunctions \
   *  --rg-id BeatAML_FlowCell1_13-00098 --rg SM:13-00098 \
   *  --rg PL:ILLUMINA --rg LB:LIB_13-00098
  */
    @Input(shortName = "i", doc = "path to reference indexes")
    var reference_base: File = _

    @Input(shortName = "r", doc = "path to first fastq file")
    var fastq1: File = _

    @Input(shortName = "R", doc = "path to second fastq file")
    var fastq2: File = _

    @Output(shortName = "o", doc = "Output BAM file name")
    var outBAM: File = _

    @Argument(shortName = "rg-id", doc = "read group ID")
    var rgheaderID: String = _

    @Argument(shortName = "rg", doc = "<tag:value> to the read group header in the mapping output - sample info")
    var rgheaderSM: String = _

    @Argument(shortName = "rg", doc = "<tag:value> to the read group header in the mapping output - platform info")
    var rgheaderPL: String = _

    @Argument(shortName = "rg", doc = "<tag:value> to the read group header in the mapping output - library info")
    var rgheaderLIB: String = _

    @Argument(shortName = "S", doc = "specify the orientation of the two reads from the same pair")
    var readsOrientation: String = _

    @Argument(shortName = "T", doc = "number of CPUs/threads to use for mapping")
    var num_threads: Int = 7

    @Argument(shortName = "I", doc = "number of Indel bases allowed in the mapping, default = 5, max = 200bp")
    var num_indels: Int = 5

    @Argument(shortName = "d", doc = "minimum fragment length; default 50")
    var minFragLength: Int = 50

    @Argument(shortName = "D", doc = "maximum fragment length; default 600")
    var maxFragLength: Int = 600

    def commandLine = required("subjunc") +
      required("-u") +
      required("-i", reference_base) +
      required("-r", fastq1) +
      required("-R", fastq2) +
      required("-o", outBAM) +
      required("-I", num_indels) +
      required("-T", num_threads) +
      required("-d", minFragLength) +
      required("-D", maxFragLength) +
      required("-S", readsOrientation) +
      required("--rg-id", rgheaderID) +
      required("--rg", rgheaderSM) +
      required("--rg", rgheaderPL) +
      required("--rg", rgheaderLIB)
  }

  /*
  * ====================================
  * Run featureCounts on subjunc output
  * =====================================
  */
  class FeatureCounts extends CommandLineFunction {

    /*
      $fc_exec -a $gtf_file -o $output -F GTF -t exon -g gene_id -s 2 -C -T 10 -p -P -d 10 -D 600 -B $cur_sams
    */

    @Input(shortName = "a", doc = "Annotation file")
    var annot_file: File = _

    @Output(shortName = "o", doc = "name of the output file")
    var outfile: File = _

    @Argument(shortName = "F", doc = "format of the annotation file")
    var annot_format: String = "GTF"

    @Argument(shortName = "t", doc = "feature type to look in GTF file; default exon")
    var featureType: String = "exon"

    @Argument(shortName = "g", doc = "attribute type used to group features")
    var attrType: String = "gene_id"

    @Argument(shortName = "s", doc = "indicate if strand specific read counting should be performed; 0 (unstranded), 1(stranded), 2(reversely stranded)")
    var strandSpecificity: Int = _

    @Argument(shortName = "T", doc = "number of CPUs/threads")
    var num_threads: Int = 7

    @Argument(shortName = "d", doc = "minimum fragment length; default 50")
    var minFragLength: Int = 50

    @Argument(shortName = "D", doc = "maximum fragment length; default 600")
    var maxFragLength: Int = 600

    @Input(doc = "BAM files")
    var bamfiles: File = _

    def commandLine = required("featureCounts") +
      required("-a", annot_file) +
      required("-o", outfile) +
      required("-F", annot_format) +
      required("-t", featureType) +
      required("-g", attrType) +
      required("-s", strandSpecificity) +
      required("-C") +
      required("-T", num_threads) +
      required("-p") +
      required("-P") +
      required("-d", minFragLength) +
      required("-D", maxFragLength) +
      required("-B") +
      required(bamfiles)
  }

  var r_func = List("args <- commandArgs(TRUE)",
    "" ,
    "fc.files <- args[2:length(args)]" ,
    "" ,
    "expr.list <- lapply(fc.files, function(x)" ,
    "   {" ,
    "      temp.mat <- read.delim(x, sep=\"\\t\", comment.char=\"#\", stringsAsFactors = F)" ,
    "      colnames(temp.mat)[7] <- sub(\"\\\\.\", \"-\", regmatches(colnames(temp.mat)[7], regexec(\"Sample_([[:alnum:][:punct:]]+)\\\\.bam\", colnames(temp.mat)[7]))[[1]][2])" ,
    "      return(temp.mat)" ,
    "    })" ,
    "" ,
    "    all.rows <- unique(sapply(expr.list, nrow))" ,
    "" ,
    "    stopifnot(length(all.rows) == 1)" ,
    "" ,
    "    expr.mat <- Reduce(function(x,y) merge(x,y, by=names(expr.list[[1]])[1:6], all=F, sort=F), expr.list)" ,
    "" ,
    "    stopifnot(nrow(expr.mat) == all.rows)" ,
    "" ,
    "    write.table(expr.mat, file=args[1], sep=\"\\t\", col.names=T, row.names=F, quote=F)")

  class SummarizeFC extends CommandLineFunction{

    /*summarize per sample FeatureCounts results*/

    @Input(doc = "summarization script")
    var summarization_script: File = _

    @Input(doc = "featureCounts files")
    var inp_files: List[File] = _

    @Output(doc="output files")
    var out_file: File = _

    def commandLine = required("Rscript") +
      required("--vanilla") +
      required(summarization_script) +
      required(out_file) +
      repeat(inp_files)

  }

  def script {

    val summarization_script = "temp_script.R"

    val scr_pw = new PrintWriter(summarization_script)

    for (i <- r_func){
      scr_pw.write(i + "\n")
    }

    scr_pw.close()


    //(sample_name zip align_files).toMap
    val sample_file_list = scala.collection.mutable.ListBuffer.empty[SampleFile]

    for (cur_samp <- inp_fastqs) {

      sample_file_list.append(SampleFile(cur_samp))

    }


    val samp_map = sample_file_list.groupBy(_.sample)

    val file_separator = System.getProperty("file.separator")

    val result_buffer = scala.collection.mutable.ListBuffer.empty[File]

    for ((sample, vals) <- samp_map) {

      val read_type = vals.toList.groupBy(_.isR1) //(if (_.isR1) "R1" else "R2")

      var fastq_map = Map(("R1", read_type(true)(0).file), ("R2", read_type(false)(0).file))

      val out_flowcell_base = qscript.out_base_dir + file_separator + read_type(true)(0).flowcell

      //could check for no false values in read_type as a way of determining paired vs unpaired reads

      if(read_type.contains(true) && qscript.should_trim){

        val sum_reads = new TrimReads
        val sum_reads_out_1 = new File(out_flowcell_base.getPath + "/processed_fastq/" + fastq_map("R1").getName().replace(".fastq.gz", ".trim.fastq.gz"))
        val sum_reads_out_2 = new File(sum_reads_out_1.getPath.replace("_R1_", "_R2_"))

        sum_reads.jobOutputFile = sum_reads_out_1.getPath + ".out"
        sum_reads.jobErrorFile = sum_reads_out_1.getPath + ".err"
        sum_reads.R1 = fastq_map("R1")
        sum_reads.R2 = fastq_map("R2")
        sum_reads.R1_out = sum_reads_out_1
        sum_reads.R2_out = sum_reads_out_2
        sum_reads.memoryLimit = 4
        sum_reads.nCoresRequest = qscript.num_threads
        sum_reads.wallTime = 600
        add(sum_reads)

        fastq_map = Map(("R1", sum_reads_out_1), ("R2", sum_reads_out_2))

      }else if (qscript.should_trim) {

        throw new QException("Trimming is only supported for paired end reads currently")

      }

      val subjunc_inst = new Subjunc
      subjunc_inst.reference_base = qscript.subjunc_indexes
      subjunc_inst.fastq1 = fastq_map("R1")
      subjunc_inst.fastq2 = fastq_map("R2")

      val alignments_dir = new File(out_flowcell_base + file_separator + "subjunc_alignments")
      subjunc_inst.outBAM = alignments_dir + file_separator + "Sample_" + sample + ".bam"
      subjunc_inst.readsOrientation = qscript.subjunc_readsorientation

      //--rg-id BeatAML_FlowCell1_13-00658 --rg SM:13-00658  --rg PL:ILLUMINA --rg LB:LIB_13-00658
      subjunc_inst.rgheaderID = "BeatAML_" +  read_type(true)(0).samplegroup + "_" + read_type(true)(0).flowcell + "_" + sample
      subjunc_inst.rgheaderSM = "SM:" + sample
      subjunc_inst.rgheaderPL = "PL:" + "ILLUMINA"
      subjunc_inst.rgheaderLIB = "LB:LIB_" + sample

      subjunc_inst.jobOutputFile = out_flowcell_base + file_separator + "logs" + file_separator + "subjunc_Sample_" + sample + ".out"
      subjunc_inst.jobErrorFile = out_flowcell_base + file_separator + "logs" + file_separator + "subjunc_Sample_" + sample + ".err"
      subjunc_inst.memoryLimit = qscript.mem_limit
      subjunc_inst.nCoresRequest = num_threads
      subjunc_inst.wallTime = 600
      add(subjunc_inst)

      //featurecounts

      var fC_inst = new FeatureCounts
      fC_inst.annot_file = qscript.annotation_file
      fC_inst.strandSpecificity = qscript.fc_orientation
      fC_inst.outfile = out_flowcell_base + file_separator + "Sample_" + sample + "_subjunc_featureCounts.txt"
      fC_inst.bamfiles = subjunc_inst.outBAM

      fC_inst.jobOutputFile = out_flowcell_base + file_separator + "logs" + file_separator + "Sample_" + sample + "_subjunc_featureCounts.out"
      fC_inst.jobErrorFile = out_flowcell_base + file_separator + "logs" + file_separator + "Sample_" + sample + "_subjunc_featureCounts.err"
      fC_inst.memoryLimit = qscript.mem_limit
      fC_inst.nCoresRequest = num_threads
      fC_inst.wallTime = 100
      add(fC_inst)

      result_buffer.append(new File(fC_inst.outfile))

    }

    val summfc = new SummarizeFC

    val summfc_out = new File(result_buffer(0).getParent + file_separator + "combined_featureCounts.txt")

    summfc.jobOutputFile = result_buffer(0).getParent + file_separator + "logs" + file_separator + summfc_out.getName.replace(".txt", ".out")
    summfc.jobErrorFile = summfc.jobOutputFile.replace(".out", ".err")

    summfc.summarization_script = summarization_script
    summfc.out_file = summfc_out
    summfc.inp_files = result_buffer.toList

    add(summfc)



    }

}