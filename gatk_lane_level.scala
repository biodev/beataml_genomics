//adapted from http://gatkforums.broadinstitute.org/discussion/3060/which-data-processing-steps-should-i-do-per-lane-vs-per-sample/

//Note, the Condor Engine, can only output to either STDOUT or STDERR.
//$Version$

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.utils.interval.IntervalSetRule
import org.broadinstitute.gatk.utils.interval.IntervalSetRule._
import org.broadinstitute.gatk.queue.util.Logging

import scala.collection.mutable.Set

class MarkDuplicatesSample extends MarkDuplicates with PicardBamFunction{
        override def commandLine = super.commandLine + "READ_NAME_REGEX=null"
}

class GATK_pipeline extends QScript with Logging {

  qscript =>

  @Argument(doc="Scatter Count", required=false)
  var scatter_count: Int = 4

  @Argument(doc="Memory Limit (GB)", required=false)
  var mem_limit: Int = 30

  @Argument(doc="Number of Threads", required=false)
  var num_threads: Int = 7

  var bam_suffix: String = ".bam"

  @Argument(doc = "Base directory", required = false)
  var base_dir: String = "/share/data/resources/gatk_v3.3/"

  @Argument(doc = "Base directory containing 'alignments' directory. Creates 'sample_alignments' directory under here.", required = true)
  var input_align_basedir: String = _

  @Argument(doc="Relative path to sequence capture coordinates (default: bundle_2_8/nextera_v1_2_for_gatk.bed)", required=false)
  var interval_coords: String = "bundle_2_8/nextera_v1_2_for_gatk.bed"

  //for nimblegen it should be bundle_2_8/Nimblegen_SeqCap_EZ_v3.bed

  trait BasicArgs extends org.broadinstitute.gatk.queue.extensions.gatk.CommandLineGATK {
    this.ip = 500
    this.isr = IntervalSetRule.UNION
    this.R = new File(qscript.base_dir + "bundle_2_8/human_g1k_v37.fasta")
    this.L = Seq(new File(qscript.base_dir + qscript.interval_coords))
    this.jarFile = new File(qscript.base_dir + "programs/GenomeAnalysisTK.jar")
  }
        
  def script(){

    val out_dir = new StringBuilder(qscript.input_align_basedir + "/sample_alignments/")

    val align_files = new File(qscript.input_align_basedir + "/alignments").listFiles

    val all_samples : Set[String] = Set()

    val real_files = scala.collection.mutable.ListBuffer.empty[File]

    val out_dir_file = new File(out_dir.toString)

    if (!out_dir_file.exists()) out_dir_file.mkdir();

    for (file <- align_files if file.getName.endsWith(qscript.bam_suffix)){

      //sort and convert sam file to bam
      // Nov 10th, 2015
      // sorting and SAM -> BAM conversion is now being done as
      // part of Picard Alignments metrics function
      // let's comment out this part that converts sorts SAM to BAM

      //val ssam_outfile = out_dir.clone()
      //ssam_outfile ++= file.getName().replace(".sam", ".bam")
      //val ssam_stdout = ssam_outfile + ".out"
      //val ssam = new SortSam
      //ssam.jobOutputFile = ssam_stdout
      //ssam.input = Seq(file)
      //ssam.output = new File(ssam_outfile.toString)
      //ssam.maxRecordsInRam=5000000
      //ssam.memoryLimit=qscript.mem_limit
      //ssam.wallTime = 420
      //add(ssam)

      //do the initial lane level mark dup

      val md_lane_outfile = out_dir.clone()

      md_lane_outfile ++= file.getName().replace(".bam", ".dedup.bam")
                
      val md_lane_stdout = md_lane_outfile + ".out"
                
      val md_lane = new MarkDuplicates

      md_lane.jobOutputFile = md_lane_stdout
      md_lane.input = Seq(file)
      md_lane.output = new File(md_lane_outfile.toString)
      md_lane.metrics = new File(md_lane_outfile.toString().replace(".bam", ".metrics"))
      md_lane.MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 1000
      md_lane.maxRecordsInRam=5000000
      md_lane.memoryLimit=qscript.mem_limit
      md_lane.wallTime = 420
      add(md_lane)
                
      //indel realignment
                
      val rtc_outfile = md_lane_outfile.toString().replace(".dedup.bam", "") + ".rtc.list"

      val rtc_stdout = rtc_outfile + ".out"

      val rtc = new RealignerTargetCreator with BasicArgs
      rtc.jobOutputFile = rtc_stdout
      //the rtc might have a bug when running with multiple threads
      rtc.nt = 1
      rtc.I = List(new File(md_lane_outfile.toString))
      rtc.known = List(new File(qscript.base_dir + "bundle_2_8/Mills_and_1000G_gold_standard.indels.b37.vcf"), new File(qscript.base_dir + "bundle_2_8/1000G_phase1.indels.b37.vcf"))
      rtc.out = new File(rtc_outfile.toString)
      rtc.memoryLimit=2*qscript.num_threads
      rtc.wallTime = 420
      add(rtc)
                
                
      val indel_r = new IndelRealigner with BasicArgs

      val indel_r_outfile = md_lane_outfile.toString().replace(".bam", ".realign.bam")

      val indel_r_stdout = indel_r_outfile + ".out"

      indel_r.jobOutputFile = indel_r_stdout
      indel_r.targetIntervals = new File(rtc_outfile.toString)
      indel_r.I = List(new File(md_lane_outfile.toString))
      indel_r.model = org.broadinstitute.gatk.tools.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_READS
      indel_r.known = List(new File(qscript.base_dir + "bundle_2_8/Mills_and_1000G_gold_standard.indels.b37.vcf"), new File(qscript.base_dir + "bundle_2_8/1000G_phase1.indels.b37.vcf"))
      indel_r.out = new File(indel_r_outfile.toString)
      indel_r.scatterCount = qscript.scatter_count
      indel_r.memoryLimit=4
      indel_r.wallTime = 1440
      add(indel_r)
                
      all_samples ++= Set(file.getName().replace(".sam", "").split("_")(2))
      real_files.append(new File(indel_r_outfile.toString))
    }
            
    for (i <- all_samples)
    {
            
      //then aggregate the samples and perform base recalibration and another round of dedupping
            
      val file_buf = scala.collection.mutable.ListBuffer.empty[File]
            
      for(j <- real_files if j.getName().indexOf(i) != -1)
      {
        file_buf.append(j)
      }
                
      //for each sample run baserecalibration and print to a sample-level bam
                
      val base_recal_initial = new BaseRecalibrator with BasicArgs
		
      val base_recal_initial_outfile = out_dir.clone()
      base_recal_initial_outfile ++= "Sample_" + i + "_recal_data.tab"

      val base_recal_initial_stdout = base_recal_initial_outfile + ".out"

      base_recal_initial.jobOutputFile = base_recal_initial_stdout
      base_recal_initial.I = file_buf.toSeq
      base_recal_initial.knownSites = List(new File(qscript.base_dir + "bundle_2_8/1000G_phase1.indels.b37.vcf"), new File(qscript.base_dir + "bundle_2_8/dbsnp_138.b37.vcf"), new File(qscript.base_dir + "bundle_2_8/Mills_and_1000G_gold_standard.indels.b37.vcf"))
      base_recal_initial.nct = qscript.num_threads
      base_recal_initial.out = new File(base_recal_initial_outfile.toString)
      base_recal_initial.scatterCount = qscript.scatter_count
      base_recal_initial.memoryLimit=(qscript.num_threads/2)+1
      base_recal_initial.wallTime = 420
      add(base_recal_initial)


      val print_reads = new PrintReads with BasicArgs

      val print_reads_out = out_dir.clone()
      print_reads_out ++= "Sample_" + i + ".bam"
      val print_reads_stdout = print_reads_out + ".out"
      print_reads.jobOutputFile = print_reads_stdout
      print_reads.nct = qscript.num_threads
      print_reads.I = file_buf.toSeq
      print_reads.BQSR = new File(base_recal_initial_outfile.toString)
      print_reads.o = new File(print_reads_out.toString)
      print_reads.memoryLimit=qscript.num_threads
      print_reads.wallTime = 420
      add(print_reads)

      //then run a sample-level deduplication

      val md_samp_out = print_reads_out.toString().replace(".bam", ".dedup.bam")

      val md_samp_stdout = md_samp_out + ".out"

      val md_samp = new MarkDuplicatesSample
      md_samp.jobOutputFile = md_samp_stdout
      md_samp.input = Seq(new File(print_reads_out.toString))
      md_samp.output = new File(md_samp_out.toString)
      md_samp.metrics = new File(print_reads_out.toString().replace(".bam", ".metrics"))
      md_samp.MAX_FILE_HANDLES_FOR_READ_ENDS_MAP = 1000
      md_samp.maxRecordsInRam=5000000
      md_samp.memoryLimit=qscript.mem_limit
      md_samp.wallTime = 420
      add(md_samp)
    }

  }

}
