import java.io.File

import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction

import scala.io._
import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.utils.interval.IntervalSetRule
import org.broadinstitute.gatk.queue.util.Logging


case class SamplePairs(line: String) {
  val data = line.split(",")
  val patient_id = data(0)
  val AML = data(1).replace("\"", "")
  val Skin = data(2).replace("\"", "")
}

case class VCFPos(line: String) {
  val data = line.split("\t")
  val chrom = data(0)
  val pos_start = data(1)
  val ref = data(3)
  val alt = data(4)
}

class somatic_pipeline extends QScript with Logging {
  qscript =>

  @Argument(doc = "Scatter Count")
  var scatter_count: Int = _

  @Argument(doc = "Memory Limit (GB)")
  var mem_limit: Int = _

  @Argument(doc = "Number of Threads")
  var num_threads: Int = _

  @Argument(doc="Base directory")
  var base_dir: String = _

  @Argument(doc="Genotyper(s) to run, one of (varscan, mutect), can be specified multiple times")
  var genotyper: List[String] = Nil

  @Argument(doc="VarScan2 tumor purity (varscan2 only)", required=false)
  var tumor_purity: Float = 1

  @Argument(doc="VarScan2 normal purity (e.g. .95 here would imply expected normal allele fraction of .05)", required=false)
  var normal_purity: Float = 1

  @Argument(doc="Tumor/normal pair CSV file (relative to specified Base directory)")
  var tn_file: String  = _
  //base_dir+"pipeline/sample_pairs_11_21_2014.csv"

  @Argument(doc="output directory under which 'paired_genotypes' directory is to be created")
  var out_base : String = _

  @Argument(doc="path to alignments", required=false)
  var alignment_path: String = "paired_alignments"

  @Argument(doc="reference path relative to base director (default: bundle_2_8/human_g1k_v37.fasta)", required=false)
  var ref_path: String = "bundle_2_8/human_g1k_v37.fasta"

  @Argument(doc="suffix expected for the bam files", required=false)
  var bam_suffix: String = ".dedup.bam"

  @Argument(doc="should indel realignment NOT be performed on the tumor/normal samples together?", required=false)
  var dont_realign: Boolean = false

  @Argument(doc="Interval file, should be in the bundle_2_8 directory of base_dir, and have available an expanded UNION version (e.g. suffixed with _UNION_500.bed)", required=false)
  var interval_file: File = "nextera_v1_2_for_gatk.bed"

  trait BasicArgs extends CommandLineGATK {
    this.ip = 500
    this.isr = IntervalSetRule.UNION
    this.R = new File(base_dir+qscript.ref_path)
    this.L = Seq(new File(base_dir+"bundle_2_8/"+qscript.interval_file))
    this.jarFile = new File(base_dir+"programs/GenomeAnalysisTK.jar")
    this.memoryLimit = qscript.mem_limit
  }


  class Samtools extends CommandLineFunction{

    //from http://varscan.sourceforge.net/somatic-calling.html
    //and
    //http://varscan.sourceforge.net/support-faq.html#read-counts-different
    //under Warnings about "resetting normal" or "resetting tumor" file
    //per the faq, try combining the tumor and normal samples into one mpileup

    @Input(doc="reference")
    var ref: File = _

    @Input(doc="intervals")
    var ints: File = _

    @Input(doc="normal bam")
    var normal_bam: File = _

    @Input(doc="tumor bam")
    var tumor_bam: File = _

    @Output(doc="output pileup")
    var out_pileup: File = _

    //samtools mpileup -f [reference sequence] [BAM file(s)] >myData.mpileup

    override def commandLine = required(base_dir+"programs/samtools-1.1/samtools") + required("mpileup") + required("-q", 1) + required("-B") +
      required("-d", 16000) + required("-f", ref) + required("-l",ints) + required(normal_bam) + required(tumor_bam) + required(">", escape=false) + required(out_pileup)
  }

  class Mutect extends CommandLineGATK with BasicArgs{

    analysisName = "MuTect"
    analysis_type = "MuTect"
    javaMainClass = "org.broadinstitute.sting.queue.extensions.gatk.CommandLineGATK"
    jarFile = base_dir+"programs/mutect-src/mutect/target/mutect-1.1.7.jar"

    @Input(doc="DBSNP or known callset to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=false)
    var dbSNP: File =new File(base_dir+"mutect_bundle/dbsnp_132_b37.leftAligned.vcf")

    //    @Input(doc="Panel Of Normals or known artifact sites to use (must be in VCF format)", fullName="panel_of_normals", shortName="pon", required=false)
    //    var pon: Seq[File] = Seq()

    @Input(doc="COSMIC sites to use (must be in VCF format)", fullName="cosmic", shortName="C", required=false)
    var cosmic: File = new File(base_dir+"mutect_bundle/b37_cosmic_v54_120711.vcf")

    @Argument(doc="tumor sample name", required = false)
    var tumor_name: String = _

    @Argument(doc="normal sample name", required = false)
    var normal_name: String = _

    @Output(doc="output", required=true)
    var vcf: File = _

    @Output(doc="output", required=true)
    var out: File = _

    @Argument(doc="max alt allele count in normal")
    var max_n_alleles: String = _

    @Argument(doc="max alt allele frequency in normal")
    var max_n_allele_freq: String = _


    override def commandLine = super.commandLine  + required("-dbsnp", dbSNP) + required("-cosmic", cosmic) +
      required("--tumor_sample_name", tumor_name) + required("--normal_sample_name", normal_name) +
      optional("--max_alt_alleles_in_normal_count", max_n_alleles) + optional("--max_alt_allele_in_normal_fraction",  max_n_allele_freq)+
      required("-vcf", vcf) + required("-o", out)

  }

  class Varscan extends JavaCommandLineFunction{

    //http://varscan.sourceforge.net/support-faq.html#read-counts-different
    //under Warnings about "resetting normal" or "resetting tumor" file
    //per the faq, try combining the tumor and normal samples into one mpileup


    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc="input pileup", required=true)
    var pileup: File = _

    @Output(doc="indel output")
    var indel_file: File = _

    @Output(doc="snp output")
    var snp_file: File = _

    override def commandLine = super.commandLine + required("somatic")  + required(pileup) +
      required("--tumor-purity", qscript.tumor_purity) + required("--normal-purity", qscript.normal_purity) +
      required("--min-coverage", 3) + required("--min-var-freq", 0.08) + required("--p-value", 0.10) + required("--somatic-p-value", 0.05) +
      required("--strand-filter", 0) + required("--output-vcf", 1) + required("--mpileup", 1) + required("--output-snp", snp_file) + required("--output-indel", indel_file)

  }

  class VarscanProcessSomatic extends JavaCommandLineFunction{

    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc="input variant file")
    var input_var: File = _

    @Output(doc="output variant file")
    var output_var: File = _

    @Output(doc="high confidence variant file")
    var output_hc_var: File = _

    override def commandLine = super.commandLine + required("processSomatic") + required(input_var) +
      required("--min-tumor-freq", .1) + // - Minimum variant allele frequency in tumor [0.10]
      optional("--max-normal-freq", .05) + // - Maximum variant allele frequency in normal [0.05]
      required("--p-value", .07) // - P-value for high-confidence calling
  }

  class VarscanGetPos extends InProcessFunction {

    @Input(doc="")
    var input_file: File = _

    @Output(doc="")
    var output_file: File = _

    def run {

      val vcf_lines = Source.fromFile(input_file) getLines() filterNot(_.startsWith("#")) map (l => VCFPos(l))

      val pw = new java.io.PrintWriter(output_file)

      for (i <- vcf_lines) {

        //fix deletion positions as in: https://sourceforge.net/p/varscan/discussion/1073559/thread/b2c0c9cb/#195a

        var pos = i.pos_start.toInt

        if (i.ref.length > i.alt.length){
          pos = pos + 1
        }

        pw.write(i.chrom + "\t" + pos + "\t" + pos + "\n")
      }

      pw.close()

    }

  }

  class SamtoolsIndex extends CommandLineFunction{

    @Input(doc="")
    var inp_bam: File = _

    @Output(doc="")
    var out_idx: File = _

    override def commandLine = required("samtools") + required("index", inp_bam)

  }

  class DeleteFiles extends CommandLineFunction{
    @Input(doc="")
    var bam_id_list: List[File] = _

    override def commandLine = required("rm") + repeat(bam_id_list)
  }

  class VarscanBamReadCount extends CommandLineFunction{

    //parameters from http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1504s44/full

    @Input(doc="")
    var input_intervals: File = _

    @Input(doc="")
    var input_bam: File = _

    @Input(doc="")
    var ref: File = _

    @Output(doc="")
    var output_readcount: File = _

    override def commandLine = required(base_dir+"programs/other_varscan/v2_4_1/bam-readcount-bin/bin/bam-readcount") + required("-q", 1) + required("-w", 1) +
      required("-l", input_intervals) + required("-f", ref) + required(input_bam) + required(">", escape=false) + required(output_readcount)
  }


  class VarscanFpfilter extends JavaCommandLineFunction{

    //parameters from: https://github.com/dkoboldt/varscan/blob/master/VarScan.v2.4.1.description.txt

    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc="")
    var input_vcf: File = _

    @Input(doc="")
    var input_readcount: File = _

    @Output(doc="")
    var outfile: File = _

    override def commandLine = super.commandLine + required("fpfilter") + required(input_vcf) + required(input_readcount) + required("--output-file", outfile) +
      required("--keep-failures") +
      required("--min-var-count", 3) +
      required("--min-var-count-lc", 1) +
      required("--min-strandedness", 0) +
      required("--min-var-basequal", 30) +
      required("--min-ref-readpos", 0.20) +
      required("--min-ref-dist3", 0.20) +
      required("--min-var-readpos", 0.15) +
      required("--min-var-dist3", 0.15) +
      required("--max-rl-diff", 0.05) +
      required("--max-mapqual-diff", 10) +
      required("--min-ref-mapqual", 20) +
      required("--min-var-mapqual", 30) +
      required("--max-var-mmqs", 100) +
      required("--max-ref-mmqs", 50)
  }

  class VarscanOldFpfilter extends CommandLineFunction{

    //note this script is under version control, it is an adapted version of a script provided with the older varscan2.

    @Input(doc="")
    var input_vcf: File = _

    @Input(doc="")
    var input_readcount: File = _

    @Argument(doc="")
    var outfile_base: String = _

    override def commandLine = required("perl") + required(base_dir+"programs/other_varscan/v2_4_1/fpfilter_vcf.pl") +
      required("--output-basename", outfile_base) +
      required(input_vcf) +
      required(input_readcount)
  }

  class FilterVarscan extends CommandLineFunction{

    @Input(doc="input vcf file")
    var vcf_input: File = _

    @Output(doc="out filtered vcf file")
    var vcf_output: File = _


    def commandLine = required("grep", "-v") +
      required(";SS=5;") +
      required(vcf_input) +
      required(">", escape=false) + required(vcf_output)
  }


  def script {

    val ref = new File(base_dir+qscript.ref_path)

    val tn_pairs = Source.fromFile(tn_file) getLines() map (l => SamplePairs(l))

    val align_files = new File(alignment_path).listFiles filter (_.endsWith(qscript.bam_suffix))

    val sample_name = align_files map (_.getName().replace(qscript.bam_suffix, "").split("_")(1))

    val file_map = (sample_name zip align_files).toMap

    val out_dir = new StringBuilder(qscript.out_base + "/paired_genotypes/")

    for (i <- tn_pairs) {

      if (file_map.contains(i.AML) && file_map.contains(i.Skin)) {

        var normal_bam_str = file_map(i.Skin)
        var tumor_bam_str = file_map(i.AML)
        val base_out = out_dir + "/Pair_"+i.AML+"_AML_"+i.Skin +"_Skin"
        val base_align_out = qscript.out_base + "/paired_alignments/Pair_"+i.AML+"_AML_"+i.Skin +"_Skin"

        if (dont_realign == false){
          //realign indels for the pair

          //create a combined bam file with both samples being realigned...
          //  targetcreation

          val rtc_outfile = base_align_out + ".rtc.list"

          val rtc_stdout = rtc_outfile + ".out"

          val rtc = new RealignerTargetCreator with BasicArgs
          rtc.jobOutputFile = rtc_stdout
          rtc.jobErrorFile = rtc_stdout.replace(".out", ".err")
          rtc.I = List(file_map(i.AML), file_map(i.Skin))
          rtc.known = List(new File(base_dir+"bundle_2_8/Mills_and_1000G_gold_standard.indels.b37.vcf"),
            new File(base_dir+"bundle_2_8/1000G_phase1.indels.b37.vcf"))
          rtc.out = new File(rtc_outfile.toString)
          rtc.wallTime = 1440
          add(rtc)

          //indelrealigner
          
          normal_bam_str = qscript.out_base + "/paired_alignments/" +  "Sample_"+i.Skin + "_Skin.realign.bam"
          tumor_bam_str = qscript.out_base + "/paired_alignments/"  +  "Sample_"+i.AML + "_AML.realign.bam"

          val indel_r = new IndelRealigner with BasicArgs

          val indel_r_out = new File(rtc_outfile.toString.replace(".rtc.list", ".realigned.bam"))

          indel_r.jobOutputFile = indel_r_out.getPath() + ".out"
          indel_r.jobErrorFile = indel_r_out.getPath() + ".err"
          indel_r.targetIntervals = new File(rtc_outfile.toString)
          indel_r.I = List(file_map(i.AML), file_map(i.Skin))
          indel_r.model = org.broadinstitute.gatk.tools.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_READS
          indel_r.known = List(new File(base_dir+"bundle_2_8/Mills_and_1000G_gold_standard.indels.b37.vcf"),
            new File(base_dir+"bundle_2_8/1000G_phase1.indels.b37.vcf"))
          //indel_r.nWayOut = fm_name.getAbsolutePath()
          indel_r.out = indel_r_out
          indel_r.wallTime = 1440
          add(indel_r)

          //printreads is used here to split the files as they need to be part of an @Output in order to work properly
          val zip_list =  List(i.AML, i.Skin) zip List(tumor_bam_str, normal_bam_str)

          for((cur_samp, samp_file) <- zip_list)
          {
            val print_r = new PrintReads with BasicArgs

            val print_r_out = new File(samp_file)

            print_r.jobOutputFile = print_r_out.getPath() + ".out"
            print_r.jobErrorFile = print_r_out.getPath() + ".err"
            print_r.I = Seq(new File(indel_r_out))
            print_r.out = new File(print_r_out)
            print_r.sample_name = Seq(cur_samp)
            print_r.nct = qscript.num_threads
            print_r.wallTime = 1440
            print_r.generate_md5 = true
            add(print_r)
          }

        }

        //then carry out the genotyping

        if (qscript.genotyper.find(_=="mutect") != None) {
          val mutect = new Mutect

          val mutect_out = base_out + ".mutect"


          val tagged_samps = Seq(
            new TaggedFile(new File(tumor_bam_str), "tumor"),
            new TaggedFile(new File(normal_bam_str), "normal")
          )

          mutect.jobOutputFile = mutect_out + ".out"
          mutect.jobErrorFile = mutect_out + ".err"
          mutect.I = tagged_samps
          mutect.tumor_name = i.AML + "_AML"
          mutect.normal_name = i.Skin + "_Skin"
          mutect.out = new File(mutect_out + ".stats")
          mutect.vcf = new File(mutect_out + ".vcf")
          mutect.max_n_alleles = "1000000"
          mutect.max_n_allele_freq = "1.0"
          mutect.wallTime = 300
          add(mutect)

          //for the somatic calls, restrict to those that pass filter
          val mut_sels = new SelectVariants with BasicArgs

          val mut_sels_out = new File(mutect_out + ".Somatic.filter.vcf")

          mut_sels.jobOutputFile = mut_sels_out.getPath.replace(".vcf", ".out")
          mut_sels.jobErrorFile = mut_sels_out.getPath.replace(".vcf", ".err")
          mut_sels.variant = new File(mutect_out + ".vcf")
          mut_sels.out = mut_sels_out
          mut_sels.excludeFiltered = true
          mut_sels.wallTime = 10

          add(mut_sels)

        }

        if (qscript.genotyper.find(_=="varscan") != None) {

          //now compute pileup on t/n samples seperately
          //note that GATK does not compute BAQ by default so should be ok computing it here.

          //pileup on the tumor sample

          val sam_tumor = new Samtools
          val sam_tumor_out = base_out + ".mpileup"

          sam_tumor.jobOutputFile = sam_tumor_out + ".out"
          sam_tumor.jobErrorFile = sam_tumor_out + ".err"

          sam_tumor.tumor_bam = tumor_bam_str
          sam_tumor.normal_bam = normal_bam_str
          sam_tumor.out_pileup = new File(sam_tumor_out)
          sam_tumor.wallTime = 300
          sam_tumor.ref = ref
          sam_tumor.memoryLimit = qscript.mem_limit
          sam_tumor.ints = new File(base_dir+"bundle_2_8/" + qscript.interval_file.replace(".bed", "_UNION_500.bed"))
          add(sam_tumor)

          //run varscan2 on the t/n pileups
          //procedure adapted from: http://onlinelibrary.wiley.com/doi/10.1002/0471250953.bi1504s44/full with respect to normal contamination and new suggestions from the DREAM challenge


          val var_geno = new Varscan
          val var_geno_out = base_out + ".varscan2.t."+(qscript.tumor_purity * 100).toInt.toString() + ".n."+(qscript.normal_purity * 100).toInt.toString()
          var_geno.jobOutputFile = var_geno_out + ".out"
          var_geno.jobErrorFile = var_geno_out + ".err"

          var_geno.pileup = new File(sam_tumor_out)
          var_geno.snp_file = new File(var_geno_out + ".snp.vcf")
          var_geno.indel_file = new File(var_geno_out + ".indel.vcf")
          var_geno.wallTime = 300
          var_geno.memoryLimit = qscript.mem_limit
          add(var_geno)

          for (var_file <- List(new File(var_geno_out + ".snp.vcf"), new File(var_geno_out + ".indel.vcf"))){

            val interval_file = new File(var_file.getPath().replace(".vcf", ".interval"))

            val vcf_pos = new VarscanGetPos

            vcf_pos.input_file = var_file
            vcf_pos.output_file = interval_file
            add(vcf_pos)

            //from here, don't use somaticFilter: https://sourceforge.net/p/varscan/discussion/1073559/thread/b2c0c9cb/#de5b

            val brc = new VarscanBamReadCount
            val brc_out =  new File (var_file.getPath().replace(".vcf", ".brc"))

            brc.jobOutputFile = brc_out.getPath() + ".out"
            brc.jobErrorFile = brc_out.getPath() + ".err"
            brc.input_bam =tumor_bam_str
            brc.input_intervals = interval_file
            brc.output_readcount = brc_out
            brc.ref = ref
            brc.wallTime = 120
            brc.memoryLimit = qscript.mem_limit
            add(brc)

            val fpfilter = new VarscanFpfilter

            val fpfilter_out = new File(var_file.getPath().replace(".vcf", ".flag.vcf"))
            fpfilter.jobOutputFile = fpfilter_out.getPath() + ".out"
            fpfilter.jobErrorFile = fpfilter_out.getPath() + ".err"
            fpfilter.input_readcount = brc_out
            fpfilter.input_vcf = var_file
            fpfilter.outfile = fpfilter_out
            fpfilter.wallTime = 10
            fpfilter.memoryLimit = 10
            add(fpfilter)

            //filter out the unknown category as it causes issues with VCF standard

            val gene_filt_vep = new FilterVarscan

            val gene_filt_vep_out = new File(fpfilter_out.getPath.replace(".flag.vcf", ".flag.fixed.vcf"))

            gene_filt_vep.jobOutputFile = gene_filt_vep_out + ".out"
            gene_filt_vep.jobErrorFile = gene_filt_vep_out + ".err"
            gene_filt_vep.vcf_input = fpfilter_out
            gene_filt_vep.vcf_output = gene_filt_vep_out
            gene_filt_vep.wallTime = 5
            add(gene_filt_vep)

            var varscan_to_filter = gene_filt_vep_out
            var rep_suffix = ".flag.fixed"

            if (var_file.getName.matches(".*indel.*")){

              //also left align and trim if indels
              val lalign = new LeftAlignAndTrimVariants

              val lalign_out = new File(gene_filt_vep_out.getPath.replace(".vcf", ".la.vcf"))

              lalign.jobOutputFile = lalign_out.getPath + ".out"
              lalign.jobErrorFile = lalign_out.getPath + ".err"
              lalign.R = new File(base_dir+"bundle_2_8/human_g1k_v37.fasta")
              lalign.variant = gene_filt_vep_out
              lalign.o = lalign_out
              lalign.splitMultiallelics = true
              lalign.wallTime=24
              lalign.memoryLimit = qscript.mem_limit
              add(lalign)

              varscan_to_filter = lalign_out
              rep_suffix = ".flag.fixed.la"
            }

            val var_filt = new VarscanProcessSomatic


            val var_filt_out = new File(varscan_to_filter.getPath().replace(".vcf", ".Somatic.vcf"))
            val var_hc_out = new File(varscan_to_filter.getPath().replace(".vcf", ".Somatic.hc.vcf"))

            var_filt.jobOutputFile = var_filt_out.getPath() + ".out"
            var_filt.jobErrorFile = var_filt_out.getPath() + ".err"
            var_filt.input_var = varscan_to_filter
            var_filt.output_var = var_filt_out
            var_filt.output_hc_var = var_hc_out
            var_filt.wallTime = 10
            var_filt.memoryLimit = 2
            add(var_filt)

            val var_sels = new SelectVariants with BasicArgs

            val var_sels_out = new File(var_hc_out.getPath().replace(rep_suffix, "").replace(".hc.vcf", ".filter.vcf"))

            var_sels.jobOutputFile = var_sels_out.getPath.replace(".vcf", ".out")
            var_sels.jobErrorFile = var_sels_out.getPath.replace(".vcf", ".err")
            var_sels.variant = var_hc_out
            var_sels.out = var_sels_out
            var_sels.excludeFiltered = true
            var_sels.wallTime = 10

            add(var_sels)


            //also run a filtering equivalent to the way it was done in v2.3.7
            val oldfp = new VarscanOldFpfilter

            val oldfp_out = new File(var_filt_out.getPath().replace(rep_suffix, "").replace("Somatic.vcf", "prev.Somatic"))

            oldfp.jobOutputFile = oldfp_out.getPath + ".out"
            oldfp.jobErrorFile = oldfp_out.getPath + ".err"
            oldfp.outfile_base = oldfp_out
            oldfp.input_readcount = brc_out
            oldfp.input_vcf = var_filt_out
            oldfp.wallTime=2
            add(oldfp)

          }

        }

      }

    }

  }

}
