import java.io.File

import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import org.broadinstitute.gatk.utils.interval.IntervalSetRule
import org.broadinstitute.gatk.queue.util.Logging

import scala.io.Source


case class VCFPos(line: String) {
  val data = line.split("\t")
  val chrom = data(0)
  val pos_start = data(1)
  val ref = data(3)
  val alt = data(4)
}

class tumor_only extends QScript with Logging {
  qscript =>

  @Argument(doc = "Memory Limit (GB)")
  var mem_limit: Int = _

  @Argument(doc="Skip sample indel realignment/genotyping")
  var skip_realignment: String = _

  @Argument(doc = "Base directory")
  var base_dir: String = _

  @Argument(doc="Output directory", required=true)
  var out_dir: File = _

  @Argument(doc="BAM files to genotype, can be specified more than once", required=true)
  var inp_bams: List[File] = _

  @Argument(doc="Interval file, should be in the bundle_2_8 directory of base_dir, and have available an expanded UNION version (e.g. suffixed with _UNION_500.bed)")
  var interval_file: File = "nextera_v1_2_for_gatk.bed"

  trait BasicArgs extends CommandLineGATK {
    this.ip = 500
    this.isr = IntervalSetRule.UNION
    this.R = new File(base_dir+"bundle_2_8/human_g1k_v37.fasta")
    this.L = Seq(new File(base_dir+"bundle_2_8/" + qscript.interval_file))
    this.jarFile = new File(base_dir+"programs/GenomeAnalysisTK.jar")
    this.memoryLimit = qscript.mem_limit
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

  class VarscanTumorOnly extends JavaCommandLineFunction{

    //http://varscan.sourceforge.net/support-faq.html#read-counts-different
    //under Warnings about "resetting normal" or "resetting tumor" file
    //per the faq, try combining the tumor and normal samples into one mpileup


    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc="input pileup", required=true)
    var pileup: File = _

    @Input(doc="a list of sample names in order, one per line")
    var sample_list: File = _

    @Argument(doc="whether to run it in indel or sn[vp] mode")
    var run_type: String = _

    @Output(doc="output vcf file")
    var out_vcf: File = _

    override def commandLine = super.commandLine + required("mpileup2"+run_type)  + required(pileup) +
      required("--variants", 1) +
      required("--output-vcf", 1) +
      required("--min-coverage", 8) +
      required("--min-var-freq", .2) +
      required("--min-freq-for-hom", .75) +
      required("--min-reads2", 2) +
      required("--strand-filter", 1) +
      required("--min-avg-qual", 15) +
      required("--p-value", .99) +
      required("--vcf-sample-list", sample_list) +
      required(">", escape=false) +
      required(out_vcf)

  }


  class Mpileup extends CommandLineFunction{

    //from http://varscan.sourceforge.net/somatic-calling.html
    //and
    //http://varscan.sourceforge.net/support-faq.html#read-counts-different
    //under Warnings about "resetting normal" or "resetting tumor" file
    //per the faq, try combining the tumor and normal samples into one mpileup

    @Input(doc="reference")
    var ref: File = _

    @Input(doc="intervals")
    var ints: File = _

    @Input(doc="bam files")
    var bam_files: List[File] = _

    @Output(doc="output pileup")
    var out_pileup: File = _

    //samtools mpileup -f [reference sequence] [BAM file(s)] >myData.mpileup

    override def commandLine = required(base_dir+"programs/samtools-1.1/samtools") + required("mpileup") + required("-q", 1) + required("-B") +
      required("-f", ref) + required("-l",ints) + required("-d", 8000) + repeat("", bam_files) +
      required(">", escape=false) + required(out_pileup)
  }

  class AddFakeNormalToMpileup extends CommandLineFunction{

    @Input(doc="")
    var orig_mpup: File = _

    @Output(doc="")
    var out_mpup: File = _

    override def commandLine = required("cat", orig_mpup) +
      required("|", escape=false) + required("perl", "-nle") + required("print \"$_\\t0\\t*\\t*\"", escape = true) +
      required(">", escape=false) + required(out_mpup)

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

  class ChangeVcfSamples extends CommandLineFunction{

    @Input(doc="")
    var inp_vcf: File = _

    @Output(doc="")
    var out_vcf: File = _

    override def commandLine = required("perl", "-ne") + required("if(/#CHROM/){s/none/No_Skin/} print") + required(inp_vcf) + required(">", escape=false) + required(out_vcf)

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

  def script() {

    if (out_dir.exists == false){
      out_dir.mkdirs()
    }

    val geno_base = new File(out_dir.getPath + "/single_genotypes")

    if (geno_base.exists == false){
      geno_base.mkdirs()
    }

    //indel realignment at the sample level

    for (file <- inp_bams) {

      val found_name = file.getName().split("[_\\.]")(1)

      var indel_r_out = file.getPath()

      if (skip_realignment != "skip"){

        val out_dir_file = new File(out_dir.getPath+"/single_alignments")

        if (! out_dir_file.exists()==false){
          out_dir_file.mkdirs()
        }

        val rtc_outfile = new File(out_dir_file.getPath() + "/" + file.getName.replace(".bam", ".rtc.list"))

        val rtc_stdout = rtc_outfile.getPath() + ".out"

        val rtc = new RealignerTargetCreator with BasicArgs
        rtc.jobOutputFile = rtc_stdout
        rtc.jobErrorFile = rtc_stdout.replace(".out", ".err")
        rtc.I = List(file)
        rtc.known = List(new File(base_dir + "bundle_2_8/Mills_and_1000G_gold_standard.indels.b37.vcf"),
          new File(base_dir + "bundle_2_8/1000G_phase1.indels.b37.vcf"))
        rtc.out = rtc_outfile
        rtc.wallTime = 1440
        add(rtc)

        //indelrealigner

        val indel_r = new IndelRealigner with BasicArgs

        indel_r_out = rtc_outfile.getPath().replace(".rtc.list", ".realign.bam")

        indel_r.jobOutputFile = indel_r_out.getPath() + ".out"
        indel_r.jobErrorFile = indel_r_out.getPath() + ".err"
        indel_r.targetIntervals = rtc_outfile
        indel_r.I = List(file)
        indel_r.model = org.broadinstitute.gatk.tools.walkers.indels.IndelRealigner.ConsensusDeterminationModel.USE_READS
        indel_r.known = List(new File(base_dir + "bundle_2_8/Mills_and_1000G_gold_standard.indels.b37.vcf"),
          new File(base_dir + "bundle_2_8/1000G_phase1.indels.b37.vcf"))
        indel_r.out = indel_r_out
        indel_r.wallTime = 1440
        indel_r.generate_md5 = true
        add(indel_r)

      }

      if (skip_realignment != "nogeno") {

        val mutect = new Mutect

        val mutect_out = new File(geno_base.getPath + "/Pair_" + found_name + "_AML_No_Skin.mutect.prev.vcf")

        val tagged_samps = Seq(
          new TaggedFile(indel_r_out, "tumor")
        )

        mutect.jobOutputFile = mutect_out + ".out"
        mutect.jobErrorFile = mutect_out + ".err"
        mutect.I = tagged_samps
        mutect.tumor_name = found_name + "_AML"
        mutect.normal_name = "No_Skin"
        mutect.out = new File(mutect_out.getPath().replace(".vcf", ".stats"))
        mutect.vcf = mutect_out
        mutect.max_n_alleles = "10000000"
        mutect.max_n_allele_freq = "1.0"
        mutect.wallTime = 1440
        add(mutect)

        //despite specifying the sample name, the output file still has 'none' as the sample name of the normal--which is true, but we want
        //it to be No_Skin

        val change_mutect = new ChangeVcfSamples

        //the 'filter' suffix doesn't denote any special filters that have been applied, just to simplify locating the final files later

        val change_mutect_out = new File(mutect_out.getPath.replace(".prev", ""))
        change_mutect.jobOutputFile = change_mutect_out.getPath + ".out"
        change_mutect.jobErrorFile = change_mutect_out.getPath + ".err"
        change_mutect.inp_vcf = mutect_out
        change_mutect.out_vcf = change_mutect_out
        change_mutect.wallTime = 10
        change_mutect.memoryLimit = 5
        add(change_mutect)

        //filter to only PASS variants

        val mut_sels = new SelectVariants with BasicArgs

        val mut_sels_out = new File(change_mutect_out.replace(".vcf", ".filter.vcf"))

        mut_sels.jobOutputFile = mut_sels_out.getPath + ".out"
        mut_sels.jobErrorFile = mut_sels_out.getPath + ".err"
        mut_sels.variant = change_mutect_out
        mut_sels.out = mut_sels_out
        mut_sels.excludeFiltered = true
        mut_sels.wallTime = 10

        add(mut_sels)

        //Start VarScan

        val sam_tumor = new Mpileup
        val sam_tumor_out = new File(geno_base.getPath + "/Sample_" + found_name + "_AML.mpileup")

        sam_tumor.jobOutputFile = sam_tumor_out.getPath + ".out"
        sam_tumor.jobErrorFile = sam_tumor_out.getPath + ".err"

        sam_tumor.bam_files = List(indel_r_out)
        sam_tumor.ref = new File(base_dir + "bundle_2_8/human_g1k_v37.fasta")
        sam_tumor.out_pileup = sam_tumor_out
        sam_tumor.ints = new File(base_dir + "bundle_2_8/" + qscript.interval_file.replace(".bed", "_UNION_500.bed"))
        sam_tumor.wallTime = 120
        sam_tumor.memoryLimit = qscript.mem_limit
        add(sam_tumor)

        //make a fake mpileup that will represent the 'normal' as well

        val sam_norm = new AddFakeNormalToMpileup
        val sam_norm_out = new File(geno_base.getPath + "/Pair_" + found_name + "_AML_No_Skin.mpileup")

        sam_norm.jobOutputFile = sam_norm_out.getPath + ".out"
        sam_norm.jobErrorFile = sam_norm_out.getPath + ".err"
        sam_norm.orig_mpup = sam_tumor_out
        sam_norm.out_mpup = sam_norm_out
        add(sam_norm)

        val v_sample_list = geno_base.getPath + "/Sample_" + found_name + "_varscan_sample_list.txt"

        val v_sample_pw = new java.io.PrintWriter(v_sample_list)

        v_sample_pw.write("TUMOR\nNORMAL\n")

        v_sample_pw.close()

        //run varscan2 on the t/n pileups


        for (v_type <- List("snp", "indel")) {

          val var_geno = new VarscanTumorOnly

          val var_geno_out = new File(sam_norm_out.getPath.replace("mpileup", "varscan2.t.100.n.100." + v_type + ".vcf"))
          var_geno.jobOutputFile = var_geno_out.getPath + ".out"
          var_geno.jobErrorFile = var_geno_out.getPath + ".err"

          var_geno.run_type = v_type
          var_geno.pileup = sam_norm_out
          var_geno.out_vcf = var_geno_out
          var_geno.sample_list = new File(v_sample_list)
          var_geno.wallTime = 120
          var_geno.memoryLimit = qscript.mem_limit
          add(var_geno)

          //run bamreadcount and fpfilter on each sample

          val interval_file = new File(var_geno_out.getPath().replace(".vcf", ".interval"))

          val vcf_pos = new VarscanGetPos

          vcf_pos.input_file = var_geno_out
          vcf_pos.output_file = interval_file
          add(vcf_pos)

          val brc = new VarscanBamReadCount
          val brc_out = new File(var_geno_out.getPath().replace(".vcf", ".brc"))

          brc.jobOutputFile = brc_out.getPath() + ".out"
          brc.jobErrorFile = brc_out.getPath() + ".err"
          brc.input_bam = indel_r_out
          brc.input_intervals = interval_file
          brc.output_readcount = brc_out
          brc.ref = new File(base_dir + "bundle_2_8/human_g1k_v37.fasta")
          brc.wallTime = 120
          brc.memoryLimit = qscript.mem_limit
          add(brc)

          val fpfilter = new VarscanFpfilter

          val fpfilter_out = new File(var_geno_out.getPath().replace(".vcf", ".flag.vcf"))
          fpfilter.jobOutputFile = fpfilter_out.getPath() + ".out"
          fpfilter.jobErrorFile = fpfilter_out.getPath() + ".err"
          fpfilter.input_readcount = brc_out
          fpfilter.input_vcf = var_geno_out
          fpfilter.outfile = fpfilter_out
          fpfilter.wallTime = 10
          fpfilter.memoryLimit = 10
          add(fpfilter)

          //fix the varscan vcf file prior to feeding it into GATK tools...

          val gene_filt_vep = new FilterVarscan

          val gene_filt_vep_out = new File(fpfilter_out.getPath.replace(".flag.vcf", ".flag.fixed.vcf"))

          gene_filt_vep.jobOutputFile = gene_filt_vep_out + ".out"
          gene_filt_vep.jobErrorFile = gene_filt_vep_out + ".err"
          gene_filt_vep.vcf_input = fpfilter_out
          gene_filt_vep.vcf_output = gene_filt_vep_out
          gene_filt_vep.wallTime = 5
          add(gene_filt_vep)

          var varscan_to_filter = gene_filt_vep_out
          var rep_suffix = ".flag.fixed.vcf"

          if (v_type == "indel") {

            //also left align and trim if indels
            val lalign = new LeftAlignAndTrimVariants

            val lalign_out = new File(gene_filt_vep_out.getPath.replace(".vcf", ".la.vcf"))

            lalign.jobOutputFile = lalign_out.getPath + ".out"
            lalign.jobErrorFile = lalign_out.getPath + ".err"
            lalign.R = new File(base_dir + "bundle_2_8/human_g1k_v37.fasta")
            lalign.variant = gene_filt_vep_out
            lalign.o = lalign_out
            lalign.splitMultiallelics = true
            lalign.wallTime = 24
            lalign.memoryLimit = qscript.mem_limit
            add(lalign)

            varscan_to_filter = lalign_out
            rep_suffix = ".flag.fixed.la.vcf"
          }

          //again limit to PASS variants

          val var_sels = new SelectVariants with BasicArgs

          val var_sels_out = new File(varscan_to_filter.replace(rep_suffix, ".filter.vcf"))

          var_sels.jobOutputFile = var_sels_out.getPath + ".out"
          var_sels.jobErrorFile = var_sels_out.getPath + ".err"
          var_sels.variant = varscan_to_filter
          var_sels.out = var_sels_out
          var_sels.excludeFiltered = true
          var_sels.wallTime = 10

          add(var_sels)

        }

      }
      }
  }

}