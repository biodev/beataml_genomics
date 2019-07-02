import java.io.File
import scala.io.Source

import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.extensions.gatk.LeftAlignAndTrimVariants
import org.broadinstitute.gatk.queue.util.Logging

import scala.io.Source

case class VCFPos(line: String) {
  val data = line.split("\t")
  val chrom = data(0)
  val pos_start = data(1)
  val ref = data(3)
  val alt = data(4)
}

case class SampleFilename (path: File){
  val split_str = path.getName.split("[_\\.]")

  val file = path
  val case_samp = split_str.slice(1,3).mkString("_")
  val control_samp = split_str.slice(3,5).mkString("_")
}

class tumor_only_varscan extends QScript with Logging {
  qscript =>

  @Argument(doc = "File of processed mpileups", required = true)
  var mpup_file: File = _

  @Argument(doc = "File of BAM files corresponding to mpileups", required = true)
  var bam_file: File = _

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

  @Argument(doc = "Use DNA or RNA fp filtering")
  var filter_type: String = _

  @Argument(doc = "reference path relative to base director (default: bundle_2_8/human_g1k_v37.fasta)", required = false)
  var ref_path: String = "bundle_2_8/human_g1k_v37.fasta"

  class FilterVarscan extends CommandLineFunction {

    @Input(doc = "input vcf file")
    var vcf_input: File = _

    @Output(doc = "out filtered vcf file")
    var vcf_output: File = _


    def commandLine = required("grep", "-v") +
      required(";SS=5;") +
      required(vcf_input) +
      required(">", escape = false) + required(vcf_output)
  }

  class VarscanFpfilterDNA extends JavaCommandLineFunction {

    //parameters from: https://github.com/dkoboldt/varscan/blob/master/VarScan.v2.4.1.description.txt

    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc = "")
    var input_vcf: File = _

    @Input(doc = "")
    var input_readcount: File = _

    @Output(doc = "")
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
      required("--max-ref-mmqs", 50) +
      required("--min-ref-avgrl", 0) +
      required("--min-var-avgrl", 0)
  }

  class VarscanFpfilterRNA extends JavaCommandLineFunction{

    //parameters from: https://github.com/dkoboldt/varscan/blob/master/VarScan.v2.4.1.description.txt

    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc="")
    var input_vcf: File = _

    @Input(doc="")
    var input_readcount: File = _

    @Output(doc="")
    var outfile: File = _

    // For RNA genotyping the min variant base quality is reduced to 21 and the read-length filters are disabled.

    override def commandLine = super.commandLine + required("fpfilter") + required(input_vcf) + required(input_readcount) + required("--output-file", outfile) +
      required("--keep-failures") +
      required("--min-var-count", 3) +
      required("--min-var-count-lc", 1) +
      required("--min-strandedness", 0) +
      required("--min-var-basequal", 21) +
      required("--min-ref-readpos", 0.20) +
      required("--min-ref-dist3", 0.20) +
      required("--min-var-readpos", 0.15) +
      required("--min-var-dist3", 0.15) +
      required("--min-ref-avgrl", 0) +
      required("--min-var-avgrl", 0) +
      required("--max-rl-diff", 1) +
      required("--max-mapqual-diff", 10) +
      required("--min-ref-mapqual", 20) +
      required("--min-var-mapqual", 30) +
      required("--max-var-mmqs", 100) +
      required("--max-ref-mmqs", 50)
  }


  class VarscanTumorOnly extends JavaCommandLineFunction {

    //http://varscan.sourceforge.net/support-faq.html#read-counts-different
    //under Warnings about "resetting normal" or "resetting tumor" file
    //per the faq, try combining the tumor and normal samples into one mpileup


    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc = "input pileup", required = true)
    var pileup: File = _

    @Input(doc = "a list of sample names in order, one per line")
    var sample_list: File = _

    @Argument(doc = "whether to run it in indel or sn[vp] mode")
    var run_type: String = _

    @Output(doc = "output vcf file")
    var out_vcf: File = _

    override def commandLine = super.commandLine + required("mpileup2" + run_type) + required(pileup) +
      required("--variants", 1) +
      required("--output-vcf", 1) +
      required("--min-coverage", 3) +
      required("--min-var-freq", .08) +
      required("--min-freq-for-hom", .75) +
      required("--min-reads2", 1) +
      required("--strand-filter", 0) +
      required("--min-avg-qual", 15) +
      required("--p-value", 0.10) +
      required("--vcf-sample-list", sample_list) +
      required(">", escape = false) +
      required(out_vcf)

  }

  class SelectVars extends CommandLineFunction {

    @Input(doc = "")
    var variant: File = _

    @Output(doc = "")
    var out: File = _

    def commandLine = required(qscript.base_dir + "/programs/bcftools/bcftools") + required("view") +
      required("-f", "PASS") +
      required("-o", out) +
      required("-O", "v") +
      required(variant)

  }

  class Vcf2MafRunner extends CommandLineFunction {

    @Input(doc = "Input VCF file") var vcf_input: File = _
    @Output(doc = "Output MAF file") var maf_output: File = _

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


  def script() {

    val out_base_file = new File(qscript.out_base)

    if (out_base_file.exists() == false) {
      out_base_file.mkdirs()
    }

    val base_vcf_dir = new File(qscript.out_base + "/vcfs")

    if (base_vcf_dir.exists() == false) {
      base_vcf_dir.mkdirs()
    }

    val base_maf_dir = new File(qscript.out_base + "/mafs")

    if (base_maf_dir.exists() == false) {
      base_maf_dir.mkdirs()
    }

    val input_mpups =  Source.fromFile(mpup_file) getLines() map (l => SampleFilename(l))

    val align_files = Source.fromFile(bam_file).getLines.toList

    val sample_name = align_files map (_.getName().split("[_\\.]")(1) + "_AML")

    val file_map = (sample_name zip align_files).toMap


    for (pup_file <- input_mpups) {

      if (file_map.contains(pup_file.case_samp) == false){
        throw new QException("Corresponding BAM file not found for " + pup_file.case_samp)
      }

      val v_sample_list = base_vcf_dir.getPath + "/Sample_" + pup_file.case_samp + "_varscan_sample_list.txt"

      val v_sample_pw = new java.io.PrintWriter(v_sample_list)

      v_sample_pw.write("TUMOR\nNORMAL\n")

      v_sample_pw.close()

      //run varscan2 on the t/n pileups


      for (v_type <- List("snp", "indel")) {

        val var_geno = new VarscanTumorOnly

        val var_geno_out = new File(base_vcf_dir.getPath + "/" + pup_file.file.getName.replace("mpileup", "varscan2.t.100.n.100." + v_type + ".vcf"))
        var_geno.jobOutputFile = var_geno_out.getPath + ".out"
        var_geno.jobErrorFile = var_geno_out.getPath + ".err"

        var_geno.run_type = v_type
        var_geno.pileup = pup_file.file
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
        brc.input_bam = file_map(pup_file.case_samp)
        brc.input_intervals = interval_file
        brc.output_readcount = brc_out
        brc.ref = new File(base_dir + "bundle_2_8/human_g1k_v37.fasta")
        brc.wallTime = 120
        brc.memoryLimit = qscript.mem_limit
        add(brc)

        val fpfilter_out = new File(var_geno_out.getPath().replace(".vcf", ".flag.vcf"))


        if (filter_type == "rna"){

          val fpfilter_r = new VarscanFpfilterRNA

          fpfilter_r.jobOutputFile = fpfilter_out.getPath() + ".out"
          fpfilter_r.jobErrorFile = fpfilter_out.getPath() + ".err"
          fpfilter_r.input_readcount = brc_out
          fpfilter_r.input_vcf = var_geno_out
          fpfilter_r.outfile = fpfilter_out
          fpfilter_r.wallTime = 10
          fpfilter_r.memoryLimit = 10
          add(fpfilter_r)

        }else{
          val fpfilter_d = new VarscanFpfilterDNA

          fpfilter_d.jobOutputFile = fpfilter_out.getPath() + ".out"
          fpfilter_d.jobErrorFile = fpfilter_out.getPath() + ".err"
          fpfilter_d.input_readcount = brc_out
          fpfilter_d.input_vcf = var_geno_out
          fpfilter_d.outfile = fpfilter_out
          fpfilter_d.wallTime = 10
          fpfilter_d.memoryLimit = 10
          add(fpfilter_d)
        }

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

        val var_sels = new SelectVars

        val var_sels_out = new File(varscan_to_filter.replace(rep_suffix, ".filter.vcf"))

        var_sels.jobOutputFile = var_sels_out.getPath + ".out"
        var_sels.jobErrorFile = var_sels_out.getPath + ".err"
        var_sels.variant = varscan_to_filter
        var_sels.out = var_sels_out
        var_sels.wallTime = 10

        add(var_sels)

        val all_maf = new VarscanVcf2MafRunner

        var all_maf_out = new File(base_maf_dir + "/" + var_sels_out.getName.replace(".vcf", ".maf"))

        all_maf.jobOutputFile = all_maf_out + ".out"
        all_maf.jobErrorFile = all_maf_out + ".err"

        all_maf.vcf_input = var_sels_out
        all_maf.maf_output = all_maf_out

        all_maf.norm_id = pup_file.control_samp
        all_maf.tumor_id = pup_file.case_samp
        all_maf.wallTime = 24
        all_maf.memoryLimit = qscript.mem_limit
        all_maf.nCoresRequest = qscript.num_threads
        add(all_maf)

      }
    }
  }
}