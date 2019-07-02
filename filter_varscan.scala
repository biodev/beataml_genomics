import java.io.{File, PrintWriter}

import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction
import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.util.Logging

import scala.io.Source


case class SampleFilename (path: File){
  val split_str = path.getName.split("[_\\.]")

  val case_samp = split_str.slice(1,3).mkString("_")
  val control_samp = split_str.slice(3,5).mkString("_")
}


class filter_varscan extends QScript with Logging {
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

  class VarscanFpfilter extends JavaCommandLineFunction {

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

  class VarscanProcessSomatic extends JavaCommandLineFunction {

    jarFile = new File(base_dir + "/programs/VarScan.v2.4.1p1.jar")

    @Input(doc = "input variant file")
    var input_var: File = _

    @Output(doc = "high confidence variant file")
    var output_hc_var: File = _

    override def commandLine = super.commandLine + required("processSomatic") + required(input_var) +
      required("--min-tumor-freq", .05) + // - Minimum variant allele frequency in tumor [0.10]
      optional("--max-normal-freq", 1) + // - Maximum variant allele frequency in normal [0.05]
      required("--p-value", .1) // - P-value for high-confidence calling
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

  def script {

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

    for (i <- input_vcfs) {

      val sample_names = SampleFilename(i)

      var gene_filt_vep_out = new File(i.getPath())

      if (i.getName.contains("indel") == false) {

        val fpfilter = new VarscanFpfilter

        val brc_out = new File(i.getPath().replace(".vcf", ".brc"))
        val fpfilter_out = new File(base_vcf_dir + "/" + i.getName().replace(".vcf", ".flag.vcf"))

        fpfilter.jobOutputFile = fpfilter_out.getPath() + ".out"
        fpfilter.jobErrorFile = fpfilter_out.getPath() + ".err"
        fpfilter.input_readcount = brc_out
        fpfilter.input_vcf = i
        fpfilter.outfile = fpfilter_out
        fpfilter.wallTime = 10
        fpfilter.memoryLimit = 10
        add(fpfilter)

        //filter out the unknown category as it causes issues with VCF standard

        val gene_filt_vep = new FilterVarscan

        gene_filt_vep_out = new File(fpfilter_out.getPath.replace(".flag.vcf", ".flag.fixed.vcf"))

        gene_filt_vep.jobOutputFile = gene_filt_vep_out + ".out"
        gene_filt_vep.jobErrorFile = gene_filt_vep_out + ".err"
        gene_filt_vep.vcf_input = fpfilter_out
        gene_filt_vep.vcf_output = gene_filt_vep_out
        gene_filt_vep.wallTime = 5
        add(gene_filt_vep)

      }

      val var_sels = new SelectVars

      val var_sels_out = new File(base_vcf_dir + "/" + gene_filt_vep_out.getName().replace(".vcf", ".filter.vcf"))

      var_sels.jobOutputFile = var_sels_out.getPath.replace(".vcf", ".out")
      var_sels.jobErrorFile = var_sels_out.getPath.replace(".vcf", ".err")
      var_sels.variant = gene_filt_vep_out
      var_sels.out = var_sels_out
      var_sels.wallTime = 10
      add(var_sels)

      val all_maf = new VarscanVcf2MafRunner

      var all_maf_out = new File(base_maf_dir + "/" + var_sels_out.getName.replace(".vcf", ".maf"))

      all_maf.jobOutputFile = all_maf_out + ".out"
      all_maf.jobErrorFile = all_maf_out + ".err"

      all_maf.vcf_input = var_sels_out
      all_maf.maf_output = all_maf_out

      all_maf.norm_id = sample_names.control_samp
      all_maf.tumor_id = sample_names.case_samp
      all_maf.wallTime = 24
      all_maf.memoryLimit = qscript.mem_limit
      all_maf.nCoresRequest = qscript.num_threads
      add(all_maf)


      if (sample_names.control_samp != "No_Skin") {

        val var_filt = new VarscanProcessSomatic

        val var_hc_out = new File(var_sels_out.getPath().replace(".vcf", ".Somatic.hc.vcf"))

        var_filt.jobOutputFile = var_hc_out.getPath() + ".out"
        var_filt.jobErrorFile = var_hc_out.getPath() + ".err"
        var_filt.input_var = var_sels_out
        var_filt.output_hc_var = var_hc_out
        var_filt.wallTime = 10
        var_filt.memoryLimit = 2
        add(var_filt)

        val som_maf = new VarscanVcf2MafRunner

        var som_maf_out = new File(base_maf_dir + "/" + var_hc_out.getName.replace(".vcf", ".maf"))

        som_maf.jobOutputFile = som_maf_out + ".out"
        som_maf.jobErrorFile = som_maf_out + ".err"

        som_maf.vcf_input = var_hc_out
        som_maf.maf_output = som_maf_out

        som_maf.norm_id = sample_names.control_samp
        som_maf.tumor_id = sample_names.case_samp
        som_maf.wallTime = 24
        som_maf.memoryLimit = qscript.mem_limit
        som_maf.nCoresRequest = qscript.num_threads
        add(som_maf)

      }

    }
  }
}
