import java.io.{File, PrintWriter}

import scala.io._
import org.broadinstitute.gatk.queue.{QException, QScript}
import org.broadinstitute.gatk.queue.util.Logging


case class SamplePairs(line: String) {
  val data = line.split(",")
  val patient_id = data(0)
  val AML = data(1).replace("\"", "")
  val Skin = data(2).replace("\"", "")
}

class pindel_pipeline extends QScript with Logging {
  qscript =>

  @Argument(doc = "Memory Limit (GB)")
  var mem_limit: Int = _

  @Argument(doc = "Number of Threads")
  var num_threads: Int = _

  @Argument(doc="Base directory")
  var base_dir: String = _

  @Argument(doc="Tumor/normal pair CSV file (relative to specified Base directory)")
  var tn_file: File  = _

  @Argument(doc="File of insert sizes")
  var ins_file: File  = _

  @Argument(doc="output directory")
  var out_dir : File = _

  @Argument(doc="chromosome")
  var pindel_chr: String = _

  @Argument(doc="bam files", required=false)
  var bam_files: List[File] = _

  @Argument(doc="reference path relative to base director (default: bundle_2_8/human_g1k_v37.fasta)", required=false)
  var ref_path: String = "bundle_2_8/human_g1k_v37.fasta"

  @Argument(doc="Path to Rscript", required=false)
  var rscript_path: String = "Rscript"

  @Argument(doc="Path to pindel executables")
  var pindel_path: String = _

  @Argument(doc="BAM suffix", required=false)
  var bam_suffix: String = ".bam"

  @Argument(doc="Ensembl tools path")
  var ens_tools_path: String =  _



  class Rscript extends CommandLineFunction{

    @Input(doc="")
    var inp_files: List[File] = _

    @Output(doc="")
    var out_files: List[File] = _

    var source_code: String = ""

    override def commandLine = required("echo", source_code) + required("|", escape=false) +
      required(qscript.rscript_path) + required("--vanilla") + required("-", escape = false) + repeat(inp_files)

  }

  class Pindel extends CommandLineFunction{

    @Input(doc="")
    var conf_file: File = _

    @Output(doc="")
    var si_out: File = _

    @Argument(doc="")
    var chr: String = _

    override def commandLine = required(pindel_path+"/pindel") + required("-T", num_threads-1) +
      required("-f", base_dir+ref_path) + required("-i", conf_file) + required("-c", chr) + required("-o", si_out.getPath().replace("_SI", ""))

  }

  class PindelVcf extends CommandLineFunction{

    @Input(doc="")
    var pindel_si: File = _

    @Output(doc="")
    var vcf_out: File = _

    override def commandLine = required(pindel_path+"/pindel2vcf") + required("-G") + required("-P", pindel_si.replace("_SI", "")) +
      required("-r", base_dir+ref_path) + required("-R", "gatk_ref") + required("-d", "20101123") + required("-v", vcf_out)


  }

  class SeperateVcf extends Rscript {

    source_code =
      """
      |    library(VariantAnnotation)
      |
      |    use.vcf <- commandArgs(trailingOnly = T)
      |
      |    temp <- readVcf(use.vcf, genome="test")
      |
      |    simple <- temp[info(temp)$SVTYPE %in% c("DEL", "INS")]
      |
      |    #remove the structural variation components so vep is more featured
      |      info(simple)$SVTYPE <- NULL
      |
      |    writeVcf(simple, file=sub(".vcf", ".simple.vcf", use.vcf))
      |
      |    complex <- temp[info(temp)$SVTYPE %in% c("DEL", "INS")==F]
      |
      |    #probably need to recode DUP:TANDEM to TDUP for VEP
      |
      |    info(complex)$SVTYPE <- ifelse(info(complex)$SVTYPE == "DUP:TANDEM", "TDUP",info(complex)$SVTYPE)
      |
      |    writeVcf(complex, file=sub(".vcf", ".complex.vcf", use.vcf))
      |
      |
    """.stripMargin

    }

    class SelectPindel extends Rscript {

      @Argument(doc="")
      var tumor_id: String = _

      @Argument(doc="")
      var var_type: String = _

      source_code =
        """
          |      library(VariantAnnotation)
          |
          |      use.vcf <- commandArgs(trailingOnly = T)
          |
          |      temp <- readVcf(use.vcf[1], genome="test")
          |
          |      sub.temp <- temp[geno(temp)$AD[,use.vcf[3],2] >= 5]
          |
          |      if (use.vcf[2] == "flt3"){
          |        sub.temp <- subsetByOverlaps(sub.temp[info(sub.temp)$SVTYPE == "TDUP"], GRanges(seqnames="13", ranges=IRanges(start=28577411, end=28674729)))
          |      }
          |
          |      writeVcf(sub.temp, file=use.vcf[4])
          |
          |
        """.stripMargin


      override def commandLine = super.commandLine + required(var_type) + required(tumor_id) + repeat(out_files)

    }

  class LeftAlign extends CommandLineFunction{

    @Input(doc="")
    var vcf_input: File = _

    @Output(doc="")
    var lalign_out: File = _

    override def commandLine = required(base_dir + "/programs/bcftools/bcftools") + required("norm") + required("-c", "e") +
      required("-f", base_dir+ref_path) + required("-m", "-any") + required("-o", lalign_out) + required("-O", "v") + required(vcf_input)

  }

  class Vcf2MafRunner extends CommandLineFunction{

    @Input(doc="Input VCF file") var vcf_input: File = _
    @Output(doc="Output MAF file") var maf_output: File = _

    @Argument(doc="tumor sample name", required = true)
    var tumor_id: String = _

    @Argument(doc="normal sample name", required = false)
    var norm_id: String = _

    //basic usage: perl vcf2maf.pl --input-vcf data/test.vcf --output-maf data/test.maf --tumor-id WD1309 --normal-id NB1308

    this.memoryLimit = qscript.mem_limit

    def commandLine = required("perl") + required(ens_tools_path.getParent() + "/vcf2maf/vcf2maf.pl") +
      required("--ref-fasta",qscript.base_dir + "vep/vep/current.fasta" ) +
      required("--vep-path", ens_tools_path + "/scripts/variant_effect_predictor/") +
      required("--vep-data", qscript.base_dir + "vep/vep") +
      required("--input-vcf", vcf_input) +
      required("--output-maf", maf_output) +
      required("--tumor-id", tumor_id) +
      optional("--normal-id", norm_id)

  }


  def script {

    val tn_pairs = Source.fromFile(tn_file) getLines() map (l => SamplePairs(l))

    val sample_name = bam_files map (_.getName().replace(qscript.bam_suffix, "").split("_")(1))
    val file_map = (sample_name zip bam_files).toMap

    val ins_size_content = Source.fromFile(ins_file).getLines.map(_.split("\t"))
    val ins_size_header = ins_size_content.next
    val insert_size_list = ins_size_content.map(ins_size_header.zip(_).toMap).toList
    val insert_size_samples = insert_size_list.map(_("SAMPLE"))
    val insert_size_map = insert_size_samples.zip(insert_size_list).toMap

    if (out_dir.exists ==false){
      out_dir.mkdirs()
    }

    for (i <- tn_pairs){
      val base_out = out_dir.getPath+"/Pair_%s_AML_%s_Skin".format(i.AML, i.Skin match {case "" => "No"; case x: String => x})
      val conf_file = new File(base_out + "_pindel.conf")

      val conf_pw = new PrintWriter(conf_file.getPath)

      if (file_map.contains(i.AML) && insert_size_map.contains(i.AML) && i.Skin==""){

        //tumor-only
        conf_pw.write("%s\t%d\t%s\n".format(file_map(i.AML).getPath,Math.max(insert_size_map(i.AML)("MEAN_INSERT_SIZE").toDouble, 100.0).round, i.AML+"_AML"))


      }else if (file_map.contains(i.AML) && insert_size_map.contains(i.AML) && file_map.contains(i.Skin) && insert_size_map.contains(i.Skin) ){

        //tumor-normal pairs
        conf_pw.write("%s\t%d\t%s\n".format(file_map(i.Skin).getPath,Math.max(insert_size_map(i.Skin)("MEAN_INSERT_SIZE").toDouble, 100.0).round, i.Skin+"_Skin"))
        conf_pw.write("%s\t%d\t%s\n".format(file_map(i.AML).getPath,Math.max(insert_size_map(i.AML)("MEAN_INSERT_SIZE").toDouble, 100.0).round, i.AML+"_AML"))

      }else{
	println(i.AML)
	println(i.Skin)
	println(file_map.contains(i.AML))
	println(insert_size_map.contains(i.AML))
	println(file_map.contains(i.Skin))
	println(insert_size_map.contains(i.Skin))
        throw new QException("Can't find all BAM files or insert size entries implied by the t/n pairing file")
      }

      conf_pw.close()

      val cur_pind = new Pindel

      val cur_pind_out = new File(conf_file.getPath.replace(".conf", "_SI"))

      cur_pind.jobOutputFile = cur_pind_out.getPath + ".out"
      cur_pind.jobErrorFile = cur_pind_out.getPath + ".err"
      cur_pind.nCoresRequest = num_threads
      cur_pind.memoryLimit = mem_limit
      cur_pind.conf_file = conf_file
      cur_pind.si_out = cur_pind_out
      cur_pind.chr = pindel_chr
      cur_pind.wallTime = 240
      add(cur_pind)

      val pind_vcf = new PindelVcf

      val pind_vcf_out = new File(cur_pind_out.getPath.replace("_SI", ".vcf"))

      pind_vcf.jobOutputFile = pind_vcf_out.getPath + ".out"
      pind_vcf.jobErrorFile = pind_vcf_out.getPath + ".err"
      pind_vcf.nCoresRequest = num_threads
      pind_vcf.memoryLimit = mem_limit
      pind_vcf.pindel_si = cur_pind_out
      pind_vcf.vcf_out  = pind_vcf_out
      pind_vcf.wallTime = 30
      add(pind_vcf)

      val sep_vcf = new SeperateVcf

      val sep_vcf_out = new File(pind_vcf_out.getPath.replace(".vcf", ".simple.vcf"))

      val sep_out_list = List(sep_vcf_out, new File(pind_vcf_out.getPath.replace(".vcf", ".complex.vcf")))

      sep_vcf.jobOutputFile = sep_vcf_out.getPath + ".out"
      sep_vcf.jobErrorFile = sep_vcf_out.getPath + ".err"
      sep_vcf.nCoresRequest = num_threads
      sep_vcf.memoryLimit = mem_limit
      sep_vcf.inp_files = List(pind_vcf_out)
      sep_vcf.out_files  = sep_out_list
      sep_vcf.wallTime = 30
      add(sep_vcf)

      for (j <- sep_out_list){

        val sel_pind = new SelectPindel

        val sel_pind_out = new File(j.getPath.replace(".vcf", ".filt.vcf"))

        sel_pind.jobOutputFile = sel_pind_out.getPath + ".out"
        sel_pind.jobErrorFile = sel_pind_out.getPath + ".err"
        sel_pind.nCoresRequest = num_threads
        sel_pind.memoryLimit = mem_limit
        sel_pind.tumor_id = i.AML + "_AML"
        sel_pind.var_type = j.getName match { case x if x.contains("simple")=>"simple"; case y if y.contains("complex")=>"flt3"}
        sel_pind.inp_files = List(j)
        sel_pind.out_files = List(sel_pind_out)
        sel_pind.wallTime = 30
        add(sel_pind)

        var lalign_out = sel_pind_out

        if (j.getName.contains("simple")){

          //left align if only simple variants

          val lalign = new LeftAlign

          lalign_out = new File(sel_pind_out.getPath.replace(".vcf", ".la.vcf"))

          lalign.jobOutputFile = lalign_out.getPath + ".out"
          lalign.jobErrorFile = lalign_out.getPath + ".err"
          lalign.nCoresRequest = num_threads
          lalign.memoryLimit = mem_limit
          lalign.vcf_input = sel_pind_out
          lalign.lalign_out = lalign_out
          lalign.wallTime = 30
          add(lalign)


        }

        //conversion to maf
        var mut_maf = new Vcf2MafRunner

        var mut_maf_out = new File(lalign_out.getPath.replace(".vcf", ".maf"))

        mut_maf.jobOutputFile = mut_maf_out.getPath + ".out"
        mut_maf.jobErrorFile = mut_maf_out.getPath + ".err"

        mut_maf.vcf_input = lalign_out
        mut_maf.maf_output = mut_maf_out

        mut_maf.norm_id = "%s_Skin".format(i.Skin match {case "" => "No"; case x: String => x})
        mut_maf.tumor_id = "%s_AML".format(i.AML)
        mut_maf.wallTime = 120
        mut_maf.memoryLimit = qscript.mem_limit
        mut_maf.nCoresRequest = qscript.num_threads
        add(mut_maf)


      }


    }


  }

}
