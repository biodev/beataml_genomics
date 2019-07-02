library(data.table)
library(stringr)

.my.short.prots <- c(Ala="A", Arg="R", Asn="N", Asp="D", Asx ="B", Cys="C", Glu="E", Gln ="Q", Glx="Z", Gly="G", His= "H", Ile ="I", Leu= "L",
                    Lys ="K", Met= "M", Phe ="F", Pro= "P", Ser= "S", Thr= "T", Trp= "W", Tyr= "Y", Val= "V", Xxx ="X", Ter ="*" )


#taken from http://dec2015.archive.ensembl.org/info/genome/variation/predicted_data.html#consequences
.ens.cons.tab <- "transcript_ablation 	A feature ablation whereby the deleted region includes a transcript feature 	SO:0001893 	Transcript ablation
splice_acceptor_variant 	A splice variant that changes the 2 base region at the 3' end of an intron 	SO:0001574 	Essential splice site
splice_donor_variant 	A splice variant that changes the 2 base region at the 5' end of an intron 	SO:0001575 	Essential splice site
stop_gained 	A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript 	SO:0001587 	Stop gained
frameshift_variant 	A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three 	SO:0001589 	Frameshift coding
stop_lost 	A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript 	SO:0001578 	Stop lost
start_lost  A codon variant that changes at least one base of the canonical start codo	SO:0002012	Start lost
transcript_amplification 	A feature amplification of a region containing a transcript 	SO:0001889 	Transcript amplification
inframe_insertion 	An inframe non synonymous variant that inserts bases into in the coding sequence 	SO:0001821 	Non synonymous coding
inframe_deletion 	An inframe non synonymous variant that deletes bases from the coding sequence 	SO:0001822 	Non synonymous coding
missense_variant 	A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved 	SO:0001583 	Non synonymous coding
protein_altering_variant	A sequence_variant which is predicted to change the protein encoded in the coding sequence	SO:0001818	Protein altering variant
splice_region_variant 	A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron 	SO:0001630 	Splice site
incomplete_terminal_codon_variant 	A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed 	SO:0001626 	Partial codon
stop_retained_variant 	A sequence variant where at least one base in the terminator codon is changed, but the terminator remains 	SO:0001567 	Synonymous coding
synonymous_variant 	A sequence variant where there is no resulting change to the encoded amino acid 	SO:0001819 	Synonymous coding
coding_sequence_variant 	A sequence variant that changes the coding sequence 	SO:0001580 	Coding unknown
mature_miRNA_variant 	A transcript variant located with the sequence of the mature miRNA 	SO:0001620 	Within mature miRNA
5_prime_UTR_variant 	A UTR variant of the 5' UTR 	SO:0001623 	5prime UTR
3_prime_UTR_variant 	A UTR variant of the 3' UTR 	SO:0001624 	3prime UTR
non_coding_transcript_exon_variant 	A sequence variant that changes non-coding exon sequence in a non-coding transcript 	SO:0001792 	Within non coding gene
intron_variant 	A transcript variant occurring within an intron 	SO:0001627 	Intronic
NMD_transcript_variant 	A variant in a transcript that is the target of NMD 	SO:0001621 	NMD transcript
non_coding_transcript_variant 	A transcript variant of a non coding RNA gene 	SO:0001619 	Within non coding gene
upstream_gene_variant 	A sequence variant located 5' of a gene 	SO:0001631 	Upstream
downstream_gene_variant 	A sequence variant located 3' of a gene 	SO:0001632 	Downstream
TFBS_ablation 	A feature ablation whereby the deleted region includes a transcription factor binding site 	SO:0001895 	Tfbs ablation
TFBS_amplification 	A feature amplification of a region containing a transcription factor binding site 	SO:0001892 	Tfbs amplification
TF_binding_site_variant 	A sequence variant located within a transcription factor binding site 	SO:0001782 	Regulatory region
regulatory_region_ablation 	A feature ablation whereby the deleted region includes a regulatory region 	SO:0001894 	Regulatory region ablation
regulatory_region_amplification 	A feature amplification of a region containing a regulatory region 	SO:0001891 	Regulatory region amplification
feature_elongation 	A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence 	SO:0001907 	Feature elongation
regulatory_region_variant 	A sequence variant located within a regulatory region 	SO:0001566 	Regulatory region
feature_truncation 	A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence 	SO:0001906 	Feature truncation
intergenic_variant 	A sequence variant located in the intergenic region, between genes 	SO:0001628 	Intergenic"

.my.consequence.order <- function(){
  
  cons.lines <- strsplit(.ens.cons.tab, "\\n")
  
  return(sapply(strsplit(cons.lines[[1]], "\\s+"), "[", 1))
}

.fix.protein.ids <- function(csq.df){
  
  #approach adapted from vcf2maf
  
  proteins <- as.character(csq.df$HGVSp)
  
  # Remove transcript ID from HGVS codon/protein changes, to make it easier on the eye
  proteins <- sub("^.*:", "", proteins, perl=T)
  
  # Remove the prefixed HGVSc code in HGVSp, if found
  proteins <- sapply(regmatches(proteins, regexec("\\(*p\\.\\S+\\)*",proteins)), "[", 1)
  
  proteins <- gsub("[\\(\\)]", "", proteins)
  
  # Create a shorter HGVS protein format using 1-letter codes
  
  for (i in names(.my.short.prots)){
    proteins <- gsub(i, .my.short.prots[i], proteins)
  }
  
  # Fix HGVSp_Short,for splice acceptor/donor variants
  
  splice.pos <- csq.df$Variant_Classification %in% c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")
  
  if (sum(splice.pos) > 0){
    
    c.pos <- as.numeric(sapply(regmatches(as.character(csq.df$HGVSc), regexec("c\\.(\\d+)", as.character(csq.df$HGVSc))), "[", 2))
    
    c.pos <- ifelse(is.na(c.pos) ==F & c.pos < 0, 1, c.pos)
    
    proteins <- ifelse((splice.pos == T) & (is.na(c.pos) ==F) , paste0("p.X", sprintf("%.0f", (c.pos + c.pos %% 3)/3), "_splice"), proteins)
  }
  
  # Fix HGVSp_Short for Silent mutations, so it mentions the amino-acid and position
  
  p.pos <- as.numeric(sapply(regmatches(as.character(csq.df$Protein_position), regexec("^(\\d+)\\/\\d+$", as.character(csq.df$Protein_position))), "[", 2))
  
  fin.prots <- ifelse(csq.df$Variant_Classification == "synonymous_variant", paste0("p.", csq.df$Amino_acids, p.pos,csq.df$Amino_acids) , proteins)
  
  return(fin.prots)
  
}

ug.geno.func <- function(x){
  x[,.(genotype=GT, split.al=lapply(strsplit(AD, ","), as.integer))][,.(genotype, allele_reads=sapply(split.al, "[", 2), total_reads=sapply(split.al, sum))]
}

varscan.gen.func <- function(x){
  x[,.(genotype=GT, allele_reads=AD, total_reads=AD+RD)]
}

csq.from.vcf <- function(vcf.file){
  
  header <- system(paste0("head -n 1000 ",vcf.file," | grep '#' "), intern=T)
  
  test <- fread(vcf.file, sep="\t", skip=length(header)-1, header=T)
  
  test$VcfRow <- as.character(seq_len(nrow(test)))
  
  info.split <- str_split(test$INFO, "[;,]", n = Inf, simplify = FALSE)
  
  test.csq <- lapply(seq_along(info.split), function(x){
    
    use.x <- info.split[[x]]
    
    which.csq <- grep("CSQ=", use.x)
    stopifnot(length(which.csq)==1)
    
    use.x[which.csq] <- sub("CSQ=", "",use.x[which.csq] )
    
    csq.end <- max(grep("\\|", use.x))
    
    #rest.split <- str_split_fixed(use.x[-(which.csq:csq.end)], "=", n=Inf)
    
    list(csq=paste(paste0(x, "|", use.x[which.csq:csq.end]), collapse="\n"))#, info=setNames(rest.split[,2], rest.split[,1]))
    
  })
  
  csq.file <- tempfile()
  
  writeLines(sapply(test.csq, "[[", "csq"), con=csq.file)
  
  csq.dt <- suppressWarnings(fread(csq.file,sep="|", header=F))
  
  csq.info <- gsub("\"", "", header[grep("CSQ", header)])
  
  stopifnot(length(csq.info)==1)
  
  csq.header <- strsplit(regmatches(csq.info, regexec("Format:\\s+(.+)>", csq.info))[[1]][2], "\\|")[[1]]
  csq.header <- append("VcfRow", csq.header)
  
  names(csq.dt) <- csq.header
  
  #pick a most deletarious consequence
  tmp <- strsplit(csq.dt$Consequence, "&")
  names(tmp) <- as.character(seq_along(tmp))
  tmp.dt <- as.data.table(stack(tmp))
  tmp.dt$values <- factor(tmp.dt$values, levels=.my.consequence.order(), ordered=T)
  var.class <- tmp.dt[,.(Variant_Classification=values[which.min(values)]),by=ind]
  var.class <- var.class[order(as.integer(as.character(var.class$ind)),decreasing = F),]
  
  stopifnot(nrow(var.class) == nrow(csq.dt))
  
  csq.dt <- cbind(csq.dt, var.class)
  
  csq.dt[,VcfRow:=as.character(VcfRow)]
  
  csq.m <- merge(test[,c("VcfRow", "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"),with=F], csq.dt, by="VcfRow")
  
  csq.m[,HGVSp_short:=.fix.protein.ids(.SD)]
  
  csq.m
}


#temp <- vcf.to.dt("renamed_vcfs/sample_14-00080_snv.mutect.vep.vcf", samples.to.include=c("14-00080"="14-00080"), geno.to.include=ug.geno.func)

vcf.to.dt <- function(vcf.file, samples.to.include=c("TUMOR"="15-00870"),geno.to.include=varscan.gen.func, output.type = c("full", "readable", "table"), parse.csq=T ){
  
  output.type <- match.arg(output.type)
  
  header <- system(paste0("head -n 1000 ",vcf.file," | grep '#' "), intern=T)
  
  test <- fread(vcf.file, sep="\t", skip=length(header)-1, header=T)
  
  test$VcfRow <- as.character(seq_len(nrow(test)))
  
  #stopifnot(nrow(test[,.N,by=FORMAT])==1)
  
  split.format <- split(test, test$FORMAT)
  
  use.test <- do.call(rbind, lapply(split.format, function(x){
    
    use.x <- x
    
    for(i in names(samples.to.include)){
      geno.names <- strsplit(use.x[["FORMAT"]][1], split=":")[[1]]
      #Some unified genotyper results with only a ./. or an incomplete record
      inp.str <- paste(ifelse(grepl("./.", use.x[[i]], fixed=T), paste("./.", paste(rep(NA, length(geno.names)-1), collapse=":"), sep=":") , use.x[[i]]), collapse="\n")
      #cause fread needs to see a '\n'
      if (grepl("\n", inp.str)==F){
        inp.str <- paste0(inp.str, "\n")
      }
      geno.info <- fread(inp.str, sep=":", header=F)
      names(geno.info) <- geno.names
      tmp.geno.info <- geno.info[,geno.to.include(.SD)]
      names(tmp.geno.info) <- paste(samples.to.include[i],names(tmp.geno.info) ,sep="_")
      
      #make a prettier version of genotype column
      which.geno <- grep("_genotype", names(tmp.geno.info))
      stopifnot(length(which.geno) == 1)
      geno.col <- names(tmp.geno.info)[which.geno]
      
      use.x <- cbind(use.x, tmp.geno.info)
      
      if (output.type == "readable"){
        use.x$old_geno <- use.x[[geno.col]]
        use.x[,c("a1", "a2"):=tstrsplit(old_geno, "/", fixed=T)]
        use.x[,new_geno:=paste(ifelse(a1=="1", ALT, REF), ifelse(a2=="1", ALT, REF), sep="/")]
        use.x[[geno.col]] <- use.x$new_geno
        use.x[,`:=`(old_geno=NULL, new_geno=NULL, a1=NULL, a2=NULL)]
      }else{
        use.x[[geno.col]] <- sapply(str_split(use.x[[geno.col]], "\\/", n=Inf), function(x) sum(x %in% c("0", ".")==F))
      }
    }
    use.x
  }))
  
  keep.test.cols <- setdiff(names(use.test), names(test))
  
  use.test <- use.test[,c("VcfRow", "#CHROM", "POS", "POS","REF", "ALT", "FILTER", keep.test.cols), with=F]
  names(use.test)[1:7] <- c("VcfRow", "Chromosome", "Start_Position", "End_Position", "REF", "ALT", "FILTER")
  use.test[,End_Position:=Start_Position+nchar(REF)-1]
  
  #get consequence info
  
  #id.vars=c("VcfRow" ,"Chromosome" ,"Start_Position" ,"End_Position" ,"REF" ,"ALT", "FILTER")
  
  if (output.type == "table"){
    
    samples <- unique(sub("_allele_reads", "", names(use.test)[grep("allele_reads", names(use.test))]))
    
    samp.list <- list(paste(samples, "genotype", sep="_"), paste(samples, "allele_reads", sep="_"), paste(samples, "total_reads", sep="_"))
    
    melt.tab <- data.table::melt( measure.vars=samp.list,  data=use.test, value.name=c("genotype", "allele_reads", "total_reads"), variable.factor=FALSE)
    
    var.tab <- data.table(variable=as.character(seq_along(samples)), sample=samples)
    
    melt.tab <- merge(melt.tab, var.tab, by="variable")
    
    melt.tab[,variable:=NULL]
  }else{
    melt.tab <- use.test
  }
  
  if (parse.csq){
    info.split <- str_split(test$INFO, "[;,]", n = Inf, simplify = FALSE)
    
    test.csq <- lapply(seq_along(info.split), function(x){
      
      use.x <- info.split[[x]]
      
      which.csq <- grep("CSQ=", use.x)
      stopifnot(length(which.csq)==1)
      
      use.x[which.csq] <- sub("CSQ=", "",use.x[which.csq] )
      
      csq.end <- max(grep("\\|", use.x))
      
      #rest.split <- str_split_fixed(use.x[-(which.csq:csq.end)], "=", n=Inf)
      
      list(csq=paste(paste0(x, "|", use.x[which.csq:csq.end]), collapse="\n"))#, info=setNames(rest.split[,2], rest.split[,1]))
      
    })
    
    csq.dt <- fread(paste(sapply(test.csq, "[[", "csq"), collapse="\n"), sep="|", header=F)
    
    csq.info <- gsub("\"", "", header[grep("CSQ", header)])
    
    stopifnot(length(csq.info)==1)
    
    csq.header <- strsplit(regmatches(csq.info, regexec("Format:\\s+(.+)>", csq.info))[[1]][2], "\\|")[[1]]
    csq.header <- append("VcfRow", csq.header)
    
    names(csq.dt) <- csq.header
    
    #get the picked consequences
    
    if (all(is.na(csq.dt$PICK))){
      message("No PICK'd consequences, keeping all")
      picked.csq <- csq.dt
    }else{
      stopifnot(all(grepl(",", test$ALT))==F && all(csq.dt$ALLELE_NUM==1))
      picked.csq <- csq.dt[is.na(PICK)==F]
    }
    
    #pick a most deletarious consequence
    tmp <- strsplit(picked.csq$Consequence, "&")
    names(tmp) <- as.character(seq_along(tmp))
    tmp.dt <- as.data.table(stack(tmp))
    tmp.dt$values <- factor(tmp.dt$values, levels=.my.consequence.order(), ordered=T)
    var.class <- tmp.dt[,.(Variant_Classification=values[which.min(values)]),by=ind]
    
    stopifnot(nrow(var.class) == nrow(picked.csq))
    
    if (output.type == "readable"){
      use.picked.csq <- picked.csq[,.(VcfRow=as.character(VcfRow), Existing_variation=gsub("&", ";", Existing_variation), in_cosmic=grepl("COSM", Existing_variation), in_dbsnp=grepl("rs\\d+", Existing_variation),
                                      Ensembl_gene=Gene, Ensembl_transcript=Feature, Protein_position, Amino_acids, HGVSp, HGVSc, Hugo_symbol=SYMBOL, SIFT, PolyPhen,  ExAC_AF=ExAC_AF)]
    }else{
      use.picked.csq <- picked.csq
      use.picked.csq[,VcfRow:=as.character(VcfRow)]
    }
    
    csq.m <- merge(use.picked.csq, var.class[,.(VcfRow=as.character(ind), Variant_Classification)], by="VcfRow")
    
    csq.m[,HGVSp_short:=.fix.protein.ids(.SD)]
    
    #merge them together
    
    ret.dt <- merge(melt.tab, csq.m, by="VcfRow", all=F)
    
    if (all(is.na(csq.dt$PICK))==F){
      stopifnot(nrow(ret.dt) == nrow(melt.tab))
    }
  }else{
    ret.dt <- melt.tab
  }
  
  all.cols <- names(ret.dt)
  
  if (output.type != "table"){
    
    geno.cols <- all.cols[apply(sapply(samples.to.include, function(x) grepl(x, all.cols, fixed=T)), 1, any)]
    
    ret.cols <- c(setdiff(all.cols, c(geno.cols, "VcfRow")), geno.cols)
    
  }else{
    
    ret.cols <- setdiff(all.cols, c("VcfRow", "Allele"))
    
  }
  
  ret.dt[,ret.cols, with=F]
}

