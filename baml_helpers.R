
#tryCatch(library(openxlsx), error=function(x) require("xlsx"))
library(reshape2)
library(xml2)
library(htmltools)

get.timestamp <- function(){
  gsub("-", "_", strsplit(as.character(Sys.time()), "\\s+")[[1]][1])
}

args <- commandArgs(TRUE)

get.ensembl.annotation <- function(){
  
  con <- dbConnect(MySQL(), groups="ensembl_78_37")
  
  #it was issuing warnings like: Unsigned INTEGER in col 0 imported as numeric
  suppressWarnings(genes <- dbGetQuery(con, "SELECT stable_id, display_label, gene.description as description, biotype FROM gene JOIN xref ON display_xref_id = xref_id"))
  
  dbDisconnect(con)
  
  return(genes)
   
}

format.thfp.results <- function(...){
  all.files <- list(...)
  files <- all.files[-1]
  out.file <- all.files[[1]]
  
  stopifnot(require("RMySQL"))
  
  top.fusion.dta <- do.call(rbind, lapply(files, function(x){
    temp.res <- read.delim(x, sep="\t", header=F, stringsAsFactors = F)
    names(temp.res) <- c("SampleName", "LeftGene","LeftChrID","LeftCoord","RightGene","RightChrID","RightCoord","Num.ofReadsSpanningFusion","Num.ofSpanningMatePairs","Num.ofMatePairsOneEndOnlySpanning", "Score")
    
    temp.res$sample.key <- with(temp.res, paste(SampleName, paste0(LeftChrID, "-", RightChrID), LeftCoord, RightCoord,  sep="."))
    
    #get the flanking sequence info from the 'potential fusions'
    pot.fusions <- readLines(file.path(dirname(x), "potential_fusion.txt"))
    num.samps <- sum(grepl("Sample_", pot.fusions))
    
    fusion.dta <- do.call(rbind, lapply(split(pot.fusions, rep(seq_len(num.samps), each=6)), function(y){
      split.y <- strsplit(y, "\\s+")
      data.frame(sample.key=paste(split.y[[1]][1:4], collapse="."), StrandInfo=split.y[[1]][5], LeftFlankingSequence=split.y[[2]][1], RightFlankingSequence=split.y[[3]][2], stringsAsFactors = F)
    }))
    
    temp.fusions <- merge(temp.res, fusion.dta, by="sample.key", all=F, incomparables=NA, sort=F)
    
    stopifnot(nrow(temp.fusions) == nrow(temp.res))
    
    temp.fusions
  }))
  
  #now annotate with biotype
  
  con <- dbConnect(MySQL(), groups="ensembl_78_37")
  
  #it was issuing warnings like: Unsigned INTEGER in col 0 imported as numeric
  suppressWarnings(genes <- dbGetQuery(con, "SELECT * FROM gene JOIN xref ON display_xref_id = xref_id"))
  
  dbDisconnect(con)
  
  #first convert the ensembl ids to something a little more useful
  #then merge in applicable biotype based on the new display_label values...if available
  
  initial.size <- nrow(top.fusion.dta)
  
  biotypes <- aggregate(biotype~display_label, function(x) paste(unique(x), collapse=";"), data=genes)
  
  for(i in c("Left", "Right")){
    new.col <- paste0(tolower(i),"_display_label")
    top.fusion.dta <- merge(top.fusion.dta, genes[,c("display_label", "stable_id")], by.x=paste0(i,"Gene"), by.y="stable_id", all.x=T, all.y=F, incomparables=NA)
    names(top.fusion.dta)[ncol(top.fusion.dta)] <- new.col
    top.fusion.dta[[new.col]] <- ifelse(is.na(top.fusion.dta[[new.col]]), top.fusion.dta[[paste0(i,"Gene")]] ,top.fusion.dta[[new.col]])
    
    top.fusion.dta <- merge(top.fusion.dta, biotypes, by.x=new.col, by="display_label", all.x=T, all.y=F, incomparables=NA)
    names(top.fusion.dta)[ncol(top.fusion.dta)] <- paste0(i,"GeneBiotype")
  }
  
  stopifnot(nrow(top.fusion.dta) == initial.size)
  
  top.fusion.dta$SampleName <- sub("Sample_", "", top.fusion.dta$SampleName)
  
  fin.top.fusions <- dplyr::select(top.fusion.dta, SampleName, LeftGene=left_display_label, LeftChrID, LeftCoord, RightGene=right_display_label, RightChrID, RightCoord, 
                                   Num.ofReadsSpanningFusion, Num.ofSpanningMatePairs, Num.ofMatePairsOneEndOnlySpanning, StrandInfo, LeftFlankingSequence, 
                                   RightFlankingSequence, LeftGeneBiotype, RightGeneBiotype)
  
  write.table(fin.top.fusions, sep="\t", file=out.file, col.names=T, row.names=F, quote=T)
  
}

#these should be the combined_featureCounts.txt files
annotate.gene.counts <- function(count.matrices, out.file.base="test", include.log=F, samples=NULL){
  
  message("Getting genes")
  
  stopifnot(require("RMySQL"))
  
  con <- dbConnect(MySQL(), groups="ensembl_78_37")
  
  #it was issuing warnings like: Unsigned INTEGER in col 0 imported as numeric
  suppressWarnings(genes <- dbGetQuery(con, "SELECT * FROM gene JOIN xref ON display_xref_id = xref_id"))
  
  dbDisconnect(con)
  
  mat.list <- lapply(count.matrices, function(x){
    read.delim(x, sep="\t", comment.char="#",check.names=F, stringsAsFactors = F)
  })
  
  cur.mat <- Reduce(function(x,y){
    
    merge(x,y, by=c("Geneid", "Chr", "Start", "End", "Strand", "Length"))
    
  }, mat.list)
  
  stopifnot(length(unique(sapply(mat.list, nrow))) == 1 && nrow(cur.mat) == nrow(mat.list[[1]]))
  
  chrs <- sapply(strsplit(cur.mat$Chr, ";"), unique)
  stopifnot(is.character(chrs))
  cur.mat$Chr <- chrs
  
  strands <- sapply(strsplit(cur.mat$Strand, ";"), unique)
  stopifnot(is.character(strands))
  cur.mat$Strand <- strands
  
  starts <- sapply(strsplit(cur.mat$Start, ";"), function(x) min(as.numeric(x)))
  names(cur.mat)[which(names(cur.mat)=="Start")] <- "Exon_Start"
  cur.mat$Start <- starts
  
  ends <- sapply(strsplit(cur.mat$End, ";"), function(x) max(as.numeric(x)))
  names(cur.mat)[which(names(cur.mat)=="End")] <- "Exon_End"
  cur.mat$End <- ends
  
  cur.mat.merge <- merge(cur.mat, genes[,c("stable_id", "display_label")], by.x="Geneid", by.y="stable_id", all.x=T, all.y=F, incomparables=NA)
  
  stopifnot(nrow(cur.mat.merge) == nrow(cur.mat) && sum(is.na(cur.mat.merge)) == 0)
  
  fin.mat <- dplyr::select(cur.mat.merge, Gene=Geneid, Symbol=display_label, Chr, Exon_Start, Exon_End, Strand, Length, GeneStart=Start, GeneEnd=End, dplyr::everything())
  
  if ((missing(samples) || is.null(samples) || all(is.na(samples)))==F){
    
    all.samples <- names(fin.mat)[grepl("\\d{2}-\\d{5}", names(fin.mat))]
    
    fin.mat <- fin.mat[,-which(names(fin.mat) %in% setdiff(all.samples, samples))]
  }
  
  #raw counts
  
  write.table(fin.mat, file=paste0(out.file.base, "_rawcounts_",get.timestamp(),".csv"), sep=",", row.names=F, col.names=T,quote=T)
  
  #cpm-counts
  
  count.portion <- dplyr::select(fin.mat, dplyr::matches("\\d{2}-\\d{5}"))
  
  cpm.portion <- edgeR::cpm(count.portion, log=FALSE)
  
  cpm.all <- cbind(dplyr::select(fin.mat, -dplyr::matches("\\d{2}-\\d{5}")), cpm.portion)
  
  write.table(cpm.all, file=paste0(out.file.base, "_cpm_",get.timestamp(),".csv"), sep=",", row.names=F, col.names=T,quote=T)
  
  if (include.log){
    
    l2.cpm.portion <- edgeR::cpm(count.portion, prior.count=2,log=TRUE)
    
    l2.cpm.all <- cbind(dplyr::select(fin.mat, -dplyr::matches("\\d{2}-\\d{5}")), l2.cpm.portion)
    
    write.table(l2.cpm.all, file=paste0(out.file.base, "_log2_cpm_",get.timestamp(),".csv"), sep=",", row.names=F, col.names=T,quote=T)
    
  }
}



setup.sgs <- function(directory=".", pair.file){
  
  if (missing(pair.file) || is.null(pair.file) || all(is.na(pair.file))){
    stop("ERROR: Please specify a valid pairing CSV file")
  }
  
  pairs <- read.csv(pair.file, stringsAsFactors = F)
  
  melt.pairs <- melt(pairs, measure.vars=c("AML", "Normal"), variable.name="type", value.name="samples")
  melt.pairs$type <- as.character(melt.pairs$type)
  
  all.files <- list.files(directory, full.names=T)
  
  cgs <- all.files[grep("CaptureGroup", all.files)]
  
  #make sure to only keep directories
  
  cgs <- cgs[dir.exists(cgs)]
  
  if (length(cgs) > 0){
    
    cg.bams <- unlist(lapply(file.path(cgs, "sample_alignments"), list.files, pattern="Sample_\\d+.*.dedup.ba[mi]$", full.names=T))
    
    if (length(cg.bams) == 0){
      stop("ERROR: No bam files found in the specified CaptureGroups, has the preprocessing been run?")
    }
    
    cg.bams.dta <- data.frame(bams=cg.bams, samples=sapply(strsplit(basename(cg.bams), "[_\\.]"), "[", 2), stringsAsFactors = F)
    
    found.cg.bams <- merge(cg.bams.dta, melt.pairs, by="samples", all=F)
    
    patient.sample.count <- aggregate(samples~patient_id, length, data=found.cg.bams)
    
    exp.count <- aggregate(type~patient_id, c, data=melt.pairs[complete.cases(melt.pairs),], simplify=F)
    exp.count$samples <- lengths(exp.count$type)*2
    
    exp.count$classification <- NA
    exp.count$classification[sapply(exp.count$type, function(x) all(c("Normal", "AML") %in% x))] <- "SomaticGenotyping"
    exp.count$classification[sapply(exp.count$type, function(x) all(x == "AML"))] <- "TumorOnly"
    exp.count$classification[sapply(exp.count$type, function(x) all(x == "Normal"))] <- "NormalOnly"
    
    if (any(is.na(exp.count$classification))){
      stop("ERROR: can't assign a classification to sample/patient relationship")
    }
    
    stopifnot(all(patient.sample.count$samples %in% c(4,2)))
    
    #this will subset to only those complete patients that have sequencing though technically a single sample from a tumor-only could be genotyped
    patient.sample.count <- merge(patient.sample.count, exp.count, by=c("patient_id", "samples"), all=F)
    
    bams.for.genotyping <- merge(found.cg.bams, patient.sample.count[,c("patient_id", "classification")], by="patient_id", all=F)
    
    split.geno.bams <- split(bams.for.genotyping, bams.for.genotyping$classification)
    
    for(i in names(split.geno.bams)){
      
      #again, make sure to only keep directories
      
      sgs <- all.files[grep(i, all.files)]
      
      sgs <- sgs[dir.exists(sgs)]
      
      use.sub.dir <- paste0(ifelse(i == "SomaticGenotyping", "paired", "single"),"_alignments")
      
      if (length(sgs) == 0){
        
        #put everyone in the same SomaticGenotyping group by creating symlinks back to the capturegroups
        
        cur.sg <- 1
        keep.bams <- split.geno.bams[[i]]
        
      }else{
        
        #determine who has already been assigned to a SomaticGenotyping group and place the rest in another group
        
        existing.files <- unlist(lapply(file.path(sgs, use.sub.dir), list.files, pattern="Sample_\\d+.*.dedup.ba[mi]$"))
        
        if (length(existing.files) > 0){
          keep.bams <- split.geno.bams[[i]][basename(split.geno.bams[[i]]$bams) %in% basename(existing.files) == F,]
        }
        
        cur.sg <- max(as.numeric(sub(i, "", basename(sgs)))) + 1
        
      }
      
      if (nrow(keep.bams) > 0){
        message("creating ", paste0(i, cur.sg))
        dir.create(file.path(paste0(i, cur.sg), use.sub.dir), recursive=T)
        
        invisible(all(file.symlink(from=normalizePath(keep.bams$bams), to=file.path(paste0(i, cur.sg), paste0(use.sub.dir, "/")))))
      }
      
    }
    
  }else{
    stop("ERROR: No CaptureGroup directories are available")
  }
  
}


setup.capturegroups <- function(directory=".", pair.file, expected.lanes=c(AML=5, Normal=3)){
  
  if (missing(pair.file) || is.null(pair.file) || all(is.na(pair.file))){
    stop("ERROR: Please specify a valid pairing CSV file")
  }
  
  pairs <- read.csv(pair.file, stringsAsFactors = F)
  
  melt.pairs <- melt(pairs, measure.vars=c("AML", "Normal"), variable.name="type", value.name="samples")
  melt.pairs$type <- as.character(melt.pairs$type)
  
  all.files <- list.files(directory, full.names=T)
  
  fcs <- all.files[grep("FlowCell", all.files)]
  
  #make sure to only keep directories
  
  fcs <- fcs[dir.exists(fcs)]
  
  if (length(fcs) > 0){
    
    fc.bams <- unlist(lapply(file.path(fcs, "alignments"), list.files, pattern=".*ba[mi]$", full.names=T))
    
    if (length(fc.bams) == 0){
      stop("ERROR: No bam files found in the specified FlowCells, has the summarization/alignment procedure been run?")
    }
    
    fc.bams.dta <- data.frame(bams=fc.bams, samples=sapply(strsplit(basename(fc.bams), "[_\\.]"), "[", 3), stringsAsFactors = F)
    
    found.fc.bams <- merge(fc.bams.dta, melt.pairs, by="samples", all=F)
    
    #all the bams should be accounted for in the pairs file
    
    stopifnot(nrow(found.fc.bams) == nrow(fc.bams.dta))
    
    sample.lane.count <- aggregate(bams~samples+type, length, data=found.fc.bams[grepl("\\.bam$", found.fc.bams$bams),])
    
    sample.lane.count <- cbind(sample.lane.count, exp.lanes=expected.lanes[sample.lane.count$type], stringsAsFactors=F)
    
    #There shouldn't be anymore lanes than expected
    stopifnot(all(sample.lane.count$bams <= sample.lane.count$exp.lanes))
    
    keep.samples <- sample.lane.count[sample.lane.count$bams == sample.lane.count$exp.lanes,]
    
    if (nrow(keep.samples) == 0){
      stop("ERROR: No samples have reached the desired lane counts, no CaptureGroups will be created")
    }
    
    bams.for.preproc <- found.fc.bams[found.fc.bams$samples %in% keep.samples$samples,]
    
    #now figure out if capturegroups already exist and then add these to new cg directories
    
    cgs <- all.files[grep("CaptureGroup", all.files)]
    
    cgs <- cgs[dir.exists(cgs)]
    
    if (length(cgs) == 0){
      
      #Determine the starting point for the capturegroups and place them in there sequentially
      
      cur.cg <- 0
      
    }else{
      
      #determine who has already been assigned to a CaptureGroup and place the rest in new groups
      
      existing.files <- unlist(lapply(file.path(cgs, "alignments"), list.files, pattern=".*ba[mi]$"))
      
      if (length(existing.files) > 0){
        bams.for.preproc <- bams.for.preproc[basename(bams.for.preproc$bams) %in% basename(existing.files) == F,]
        keep.samples <- keep.samples[keep.samples$samples %in% bams.for.preproc$samples,]
      }
      
      if (nrow(bams.for.preproc) == 0){
        stop("ERROR: No BAMs left to preprocess, nothing will be done")
      }
      
      cur.cg <- max(as.numeric(sub("CaptureGroup", "", basename(cgs))))
      
    }
    
    num.cgs <- ceiling(nrow(keep.samples) / 12)
    
    #assume that ordering roughly follows cg assignment
    
    keep.samples <- keep.samples[order(keep.samples$samples),]
    
    keep.samples$cg <- rep(seq_len(num.cgs), each=12)[seq_len(nrow(keep.samples))]
    
    bam.assignments <- merge(bams.for.preproc, keep.samples[,c("samples", "type", "cg")], by=c("samples", "type"), all=F)
    
    split.bams <- split(bam.assignments, bam.assignments$cg)
    
    for(i in names(split.bams)){
      
      use.cg <- cur.cg + as.numeric(i)
      
      dir.create(file.path(paste0("CaptureGroup", use.cg), "alignments"), recursive=T)
      
      stopifnot(all(file.symlink(from=normalizePath(split.bams[[i]]$bams), to=file.path(paste0("CaptureGroup", use.cg), "alignments/"))))
    }
    
  }else{
    stop("ERROR: No FlowCell directories are available")
  }

}


make.fastqc.report <- function(fastqc.dir, report.name="test.html", num.per.line=4, img.attrs=c(width=500, height=350), sample.name.func=function(x) x[2:5]){
  
  if (missing(fastqc.dir) || is.null(fastqc.dir) || all(is.na(fastqc.dir))){
    stop("ERROR: Please provide a path to sample-level FastQC results")
  }
  
  all.files <- list.files(fastqc.dir, pattern="[RD]NA.+html$", full.names=T)
  
  if (length(all.files) == 0){
    stop("ERROR: Specified directory does not have FastQC html files")
  }
  
  status.levs <- c("OK", "WARN", "FAIL")
  
  all.images <- lapply(all.files, function(x){
    
    temp.html <- read_html(x)
    modules <- xml_find_all(temp.html, ".//div[@class='module']")
    
    module.names <- xml_text(xml_find_all(modules, ".//h2"))
    
    names(modules) <- module.names
    
    module.status <- gsub("\\[|\\]", "", xml_attr(xml_find_all(modules, ".//h2/img"), "alt"))
    
    names(module.status) <- module.names
    
    module.images <- xml_find_all(modules, ".//p/img")
    
    for(i in names(img.attrs)){
      xml_attr(module.images, i) <- as.character(img.attrs[i])
    }
    
    return(list(module.images=module.images, module.status=factor(module.status[names(module.images)], levels=status.levs, ordered=T)))
    
  })
  
  names(all.images) <- basename(all.files)
  
  samp.list <- strsplit(basename(all.files), "_")
  samp.dta <- setNames(data.frame(t(sapply(samp.list, sample.name.func)), stringsAsFactors = F), c("lab_id", "core_sample_id", "lane", "pair"))
  samp.dta$files <- basename(all.files)
  
  samp.dta.ord <- samp.dta[do.call(order, samp.dta[,c("lab_id", "lane", "pair")]),]
  
  num.lines <- ceiling(nrow(samp.dta.ord)/num.per.line)
  
  samp.dta.ord$lines <- rep(seq_len(num.lines), each=num.per.line, length.out=nrow(samp.dta.ord))
  
  split.samps <- split(samp.dta.ord, samp.dta.ord$lines)
  
  image.names <- unique(unlist(lapply(all.images, function(x) names(x$module.images))))
  
  cur.sample.vals <- character()
  
  all.dts <- lapply(image.names, function(x){
    
    all.samps <- lapply(split.samps, function(y){
      
      all.lines <- lapply(seq_len(nrow(y)), function(z){
        
        stopifnot(y$files[z] %in% names(all.images))
        stopifnot("module.images" %in% names(all.images[[basename(y$files[z])]]))
        
        if (x %in% names(all.images[[basename(y$files[z])]]$module.images)){
          
          sample.dt.status <- as.character(all.images[[y$files[z]]]$module.status[x])
          html.el <- HTML(as.character(all.images[[basename(y$files[z])]]$module.images[x]))
        }
        else{
          #if an image isn't present, assume that it is ok and default to a blank image
          sample.dt.status <- "OK"
          html.el <- img(class="indented", src="//:0", alt=x, width="500", height="350")
        }
        
        sample.status <- paste(y$lab_id[z], gsub(" ", "_", x), sample.dt.status , sep="_")
          
        if (sample.status %in% cur.sample.vals == F){
          cur.sample.vals <<- append(cur.sample.vals, sample.status)
          
          tl <- tagList(a(name=sample.status), p(paste(y$lab_id[z], y$lane[z], y$pair[z], sample.dt.status), style="text-align: center"),
                        html.el)
        }else{
          tl <- tagList(p(paste(y$lab_id[z], y$lane[z], y$pair[z], sample.dt.status), style="text-align: center"),
                        html.el)
        }
        
        tags$td(tl)
        
      })
      
      tags$tr(do.call(tagList, all.lines))
      
    })
    
    tags$table(
      
      tags$thead(tags$tr(tags$a(name=gsub(" ", "_", x)), tags$th(h2(x), colspan=num.per.line))),
      tags$tbody(do.call(tagList, all.samps))
    )
    
  })
  
  #summarize the QC statuses for all the lanes/read pairs of a sample
  
  split.files <- split(samp.dta.ord$files, samp.dta.ord$lab_id)
  
  sample.qc.summary <- lapply(names(split.files), function(x){
    
    temp.tds <- lapply(image.names, function(y){
      
      status.vals <- sapply(all.images[split.files[[x]]], function(z){
        
        if (y %in% names(z$module.status)){
          z$module.status[y]
        }else{
          factor("OK", levels=c("OK", "WARN", "FAIL"), ordered=T)
        }
        
      })
      
      tags$td(a(href=paste0("#", x, "_", gsub(" ", "_", y), "_", 
                            as.character(status.vals[which.max(status.vals)])), as.character(status.vals[which.max(status.vals)])), style="text-align: center")
      
    })
    
    temp.tds <- append(list(tags$td(x)), temp.tds)
    tags$tr(do.call(tagList, temp.tds))
  })
  
  tab.head <- lapply(c("Sample", image.names), function(x){
    if (x != "Sample"){
      use.x <- tags$a(x, href=paste0("#", gsub(" ", "_", x)))
    }else{
      use.x <- x
    }
    
    tags$th(use.x, style="text-align: center")
  })
  
  #end sample summary
  
  #make the final html page
  
  page <- tags$html(
    
    tags$body(
      tags$table(tags$thead(
        tags$tr(tags$th(h2("Sample QC Summary"), colspan=length(tab.head))),
        tags$tr(do.call(tagList, tab.head))
      ),
      tags$tbody(do.call(tagList, sample.qc.summary)), style="border: 1px solid #000000;", border=1),
      do.call(tagList, all.dts)
    )
  )
  
  save_html(page, file=report.name)
}

if (length(args) > 0){
  
  if (exists(args[1])){
    do.call(args[1], as.list(args[-1]))
  }
  
}
