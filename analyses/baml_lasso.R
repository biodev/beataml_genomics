library(data.table)
library(glmnet)

##Note also need openxlsx, reshape2 which are called via ::

preprocess <- function(){

  ##sample clinical info
  
  clin.dt <- fread("All_sample_data.csv", sep=",", header=T, skip=3)
  
  ##process drugs (from paper supplement/sage)
  
  final.drugs <- data.table(openxlsx::read.xlsx("TableS10_drug_response.xlsx"))
  
  split.drugs <- split(final.drugs, by="inhibitor")
  
  ##process variants (retrieved from vizome)
  
  var.dt <- fread("All_variants.csv", sep=",", header=T, skip=2)
  
  #complete set of samples that have wes
  wes.clin <- clin.dt[Included_2018_DNAseqAnalysis == "y"]
  
  ##add in the fusion data as well
  
  wes.clin[,new_value:=make.names(WHO_Fusion)]
  
  ##add in fusion categories as symbol
  cat.dt <- var.dt[NA][seq_len(wes.clin[,.N])]
  cat.dt$sample_id <- wes.clin$SampleID
  cat.dt$gene <- wes.clin$new_value
  
  cat.dt <- cat.dt[gene %in% c("None", "Unknown")==F]
  
  var.dt <- rbind(var.dt, cat.dt)
  
  ##end variant
  
  ##expression,focusing on the WGCNA modules here.
  
  ##remove the grey module
  use.mes <- as.matrix(read.delim("wgcna_modules.txt", sep="\t", header=T))
  
  rna.clin <- clin.dt[Included_2018_RNAseqAnalysis == "y"]
  
  ##end expression

  save(split.drugs, rna.clin, wes.clin, use.mes, var.dt, clin.dt, file="preproc_data.RData")
  
}

.get.boot.list <- function(use.mat, auc.vec, drug, stand){
  
  #only doing 10 75/25 training/test splits here for now but used 100 for the paper
  lapply(1:10, function(y){
    train.samp <- sample.int(length(auc.vec), size=round(length(auc.vec)*.75), replace=F)
    test.samp <- setdiff(seq_along(auc.vec), train.samp)
    
    #again, only doing 100 bootstraps here for demonstration
    boot.vals <- lapply(1:100, function(y){
      
      boot.samps <- sample.int(n=length(train.samp), size=length(train.samp), replace=T)
      
      #divide the bootstrapped samples into 10-folds as in page 13 here https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet_beta.pdf
      foldid <- sample(1:10,size=length(boot.samps),replace=TRUE)
      
      list(iter=y, boot=boot.samps,foldid=foldid)
    })
    
    list(rep=y, mat=use.mat, auc=auc.vec, train=train.samp, test=test.samp, boots=boot.vals, drug=drug, stand=stand)
  })
}

.get.earliest.samples <- function(pat.info, ...){
  
  samp.lists <- list(...)
  
  use.samps <- copy(pat.info)
  
  for(i in seq_along(samp.lists)){
    use.samps <- use.samps[SampleID %in% samp.lists[[i]]]
  }
  
  if (nrow(use.samps) > 0){
    
    use.samps[,collection_diff:=suppressWarnings(as.numeric(TimeOfCollectionRelativeToInclusion))]
    
    use.samps <- use.samps[is.na(collection_diff)==F]
    
    min.samps <- use.samps[,.(collection_diff=min(collection_diff), all_dates=paste(collection_diff, collapse=",")),by=PatientID]
    samps.m <- merge(use.samps, min.samps[,.(PatientID, collection_diff)], by=c("PatientID", "collection_diff"))
    
    #to get rid of duplicate patients (which are usual PB and BM), sort by tissue type such that BM is first and then !duplicated
    
    fin.samps <- samps.m[order(SpecimenType)][!duplicated(PatientID)]
    
    stopifnot(nrow(min.samps) == nrow(fin.samps))
  
  }else{
    
    fin.samps <- copy(use.samps)  
  
  }
  
  fin.samps
}

#Since the main random part of the code has been carried out sequentially this can also be run in parallel in any desired fashion
.run.lasso <- function(auc.vec, tmp.var.mat, train.samp, test.samp, boot.vals, stand, drug){
  
  stopifnot(length(auc.vec) == nrow(tmp.var.mat))
  
  train.mat <- tmp.var.mat[train.samp,]
  train.auc <- auc.vec[rownames(tmp.var.mat)][train.samp]
  
  #train models on the bootrap samples using 10-fold cv (The 10 folds provided)
  boot.res <- lapply(boot.vals, function(y){
    
    test.cv <- cv.glmnet(x=train.mat[y$boot,],y=as.numeric(train.auc[y$boot]),alpha=1, foldid=y$foldid, family="gaussian", standardize=stand, intercept=T)
    
    list(boot=y$iter, model=test.cv)
    
  })
  
  #predict using the test data
  test.res <- do.call(rbind, lapply(boot.res, function(y){
    
    preds.test <- predict(y$model, newx=tmp.var.mat[test.samp,], s="lambda.1se")
    
    data.table(boot=y$boot, test_samp=test.samp, pred_auc=preds.test[,1])
    
  }))
  
  #average the predictions
  test.sum <- test.res[,.(pred_auc=mean(pred_auc)),by=test_samp]
  
  tmp.dt <- merge(data.table(test_samp=test.samp, test_auc=auc.vec[rownames(tmp.var.mat)][test.samp]), test.sum, by="test_samp")
  
  tmp.dt$inhibitor <- drug
  
  #also capture the resulting coefficients
  coef.mat <- do.call(cbind, lapply(boot.res, function(y){
    as.matrix(coef(y$model, s="lambda.1se"))
  } ))
  
  list(preds=tmp.dt, coefs=coef.mat, models=boot.res, drug=drug)
}



run.combined.datatypes <- function(){
  
  load("preproc_data.RData")
  
  run.comb <- sapply(split.drugs, function(x){
    
    cur.drug <- copy(x)
    
    message(cur.drug$inhibitor[1])
    
    use.samps <- .get.earliest.samples(clin.dt, rna.clin$SampleID , wes.clin$SampleID, cur.drug$lab_id)
    
    if (use.samps[,.N] >= 200){
      
      use.exprs <- use.mes[use.samps$SampleID,]
      fin.drugs <- cur.drug[lab_id %in% use.samps$SampleID]
      
      use.vars <- var.dt[sample_id %in% use.samps$SampleID]
      use.vars <- use.vars[,samp_count:=length(unique(sample_id)), by=gene][samp_count >= 10]
      
      use.vars[,samp_fac:=factor(sample_id, levels=use.samps$SampleID)]
      
      use.vars[,value:=1]
      
      tmp.var.mat <- reshape2::acast(samp_fac~gene,data=use.vars, value.var="value", fun.aggregate = function(x) as.integer(length(x) > 0), drop=F)
      
      stopifnot(all(rownames(tmp.var.mat) %in% rownames(tmp.var.mat)) && nrow(tmp.var.mat) == nrow(use.exprs))
      
      tmp.var.mat <- cbind(tmp.var.mat, use.exprs[rownames(tmp.var.mat),])
      
      auc.vec <- setNames(fin.drugs$auc, fin.drugs$lab_id)[rownames(tmp.var.mat)]
      
      #note the stand=T here means that we are having glmnet standardize the Xs here
      drug.run <- .get.boot.list(tmp.var.mat, auc.vec, cur.drug$inhibitor[1], T)
      
      drug.run
      
    }else{
      
      NULL
    
    }
    
  })
  
  #only 5 drugs for now
  drug.res.list <- lapply(run.comb[1:5], function(run.list){
    
    res.list <- lapply(run.list, function(x){
      
      if (is.null(x)){
        NULL
      }else{
        tmp.list <- .run.lasso(x$auc, x$mat, x$train, x$test, x$boots, x$stand, x$drug)
        
        tmp.list$index <- x$rep
        
        tmp.list$preds$index <- x$rep
        
        tmp.list
      }
     
    })
    
  })
  
  #summarize
  
  res.list <- lapply(drug.res.list, function(tmp.list){
    
    pred.dt <- do.call(rbind, lapply(tmp.list, function(y){
      tmp.dt <- y$preds
      tmp.dt$inhibitor <- y$drug
      tmp.dt
    }))
    
    coef.mat <- do.call(cbind, lapply(tmp.list, function(y){
      tmp.mat <- y$coef
      colnames(tmp.mat) <- paste(y$preds$index[1], seq_len(ncol(tmp.mat)), sep="_")
      tmp.mat
    }))
    
    list(preds=pred.dt, coefs=coef.mat, drug=pred.dt$inhibitor[1], datatype="combined")
    
  })
  
  #coerce into prediction results
  
  pred.sum <- do.call(rbind, lapply(res.list, "[[", "preds"))
  
  #pull out list of coefficient matrices
  
  coef.list <- lapply(res.list, "[[" ,"coefs")
  
  save(pred.sum, coef.list, file="lasso_results.Rdata")
}

