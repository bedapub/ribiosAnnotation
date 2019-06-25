summarizeProbesets <- function(eset,
                               index.name,
                               fun=mean,
                               keep.nonindex=FALSE,
                               keep.featureNames=FALSE, ...) {
  
  if(missing(index.name) || !index.name %in% colnames(fData(eset))) {
    stop("'index.name' must be a valid column name in fData(", as.character(match.call()$eset), ")")
  }
  fun <- match.fun(fun)
  
  probe.index <- as.character(fData(eset)[,index.name])
  probe.has.index <- !is.na(probe.index) & probe.index != "" 
  eset.indexed <- eset[probe.has.index,]
  
  probe.indexed.fac <- factor(fData(eset.indexed)[,index.name])
  eset.fun <- summarizeRows(exprs(eset.indexed),
                            probe.indexed.fac,
                            fun=fun, ...)

  fun.fd.match <- match(levels(probe.indexed.fac), probe.index)
  eset.remain <- eset[fun.fd.match,]
  
  if(keep.featureNames) {
    newFeatureNames <- rownames(fData(eset.remain))
    ## featureNames(eset.remain) <- rownames(fData(eset.remain))
  } else {
    newFeatureNames <- rownames(eset.fun)
    ## rownames(fData(eset.remain)) <- featureNames(eset.remain)
  }
  rownames(eset.fun) <- featureNames(eset.remain) <- newFeatureNames
  exprs(eset.remain) <- eset.fun
  
  if(keep.nonindex) {
    eset.inval <- eset[!probe.has.index,]
    resExprs <- rbind(exprs(eset.remain),
                      exprs(eset.inval))
    resFdata <- rbind(fData(eset.remain),
                      fData(eset.inval))
    eset.remain <- new("ExpressionSet",
                       exprs=resExprs,
                       featureData=new("AnnotatedDataFrame", resFdata),
                       phenoData=phenoData(eset))
  }
  return(eset.remain)
}
