## transformations
matrixToLongTable <- function(x, valueLabel="value", rowLabel="row", colLabel="col") {
  val <- as.vector(x)
  rowLabels <- rep(rownames(x), ncol(x))
  colLabels <- rep(colnames(x), each=nrow(x))
  res <- data.frame(rowLabl=rowLabels,
                    colLabels=colLabels,
                    value=val)
  colnames(res) <- c(rowLabel, colLabel, valueLabel)
  return(res)
}

#' Detect if any column has an empty string as name and fix
#' @param df
#' 
#' If any column has an empty string as name, its replaced by the prefix appended by an index starting from 1
#' @examples 
#' testDf <- data.frame("Col1"=LETTERS[1:3], "Col2"=letters[2:4])
#' colnames(testDf) <- c("", "")
#' testDf
#' fixEmptyColumnName(testDf)
#' fixEmptyColumnName(testDf, prefix="fData")
fixEmptyColumnName <- function(df, prefix="X") {
  isEmptyColName <- colnames(df)==""
  if(any(isEmptyColName)) {
    colnames(df)[isEmptyColName] <- paste0(prefix,
                                          seq(along=which(isEmptyColName)))
  }
  return(df)
}

vectorizeExprs <- function(exp) {
  if(is.data.frame(exp)) {
    expVec <- unlist(exp)
    rownames(expVec) <- rownames(exp)
    exprsLong <- as.data.frame(expVec)
  } else {
    exprsLong <- as.data.frame(as.vector(exp))    
  }
  return(exprsLong)
}

#' Transform eSet to long data.frame
#' 
#' @param x An \code{eSet} object
#' @param exprsFun A function to extract expression values, by default \code{exprs}
#' @param includeOtherAssayData Logical, whether other elements in the \code{assayData} environment (if present) should be returned.
#' 
#' The function extracts exprs (and other values in the \code{assayData} environment), and return it in a long data.frame format with phenotypic data
#' 
#' @examples 
#' data(ribios.ExpressionSet, package="ribiosExpression")
#' exprsLongTbl <- eSetToLongTable(ribios.ExpressionSet)
#' seLongTbl <- eSetToLongTable(ribios.ExpressionSet, exprsFun=function(eset) Biobase::assayData(eset)$se.exprs)
eSetToLongTable <- function(x, 
                            exprsFun=function(eset) Biobase::exprs(eset),
                            includeOtherAssayData=FALSE) {
  exp <- do.call(exprsFun, list(x))
  exprsLong <- vectorizeExprs(exp)
  colnames(exprsLong) <- "exprs"
  
  if(includeOtherAssayData) {
    otherAssayDataNames <- setdiff(ls(assayData(x)), "exprs")
    if(length(otherAssayDataNames)>0) {
      for(assay in otherAssayDataNames) {
        ad <- get(assay, assayData(x))
        adVec <- vectorizeExprs(ad)
        exprsLong <- cbind(exprsLong, adVec)
        colnames(exprsLong)[ncol(exprsLong)] <- assay
      }
    }
  }
  
  fData(x) <- fixEmptyColumnName(fData(x), prefix="fData.")
  pData(x) <- fixEmptyColumnName(pData(x), prefix="pData.")
  
  fDataCol <- colnames(fData(x))
  pDataCol <- colnames(pData(x))
  pfCommon <- intersect(fDataCol, pDataCol)
  
  for(i in colnames(fData(x))) {
    if(i %in% pfCommon) {
      inew <- sprintf("fData.%s", i)
    } else {
      inew <- i
    }
    exprsLong[, inew] <- rep(fData(x)[,i], dim(x)[2])
  }

  for(j in colnames(pData(x))) {
    if(j %in% pfCommon) {
      jnew <- sprintf("pData.%s", j)
    } else {
      jnew <- j
    }
    exprsLong[, jnew] <- rep(pData(x)[,j], each=dim(x)[1])
  }
  return(exprsLong)
}

