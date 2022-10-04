#' Sort the annotation table by query IDs
#' @param anno A \code{data.frame} containing annotations
#' @param ids A vector of character or integer, identifiers used to query
#'  the annotation
#' @param id_column Character, column of the data frame where ids can be found. 
#' @param multi In case that an identifier appears more than once in 
#' \code{id_column}, should all rows be returned or only the first row? Default:
#' \code{FALSE}, namely only the first row is returned.
#' 
#' @return A \code{data.frame} sorted by the query identifiers, with the
#' column \code{id_column} containing exactly the same value as \code{ids}.
#' If the identifiers are unique and if they do not contain NA, they are used as
#' the row names of the data.frame; otherwise, \code{NULL} will be used.
#' 
#' @examples 
#' myAnno <- data.frame(GeneID=c(4,6,5), GeneName=c("Gene4", "Gene6", "Gene5"))
#' inputIds <- c("6", "5", "6", "4", "NotAGeneID")
#' sortAnnotationByQuery(myAnno, inputIds, "GeneID")
#' 
#' myAnno2 <- data.frame(GeneID=c(4,6,5, 5), 
#'                       GeneName=c("Gene4", "Gene6", "Gene5", "Gene5.V2"))
#' inputIds <- c("6", "5", "6", "4", "NotAGeneID")
#' sortAnnotationByQuery(myAnno2, inputIds, "GeneID")
#' sortAnnotationByQuery(myAnno2, inputIds, "GeneID", multi=TRUE)
#' @export
sortAnnotationByQuery <- function(anno, ids, id_column="GeneID", multi=FALSE) {
  stopifnot(id_column %in% colnames(anno))
  anno <- unique(anno)
  resInd <- ribiosUtils::matchColumnIndex(ids, anno, id_column, multi=multi)
  if(is.list(resInd)) {
    resIndLen <- sapply(resInd, length)
    ids <- rep(ids, resIndLen)
    resInd <- unlist(resInd)
  }
  res <- anno[resInd, , drop=FALSE]
  rownames(res) <- id2rownames(ids)
  res[, id_column] <- ids
  return(res)
}
