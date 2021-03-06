#' @include annotateAnyIDs.R
NULL

#' Annotate Ensembl gene/transcript/protein identifiers
#' 
#' Annotate Ensembl identifiers, optionally with human orthologue mappings
#' 
#' @param ids Character vector, Ensembl gene (ENSG)/transcript (ENST)/protein
#' (ENSG) identifiers . It can contain \code{NA} or \code{NULL}
#' @param orthologue Logical, whether human orthologues should be returnd
#' @param multiOrth Logical, whether multiple mapped orthologues should be
#' returned or not. Only useful when \code{orthologue} is \code{TRUE}, By
#' default \code{FALSE}. Be cautious with the \code{TRUE} option: in this case
#' the returning \code{data.frame} may have different row numbers as the input
#' vector.
#' @return A \code{data.frame} object containing the annotations.
#' 
#' If \code{orthologue} is set to \code{FALSE}, the data frame contains
#' following columns: \code{ENSEMBL_ID, GeneID, GeneSymbol, GeneName and TaxID}
#' 
#' If \code{orthologue} is \code{TRUE}, the data frame contains
#' \code{ENSEMBL_ID, GeneID, GeneSymbol,TaxID, OrigTaxID, OrigGeneID,
#' OrigGeneSymbol, OrigGeneName}. Note that \code{GeneID, GeneSymbol, TaxID}
#' contains the information of mapped orthologues, while \code{OrigXXX}
#' contains the information of queried genes.
#' 
#' If \code{multiOrth} is \code{TRUE} (only valid when \code{orthologue} is
#' \code{TRUE}), multiples orthologues mapped to the same gene are all
#' returned. This means that the result data frame may contain more rows than
#' the length of input vector. If set to \code{FALSE}, the resulting data frame
#' contains exact number of rows as the input vector length.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso See \code{\link{gtiChipAnnotation}} to get annotation for all
#' probesets in a chip.
#' 
#' See \code{\link{annotateGeneIDs}} to get annotation for Entrez GeneIDs.
#' 
#' See \code{\link{annotateProbesets}} to get annotation for probesets.
#' @importFrom ribiosUtils matchColumnIndex
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' ## normal use
#' annotateEnsembl(ids=c("ENSG00000197535", "ENSG00000105221", NA))
#' 
#' annotateEnsembl(ids=c("ENSG00000112062", "ENST00000229795", "ENSP00000229795", NULL))
#' 
#' annotateEnsembl(ids="ENSMUSG00000031710", orthologue=TRUE)
#' 
#' ## with version numbers
#' annotateEnsembl(ids=c("ENSG00000197535.1", "ENSG00000105221.2", "ENSG00000112062"))
#' 
#' annotateEnsembl(ids=c("ENSMUSG00000031710.1", 
#'   "ENSMUSG00000031710",
#'   "ENSMUSG00000041272.1"), orthologue=TRUE, multiOrth=TRUE)
#'   
#' options(error=NULL)
#' 
#' @export annotateEnsembl
annotateEnsembl <- function (ids, orthologue = FALSE, multiOrth = FALSE) {
  idsWoVersion <- removeEnsemblVersion(ids)
  comm <- paste("SELECT e.ENSEMBL_ID, c.RO_GENE_ID,c.GENE_SYMBOL, c.DESCRIPTION, c.TAX_ID", 
                " FROM GTI_GENES c INNER JOIN ENSEMBL_GENE e ON c.RO_GENE_ID=e.GENE_ID ", sep = "")
  ann <- querydbTmpTbl(comm, "e.ENSEMBL_ID", idsWoVersion, dbName(),
                       binUser(), 
                       binPwd())
  keyCol <- "EnsemblID"
  cnames <- c(keyCol, "GeneID", "GeneSymbol", "GeneName", "TaxID")
  conames <- c(keyCol, "OrigGeneID", "OrigGeneSymbol", "OrigGeneName", 
               "OrigTaxID")
  if (!orthologue) {
    colnames(ann) <- cnames
    res <- ann
  } else {
    colnames(ann) <- conames
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, 
                                           multiOrth = multiOrth)
    if (multiOrth) {
      res <- merge(ann, ort, by = "OrigGeneID", all.x = TRUE)
    }
    else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", 
                            multi = FALSE)
      res <- cbind(ann, ort.re[, -1L])
    }
    res <- putColsFirst(res, c(keyCol, "GeneID", "GeneSymbol", "TaxID", 
                               "OrigTaxID", "OrigGeneID", "OrigGeneSymbol", "OrigGeneName"))
  }
  resInd <- matchColumnIndex(idsWoVersion, res, keyCol, multi = orthologue && multiOrth)
  res <- res[unlist(resInd),,drop=FALSE]
  res[,keyCol] <- ids
  rownames(res) <- id2rownames(res[,keyCol])
  return(res)
}

#' Annotate EntrezGeneID with Ensemble identifiers
#' 
#' 
#' @param ids Character or integer vector, Entrez Gene IDs. It can contain
#'  \code{NA} or \code{NULL}.
#' @param orthologue Logical, whether human orthologues should be returnd
#' @param multiOrth Logical, whether multiple mapped orthologues should be
#' returned or not. Only useful when \code{orthologue} is \code{TRUE}, By
#' default \code{FALSE}. Be cautious with the \code{TRUE} option: in this case
#' the returning \code{data.frame} may have different row numbers as the input
#' vector.
#' @param type Character string, one of the following options: \code{gene}, 
#' \code{transcript}, and \code{protein}.
#' @return A \code{data.frame} object containing the annotations.
#' 
#' If \code{orthologue} is set to \code{FALSE}, the data frame contains
#' following columns: \code{GeneID, ENSEMBL_ID, GeneSymbol, GeneName and TaxID}
#' 
#' If \code{orthologue} is \code{TRUE}, the data frame contains
#' \code{GeneID, ENSEMBL_ID, GeneSymbol,TaxID, OrigTaxID, OrigGeneID,
#' OrigGeneSymbol, OrigGeneName}. Note that \code{GeneID, GeneSymbol, TaxID}
#' contains the information of mapped orthologues, while \code{OrigXXX}
#' contains the information of queried genes.
#' 
#' If \code{multiOrth} is \code{TRUE} (only valid when \code{orthologue} is
#' \code{TRUE}), multiples orthologues mapped to the same gene are all
#' returned. This means that the result data frame may contain more rows than
#' the length of input vector. If set to \code{FALSE}, the resulting data frame
#' contains exact number of rows as the input vector length.
#'
#' If either \code{transcript} or \code{protein} type is used, more than one 
#' items will be returned for each input ID.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso See \code{\link{annotateEnsembl}} to annotate Ensembl IDs.
#' 
#' See \code{\link{annotateGeneIDs}} to get annotation for Entrez GeneIDs.
#' 
#' See \code{\link{annotateProbesets}} to get annotation for probesets.
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' ## normal use
#' annotateGeneIDsWithEnsembl(c(1234, 1432, NULL))
#' annotateGeneIDsWithEnsembl(c(1234, 1432, NULL), type="transcript")
#' annotateGeneIDsWithEnsembl(c(1234, 1432, NULL), type="protein")
#' 
#' annotateGeneIDsWithEnsembl(ids="22227", orthologue=TRUE)
#' ## it is possible to mix ids of different species
#' annotateGeneIDsWithEnsembl(ids=c("1234", "22227"), orthologue=TRUE)
#' 
#' options(error=NULL)
#' 
#' @export annotateGeneIDsWithEnsembl
annotateGeneIDsWithEnsembl <- function (ids, 
                                        orthologue = FALSE, 
                                        multiOrth = FALSE,
                                        type=c("gene", "transcript", "protein")) {
  type <- match.arg(type)
  typeID <- c("gene"="26", "transcript"="30", "protein"="31")[type]
  comm <- paste("SELECT c.RO_GENE_ID, e.ENSEMBL_ID, c.GENE_SYMBOL, c.DESCRIPTION, c.TAX_ID", 
                " FROM GTI_GENES c INNER JOIN ENSEMBL_GENE e ON c.RO_GENE_ID=e.GENE_ID ", 
                paste0("WHERE e.TYPE_ID=", typeID), " ",
                sep = "")
  ann <- querydbTmpTbl(comm, "c.RO_GENE_ID", ids, dbName(),
                       binUser(), 
                       binPwd())
  keyCol <- ifelse(orthologue, "OrigGeneID", "GeneID")
  cnames <- c(keyCol, "EnsemblID", "GeneSymbol", "GeneName", "TaxID")
  conames <- c(keyCol, "EnsemblID", "OrigGeneSymbol", "OrigGeneName", 
               "OrigTaxID")
  if (!orthologue) {
    colnames(ann) <- cnames
    res <- ann
  } else {
    colnames(ann) <- conames
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, 
                                           multiOrth = multiOrth)
    if (multiOrth) {
      res <- merge(ann, ort, by = "OrigGeneID", all.x = TRUE)
    }  else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", 
                            multi = FALSE)
      res <- cbind(ann, ort.re[, -1L])
    }
    res <- putColsFirst(res, c("OrigGeneID", "GeneID", "GeneSymbol", "TaxID", 
                               "OrigTaxID", "OrigGeneSymbol", "OrigGeneName"))
  }
  isTranscriptOrProtein <- type!="gene"
  res <- matchColumn(ids, res, keyCol, multi = (orthologue && multiOrth) || isTranscriptOrProtein)
  rownames(res) <- id2rownames(res[,keyCol])
  return(res)
}
