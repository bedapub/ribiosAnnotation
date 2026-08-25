#' @include utils.R annotateAnyIDs.R sortAnnotationByQuery.R backend.R
#' @include appendHumanOrthologsWithNCBI.R
NULL

#' Annotate Entrez GeneIDs
#' 
#' @param ids A vector of integers or characters, encoding NCBI Entrez GeneIDs.
#' It can contain \code{NA} or \code{NULL}.
#' @param orthologue Logical, whether human orthologues should be returned. 
#' Default: \code{FALSE}
#' @param multiOrth Logical, whether mutliple orthologues should be returned if
#' exist. Deafult: \code{FALSE}
#' @param backend Character string, data backend to use. One of
#' \sQuote{\code{mongodb}} (default) or \sQuote{\code{bioconductor}}.
#' @return A \code{data.frame} object containing the annotations:
#' * GeneID EntrezGeneID
#' * GeneSymbol Official gene symbol
#' * Description Gene description
#' * TaxID Taxonomy ID
#' * Type Gene type
#' 
#' If \code{orthologue} is \code{TRUE}, following columns are appended:
#' * HumanGeneID
#' * HumanGeneSymbol
#' * HumanDescription
#' * HumanType
#' @seealso \code{\link{annotateGeneIDsWithoutHumanOrtholog}} and
#' \code{\link{annotateGeneIDsWithHumanOrtholog}}
#' @examples 
#' \dontrun{
#'   annotateGeneIDs(ids=c(780, 5982, 3310))
#'   annotateGeneIDs(ids=c(780, 5982, 3310, NA), orthologue=TRUE)
#'   annotateGeneIDs(ids=c(780, 1506, 1418,
#'                         114483548, 57300, 
#'                         20, 1506, 102129055),
#'                         orthologue=TRUE)
#' }
#' @export
annotateGeneIDs <- function(ids, orthologue=FALSE, multiOrth=FALSE,
                            backend = NULL) {
  backend <- normalizeAnnotationBackend(backend)
  if(!orthologue) {
    res <- annotateGeneIDsWithoutHumanOrtholog(ids, backend = backend)
  } else {
    res <- annotateGeneIDsWithHumanOrtholog(ids, multiOrth = multiOrth,
                                            backend = backend)
  }
  return(res)
}

#' Annotate Entrez GeneIDs without querying human orthologs
#' 
#' @param ids A vector of integers or characters, encoding NCBI Entrez GeneIDs.
#' It can contain \code{NA} or \code{NULL}.
#' @param backend Character string, data backend to use. One of
#' \sQuote{\code{mongodb}} (default) or \sQuote{\code{bioconductor}}.
#' @return A \code{data.frame} object containing the annotations:
#' * GeneID EntrezGeneID
#' * GeneSymbol Official gene symbol
#' * Description Gene description
#' * TaxID Taxonomy ID
#' * Type Gene type
#'
#' @note \code{annotatemRNAs} is an alias of \code{annotateRefSeqs}
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#'
#' @details The collection \code{ncbi_gene_info} is used.
#' 
#' @examples
#' 
#' \dontrun{
#'   annotateGeneIDsWithoutHumanOrtholog(ids=c(780, 5982, 3310))
#'   annotateGeneIDsWithoutHumanOrtholog(ids=c(780, 5982, 3310, NA))
#'   annotateGeneIDsWithoutHumanOrtholog(ids=c(780, 5982, 3310, NULL))
#'   annotateGeneIDsWithoutHumanOrtholog(ids=c(780, 5982, 3310, "1418"))
#'   annotateGeneIDsWithoutHumanOrtholog(ids=c(780, 5982, 3310, "NotValidGeneID"))
#'   annotateGeneIDsWithoutHumanOrtholog(ids=c(780, 5982, 3310, 1418, 5982))
#'   annotateGeneIDsWithoutHumanOrtholog(ids=c(780, 5982, 1418, 5982, 25120,
#'                         114483548, 57300, 20, 1506,
#'                         1545, 102129055))
#' }
#' 
#' @export
annotateGeneIDsWithoutHumanOrtholog <- function(ids, backend = NULL) {
  backend <- normalizeAnnotationBackend(backend)
  GeneID <- GeneSymbol <- Description <- TaxID <- Type <- NULL
  
  intID <- suppressWarnings(as.integer(ids))
  validID <- intID[!is.na(intID)]

  if (backend == "mongodb") {
    giCon <- connectMongoDB(instance="bioinfo_read",
                            collection="ncbi_gene_info")

    speciesFieldsJson <- returnFieldsJson(c("Symbol", "description", "geneId",
                                            "taxId", "type_of_gene"))
    query <- paste0('{"geneId":{"$in":[', paste(as.character(validID), collapse=","),']}}')
    genes <- giCon$find(query, fields=speciesFieldsJson)
  } else {
    genes <- biocAnnotateGeneIDs(validID) %>%
      dplyr::rename(Symbol = GeneSymbol,
                    description = Description,
                    geneId = GeneID,
                    taxId = TaxID,
                    type_of_gene = Type)
  }

  if(nrow(genes)==0) {
    res <- data.frame(GeneID=ids,
                      GeneSymbol=NA,
                      Description=NA,
                      TaxID=NA,
                      Type=NA)
  } else {
    res <- genes %>%
      dplyr::rename('GeneID'='geneId',
                    'GeneSymbol'='Symbol',
                    'Description'='description',
                    "TaxID"="taxId",
                    "Type"="type_of_gene") %>%
      dplyr::select(GeneID, GeneSymbol, Description, TaxID, Type)
  }
  res <- sortAnnotationByQuery(res, ids, "GeneID", multi=FALSE)
  return(res)
}

#' Annotate Entrez GeneIDs with the query of human orthologs
#' @param ids Vector of integer or character strings, EntrezIDs to be annotated
#' @param multiOrth Logical, whether mutliple orthologues should be returned if
#' exist. Deafult: \code{FALSE}
#' @param backend Character string, data backend to use. One of
#' \sQuote{\code{mongodb}} (default) or \sQuote{\code{bioconductor}}.
#' @return A \code{data.frame} object containing the annotations:
#' * GeneID EntrezGeneID
#' * GeneSymbol Official gene symbol
#' * Description Gene description
#' * TaxID Taxonomy ID
#' * Type Gene type
#' * HumanGeneID
#' * HumanGeneSymbol
#' * HumanDescription
#' * HumanType
#' 
#' @examples
#' 
#' \dontrun{
#'   annotateGeneIDsWithHumanOrtholog(ids=c(780, 5982, 3310))
#'   annotateGeneIDsWithHumanOrtholog(ids=c(780, 5982, 3310, NA))
#'   annotateGeneIDsWithHumanOrtholog(ids=c(780, 5982, 3310, NULL))
#'   annotateGeneIDsWithHumanOrtholog(ids=c(780, 5982, 3310, "1418"))
#'   annotateGeneIDsWithHumanOrtholog(ids=c(780, 5982, 3310, "NotValidGeneID"))
#'   annotateGeneIDsWithHumanOrtholog(ids=c(780, 5982, 3310, 1418, 5982))
#'   annotateGeneIDsWithHumanOrtholog(ids=c(780, 5982, 1418, 5982, 25120,
#'                         114483548, 57300, 20, 1506,
#'                         1545, 102129055))
#' }
#' 
#' @export
annotateGeneIDsWithHumanOrtholog <- function(ids, multiOrth=FALSE,
                                             backend = NULL) {
  backend <- normalizeAnnotationBackend(backend)
  HumanGeneID <- HumanGeneSymbol <- HumanDescription <- Type <- NULL
  GeneSymbol <- GeneID <- Description <- NULL
  
  geneIdAnno <- annotateGeneIDsWithoutHumanOrtholog(ids, backend = backend)
  appGeneIDanno <- appendHumanOrthologsWithNCBI(geneIdAnno,
                                                multiOrth = multiOrth,
                                                backend = backend)
  res <- sortAnnotationByQuery(appGeneIDanno, ids, "GeneID", multi = multiOrth)
  
  return(res)
}
