#' @include orderAnnotationByQuery.R
NULL


#' Annotate human orthologs with data from NCBI
#' @param geneids Integer GeneIDs, can contain human GeneIDs.
#' @param multiOrth Logical, whether one gene is allowed to map to multiple
#' human orthologs? Default value is \code{FALSE}, i.e. only the first 
#' human ortholog (random choice) is returned.
#' 
#' @return A \code{data.frame} containing following columns:
#' * GeneID: Input GeneID
#' * TaxID: Taxonomy ID of the input gene
#' * HumanGeneID: Human Entrez GeneID
#' 
#' This function annotates human orthologs for any non-human genes. Use 
#' \code{\link{annotateHumanOrthologsWithNCBI}} if you are not sure
#' whether all input GeneIDs are non-human.
#' 
#' @seealso \code{\link{annotateHumanOrthologsWithNCBI}}
#' 
#' @details \code{ncbi_gene_info} and \code{ncbi_gene_orthologs} collections 
#' are used.
#' 
#' @importFrom dplyr mutate
#' @examples
#' \dontrun{
#' annotateNonHumanGenesHumanOrthologsWithNCBI(c(25120, 114483548, 57300, 20))
#' annotateNonHumanGenesHumanOrthologsWithNCBI(c(25120, 114483548, 57300, 
#'                                               20, 1506, 1545, 102129055))
#' }
#' @export
annotateNonHumanGenesHumanOrthologsWithNCBI <- function(geneids,
                                                        multiOrth=FALSE) {
  
  validIDs <- suppressWarnings(unique(as.integer(as.character(geneids))))
  validIDs <- validIDs[!is.na(validIDs)]
  
  GeneID <- TaxID <- HumanGeneID <- NULL
  
  giCon <- connectMongoDB(instance="bioinfo_read",
                          collection="ncbi_gene_orthologs")
  
  ## in ortholog table, taxIs and geneId are human, while other_* are 
  ## other species
  orthologFields <- returnFieldsJson(c("geneId",
                                       "other_taxId", "other_geneId"))
  ortQuery <- paste0('{"other_geneId":{"$in":[', 
                     paste(as.character(validIDs), collapse=","),
                     ']}, "taxId":9606}')
  ortDf <- giCon$find(ortQuery,
                      fields=orthologFields) %>%
    dplyr::rename('GeneID'='other_geneId',
                  'TaxID'='other_taxId',
                  'HumanGeneID'='geneId') %>%
    dplyr::select(GeneID, TaxID, HumanGeneID)
  
  res <- sortAnnotationByQuery(ortDf, 
                               geneids, 
                               id_column="GeneID",
                               multi=multiOrth)
  return(res)
}



#' Annotate human orthologs with data from NCBI
#' @param geneids Integer GeneIDs, can contain human GeneIDs.
#' @param multiOrth Logical, whether one gene is allowed to map to multiple
#' human orthologs? Default value is \code{FALSE}, i.e. only the first 
#' human ortholog (random choice) is returned.
#' 
#' @return A \code{data.frame} containing following columns:
#' * GeneID: Input GeneID
#' * TaxID: Taxonomy ID of the input gene
#' * HumanGeneID: Human Entrez GeneID
#' 
#' This function annotates human orthologs for any GeneID, including human
#' genes, in which case the ortholog will be itself. Use 
#' \code{\link{annotateNonHumanGenesHumanOrthologsWithNCBI}} if you are sure
#' that input GeneIDs do not come from human.
#' 
#' \seealso \code{\link{annotateNonHumanGenesHumanOrthologsWithNCBI}} 
#' 
#' @details \code{ncbi_gene_info} and \code{ncbi_gene_orthologs} collections 
#' are used.
#' @importFrom dplyr mutate
#' @examples
#' \dontrun{
#' annotateHumanOrthologsWithNCBI(c(25120, 114483548, 57300, 20))
#' annotateHumanOrthologsWithNCBI(c(25120, 114483548, 57300, 20, 1506, 1545,
#'                                  102129055))
#' }
#' @export
annotateHumanOrthologsWithNCBI <- function(geneids, multiOrth=FALSE) {
  
  validIDs <- suppressWarnings(unique(as.integer(as.character(geneids))))
  validIDs <- validIDs[!is.na(validIDs)]
  
  GeneID <- TaxID <- HumanGeneID <- NULL
  
  ## we first find out which GeneIDs are human GeneIDs
  humanCon <- connectMongoDB(instance="bioinfo_read",
                             collection="ncbi_gene_info")
  
  humanQuery <- paste0('{"geneId":{"$in":[', 
                       paste(as.character(validIDs), collapse=","),
                       ']}, "taxId":9606}')
  humanDf <- humanCon$find(humanQuery,
                           fields=returnFieldsJson("geneId"))
  if(nrow(humanDf)==0) {
    humanDf <- data.frame(GeneID=NA, TaxID=NA, HumanGeneID=NA)
  } else {
    humanDf <- humanDf %>%
      dplyr::rename('GeneID'='geneId') %>%
      dplyr::mutate(TaxID=9606, HumanGeneID=GeneID) %>%
      dplyr::select(GeneID, TaxID, HumanGeneID)
  }
  
  humanValidGeneIDs <- unique(humanDf$GeneID)
  nonhumanValidGeneIDs <- setdiff(validIDs, humanValidGeneIDs)
  
  if(length(nonhumanValidGeneIDs)>0) {
    ortDf <- annotateNonHumanGenesHumanOrthologsWithNCBI(nonhumanValidGeneIDs,
                                                         multiOrth = multiOrth)
    combDf <- rbind(humanDf, ortDf)
  } else {
    combDf <- humanDf
  }
  
  res <- sortAnnotationByQuery(combDf, 
                               geneids, 
                               id_column="GeneID",
                               multi=multiOrth)
  return(res)
}
