#' @include utils.R
NULL


#' Annotate gene symbols without human ortholog
#' @param ids Character vector, gene symbols to be queried
#' @param taxId Integer, NCBI taxonomy ID of the species. Default value: 9606 (human). See \code{commonSpecies} for tax id of common species.
#' @return A data.frame of following columns
#' \describe{
#'   \item{GeneID}{Entrez Gene ID}
#'   \item{GeneSymbol}{Official gene symbols}
#'   \item{Description}{Description}
#'   \item{TaxID}{NCBI Taxonomy ID}
#'   \item{Type}{Gene type}
#' }
#' @examples
#' \dontrun{
#'     annotateGeneSymbolsWithoutHumanOrtholog(c("AKT1", "ERBB2", 
#'                                             "NoSuchAGene", "TGFBR1"), 9606)
#'     annotateGeneSymbolsWithoutHumanOrtholog(c("Akt1", "Erbb2", 
#'                                             "NoSuchAGene", "Tgfbr1"), 10090)
#' }
#' @export
annotateGeneSymbolsWithoutHumanOrtholog <- function(ids,
                                                    taxId=9606) {
  GeneID <- GeneSymbol <- Description <- TaxID <- Type <- NULL
  
  giCon <- connectMongoDB(instance="bioinfo_read",
                          collection="ncbi_gene_info")
  
  qids <- paste0('"', ids, '"')
  speciesFieldsJson <- returnFieldsJson(c("Symbol", "description", "geneId",
                                          "taxId", "type_of_gene"))
  query <- paste0('{"Symbol":{"$in":[', paste(qids, collapse=","),']}, "taxId":', as.character(taxId), '}')
  genes <- giCon$find(query, fields=speciesFieldsJson) 
  if(nrow(genes)==0) {
    res <- data.frame(GeneSymbol=ids,
                      GeneID=NA,
                      Description=NA,
                      TaxID=taxId,
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
  res <- sortAnnotationByQuery(res, ids, "GeneSymbol", multi=FALSE) %>%
    dplyr::mutate(TaxID=taxId)
  return(res)
}

#' Annotate GeneSymbol with human ortholog
#' @param ids Character strings, gene symbols
#' @param taxId Integer, NCBI taxonomy ID. Default value: 9606 (human). See \code{commonSpecies} for tax id of common species.
#' @param multiOrth Logical, only valid when orthologue is set to TRUE,  whether
#' multiple orthologues are returned
#' @return A data.frame containing following columns:
#' \describe{
#'   \item{GeneID}{Entrez Gene ID}
#'   \item{GeneSymbol}{Official gene symbols}
#'   \item{Description}{Description}
#'   \item{TaxID}{NCBI Taxonomy ID}
#'   \item{Type}{Gene type}
#'   \item{HumanGeneID}{Human orthologue Entrez GeneID}
#'   \item{HumanGeneSymbol}{Human orthologue official gene symbol}
#'   \item{HumanDescription}{Human orthologue gene description}
#'   \item{HumanType}{Human orthologue gene type}
#' }
#' @examples
#' \dontrun{
#'     annotateGeneSymbolsWithHumanOrtholog(c("Akt1", "Erbb2",
#'                                          "NoSuchAGene", "Tgfbr1"), 
#'                                          taxId=10090, multiOrth=FALSE)
#' }
#' @export
annotateGeneSymbolsWithHumanOrtholog <- function(ids, taxId, multiOrth=FALSE) {
  HumanGeneID <- HumanGeneSymbol <- HumanDescription <- Type <- NULL
  GeneSymbol <- GeneID <- Description <- NULL
  
  geneIdAnno <- annotateGeneSymbolsWithoutHumanOrtholog(ids, taxId=taxId)
  appGeneIDanno <- appendHumanOrthologsWithNCBI(geneIdAnno, multiOrth=multiOrth)
  res <- sortAnnotationByQuery(appGeneIDanno, ids, "GeneSymbol", multi = multiOrth)
  
  return(res)
}


#' Annotate GeneSymbols
#' @param ids Character strings, gene symbols
#' @param taxId Integer, NCBI taxonomy ID. Default value: 9606 (human). See \code{commonSpecies} for tax id of common species.
#' @param orthologue Logical, whether orthologues are to be returned
#' @param multiOrth Logical, only valid when orthologue is set to TRUE, whether
#' multiple orthologues are returned
#' 
#' @seealso The function is a convenient wrapper of two functions: \code{annotateGeneSymbolsWithoutHumanOrtholog} and \code{annotateGeneSymbolsWithHumanOrtholog}.
#' @return A data.frame containing following columns
#' \describe{
#'   \item{GeneID}{Entrez Gene ID}
#'   \item{GeneSymbol}{Official gene symbols}
#'   \item{Description}{Description}
#'   \item{TaxID}{NCBI Taxonomy ID}
#'   \item{Type}{Gene type}
#' }
#' If \code{orthologue} is \code{TRUE}, then additional columns are appended:
#' \describe{
#'   \item{HumanGeneID}{Human orthologue Entrez GeneID}
#'   \item{HumanGeneSymbol}{Human orthologue official gene symbol}
#'   \item{HumanDescription}{Human orthologue gene description}
#'   \item{HumanType}{Human orthologue gene type}
#' }
#' @examples
#' \dontrun{
#'     annotateGeneSymbols(c("AKT1", "ERBB2", "NoSuchAGene", "TGFBR1"), 9606)
#'     annotateGeneSymbols(c("Akt1", "Erbb2", "NoSuchAGene", "Tlr7"), 
#'                         taxId=10116, orthologue=FALSE)
#'     annotateGeneSymbols(c("Akt1", "Erbb2", "NoSuchAGene", "Tlr7"), 
#'                         taxId=10116, orthologue=TRUE)
#'     annotateGeneSymbols(c("Akt1", "Erbb2", "NoSuchAGene", "Tlr7"), 
#'                         taxId=10116, orthologue=TRUE, multiOrth=TRUE)
#' }
#' @export
annotateGeneSymbols <- function(ids, taxId=9606, orthologue=FALSE, multiOrth=FALSE) {
  if(!orthologue) {
    res <- annotateGeneSymbolsWithoutHumanOrtholog(ids, taxId=taxId)
  } else {
    res <- annotateGeneSymbolsWithHumanOrtholog(ids, taxId=taxId, multiOrth=TRUE)
  }
  return(res)
	return(NULL)
}
