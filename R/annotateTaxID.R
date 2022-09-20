#' @include utils.R
NULL

#' Annotations of all genes associated with the given TaxID
#' 
#' The function returns annotations (see details below) of all features
#' (probably probesets) associated with the given taxon.
#' 
#' The function reads from the backend, the MongoDB bioinfo database.
#' 
#' @param taxid Integer, the TaxID of the species in interest. For
#' instance \sQuote{9606} for Homo sapiens.
#' @return A \code{data.frame} object with very similar structure as the
#' \code{EG_GENE_INFO} table in the database.
#' 
#' Rownames of the \code{data.frame} are set to \code{NULL}.
#' 
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{gtiChipAnnotation}} for annotating chips.
#' @references \url{http://bioinfo.bas.roche.com:8080/bicgi/gti_ui.cgi}
#' @examples
#' 
#' \dontrun{
#' hsAnno <- annotateTaxID("9606")
#' dim(hsAnno)
#' head(hsAnno)
#' 
#' hsMtAnno <- annotateTaxID("10092")
#' dim(hsMtAnno)
#' head(hsMtAnno)
#' }
#'
#' @importFrom dplyr any_of rename
#' @importFrom mongolite mongo
#' @importFrom rjson toJSON
#' @importFrom magrittr '%>%'
#' @export annotateTaxID
annotateTaxID <- function(taxid) {
  if(missing(taxid)) {
    stop("'taxid' cannot be missing.")
  }
  if(!is.integer(taxid)) {
    taxid <- as.integer(taxid)
  }
  if(is.na(taxid) || !is.integer(taxid)) {
    stop("'taxid' should be an integer.")
  }
  
  GeneID <- GeneSymbol <- Description <- NULL

  giCon <- connectMongoDB(instance="bioinfo_read",
                          collection="ncbi_gene_info")
  
  speciesFields <- c("Symbol", "description", "geneId",
                     "chromosome", "map_location", "type_of_gene")
  speciesFieldsJson <- returnFieldsJson(speciesFields)
  renames <- c('GeneID'='geneId',
               'GeneSymbol'='Symbol',
               'Description'='description',
               'Chromosome'='chromosome',
               'MapLocation'='map_location',
               'TypeOfGene'='type_of_gene')
  speciesDf <- giCon$find(rjson::toJSON(list(taxId=taxid)),
                          fields=speciesFieldsJson) %>%
    dplyr::rename(dplyr::any_of(renames)) %>%
    dplyr::select(dplyr::any_of(names(renames)))
  rownames(speciesDf) <- speciesDf$GeneID
  return(speciesDf)
}
