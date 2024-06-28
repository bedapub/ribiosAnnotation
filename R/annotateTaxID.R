#' @include utils.R
NULL

#' Annotations of all genes associated with the given TaxID
#' 
#' The function returns annotations (see details below) of all features
#' (probably probesets) associated with the given taxon.
#' 
#' The function reads from the backend, the MongoDB bioinfo database.
#' 
#' @param taxId Integer, the TaxID of the species in interest. For
#' instance \sQuote{9606} for Homo sapiens.
#' @param orthologue Logical, whether human orthologues should be returned
#' @param multiOrth Logical, in case \code{orthologue} is set to \code{TRUE}, whether
#' to return multiple orthologues for one gene. Default: FALSE
#' @return A \code{data.frame} object with very similar structure as the
#' \code{EG_GENE_INFO} table in the database. In case \code{orthologue} is \code{TRUE}, additional
#' columns containing human orthologue information are returned.
#' 
#' Rownames of the \code{data.frame} are set to \code{NULL}.
#' 
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
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
#' 
#' mtOrthAnno <- annotateTaxID(10090, orthologue=TRUE)
#' dim(mtOrthAnno)
#' head(mtOrthAnno)
#' 
#' pigMultiOrthAnno <- annotateTaxID(9823, orthologue=TRUE, multiOrth=TRUE)
#' dim(pigMultiOrthAnno)
#' head(pigMultiOrthAnno)
#' }
#'
#' @importFrom dplyr any_of rename mutate
#' @importFrom mongolite mongo
#' @importFrom rjson toJSON
#' @importFrom magrittr '%>%'
#' @export annotateTaxID
annotateTaxID <- function(taxId, orthologue=FALSE, multiOrth=FALSE) {
  taxId <- checkSingleIntegerTaxId(taxId)

  
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
               'Type'='type_of_gene')
  speciesDf <- giCon$find(rjson::toJSON(list(taxId=taxId)),
                          fields=speciesFieldsJson) %>%
    dplyr::rename(dplyr::any_of(renames)) %>%
    dplyr::select(dplyr::any_of(names(renames))) %>%
    dplyr::mutate(TaxID=taxId)
  rownames(speciesDf) <- speciesDf$GeneID

  if(orthologue) {
    res <- appendHumanOrthologsWithNCBI(speciesDf, multiOrth=multiOrth)
  } else {
    res <- speciesDf
  }
  
  return(res)
}
