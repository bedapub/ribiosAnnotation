#' Common species taxonomy IDs
#' @format ## `commonSpecies`
#' The data.frame contains three columns: 
#' \describe{
#'   \item{TaxID}{NCBI taxonomy ID}
#'   \item{ScientificName}{Scientific name}
#'   \item{CommonName}{Common name
#' }
"commonSpecies"

commonTaxIDs <- c("Homo sapiens"=9606,
                  "Pan troglodytes"=9598,
                  "Macaca fascicularis"=9541,
                  "Mus musculus"=10090,
                  "Rattus norvegicus"=10116,
                  "Oryctolagus cuniculus"=9986,
                  "Cricetinae gen. sp."=36483,
                  "Canis lupus familiaris"=9615,
                  "Sus scrofa"=9823,
                  "Gallus gallus"=9031,
                  "Naio rerio"=7955)
#' Get all taxonomy ID and scientific names offered by NCBI
#' @return A data.frame containing two columns, \code{TaxID} and \code{ScientificName}.
#' @examples
#' \dontrun{
#'     all_tax_ids <- getAllTaxIDs()
#' }
#' @export
getAllTaxIDs <- function() {
  taxId <- scientific_name <- NULL
  giCon <- connectMongoDB(instance = "bioinfo_read", 
                          collection = "ncbi_taxonomy")
  taxFieldsJson <- returnFieldsJson(c("taxId", "scientific_name"))
  renames <- c("TaxID"="taxId", "ScientificName"="scientific_name")
  speciesDf <- giCon$find(fields=taxFieldsJson) %>%
    dplyr::rename(dplyr::any_of(renames)) %>% 
    dplyr::select(dplyr::any_of(names(renames)))
  rownames(speciesDf) <- speciesDf$TaxID
  return(speciesDf)
}
