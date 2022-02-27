#' @include utils.R
NULL

#' Retrieve human orthologs of genes of another species with its Taxonomy ID
#' @param taxid An integer, a NCBI taxonomy ID to identify a species, 
#'    for instance \code{10116} for rat, \code{10090} for mouse, and \code{9541}
#'    for cyno (crab-eating macaque).
#' @return A \code{data.frame} contains following columns:
#'  \itemize{
#'     \item{\code{GeneID}}{NCBI Gene ID of the query species}
#'     \item{\code{GeneSymbol}}{NCBI Gene symbol of the query species}
#'     \item{\code{Description}}{Gene description of the query species}
#'     \item{\code{HumanGeneID}}{NCBI Gene ID of the human homolog}
#'     \item{\code{HumanGeneSymbol}}{NCBI Gene symbol of the human homolog}
#'     \item{\code{HumanDescription}}{Gene description of the human homolog}
#'  }
#' @note
#' To query NCBI taxonomy IDs from free-text search, visit [NCBI Taxonomy 
#'   Browser](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi).
#' @importFrom magrittr %>%
#' @importFrom rjson fromJSON toJSON
#' @importFrom mongolite mongo
#' @importFrom dplyr right_join left_join select
#' @examples 
#' \dontrun{
#' ## human orthologs of rat genes
#' ratOrths <- humanOrthologsByTaxID(10116)
#' ## human orthologs of mouse genes
#' mouseOrths <- humanOrthologsByTaxID(10090)
#' ## human orthologs of cyno genes (crab-eating macaque, Macaca fascicularis)
#' cynoOrths <- humanOrthologsByTaxID(9541)
#' }
#' @export
humanOrthologsByTaxID <- function(taxid) {
  GeneID <- GeneSymbol <- Description <- NULL
  HumanGeneID <- HumanGeneSymbol <- HumanDescription <- NULL
  
  bioinfoReadSecrets <- loadMongodbSecrets(instance="bioinfo_read")
  bioinfoReadURL <- sprintf("mongodb://%s:%s@%s:%s/%s?authSource=%s", 
                            bioinfoReadSecrets$username,
                            bioinfoReadSecrets$password,
                            bioinfoReadSecrets$hostname, 
                            bioinfoReadSecrets$port,
                            bioinfoReadSecrets$dbname,
                            bioinfoReadSecrets$dbname)
  giCon <- mongolite::mongo(collection='ncbi_gene_info',
                            db=bioinfoReadSecrets$dbname,
                            url=bioinfoReadURL)
  goCon <- mongolite::mongo(collection='ncbi_gene_orthologs',
                            db=bioinfoReadSecrets$dbname,
                            url=bioinfoReadURL)
  
  speciesFields <- c("Symbol", "description", "geneId")
  speciesFieldsJson <- returnFieldsJson(speciesFields)
  speciesDf <- giCon$find(rjson::toJSON(list(taxId=taxid)),
                          fields=speciesFieldsJson) %>%
    dplyr::rename('GeneID'='geneId',
                  'GeneSymbol'='Symbol',
                  'Description'='description')
  humanFields <- c("Symbol", "description", "geneId")
  humanFieldsJson <- returnFieldsJson(humanFields)
  humanDf <- giCon$find(rjson::toJSON(list(taxId=9606)),
                        fields=humanFieldsJson) %>%
    dplyr::rename('HumanGeneID'='geneId',
                  'HumanGeneSymbol'='Symbol',
                  'HumanDescription'='description')
  ## in ortholog table, taxIs and geneId are human, while other_* are other species
  orthoFields <- c("geneId", "other_geneId")
  speciesHumanDf <- goCon$find(rjson::toJSON(list(taxId=9606,
                                                  other_taxId=taxid)),
                               fields=returnFieldsJson(orthoFields)) %>%
    dplyr::rename('HumanGeneID'='geneId',
                  'GeneID'='other_geneId')
  speciesHumanOrt <- speciesHumanDf %>%
    right_join(speciesDf, by="GeneID") %>%
    left_join(humanDf, by="HumanGeneID") %>%
    dplyr::select(GeneID, GeneSymbol, Description,
                  HumanGeneID, HumanGeneSymbol, HumanDescription)
  return(speciesHumanOrt)
}
