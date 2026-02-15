#' @include utils.R removeEnsemblVersion.R sortAnnotationByQuery.R
NULL

#' Annotate any identifiers
#' 
#' This annotates any identifies that can be recognized by GTI.
#' 
#' 
#' @param ids A vector of identifiers. They must be of the same type. Supported 
#' types include Entrez GeneID, GeneSymbol, Probesets, UniProt identifiers, 
#' NCBI RefSeq mRNA identifiers, and Ensembl gene identifiers (with possible 
#' version suffixes).
#' 
#' @param orthologue Logical, is orthologous mapping needed?
#' @param multiOrth Logical, is more than one orthologs allowed
#' @return A \code{data.frame} containing annotation information. Following
#'    columns exist at least: 
#' \enumerate{
#'   \item \code{Input} Input string, it will be in the first column.
#'   \item \code{IDType} Input ID type
#'   \item \code{GeneID} (Human) Entrez GeneID
#'   \item \code{GeneSymbol} (Human) official gene symbol
#'   \item \code{GeneName} (Human) gene name
#'   \item \code{TaxID} NCBI taxonomy ID
#' }
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{annotateGeneIDs}}, \code{\link{annotateGeneSymbols}}
#' @importFrom ribiosUtils putColsFirst matchColumn
#' @examples
#' \dontrun{
#' # GeneID
#' annotateAnyIDs(ids=c(780, 5982, 3310, NA))
#'
#' # GeneSymbol
#' annotateAnyIDs(ids=c("DDR1", "RFC2", "HSPA6", "HSAP6"))
#'
#' # Probesets
#' myprobes <- c("1000_at", "1004_at", "1002_f_at", "nonsense_at")
#' annotateAnyIDs(myprobes)
#'
#' # UniProt
#' annotateAnyIDs(ids=c("P38398", "Q8NDF8"))
#'
#' # EnsEMBL
#' ensemblIDs <- c("ENSG00000197535", "ENST00000399231.7", "ENSP00000418960.2")
#' annotateAnyIDs(ensemblIDs)
#'
#' # RefSeq
#' annotateAnyIDs(c("NM_000235", "NM_000498"))
#' }
#'
#' @export annotateAnyIDs
annotateAnyIDs <- function(ids, orthologue = FALSE, multiOrth = FALSE) {
  validIDs <- removeEnsemblVersion(ids)
  validIDs <- validIDs[!is.na(validIDs)]
  
  Input <- GeneID <- TaxID <- IDType <- NULL
  
  giCon <- connectMongoDB(instance="bioinfo_read",
                          collection="featureanno")
  
  fieldsJson <- returnFieldsJson(c("feature_id", "tax_id", "gene_id", "id_type"))
  query <- paste0('{"feature_id":{"$in":[', 
                  paste(paste0("\"", as.character(validIDs), "\""), 
                        collapse=","),']}}')
  genes <- giCon$find(query, fields=fieldsJson)
  if(nrow(genes)==0) {
    featAnno <- data.frame(Input=ids, 
                           GeneID=NA,
                           TaxID=NA,
                           IDType=NA)
  } else {
    featAnno <- genes %>%
      dplyr::rename("Input"="feature_id",
                    'GeneID'='gene_id',
                    "TaxID"="tax_id",
                    "IDType"="id_type") %>%
      dplyr::select(Input, GeneID, TaxID, IDType)
  }
  
  if(orthologue) {
    res <- appendHumanOrthologsWithNCBI(featAnno, multiOrth = multiOrth)
  } else {
    geneanno <- annotateGeneIDsWithoutHumanOrtholog(featAnno$GeneID) %>%
      dplyr::select(-TaxID)
    res <- left_join(featAnno, geneanno, by="GeneID")
  }
  res <- sortAnnotationByQuery(res, ids, "Input", multi=FALSE)
  return(res)
}
