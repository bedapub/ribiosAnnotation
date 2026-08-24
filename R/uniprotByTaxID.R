#' @include utils.R sortAnnotationByQuery.R
NULL

#' Get Uniprot annotation with NCBI Taxonomy ID
#' @param taxid NCBI Taxonomy ID
#' @param orthologue Logical, whether human orthologues should be appended to
#' the annotation
#' @param multiOrth Logical, whether to return all orthologues or the (randomly)
#' selected top one if multiple exist. Only valid when \code{orthologue} is set 
#' as \code{TRUE}.
#' 
#' @return A \code{data.frame} with UniProt accessions and gene annotations.
#' @seealso
#' * \code{\link{annotateUniprotAccession}}, which annotates Uniprot accessions
#' * \code{\link{annotateTaxID}}, which annotates genes given TaxID.
#' @examples 
#' \dontrun{
#' humanUniprot <- uniprotByTaxID(9606)
#' }
#' @importFrom tidyr unnest
#' @export
uniprotByTaxID <- function(taxid, orthologue=FALSE, multiOrth=FALSE) {
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
  accessions <- Accession <- EntryName <- TaxID <- EnsemblGeneID <- NULL
  
  giCon <- connectMongoDB(instance="bioinfo_read",
                          collection="uniprot")
  
  fieldsJson <- returnFieldsJson(c("UniProtKB-AC",
                                   "UniProtKB-ID", "geneID", 
                                   "NCBI-taxon", "Ensembl"))
  query <- sprintf('{"NCBI-taxon":%d, "Ensembl":{"$exists": true}}', taxid)
  uniprots <- giCon$find(query, fields=fieldsJson) %>%
    tidyr::unnest("geneID") %>%
    tidyr::unnest("Ensembl")
  if(nrow(uniprots)==0) {
    featAnno <- data.frame(Accession=accessions, 
                           EntryName=NA,
                           GeneID=NA,
                           TaxID=NA,
                           EnsemblGeneID=NA)
  } else {
    featAnno <- uniprots %>%
      dplyr::rename("Accession"="UniProtKB-AC",
                    "EntryName"="UniProtKB-ID",
                    'GeneID'='geneID',
                    "TaxID"="NCBI-taxon",
                    "EnsemblGeneID"="Ensembl") %>%
      dplyr::mutate(GeneID=replace(GeneID, is.null(GeneID), NA)) %>%
      dplyr::select(Accession, EntryName, GeneID, TaxID, EnsemblGeneID) %>%
      unique
  }
  
  if(orthologue) {
    res <- appendHumanOrthologsWithNCBI(featAnno, multiOrth = multiOrth)
  } else {
    geneanno <- annotateGeneIDsWithoutHumanOrtholog(featAnno$GeneID) %>%
      dplyr::select(-TaxID)
    res <- left_join(featAnno, geneanno, by="GeneID")
  }
  res <- unique(res)
  return(res)
}
