#' @include utils.R sortAnnotationByQuery.R
NULL

#' Annotate UniProt accessions or names
#' @param accessions Character strings, UniProt accessions or names
#' @param orthologue Logical, whether orthologues are returned
#' @param multiOrth Logical, only valid if \code{orthologue} is \code{TRUE},
#' whether multiple orthologues are returned instead of only one.
#' 
#' @examples 
#' \dontrun{
#' annotateUniprotAccession(c("B4E0K5"))
#' }
annotateUniprotAccession <- function(accessions,
                                     orthologue=FALSE,
                                     multiOrth=FALSE) {
  validIDs <- accessions[!is.na(accessions)]
  
  Accession <- EntryName <- GeneID <- TaxID <- Ensembl <- NULL
  EnsemblGeneID <- NULL
  
  giCon <- connectMongoDB(instance="bioinfo_read",
                          collection="uniprot")
  
  fieldsJson <- returnFieldsJson(c("UniProtKB-AC",
                                   "UniProtKB-ID", "geneID", "NCBI-taxon", "Ensembl"))
  query <- paste0('{"UniProtKB-AC":{"$in":[', 
                  paste(paste0("\"", as.character(validIDs), "\""), 
                        collapse=","),']}}')
  uniprots <- giCon$find(query, fields=fieldsJson)
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
      dplyr::select(Accession, EntryName, GeneID, TaxID, EnsemblGeneID)
  }
  
  if(orthologue) {
    res <- appendHumanOrthologsWithNCBI(featAnno, multiOrth = multiOrth)
  } else {
    geneanno <- annotateGeneIDsWithoutHumanOrtholog(featAnno$GeneID) %>%
      dplyr::select(-TaxID)
    res <- left_join(featAnno, geneanno, by="GeneID")
  }
  res <- sortAnnotationByQuery(res, accessions, "Accession", multi=FALSE)
  return(res)
}
