#' @include utils.R sortAnnotationByQuery.R backend.R
NULL

#' Annotate UniProt accessions or names
#' @param accessions Character strings, UniProt accessions or names
#' @param orthologue Logical, whether orthologues are returned
#' @param multiOrth Logical, only valid if \code{orthologue} is \code{TRUE},
#' whether multiple orthologues are returned instead of only one.
#' @param backend Character string, data backend to use. One of
#' \code{"mongodb"} (default) or \code{"bioconductor"}.
#' 
#' @examples 
#' \dontrun{
#' annotateUniprotAccession(c("B4E0K5"))
#' }
annotateUniprotAccession <- function(accessions,
                                     orthologue=FALSE,
                                     multiOrth=FALSE,
                                     backend = NULL) {
  backend <- normalizeAnnotationBackend(backend)
  validIDs <- accessions[!is.na(accessions)]
  
  Accession <- EntryName <- GeneID <- TaxID <- Ensembl <- NULL
  EnsemblGeneID <- NULL
  
  if (backend == "mongodb") {
    giCon <- connectMongoDB(instance="bioinfo_read",
                            collection="uniprot")

    fieldsJson <- returnFieldsJson(c("UniProtKB-AC",
                                     "UniProtKB-ID", "geneID", "NCBI-taxon", "Ensembl"))
    query <- paste0('{"UniProtKB-AC":{"$in":[',
                    paste(paste0("\"", as.character(validIDs), "\""),
                          collapse=","),']}}')
    uniprots <- giCon$find(query, fields=fieldsJson)
  } else {
    uniprots <- biocAnnotateUniProt(validIDs) %>%
      dplyr::rename("UniProtKB-AC" = "Accession",
                    "UniProtKB-ID" = "EntryName",
                    "geneID" = "GeneID",
                    "NCBI-taxon" = "TaxID",
                    "Ensembl" = "EnsemblGeneID")
  }

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
    res <- appendHumanOrthologsWithNCBI(featAnno,
                                        multiOrth = multiOrth,
                                        backend = backend)
  } else {
    geneanno <- annotateGeneIDsWithoutHumanOrtholog(featAnno$GeneID,
                                                    backend = backend) %>%
      dplyr::select(-TaxID)
    res <- left_join(featAnno, geneanno, by="GeneID")
  }
  res <- sortAnnotationByQuery(res, accessions, "Accession", multi=FALSE)
  return(res)
}
