#' @include removeEnsemblVersion.R
NULL

#' Annotate Enesembl GeneIDs
#' 
#' @param ids A vector of EnsemblGeneIDs in form of 
#' \code{ENS(species)(object type)(identifier).(version)}. The version is
#' optional.
#' @param orthologue Logical, whether human orthologues should be returned. 
#' Default: \code{FALSE}
#' @param multiOrth Logical, whether mutliple orthologues should be returned if
#' exist. Deafult: \code{FALSE}
#' @return A \code{data.frame} object containing the annotations:
#' * GeneID EntrezGeneID
#' * GeneSymbol Official gene symbol
#' * Description Gene description
#' * TaxID Taxonomy ID
#' * Type Gene type
#' 
#' If \code{orthologue} is \code{TRUE}, following columns are appended:
#' * HumanGeneID
#' * HumanGeneSymbol
#' * HumanDescription
#' * HumanType
#' @seealso \code{\link{annotateEnsemblGeneIDsWithoutHumanOrtholog}} and
#' \code{\link{annotateEnsemblGeneIDsWithHumanOrtholog}}
#' @examples 
#' \dontrun{
#'   annotateEnsemblGeneIDs(ids=c("ENSG00000236453", "ENSG00000170782",
#'                                "ENSG00000187867"))
#'   annotateEnsemblGeneIDs(ids=c("ENSG00000236453", "ENSG00000170782",
#'                                "ENSG00000187867", NA), orthologue=TRUE)
#'   annotateEnsemblGeneIDs(ids=c("ENSG00000174827", "ENSMUSG00000038298",
#'                                "ENSG00000198483", "ENSMUSG00000038354",
#'                                "ENSRNOG00000054947", "ENSG00000278099"),
#'                         orthologue=TRUE)
#' }
#' @export
annotateEnsemblGeneIDs <- function(ids, orthologue=FALSE, multiOrth=FALSE) {
  if(!orthologue) {
    res <- annotateEnsemblGeneIDsWithoutHumanOrtholog(ids)
  } else {
    res <- annotateEnsemblGeneIDsWithHumanOrtholog(ids, multiOrth=TRUE)
  }
  return(res)
}

#' Annotate EnsEMBL GeneID with data from EnsEMBL
#' @param ids Character strings, Ensembl GeneIDs in form of 
#' \code{ENS(species)(object type)(identifier).(version)}. The version is
#' optional.
#' @return A \code{data.frame} containing following columns:
#' \itemize{
#'   \item{EnsemblID}: The input EnsemblID
#'   \item{GeneID}: NCBI GeneID
#'   \item{GeneSymbol}: Official gene symbol
#'   \item{Description}: Gene description
#'   \item{TaxID}: Taxonomy ID
#' }
#' This function uses data from EnsEMBL to annotate EnsEMBL GeneIDs. For most 
#' users, it is recommended to use \code{\link{annotateEnsemblGeneIDs}}, 
#' because it uses both data from EnsEMBL and data from NCBI to perform the
#' task.
#' 
#' @details The \code{ensembl_genes} collection is used. Note that Ensembl 
#' IDs often refer to novel transcripts which do not have identifiers in other
#' databases like NCBI Genes. If an EnsemblID is invalid or obsolete, the fields
#' \code{GeneName} and \code{TaxID} will be NA.
#' 
#' @seealso Function \code{\link{annotateEnsemblGeneIDsWithNCBI}} annotates
#' EnsEMBL GeneIDs with data from NCBI, and \code{\link{annotateEnsemblGeneIDs}}
#' annotates EnsEMBL GeneIDs with both data from EnsEMBL and data from NCBI.
#' @importFrom magrittr %>%
#' @importFrom ribiosUtils matchColumnIndex
#' @importFrom dplyr left_join
#' @examples
#' \dontrun{
#' ensIDs <- readLines(system.file(file.path("extdata/ribios_annotate_testdata",
#'                     "ensemble_geneids.txt"), package="ribiosAnnotation"))
#' ensAnno <- annotateEnsemblGeneIDsWithEnsembl(ensIDs)
#' }
#' @export
annotateEnsemblGeneIDsWithEnsembl <- function(ids) {
  ids <- as.character(ids)
  uvids <- removeEnsemblVersion(ids)
  
  input <- data.frame(EnsemblID=ids,
                      UVID=uvids)
  EnsemblID <- GeneSymbol <-  GeneID <- Description <- TaxID <- NULL
  
  giCon <- giCon <- connectMongoDB(instance="bioinfo_read",
                                   collection='ensembl_genes')
  
  speciesFields <- c("symbol", "description", "geneId", "entrezgeneId",
                     "taxId")
  speciesFieldsJson <- returnFieldsJson(speciesFields)
  query <- paste0('{"geneId":{"$in":[', 
                  paste("\"", 
                        as.character(uvids), 
                        "\"", collapse=",", sep=""),']}}')
  genes <- giCon$find(query, fields=speciesFieldsJson)
  if(nrow(genes)==0) {
    warning("No annotation was found")
    res <- data.frame(EnsemblID=ids,
                      GeneID=NA,
                      GeneSymbol=NA,
                      Description=NA,
                      TaxID=NA)
  } else {
    res <- genes %>%
      dplyr::rename('UVID'='geneId',
                    'GeneID'='entrezgeneId',
                    'GeneSymbol'='symbol',
                    'Description'='description',
                    "TaxID"="taxId") %>%
      dplyr::left_join(input, by="UVID") %>%
      dplyr::select(EnsemblID, GeneID, GeneSymbol, Description, TaxID)
      resInd <- ribiosUtils::matchColumnIndex(ids, res, "EnsemblID")
      res <- res[resInd, , drop=FALSE]
      res$EnsemblID <- ids
  }
  rownames(res) <- id2rownames(ids)
  return(res)
}

#' Annotate EnsEMBL GeneID with data from NCBI
#' @param ids Character strings, Ensembl GeneIDs in form of 
#' \code{ENS(species)(object type)(identifier).(version)}. The version is
#' optional.
#' @return A \code{data.frame} containing following columns:
#' \itemize{
#'   \item{EnsemblID}: The input EnsemblID
#'   \item{GeneID}: NCBI GeneID
#'   \item{GeneSymbol}: Official gene symbol
#'   \item{Description}: Gene description
#'   \item{TaxID}: Taxonomy ID
#'   \item{Type}: Gene type
#' }
#' This function uses data from NCBI to annotate EnsEMBL GeneIDs. For most 
#' users, it is recommended to use \code{\link{annotateEnsemblGeneIDs}}, 
#' because it uses both data from EnsEMBL and data from NCBI to perform the
#' task.
#' 
#' @details The \code{ncbi_gene2ensembl} collection is used.
#' @seealso Function \code{\link{annotateEnsemblGeneIDsWithEnsembl}} annotates
#' EnsEMBL GeneIDs with data from Ensembl, and 
#' \code{\link{annotateEnsemblGeneIDs}} annotates EnsEMBL GeneIDs with both 
#' data from EnsEMBL and data from NCBI.
#' @importFrom magrittr %>%
#' @importFrom ribiosUtils matchColumnIndex
#' @importFrom dplyr left_join
#' @examples
#' \dontrun{
#' ensIDs <- readLines(system.file(file.path("extdata/ribios_annotate_testdata",
#'                     "ensemble_geneids.txt"), package="ribiosAnnotation"))
#' ncbiAnno <- annotateEnsemblGeneIDsWithNCBI(ensIDs)
#' }
#' @export
annotateEnsemblGeneIDsWithNCBI <- function(ids) {
  ids <- as.character(ids)
  uvids <- removeEnsemblVersion(ids)
  input <- data.frame(EnsemblID=ids,
                      UVID=uvids)
  
  EnsemblID <- GeneID <- NULL
  
  giCon <- connectMongoDB(instance="bioinfo_read",
                          collection='ncbi_gene2ensembl')
  
  queryFields <- c("geneId", "Ensembl_geneId")
  queryFieldsJson <- returnFieldsJson(queryFields)
  query <- paste0('{"Ensembl_geneId":{"$in":[', 
                  paste("\"", as.character(uvids), "\"", collapse=",", sep=""),']}}')
  genes <- giCon$find(query, fields=queryFieldsJson) 
  if(nrow(genes)==0) {
    warning("No annotation was found")
    res <- data.frame(EnsemblID=ids,
                      GeneID=NA,
                      GeneSymbol=NA,
                      Description=NA,
                      TaxID=NA,
                      Type=NA)
  } else {
    resE2N <- genes %>%
      dplyr::rename('UVID'='Ensembl_geneId',
                    'GeneID'='geneId') %>%
      dplyr::left_join(input, by="UVID") %>%
      dplyr::select(EnsemblID, GeneID)
    resAnno <- annotateGeneIDs(ids=resE2N$GeneID)
    res <- dplyr::left_join(resE2N, resAnno, by="GeneID")
    resInd <- ribiosUtils::matchColumnIndex(ids, res, "EnsemblID")
    res <- res[resInd, , drop=FALSE]
    res$EnsemblID <- ids
  }
  rownames(res) <- id2rownames(ids)
  return(res)
}

#' Annotate Ensembl GeneIDs with data from both EnsEMBL and NCBI
#' @param ids A vector of character strings, Ensembl GeneIDs in form of 
#' \code{ENS(species)(object type)(identifier).(version)}. The version is
#' optional.
#' @return A \code{data.frame} containing following columns:
#' \itemize{
#'   \item{EnsemblID}: The input EnsemblID
#'   \item{GeneID}: NCBI GeneID
#'   \item{GeneSymbol}: Official gene symbol
#'   \item{Description}: Gene description
#'   \item{TaxID}: Taxonomy ID
#'   \item{Type}: Gene type
#' }
#' @details First, both EnsEMBL and NCBI annotation is queried. Next, we use
#' the NCBI annotation as the template. Finally, we take the EnsEMBL annotation 
#' for those genes that are annotated by EnsEMBL but not by NCBI, merging the
#' information from both sources.
#' @examples
#' \dontrun{
#' ensIDs <- readLines(system.file(file.path("extdata/ribios_annotate_testdata",
#'                     "ensemble_geneids.txt"), package="ribiosAnnotation"))
#' enAnno <- annotateEnsemblGeneIDsWithoutHumanOrtholog(ensIDs)
#' }
#' @export
annotateEnsemblGeneIDsWithoutHumanOrtholog <- function(ids) {
  ensAnno <- annotateEnsemblGeneIDsWithEnsembl(ids)
  ncbiAnno <- annotateEnsemblGeneIDsWithNCBI(ids)

  res <- ncbiAnno
  isEnsSuccess <- !is.na(ensAnno$TaxID)
  isNcbiFailure <- is.na(ncbiAnno$GeneID)
  isToReplace <- isEnsSuccess & isNcbiFailure
  replaceCols <- c("GeneID", "GeneSymbol", "Description", "TaxID")
  for (column in replaceCols) {
    res[isToReplace, column] <- ensAnno[isToReplace, column]
  }
  return(res)
}

#' Annotate Ensembl GeneIDs while appending human orthologs
#' @param ids A vector of character strings, Ensembl GeneIDs in form of 
#' \code{ENS(species)(object type)(identifier).(version)}. The version is
#' optional.
#' @return A \code{data.frame} containing following columns:
#' \itemize{
#'   \item{EnsemblID}: The input EnsemblID
#'   \item{GeneID}: NCBI GeneID
#'   \item{GeneSymbol}: Official gene symbol
#'   \item{Description}: Gene description
#'   \item{TaxID}: Taxonomy ID
#'   \item{Type}: Gene type
#'   \item{HumanGeneID}: NCBI GeneID of the human orthologue
#'   \item{HumanGeneSymbol}: Official gene symbol of the human orthologue
#'   \item{HumanDescription}: Gene description of the human orthologue
#'   \item{HumanType}: Gene type of the human orthologue
#' }
#' @note Currently the human orthologs are looked up in NCBI. It remains to
#' be changed to EnsEMBL
#' @examples
#' \dontrun{
#' ensIDs <- readLines(system.file(file.path("extdata/ribios_annotate_testdata",
#'                     "ensemble_geneids.txt"), package="ribiosAnnotation"))
#' enAnnoHumanOrt <- annotateEnsemblGeneIDsWithHumanOrtholog(ensIDs)
#' }
#' @export
annotateEnsemblGeneIDsWithHumanOrtholog <- function(ids, multiOrth=FALSE) {
  HumanGeneID <- HumanGeneSymbol <- HumanDescription <- Type <- NULL
  GeneSymbol <- GeneID <- Description <- NULL
  
  ensAnno <- annotateEnsemblGeneIDsWithoutHumanOrtholog(ids)
  appEnsAnno <- appendHumanOrthologsWithNCBI(ensAnno, multiOrth = multiOrth)
  res <- sortAnnotationByQuery(appEnsAnno, ids, id_column="EnsemblID")
  return(res)
}