#' @include removeEnsemblVersion.R
NULL

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
#' enAnno <- annotateEnsemblGeneIDs(ensIDs)
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
  res$Source <- ifelse(isToReplace, "EnsEMBL", "NCBI")
  return(res)
}

annotateEnsemblGeneIDsWithHumanOrtholog <- function(ids) {
  HumanGeneID <- HumanGeneSymbol <- HumanDescription <- Type <- NULL
  GeneSymbol <- GeneID <- Description <- NULL
  
  ensAnno <- annotateEnsemblGeneIDsWithoutHumanOrtholog(ids)
  if(!any(is.na(geneIdAnno$TaxID)) && all(geneIdAnno$TaxID==9606)) {
    res <- geneIdAnno %>%
      dplyr::mutate(HumanGeneID=GeneID,
                    HumanGeneSymbol=GeneSymbol,
                    HumanDescription=Description,
                    HumanType=Type)
  } else {
    orthologs <- annotateHumanOrthologsWithNCBI(geneIdAnno$GeneID,
                                                multiOrth=multiOrth) %>%
      dplyr::select(GeneID, HumanGeneID)
    orthologAnno <- annotateGeneIDsWithoutHumanOrtholog(orthologs$HumanGeneID) %>%
      dplyr::select(GeneID, HumanGeneSymbol=GeneSymbol, HumanDescription=Description,
                    HumanType=Type)
    if(hasCharOrIsNa) {
      orthologAnno$GeneID <- as.character(orthologAnno$GeneID)
    }
    res <- left_join(geneIdAnno, orthologs, by="GeneID") %>%
      left_join(orthologAnno, by="GeneID") %>% unique
  }
  res <- sortAnnotationByQuery(res, ids, "GeneID", multi = multiOrth)
  return(res)
}

#' Annotate EntrezGeneID with Ensemble identifiers
#' 
#' 
#' @param ids Character or integer vector, Entrez Gene IDs. It can contain
#'  \code{NA} or \code{NULL}.
#' @param orthologue Logical, whether human orthologues should be returnd
#' @param multiOrth Logical, whether multiple mapped orthologues should be
#' returned or not. Only useful when \code{orthologue} is \code{TRUE}, By
#' default \code{FALSE}. Be cautious with the \code{TRUE} option: in this case
#' the returning \code{data.frame} may have different row numbers as the input
#' vector.
#' @param type Character string, one of the following options: \code{gene}, 
#' \code{transcript}, and \code{protein}.
#' @return A \code{data.frame} object containing the annotations.
#' 
#' If \code{orthologue} is set to \code{FALSE}, the data frame contains
#' following columns: \code{GeneID, ENSEMBL_ID, GeneSymbol, GeneName and TaxID}
#' 
#' If \code{orthologue} is \code{TRUE}, the data frame contains
#' \code{GeneID, ENSEMBL_ID, GeneSymbol,TaxID, OrigTaxID, OrigGeneID,
#' OrigGeneSymbol, OrigGeneName}. Note that \code{GeneID, GeneSymbol, TaxID}
#' contains the information of mapped orthologues, while \code{OrigXXX}
#' contains the information of queried genes.
#' 
#' If \code{multiOrth} is \code{TRUE} (only valid when \code{orthologue} is
#' \code{TRUE}), multiples orthologues mapped to the same gene are all
#' returned. This means that the result data frame may contain more rows than
#' the length of input vector. If set to \code{FALSE}, the resulting data frame
#' contains exact number of rows as the input vector length.
#'
#' If either \code{transcript} or \code{protein} type is used, more than one 
#' items will be returned for each input ID.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso See \code{\link{annotateEnsembl}} to annotate Ensembl IDs.
#' 
#' See \code{\link{annotateGeneIDs}} to get annotation for Entrez GeneIDs.
#' 
#' See \code{\link{annotateProbesets}} to get annotation for probesets.
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' ## normal use
#' annotateGeneIDsWithEnsembl(c(1234, 1432, NULL))
#' annotateGeneIDsWithEnsembl(c(1234, 1432, NULL), type="transcript")
#' annotateGeneIDsWithEnsembl(c(1234, 1432, NULL), type="protein")
#' 
#' annotateGeneIDsWithEnsembl(ids="22227", orthologue=TRUE)
#' ## it is possible to mix ids of different species
#' annotateGeneIDsWithEnsembl(ids=c("1234", "22227"), orthologue=TRUE)
#' 
#' options(error=NULL)
#' 
#' @export annotateGeneIDsWithEnsembl
annotateGeneIDsWithEnsembl <- function (ids, 
                                        orthologue = FALSE, 
                                        multiOrth = FALSE,
                                        type=c("gene", "transcript", "protein")) {
  .Deprecated()
  type <- match.arg(type)
  typeID <- c("gene"="26", "transcript"="30", "protein"="31")[type]
  comm <- paste("SELECT c.RO_GENE_ID, e.ENSEMBL_ID, c.GENE_SYMBOL, c.DESCRIPTION, c.TAX_ID", 
                " FROM GTI_GENES c INNER JOIN ENSEMBL_GENE e ON c.RO_GENE_ID=e.GENE_ID ", 
                paste0("WHERE e.TYPE_ID=", typeID), " ",
                sep = "")
  ann <- querydbTmpTbl(comm, "c.RO_GENE_ID", ids, dbName(),
                       binUser(), 
                       binPwd())
  keyCol <- ifelse(orthologue, "OrigGeneID", "GeneID")
  cnames <- c(keyCol, "EnsemblID", "GeneSymbol", "GeneName", "TaxID")
  conames <- c(keyCol, "EnsemblID", "OrigGeneSymbol", "OrigGeneName", 
               "OrigTaxID")
  if (!orthologue) {
    colnames(ann) <- cnames
    res <- ann
  } else {
    colnames(ann) <- conames
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, 
                                           multiOrth = multiOrth)
    if (multiOrth) {
      res <- merge(ann, ort, by = "OrigGeneID", all.x = TRUE)
    }  else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", 
                            multi = FALSE)
      res <- cbind(ann, ort.re[, -1L])
    }
    res <- putColsFirst(res, c("OrigGeneID", "GeneID", "GeneSymbol", "TaxID", 
                               "OrigTaxID", "OrigGeneSymbol", "OrigGeneName"))
  }
  isTranscriptOrProtein <- type!="gene"
  res <- matchColumn(ids, res, keyCol, multi = (orthologue && multiOrth) || isTranscriptOrProtein)
  rownames(res) <- id2rownames(res[,keyCol])
  return(res)
}
