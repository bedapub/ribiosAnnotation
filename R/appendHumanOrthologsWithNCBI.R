#' @include annotateHumanOrthologsWithNCBI.R
NULL

#' Append human orthologs to an existing annotation dataframe
#' 
#' @param anno A \code{data.frame} containing at least two columns \code{GeneID}
#'  and \code{TaxID}. The column \code{GeneID} stores Entrez GeneIDs, which are
#'  integers (NA values and characters are tolerated). The column \code{TaxID}
#'  stores taxonomy ID, which are integers, too (again, NA values and characters
#'  are tolerated).
#' @param multiOrth Logical, whether one row is allowed to map to multiple 
#' orthologues
#' 
#' The function appends human orthologs to an existing annotation data.frame. It
#' is usually called by another function. Please make sure of what you are doing 
#' if you call it directly.
#' 
#' @note The function does not sort the rows by GeneID. It is the responsibility
#'  of the calling function to do so.
#' 
#' @return A \code{data.frame} with annotation and human orthologs appended.
#' @examples
#' \dontrun{
#' anno <- data.frame(GeneID=c(780, 1506, 114483548, 102129055, NA),
#'                    TaxID=c(9606, 9606, 10116, 9541, NA))
#' appendHumanOrthologsWithNCBI(anno)
#'
#' tol_anno <- data.frame(GeneID=c(780, 1506, 114483548, 102129055, NA, "NotV"),
#'                    TaxID=c(9606, 9606, 10116, 9541, NA, NA))
#' appendHumanOrthologsWithNCBI(tol_anno)
#' }
#' @importFrom ribiosUtils haltifnot
#' @export
appendHumanOrthologsWithNCBI <- function(anno,
                                         multiOrth=FALSE) {
  Description <- GeneID <- GeneSymbol <- Type <- HumanGeneID <- NULL
  chrGeneID <- NULL

  ribiosUtils::haltifnot(all(c("GeneID", "TaxID") %in% colnames(anno)),
                         msg=paste0("The annotation dataframe must contain ",
                         "'GeneID' and 'TaxID' columns."))
  
  anno$chrGeneID <- as.character(anno$GeneID)
  taxIsAllHumanOrNA <- with(anno,
                            setequal(as.character(TaxID[!is.na(TaxID)]),
                                     "9606"))
  if(taxIsAllHumanOrNA) {
    humanAnno <- annotateGeneIDsWithoutHumanOrtholog(unique(anno$GeneID)) %>%
      dplyr::select(HumanGeneID=GeneID,
                    HumanGeneSymbol=GeneSymbol,
                    HumanDescription=Description,
                    HumanType=Type)
    res <- left_join(anno, humanAnno,
                     by=c("GeneID"="HumanGeneID"), na_matches="never") %>%
      mutate(HumanGeneID=GeneID)
  } else {
    orthologs <- annotateHumanOrthologsWithNCBI(unique(anno$GeneID),
                                                multiOrth=multiOrth) %>%
      dplyr::mutate(chrGeneID=as.character(GeneID)) %>%
      dplyr::select(chrGeneID, HumanGeneID)
    orthologAnno <- annotateGeneIDsWithoutHumanOrtholog(unique(orthologs$HumanGeneID)) %>%
      dplyr::select(HumanGeneID=GeneID,
                    HumanGeneSymbol=GeneSymbol,
                    HumanDescription=Description,
                    HumanType=Type)
    res <- left_join(anno, orthologs, by="chrGeneID", na_matches="never") %>%
      left_join(orthologAnno, by="HumanGeneID", na_matches="never")
  }
  
  res <- unique(res) %>% dplyr::select(-chrGeneID)
  return(res)
}