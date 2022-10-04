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
  ribiosUtils::haltifnot(all(c("GeneID", "TaxID") %in% colnames(anno)),
                         msg=paste0("The annotation dataframe must contain ",
                         "'GeneID' and 'TaxID' columns."))
  
  anno$chrGeneID <- as.character(anno$GeneID)
  if(!any(is.na(anno$TaxID)) && all(anno$TaxID==9606)) {
    res <- anno %>%
      dplyr::mutate(HumanGeneID=GeneID,
                    HumanGeneSymbol=GeneSymbol,
                    HumanDescription=Description,
                    HumanType=Type)
  } else {
    orthologs <- annotateHumanOrthologsWithNCBI(anno$GeneID,
                                                multiOrth=multiOrth) %>%
      dplyr::mutate(chrGeneID=as.character(GeneID)) %>%
      dplyr::select(chrGeneID, HumanGeneID)
    orthologAnno <- annotateGeneIDsWithoutHumanOrtholog(orthologs$HumanGeneID) %>%
      dplyr::mutate(chrGeneID=as.character(GeneID)) %>%
      dplyr::select(chrGeneID, 
                    HumanGeneSymbol=GeneSymbol,
                    HumanDescription=Description,
                    HumanType=Type)
    res <- left_join(anno, orthologs, by="chrGeneID") %>%
      left_join(orthologAnno, by="chrGeneID") %>% 
      unique %>%
      dplyr::select(-chrGeneID)
  }
  
  res <- sortAnnotationByQuery(res, anno$GeneID, "GeneID", multi = multiOrth)
  return(res)
}
