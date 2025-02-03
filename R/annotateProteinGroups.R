#' @include annotateAnyIDs.R
NULL

getGenesPerIndex <- function(x) {
  if(all(is.na(x$GeneID))) {
    idx <- rep(TRUE, nrow(x))
  } else {
    idx <- which(!is.na(x$GeneID))
  }
  return(x[idx,,drop=FALSE])
}

#' Annotate protein groups for proteomics studies
#' @param ids Character, Protein groups with identifiers 
#'            (e.g. UniProt/SwissProt IDs) by the delimiter
#' @param delimiter Character, delimiter, default semicolon
#' @param orthologue Logical, whether human orthologues should be queried.
#' @param multiOrth Logical, in case of multiple orthologues, whether all of 
#' them should be returned
#' The function queries proteins in protein groups, and annotate
#' all proteins that cannot be annotated. For protein groups in which
#' no protein can be annotated, all proteins will be returned as they are,
#' without annotation
#' @returns A \code{data.frame} with following columns:
#' * ProteinGroup
#' * Protein
#' * GeneID
#' * GeneSymbol
#' * GeneName
#' * TaxID
#' In case \code{orthologue} is \code{TRUE}, human orthologue information is 
#' returned as well.
#' @examples
#' options(error=utils::recover)
#' \dontrun{
#' annotateProteinGroups(c("A0A024RBG1;Q9NZJ9", "A0A0B4J2D5;P0DPI2",
#'                            "A0A0B4J2F0;A0A0U1RRL7"))
#' }
#' options(error=NULL)
#' @importFrom dplyr group_by arrange group_modify ungroup
#' @export
annotateProteinGroups <- function(ids, delimiter=";",
                                  orthologue=FALSE, multiOrth=FALSE) {
  
  GeneID <- NULL
  GeneSymbol <- InputIDType <- index <- NULL
  idsplit <- strsplit(ids, delimiter)
  idcount <- sapply(idsplit, length)
  pg <- data.frame(index=rep(seq(along=ids), idcount),
                   ProteinGroup=rep(ids, idcount),
                   Protein=unlist(idsplit))
  pgAnno <- ribiosAnnotation::annotateAnyIDs(pg$Protein,
                                             orthologue=orthologue,
                                             multiOrth=multiOrth)
  featAnno <- left_join(pg, pgAnno, by=c("Protein"="Input"),
                        relationship="many-to-many")
  res <- featAnno %>% group_by(index) %>% arrange(GeneID) %>%
    group_modify(~getGenesPerIndex(.x)) %>%
    arrange(index)
  res <- ungroup(res) %>% 
    dplyr::select(-index) %>% 
    unique
  return(res)
}

## assignInNamespace("getGenesPerIndex", getGenesPerIndex, "ribiosAnnotation")
## assignInNamespace("annotateProteinGroups", annotateProteinGroups, "ribiosAnnotation")
