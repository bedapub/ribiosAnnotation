#' Remove version suffix from Ensembl IDs
#' @param ensemblIDs A vector of character strings. Other types of inputs are 
#'   converted.
#' @return A character vector of the same length as input
#' 
#' @examples 
#' ensemblIDs <- c("ENSG00000197535", "ENST00000399231.7", "ENSP00000418960.2")
#' removeEnsemblVersion(ensemblIDs)
#' @export
removeEnsemblVersion <- function(ensemblIDs) {
  return(gsub("\\.[0-9]+$", "", as.character(ensemblIDs)))
}