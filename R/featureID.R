#' Whether input character strings are valid feature IDs
#' @param featureIDs A vector of character strings
#' @return Logical vector of the same lenght as the input
#' 
#' Invalid feature IDs include \code{NA}, \code{"-"}, and empty string. 
#' Other features are deedm as valid.
#' 
#' @seealso \code{\link{validFeatureIDs}}
#' @examples 
#' featureIDs <- c("AMPK", "", "ACTB", "-")
#' isValidFeatureID(featureIDs)
#' @export
isValidFeatureID <- function(featureIDs) {
  if(!is.character(featureIDs))
    featureIDs <- as.character(featureIDs)
  invalid <- is.na(featureIDs) | featureIDs=="" | featureIDs=="-"
  return(!invalid)
}

#' Return valid features in a vector
#' @param featureIDs A vector of character strings
#' @return A filtered vector containing only valid feature IDs.
#' 
#' Factor input will remain factors as output, but with invalid levels dropped.
#' The output class will remain the same in case of integer or character input.
#' 
#' @seealso \code{\link{isValidFeatureID}}
#' @examples 
#' featureIDs <- c("AMPK", "", "ACTB", "-")
#' validFeatureIDs(featureIDs)
#' @export
validFeatureIDs <- function(featureIDs) {
  res <- featureIDs[isValidFeatureID(featureIDs)]
  if(is.factor(res))
    res <- droplevels(res)
  return(res)
}

#' Whether input strings look like Entrez GeneIDs
#' @param featureIDs Character strings. Input of other types are converted to them.
#' @return A logical vector of the same length as input
#' @examples
#' feats <- c("1234", "LOX", "345", "-", "", "NKX-1", "CXorf21", "Snail",
#'    "A2BC19", "P12345", "A0A023GPI8",
#'    "NM_000259", "NM_000259.3", "ENSG00000197535", "ENST00000399231.7")
#' likeGeneID(feats) 
#' likeGeneSymbol(feats)
#' likeRefSeq(feats)
#' likeEnsembl(feats)
#' likeUniProt(feats)
#' likeHumanGeneSymbol(feats)
#' @export
likeGeneID <- function(featureIDs) {
  return(grepl("^[0-9]+$", as.character(featureIDs)))
}

#' @describeIn likeGeneID tests whether input strings look like NCBI RefSeq IDs
#' @export
likeRefSeq <- function(featureIDs) {
  res <- grepl("^[N|X][M|R|G|P]_[0-9]+\\.?[0-9]*$",
               as.character(featureIDs))
  return(res)
}

#' @describeIn likeGeneID tests whether input strings look like Ensembl IDs
#' @export
likeEnsembl <- function(featureIDs) {
  res <- grepl("^ENS[T|G|P][0-9]+\\.?[0-9]*$",
               as.character(featureIDs))
  return(res)
}

#' @describeIn likeGeneID tests whether input strings look like UniProt IDs
#' @references Regular expression of UniProt accesion numbers is available at \url{https://www.uniprot.org/help/accession_numbers}. We requirea a whole-string match additionally
#' @export
likeUniProt <- function(featureIDs) {
  featureIDs <- as.character(featureIDs)
  res <- grepl("^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$",
                     featureIDs)
  return(res)
}

#' @describeIn likeGeneID tests whether input strings look like gene symbols
#' @references The HGNC guideline is available at \url{https://www.genenames.org/about/guidelines/}
#' @export
likeGeneSymbol <- function(featureIDs) {
  likeGS <- grepl("^[A-Za-z][A-Za-z0-9-]+$|^C[0-9XY]+orf[0-9]+$", as.character(featureIDs))
  isEnsembl <- likeEnsembl(featureIDs)
  isUniProt <- likeUniProt(featureIDs)
  res <- likeGS & !isEnsembl & !isUniProt
  return(res)
}

#' @describeIn likeGeneID tests whether input strings look like human 
#'   gene symbols
#' @export
likeHumanGeneSymbol <- function(featureIDs) {
  likeGS <- grepl("^[A-Z][A-Z0-9-]+$|^C[0-9XY]+orf[0-9]+$", as.character(featureIDs))
  isEnsembl <- likeEnsembl(featureIDs)
  isUniProt <- likeUniProt(featureIDs)
  res <- likeGS & !isEnsembl & !isUniProt
  return(res)
}


#' Guess the majority members of a character string look like human gene symbols
#' @param x A vector of character strings
#' @param majority A numeric value between 0 and 1, the threshold of majority 
#'   voting
#' @return A logical value
#' \code{TRUE} is only returned if at least a proportion of \code{majority} 
#' members look like human gene symbols
#' 
#' @examples
#' majorityLikeHumanGeneSymbol(c("AKT1", "AKT2", "MYOA")) # TRUE
#' majorityLikeHumanGeneSymbol(c("Akt1", "Akt2", "Myoa")) # FALSE
#' majorityLikeHumanGeneSymbol(c("AKT1", "Akt2", "MYOA"), majority=0.5) # TRUE
#' @export
majorityLikeHumanGeneSymbol <- function(x, majority=0.8) {
  isHumanGS <- likeHumanGeneSymbol(x)
  res <- mean(isHumanGS, na.rm=TRUE) >= majority
  return(res)
}

#' Guess feature ID type by majority voting
#' 
#' @param featureIDs A vector of character strings. Other input types will be 
#'   converted to character strings.
#' @param majority Numeric value between 0 and 1. If the proportion of valid
#'   feature IDs in the input matching the pattern of a certain feature type 
#'   exceeds this value, the function returns a character string representing
#'   the feature ID type.
#' @return A character string, one of the following values: 
#'   \itemize{
#'     \item \code{GeneID}
#'     \item \code{GeneSymbol}
#'     \item \code{RefSeq}
#'     \item \code{Ensembl}
#'     \item \code{UniProt}
#'     \item \code{Unknown}
#'   }. The majority voting is done in the same order
#'
#' @examples 
#' guessFeatureType(c("AKT1", "AKT2", "MAPK14"))
#' guessFeatureType(c(1,2,14,149))
#' guessFeatureType(c("NM_000259", "NM_000331"))
#' guessFeatureType(c("ENST00000613858.4", "ENST00000553916.5",
#'     "ENST00000399229.6"))
#' guessFeatureType(c("A2BC19", "P12345", "A0A023GPI8"))
#' guessFeatureType(c("CM000677.2"))
#' @export
guessFeatureType <- function(featureIDs, majority=0.5) {
  vnames <- as.character(validFeatureIDs(featureIDs))
  geneIdsLike <- likeGeneID(vnames)
  ensemblLike <- likeEnsembl(vnames)
  geneSymbolsLike <- likeGeneSymbol(vnames)
  refseqLike <- likeRefSeq(vnames)
  uniprotLike <- likeUniProt(vnames)
  
  if(mean(geneIdsLike, na.rm=TRUE)>=majority) {
    return("GeneID")
  } else if(mean(geneSymbolsLike, na.rm=TRUE)>=majority) {
    return("GeneSymbol")
  } else if(mean(refseqLike, na.rm=TRUE)>=majority) {
    return("RefSeq")
  } else if(mean(ensemblLike, na.rm=TRUE)>=majority) {
    return("Ensembl")
  } else if(mean(uniprotLike, na.rm=TRUE)>=majority) {
    return("UniProt")
  } else {
    return("Unknown")
  }
}
