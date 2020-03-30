## translate GTI chiptype to bioc chiptype
#' @export bioc2gti
bioc2gti <- function (chipname) {
  data("gtibioc", package="ribiosAnnotation")
  if(missing(chipname)) {

    supp <- !is.na(gtibioc$Bioconductor)
    vec <- gtibioc$GTI[supp]
    names(vec) <- gtibioc$Bioconductor[supp]
    return(vec)
  } else {
    matchColumn(chipname, gtibioc, "Bioconductor")$GTI
  }
}



#' Translate chiptypes between GTI Bioconductor naming conventions
#' 
#' \code{gti2bioc} converts chip types from GTI array names into Bioconductor
#' names, and \code{bioc2gti} converts Bioconductor array names to GTI names.
#' If the array name is not valid or not found, NA will be returned.
#' 
#' The translation table \code{gtibioc} was compiled manually in December 2011.
#' 
#' When the parameter \sQuote{chipname} is missing, chip types supported by
#' both GTI and Bioconductor will be printed: \code{gti2bioc} returns a
#' character vector of the Bioconductor names, and \code{bioc2gti} returns such
#' a vector of the GTI names. Both vectors have the chip types in the other
#' system as names. See examples.
#' 
#' @aliases gtibioc gti2bioc bioc2gti
#' @param chipname Character vector, chip names (types). If missing, chip types
#' supported by both GTI and Bioconductor will be printed, see details.
#' @return Chracter vector of the same length as the input
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @examples
#' 
#' bioc2gti("hgu133plus2")
#' bioc2gti(c("hgu133plus2", "hgu95av2", "bad_array"))
#' gti2bioc("HG_U95AV2")
#' gti2bioc(c("HG_U95AV2", "CANINE", "HG_U95A"))
#' 
#' ## supporting empty option
#' bioc2gti()
#' gti2bioc()
#' 
#' @export gti2bioc
gti2bioc <- function(chipname) {
  data("gtibioc", package="ribiosAnnotation")
  if(missing(chipname)) {
    supp <- !is.na(gtibioc$Bioconductor)
    vec <- gtibioc$Bioconductor[supp]
    names(vec) <- gtibioc$GTI[supp]
    return(vec)
  } else {
    matchColumn(chipname, gtibioc, "GTI")$Bioconductor
  }
}
