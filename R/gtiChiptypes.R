#' List chip names supported by the GTI system
#' 
#' The connects to the backend Oracle database, and reads a list of chip names
#' that are supported by the GTI system. These names are returned as a vector
#' of characters by default, or in a data.frame with additional information of
#' description.
#' 
#' 
#' @aliases gtiChiptypes gtiChipnames gtiArraytypes raceChiptypes raceChipnames
#' affychipNames
#' @param include.desc Whether description should be returned as well
#' @return If \code{include.desc} set to \code{FALSE}, a character vector of
#' chip names that are supported by the GTI system.
#' 
#' Otherwise, a data.frame with chip names and several columns of descriptions,
#' including species, technology and descriptions.
#' @note From package version 1.0-15, gtiChiptypes uses
#' genome.CHIP_ARRAY_TYPES@bin table to report GTI-supported chip types.
#' Comparing previous versions, the default behavior is not changing. By
#' settign 'include.desc=TRUE', more details are returned, includeing
#' Technology and Species.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com> with inputs from
#' Laura Badi <l.badi@@roche.com>
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' gtiChiptypes()
#' gtiChiptypes(include.desc=TRUE)
#' 
#' options(error=NULL)
#' 
#' @export gtiChiptypes
gtiChiptypes <- function(include.desc=FALSE) {

  state <- "SELECT ARRAY_TYPE, TECHNOLOGY, SPECIES, DESCRIPTION FROM genome.CHIP_ARRAY_TYPES"
  df <- querydb(state, db=dbName(), user=binUser(), password=binPwd())
  if(include.desc) {
    res <- df
    colnames(res) <- c("Chiptype", "Technology", "Species", "Description")
  } else {
    res <- df[,1L]
  }
  return(res)
}



#' Logical function of whether input chiptype is supported by GTI
#' 
#' Given a character vector of chiptypes, the function returns a logical
#' (boolean ) vector of the same length, each value indicating whether the
#' input chiptype is supported by GTI. Not rarely the vector is of length one.
#' 
#' It is possible to specify exceptions and case-insensitive comparisons. See
#' details below
#' 
#' The function can be used to test whether user input is supported by GTI,
#' ensuring the probesets or features of array can be properly annotated.
#' 
#' \code{exceptions} allow an extended list to be checked against. For
#' instance, currently (Feb 2012) \code{GeneID} is not one of the official GTI
#' chiptypes. However, if the software should think \code{GeneID} is also a
#' valid input, \code{exceptions} can be set as \dQuote{GeneID}.
#' 
#' @param x Mandatory, a character vector of chip types (names) to be checked
#' @param exceptions Optional, a character vector of extended list. When an
#' input chip type is included in \code{exceptions}, the corresponding
#' returning value will be set as \code{TRUE}
#' @param ignore.case Optional,logical. If set to \code{TRUE}, the comparison
#' is performed in a case-insensitive way.
#' @return A logical vector of the same length as the input, indicating whether
#' the input types are supported by GTI (\code{TRUE}) or not (\code{FALSE}).
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso The function calls internally the \code{\link{gtiChiptypes}}
#' function to get a list of GTI-supported chiptypes.
#' @examples
#' 
#' \dontrun{
#' user.input <- c("HG-U133A", "HG-U133_PLUS_2", "HUMAN_REFSEQ-8")
#' isGtiChiptype(user.input)
#' 
#' ## exceptions
#' user.extinput <- c(user.input, "GeneID", "GeneSymbol", "RefSeq")
#' isGtiChiptype(user.extinput)
#' isGtiChiptype(user.extinput, exceptions=c("GeneID", "RefSeq"))
#' 
#' ## case-insensitive comparisons
#' isGtiChiptype(tolower(user.input))
#' isGtiChiptype(tolower(user.input), ignore.case=TRUE)
#' }
#' 
#' @export isGtiChiptype
isGtiChiptype <- function(x, exceptions, ignore.case=FALSE) {
  gct <- gtiChiptypes(include.desc=FALSE)
  if(!missing(exceptions))
    gct <- c(gct, as.character(exceptions))
  if(ignore.case) {
    x <- tolower(x);
    gct <- tolower(gct)
  }
  x %in% gct
}



#' Rscript supported chiptypes
#' 
#' A wrapper function for Rscripts to prompt supported chip types for
#' annotation
#' 
#' The functions returns a vector of \code{gtiChiptypes()}, as well as
#' \code{GeneID} and \code{GeneSymbol}, as the result. These are commonly used
#' chiptypes in Rscripts. Depending on the user input, \code{annotateAnyID} can
#' be called subsequently to annotate user-input features.
#' 
#' @return A character vector of supported chip types.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{gtiChiptypes}}, and \code{\link{annotateAnyIDs}}
#' @examples
#' 
#' options(error=utils::recover)
#' scriptChiptypes()
#' options(error=NULL)
#' 
#' @export scriptChiptypes
scriptChiptypes <- function() {
  return(c("GeneID", "GeneSymbol",gtiChiptypes()))
}

#' @export gtiChipnames
gtiChipnames <- function(...) {gtiChiptypes(...)}
#' @export gtiArraytypes
gtiArraytypes <- function(...) {gtiChiptypes(...)}

##----------------------------------------##
## decrecated
##----------------------------------------##
affychipNames <- function(...) {
  .Deprecated("gtiChiptypes",package="ribiosAnnotation")
  gtiChipnames(...) 
}
raceChiptypes <- function(...) {
  .Deprecated("gtiChiptypes",package="ribiosAnnotation")
  gtiChiptypes(...) 
}
raceChipnames <- function(...) {
  .Deprecated("gtiChipnames",package="ribiosAnnotation")
  gtiChipnames(...) 
}
