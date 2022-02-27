#' @include utils.R
NULL

#' Get chromosome lengths in basepairs
#' 
#' 
#' @param organism Organism to query, supporting human, mouse and rat
#' @return A \code{data.frame} containing chromosome length information.
#' \item{Chromosome}{Chromosome names (without the \code{CHR} prefix as in the
#' database} \item{Description}{Description} \item{Length}{Length of the
#' chromosome, in basepairs (bp)} \item{Rank}{Rank of chromosomes by lengths}
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>, Isabelle Wells
#' <isabelle.wells@@roche.com>
#' @examples
#' 
#' \dontrun{
#' chrLen("human")
#' chrLen("mouse")
#' }
#' 
#' @export chrLen
chrLen <- function(organism=c("human", "mouse", "rat")) {
  organism <- match.arg(organism)
  if(organism=="human") {
    db <- "HUMANN"
  } else if(organism=="mouse") {
    db <- "MOUSEN"
  } else if (organism=="rat") {
    db <- "RATN"
  }
  
  state <- paste("SELECT SEQ, DESCR, SEQLEN, RANK ",
                 "FROM genome_sequence ",
                 "WHERE DB='", db, "'", sep="")
  ann <- querydb(state, db=dbName(), user=binUser(), password=binPwd())

  colnames(ann) <- c("Chromosome", "Description", "Length", "Rank")
  ann$Chromosome <- gsub("^CHR","", ann$Chromosome)
  rownames(ann) <- NULL
  ann
}
