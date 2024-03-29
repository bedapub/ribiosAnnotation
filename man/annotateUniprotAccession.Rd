% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateUniprotAccession.R
\name{annotateUniprotAccession}
\alias{annotateUniprotAccession}
\title{Annotate UniProt accessions or names}
\usage{
annotateUniprotAccession(accessions, orthologue = FALSE, multiOrth = FALSE)
}
\arguments{
\item{accessions}{Character strings, UniProt accessions or names}

\item{orthologue}{Logical, whether orthologues are returned}

\item{multiOrth}{Logical, only valid if \code{orthologue} is \code{TRUE},
whether multiple orthologues are returned instead of only one.}
}
\description{
Annotate UniProt accessions or names
}
\examples{
\dontrun{
annotateUniprotAccession(c("B4E0K5"))
}
}
