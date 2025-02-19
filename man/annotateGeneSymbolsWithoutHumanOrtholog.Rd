% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateGeneSymbols.R
\name{annotateGeneSymbolsWithoutHumanOrtholog}
\alias{annotateGeneSymbolsWithoutHumanOrtholog}
\title{Annotate gene symbols without human ortholog}
\usage{
annotateGeneSymbolsWithoutHumanOrtholog(ids, taxId = 9606)
}
\arguments{
\item{ids}{Character vector, gene symbols to be queried}

\item{taxId}{Integer, NCBI taxonomy ID of the species. Default value: 9606 (human). See \code{commonSpecies} for tax id of common species.}
}
\value{
A data.frame of following columns
\describe{
  \item{GeneID}{Entrez Gene ID}
  \item{GeneSymbol}{Official gene symbols}
  \item{Description}{Description}
  \item{TaxID}{NCBI Taxonomy ID}
  \item{Type}{Gene type}
}
}
\description{
Annotate gene symbols without human ortholog
}
\examples{
\dontrun{
    annotateGeneSymbolsWithoutHumanOrtholog(c("AKT1", "ERBB2", 
                                            "NoSuchAGene", "TGFBR1"), 9606)
    annotateGeneSymbolsWithoutHumanOrtholog(c("Akt1", "Erbb2", 
                                            "NoSuchAGene", "Tgfbr1"), 10090)
}
}
