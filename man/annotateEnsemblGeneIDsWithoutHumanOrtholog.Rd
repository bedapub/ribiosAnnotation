% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensembl.R
\name{annotateEnsemblGeneIDsWithoutHumanOrtholog}
\alias{annotateEnsemblGeneIDsWithoutHumanOrtholog}
\title{Annotate Ensembl GeneIDs with data from both EnsEMBL and NCBI}
\usage{
annotateEnsemblGeneIDsWithoutHumanOrtholog(ids)
}
\arguments{
\item{ids}{A vector of character strings, Ensembl GeneIDs in form of 
\code{ENS(species)(object type)(identifier).(version)}. The version is
optional.}
}
\value{
A \code{data.frame} containing following columns:
\itemize{
  \item{EnsemblID}: The input EnsemblID
  \item{GeneID}: NCBI GeneID
  \item{GeneSymbol}: Official gene symbol
  \item{Description}: Gene description
  \item{TaxID}: Taxonomy ID
  \item{Type}: Gene type
}
}
\description{
Annotate Ensembl GeneIDs with data from both EnsEMBL and NCBI
}
\details{
First, both EnsEMBL and NCBI annotation is queried. Next, we use
the NCBI annotation as the template. Finally, we take the EnsEMBL annotation 
for those genes that are annotated by EnsEMBL but not by NCBI, merging the
information from both sources.
}
\examples{
\dontrun{
ensIDs <- readLines(system.file(file.path("extdata/ribios_annotate_testdata",
                    "ensemble_geneids.txt"), package="ribiosAnnotation"))
enAnno <- annotateEnsemblGeneIDsWithoutHumanOrtholog(ensIDs)
}
}
