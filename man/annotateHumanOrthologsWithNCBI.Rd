% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateHumanOrthologsWithNCBI.R
\name{annotateHumanOrthologsWithNCBI}
\alias{annotateHumanOrthologsWithNCBI}
\title{Annotate human orthologs with data from NCBI}
\usage{
annotateHumanOrthologsWithNCBI(geneids, multiOrth = FALSE)
}
\arguments{
\item{geneids}{Integer GeneIDs, can contain human GeneIDs.}

\item{multiOrth}{Logical, whether one gene is allowed to map to multiple
human orthologs? Default value is \code{FALSE}, i.e. only the first 
human ortholog (random choice) is returned.}
}
\value{
A \code{data.frame} containing following columns:
* GeneID: Input GeneID
* TaxID: Taxonomy ID of the input gene
* HumanGeneID: Human Entrez GeneID

This function annotates human orthologs for any GeneID, including human
genes, in which case the ortholog will be itself. Use 
\code{\link{annotateNonHumanGenesHumanOrthologsWithNCBI}} if you are sure
that input GeneIDs do not come from human.
}
\description{
Annotate human orthologs with data from NCBI
}
\details{
\code{ncbi_gene_info} and \code{ncbi_gene_orthologs} collections 
are used.
}
\examples{
\dontrun{
annotateHumanOrthologsWithNCBI(c(25120, 114483548, 57300, 20))
annotateHumanOrthologsWithNCBI(c(25120, 114483548, 57300, 20, 1506, 1545,
                                 102129055))
}
}
\seealso{
\code{\link{annotateNonHumanGenesHumanOrthologsWithNCBI}}
}
