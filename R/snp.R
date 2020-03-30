#' Get HAPMAP SNP annotation and chromosomal location
#' 
#' \code{hapmapSnp} queries database with a vector of input HAPMAP SNP
#' identifiers, in the form of \code{rsXXXXX}, where \code{XXXXX} are integer
#' numbers, and returns a data.frame containing
#' 
#' The \code{hapmapSnp} function returns SNP annotations from the HapMap
#' project. When \code{genes} is set to \code{TRUE}, the function returns the
#' gene in which SNP resides if it exists, or the closest genes and the
#' distances from SNP to them. A side effect of using this information is that
#' the output \code{data.frame} usually does not have the same row number as
#' the length of the input vector, since it is quite common that a SNP is
#' annotated with two nearest genes (up- and down-stream). User may want to
#' check their results for these cases manually.
#' 
#' @param ids A character vector, containing SNP identifiers
#' @param genes Logical, whether the information of genes residing on or around
#' the SNP should be returned
#' @param flanking Logical, whether the flanking DNA sequences should be
#' returned
#' @return A \code{data.frame} object containing at least following columns:
#' \item{SNP_ID}{SNP ID in the form of rsXXXXX} \item{Chromosome}{On which
#' chromosome (currently, Feb 2012, we are using the hg19 coordinates)}
#' \item{Position}{Chromosome coordinate of the SNP} \item{Allele1}{Most
#' frequent Allele} \item{Allele2}{Secondly frequent Allele}
#' 
#' An extra column named \code{FlankingSeq} is available when \code{flanking}
#' is set to \code{TRUE}, providing flanking sequences.
#' 
#' When \code{genes=TRUE}, following columns are available
#' \item{ClosestGeneID}{Entrez GeneID of the closest gene. Note that when
#' DistanceToGene is 0, this is the GeneID of the gene where the SNP resides}
#' \item{DistanceToGene}{Distance to the closest gene (in bp). SNP in a gene is
#' annotated as 0.} \item{GeneStart}{Where does the gene start}
#' \item{GeneStop}{Where does the gene stops}
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @references \url{http://bioinfo.bas.roche.com:8080/apps/hapmap/}
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' snpIds <- c("rs884080", "rs557477", "rs531099","rs763318")
#' hapmapSnp(snpIds)
#' hapmapSnp(snpIds, flanking=TRUE, genes=TRUE)
#' hapmapSnp(snpIds, flanking=TRUE, genes=FALSE)
#' hapmapSnp(snpIds, flanking=FALSE, genes=TRUE)
#' 
#' options(error=NULL)
#' 
#' @export hapmapSnp
hapmapSnp <- function(ids,
                      genes=FALSE,
                      flanking=FALSE) {
  ## note that the table should be SNP_HAPMAP2
  state.format <- "SELECT %s FROM genome.SNP_HAPMAP2 a "
  state.sel <- c("a.SNP_ID", "a.SEQ", "a.POSITION", "a.ALLELE1", "a.ALLELE2")
  cnames <- c("SNP_ID","Chromosome",  "Position", "Allele1", "Allele2")
  if(genes) {
    state.format <- paste(state.format,
                          ", genome.snp_hapmap_gene b WHERE a.SNP_ID(+)=b.SNP_ID AND ")
    state.sel <- c(state.sel,
                   "b.RO_GENE_ID", "b.DISTANCE", "b.GENE_START", "b.GENE_STOP")
    cnames <- c(cnames, c("ClosestGeneID", "DistanceToGene", "GeneStart", "GeneStop"))
  } else {
    state.format <- paste(state.format, "WHERE ")
  }

  if(flanking) {
    state.sel <- c(state.sel, "FLANKING")
    cnames <- c(cnames, "FlankingSeq")
  }
  state <- sprintf(state.format,
                   paste(state.sel, collapse=","))
  ann <- querydbSelectIn(state,
                         "a.SNP_ID", ids,
                         db=dbName(), user=binUser(), password=binPwd())
  colnames(ann) <- cnames
  ann$Chromosome <- gsub("^CHR", "", ann$Chromosome)
  
  res <- matchColumn(ids, ann, "SNP_ID", multi=genes)
  rownames(res) <- NULL
  return(res)
}

