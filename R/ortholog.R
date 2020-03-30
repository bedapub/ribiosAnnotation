safe.as.integer <- function(x) as.integer(as.character(x))



#' Mapping between orthologous genes in human and in other taxa
#' 
#' \code{humanOrthologs} and \code{nonhumanOrthologs} are two basic functions
#' to map between orthologous genes in human and in other taxa, mainly in rat
#' and/or in mouse, using information stored in the GTI database.
#' \code{humanOrthologs} get all orthologous genes in human of a vector of
#' input genes in other taxon. \code{nonhumanOrthologs} makes the revserve
#' mapping, namely getting orthologs of human genes in other taxon.
#' 
#' \code{humanUniqOrtholog} and \code{nonhumanUniqOrtholog} are wrappers to
#' return exact one ortholog for each gene. No complex logic is implemented
#' trying to guess which ortholog is the best. Users in want of such behaviour
#' should start from \code{humanOrthologs} and \code{nonhumanOrthologs}. See
#' details below.
#' 
#' For a vector of GeneID in query, \code{humanOrthologs} and
#' \code{nonhumanOrthologs} returns a list in the same order as the input
#' vector, each list item containing a vector of Entrez GeneIDs that are
#' ortholog(s) of the input gene, or \code{NULL} if no ortholog was found.
#' 
#' \code{humanUniqOrtholog} and \code{nonhumanUniqOrtholog} are wrappers that
#' ensure each gene has has one ortholog returned. This is done by only keeping
#' the gene with the minimum GeneID. This is only a pragmatic solution, since
#' for genes with multiple Orthologs this causes loss of information. However,
#' we note that only a few genes have more than one orthologs (for example <5\%
#' in mouse), and we hypothesize using any one of them will not affect
#' high-throughput analysis, for example gene set enrichment analysis, much.
#' Choosing the minimum GeneID garantees the function runs determinsitically,
#' that means it returns the same result every time so far the underlying
#' database remains unchanged.
#' 
#' @aliases humanOrthologs humanUniqOrtholog nonhumanOrthologs
#' nonhumanUniqOrtholog
#' @param geneids A vector of integers or charaters, giving Entrez GeneIDs.
#' @param taxid NCBI Taxon ID, vector of integer or character. For
#' \code{nonhumanOrthologs}, if set to \code{NULL} (default), all orthologs
#' available will be returned. If a vector is given, only orthologs in these
#' specified taxa are returned. For \code{nonhumanUniqOrtholog}, onle one Taxon
#' ID is accepted.
#' @return \code{humanOrthologs} and \code{nonhumanOrthologs} return a list in
#' the same length and exactly the same order as the input vector of GeneIDs.
#' Each list item contains a vector of ortholog GeneIDs, or \code{NULL}. For
#' \code{nonhumanOrthologs}, if genes from more than one taxon is returned,
#' Taxon IDs are used as name of list items.
#' 
#' \code{humanUniqOrtholog} and \code{nonhumanUniqOrtholog} return a vector in
#' the same length and exactly the same order as the input vector.It contains
#' Entrez GeneIDs of orthologous genes of queried GeneIDs. For genes without
#' available orthologs, \code{NA}s are returned.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com> and Laura Badi
#' <laura.badi@@roche.com>
#' @examples
#' 
#' options(error=utils::recover)
#' humanOrthologs(c(100017, 100019, 100012, 100037258))
#' humanUniqOrtholog(c(100017, 100019, 100012, 100037258))
#' 
#' nonhumanOrthologs(c(1432, 5611, 26119))
#' nonhumanOrthologs(c(1432, 5611, 26119), taxid=10116)
#' nonhumanOrthologs(c(1432, 5611, 26119), taxid=c(10116,10090))
#' nonhumanUniqOrtholog(c(1432,5611,26119), taxid=10116)
#' options(error=NULL)
#' 
#' @export humanOrthologs
humanOrthologs <- function(geneids) {
  otbl <- querydbSelectIn("SELECT RO_GENE_ID1 as HG, RO_GENE_ID2 as NHG FROM GTI_ORTHOLOGS WHERE TAX_ID1='9606' AND ",
                          inCol="RO_GENE_ID2",
                          inValues=geneids,
                          db=dbName(), user=binUser(), password=binPwd())
  otbl <- otbl[grepl("^[0-9]*$", otbl[,"HG"]),,drop=FALSE]
  olist <- split(safe.as.integer(otbl[,"HG"]),
                 otbl[,"NHG"])
  ind <- match(geneids, names(olist))
  res <- olist[ind]
  names(res) <- geneids
  return(res)
}

#' @export humanUniqOrtholog
humanUniqOrtholog <- function(geneids) {
  allOrths <- humanOrthologs(geneids)
  sapply(allOrths, function(x) {
    if(length(x)==0) return(NA)
    min(x, na.rm=TRUE)
  })  
}

#' @export nonhumanOrthologs
nonhumanOrthologs <- function(geneids, taxid=NULL) {
  pre.state <- ifelse(is.null(taxid),
                      "SELECT RO_GENE_ID1 as HG, TAX_ID2 as TAXID, RO_GENE_ID2 as NHG FROM GTI_ORTHOLOGS WHERE TAX_ID2!='9606' AND ",
                      paste("SELECT RO_GENE_ID1 as HG, TAX_ID2 as TAXID,",
                            "RO_GENE_ID2 as NHG FROM GTI_ORTHOLOGS WHERE TAX_ID2!='9606' AND ",
                            "TAX_ID2 IN", formatIn(taxid),
                            "AND ",sep=" "))
  
  ntbl <- querydbSelectIn(pre.state,
                          inCol="RO_GENE_ID1",
                          inValues=geneids,
                          db=dbName(),user=binUser(), password=binPwd())
  ntbl <- ntbl[grepl("^[0-9]*$", ntbl[,"NHG"]),,drop=FALSE]
  nhgs <- safe.as.integer(ntbl[,"NHG"])
  if(is.null(taxid) || length(taxid)>1)
    names(nhgs) <- ntbl[,"TAXID"]
  nlist <- split(nhgs, ntbl[,"HG"])
  ind <- match(geneids, names(nlist))
  res <- nlist[ind]
  names(res) <- geneids
  return(res)
}

#' @export nonhumanUniqOrtholog
nonhumanUniqOrtholog <- function(geneids, taxid) {
  if(missing(taxid) || length(taxid)>1) {
    stop("Exact one TaxID is required to find unique orthologous genes for human genes")
  }
  allOrths <- nonhumanOrthologs(geneids, taxid=taxid)
  sapply(allOrths, function(x) {
    if(length(x)==0) return(NA)
    unname(min(x, na.rm=TRUE))
  })  
}
