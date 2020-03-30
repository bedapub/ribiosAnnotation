#' Annotations of all genes associated with the given TaxID
#' 
#' The function returns annotations (see details below) of all features
#' (probably probesets) associated with the given taxon.
#' 
#' The function reads from the backend Oracle database.
#' 
#' @param taxid Character string, the TaxID of the species in interest. For
#' instance \sQuote{9606} for Homo sapiens.
#' @return A \code{data.frame} object with very similar structure as the
#' \code{EG_GENE_INFO} table in the database.
#' 
#' Rownames of the \code{data.frame} are set to \code{NULL}.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{gtiChipAnnotation}} for annotating chips.
#' @references \url{http://bioinfo.bas.roche.com:8080/bicgi/gti_ui.cgi}
#' @examples
#' 
#' \dontrun{
#' hsAnno <- gtiTaxAnnotation("9606")
#' dim(hsAnno)
#' head(hsAnno)
#' }
#' 
#' @export gtiTaxAnnotation
gtiTaxAnnotation <- function(taxid) {
  if(missing(taxid))
    stop("'taxid' cannot be missing.")
  state <- paste("SELECT a.*, b.SEQ, b.CLEFT, b.CRIGHT, b.REVCOMP FROM (SELECT TO_CHAR(g.GENEID) as RO_GENE_ID, g.OFFICIAL_SYMBOL, g.OFFICIAL_NAME, g.SYNONYMS, g.DBXREFS, ",
                 "g.CHROMOSOME, g.MAP_LOCATION, g.GENE_TYPE, ",
                 "g.TAX_ID ",
                 "FROM bi.EG_GENE_INFO g ",
                 "where g.TAX_ID = '",taxid, "') a, genome.gti_gene_map_a@genome b ",
                 "WHERE a.RO_GENE_ID=b.RO_GENE_ID(+)", sep="")
  ann <- querydb(state, db=dbName(), user=biaUser(), password=biaPwd())
  colnames(ann) <- c("GeneID", "GeneSymbol", "GeneName",
                     "Synonyms", "Xrefs",
                     "Chromosome", "MapLocation",
                     "GeneType","TaxID",
                     "MappedChr", "CoordLeft", "CoordRight", "RevComp")
  ##post processing
  ann$RevComp <- ifelse(ann$RevComp==1L, TRUE, FALSE)
  ann$MappedChr <- gsub("^CHR", "", as.character(ann$MappedChr))
  rownames(ann) <- NULL
  return(ann)
}
