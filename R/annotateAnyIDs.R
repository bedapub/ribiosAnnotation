## if only looking after proteins, it is also possible to check the GA_PROTEIN_GENE table in the protein component of the Genome Analysis Pipeline

#' Remove version suffix from Ensembl IDs
#' @param ensemblIDs A vector of character strings. Other types of inputs are 
#'   converted.
#' @return A character vector of the same length as input
#' 
#' @examples 
#' ensemblIDs <- c("ENSG00000197535", "ENST00000399231.7", "ENSP00000418960.2")
#' removeEnsemblVersion(ensemblIDs)
#' @export
removeEnsemblVersion <- function(ensemblIDs) {
  return(gsub("\\.[0-9]+$", "", as.character(ensemblIDs)))
}

#' Annotate any identifiers
#' 
#' This annotates any identifies that can be recognized by GTI.
#' 
#' 
#' @param ids A vector of identifiers. They must be of the same type. Supported 
#' types include Entrez GeneID, GeneSymbol, Probesets, UniProt identifiers, 
#' NCBI RefSeq mRNA identifiers, and Ensembl gene identifiers (with possible 
#' version suffixes).
#' 
#' @param orthologue Logical, is orthologous mapping needed?
#' @param multiOrth Logical, is more t han one orthologs allowed
#' @return A \code{data.frame} containing annotation information.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{annotateGeneIDs}}, \code{\link{annotateGeneSymbols}}
#' and \code{\link{annotateProbesets}}
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' # GeneID
#' annotateAnyIDs(ids=c(780, 5982, 3310, NA))
#' 
#' # GeneSymbol
#' annotateAnyIDs(ids=c("DDR1", "RFC2", "HSPA6", "HSAP6"))
#' 
#' # Probesets
#' myprobes <- c("1000_at", "1004_at", "1002_f_at", "nonsense_at")
#' annotateProbesets(myprobes, chiptype="HG_U95A")
#' 
#' # UniProt
#' annotateAnyIDs(ids=c("P38398", "Q8NDF8"))
#' 
#' # EnsEMBL 
#' ensemblIDs <- c("ENSG00000197535", "ENST00000399231.7", "ENSP00000418960.2")
#' annotateAnyIDs(ensemblIDs)
#' 
#' # RefSeq
#' annotateAnyIDs(c("NM_000235", "NM_000498"))
#' 
#' options(error=NULL)
#' 
#' @export annotateAnyIDs
annotateAnyIDs <- function(ids, orthologue = FALSE, multiOrth = FALSE) {
  inputIds <- ids
  ids <- removeEnsemblVersion(inputIds)
  comm <- paste("SELECT m.ANY_ID, m.ID_TYPE, c.RO_GENE_ID,c.GENE_SYMBOL, c.DESCRIPTION, c.TAX_ID ", 
                " FROM GTI_GENES c JOIN GTI_IDMAP m ON c.RO_GENE_ID=m.RO_GENE_ID ",sep="")
  ann <- querydbTmpTbl(comm, "m.ANY_ID", ids, dbName(),
                       binUser(), 
                       binPwd())
  cnames <- c("Input", "InputIDType", "GeneID", "GeneSymbol", "GeneName", "TaxID")
  conames <- c("Input", "InputIDType", "OrigGeneID", "OrigGeneSymbol",
               "OrigGeneName", "OrigTaxID")
  cn <- "Input"
  if (!orthologue) {
    colnames(ann) <- cnames
    res <- ann
  } else {
    colnames(ann) <- conames
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, 
                                           multiOrth = multiOrth)
    if (multiOrth) {
      res <- merge(ann, ort, by = "OrigGeneID", all.x = TRUE)
    }  else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", 
                            multi = FALSE)
      res <- cbind(ann, ort.re[, -1L])
    }
    res <- putColsFirst(res, c("Input", "InputIDType", "GeneID", "GeneSymbol", "TaxID", 
                               "OrigTaxID", "OrigGeneID", "OrigGeneSymbol", "OrigGeneName"))
  }
  res <- matchColumn(ids, res, cn, multi = orthologue && multiOrth)
  if(!identical(inputIds, ids)) {
    res[, cn] <- inputIds[match(res[, cn], ids)]
  }
  rownames(res) <- id2rownames(res[, cn])
  return(res)
}
