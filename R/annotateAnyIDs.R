## if only looking after proteins, it is also possible to check the GA_PROTEIN_GENE table in the protein component of the Genome Analysis Pipeline

## note by zhangj83: this file is not yet part of the package. needs to be tested. On 04.05.2015



#' Annotate any identifiers
#' 
#' This annotates any identifies that can be recognized by GTI.
#' 
#' 
#' @param ids A vector of identifiers. They must be of the same type: GeneID,
#' GeneSymbol or Probesets.
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
#' annotateAnyIDs(ids=c(780, 5982, 3310, NA))
#' 
#' annotateAnyIDs(ids=c("DDR1", "RFC2", "HSPA6", "HSAP6"))
#' 
#' myprobes <- c("1000_at", "1004_at", "1002_f_at", "nonsense_at")
#' annotateProbesets(myprobes, chiptype="HG_U95A")
#' 
#' annotateAnyIDs(ids=c("P38398", "Q8NDF8"))
#' options(error=NULL)
#' 
#' @export annotateAnyIDs
annotateAnyIDs <- function(ids, orthologue = FALSE, multiOrth = FALSE) {
  comm <- paste("SELECT m.ANY_ID, m.ID_TYPE, c.RO_GENE_ID,c.GENE_SYMBOL, c.DESCRIPTION, c.TAX_ID ", 
                " FROM GTI_GENES c JOIN GTI_IDMAP m ON c.RO_GENE_ID=m.RO_GENE_ID ",sep="")
  ann <- querydbTmpTbl(comm, "m.ANY_ID", ids, dbName(),
                       binUser(), 
                       binPwd())
  cnames <- c("Input", "InputIDType", "GeneID", "GeneSymbol", "GeneName", "TaxID")
  conames <- c("Input", "InputIDType", "OrigGeneID", "OrigGeneSymbol", "OrigGeneName", 
               "OrigTaxID")
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
  rownames(res) <- id2rownames(res[, cn])
  return(res)
}

## testIDs <- c("O43684", "O43670")
## testAnno <- annotateAnyIDs(testIDs)
