#' @include utils.R annotateAnyIDs.R querydb.R sortAnnotationByQuery.R
NULL

id2rownames <- function(ids) {
  if(identical(anyDuplicated(ids), 0L) & !any(is.na(ids))) {
    return(ids)
  } else {
    return(NULL)
  }
}
notValid <- function(x)  is.null(x) || is.na(x) || tolower(x)=="any" || x==""

##------------------------------------------------------------##
## Annotate Probesets
## Scenario 1: annotate probesets with known chip type
## Scenario 2: annotate probesets with unknown (single) chip type
## Scenario 3: annotate probesets with mixed chip types
##------------------------------------------------------------##

## Scenario 1: annotate probesets with known chip type
## TODO: gtiChipAnnotation adds support for othologue mapping


#' Annotate human orthologues of genes with given GeneID
#' 
#' The function maps any given GeneID to their human orthologues, and annotate
#' them with essential information.
#' 
#' The function maps any GeneID to human orthologues by querying Oracle
#' database.
#' 
#' @param geneids Character or numeric vector, Entrez GeneIDs to be queried
#' @param multiOrth Logical, whether multiple orthologues to the same gene
#' should be all returned. By default \code{FALSE}
#' @return A \code{data.frame} containing following columns:
#' \item{OrigGeneID}{Entrez Gene ID in query} \item{OrigTaxID}{Taxonomy ID of
#' input Entrez Gene ID} \item{TaxID}{Taxonomy ID of orthologue Entrez Gene ID.
#' Currently it is fixed as 9606, the human TaxID.} \item{GeneID}{Entrez Gene
#' ID of human orthologue gene} \item{GeneSymbol}{GeneSymbol of human
#' orthologue gene}
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{humanOrthologs}} and \code{\link{humanUniqOrtholog}}
#' @examples
#' 
#' options(error=utils::recover)
#' annotateHumanOrthologs(c(100034253, 100036582, 1000125361, 1432))
#' options(error=NULL)
#' 
#' @importFrom ribiosUtils matchColumn putColsFirst
#' @export annotateHumanOrthologs
annotateHumanOrthologs <- function(geneids, multiOrth=FALSE) {
  geneids <- unique(geneids)
  comm <- paste("SELECT a.RO_GENE_ID2, a.TAX_ID2, 9606 as TaxID, a.RO_GENE_ID1, b.GENE_SYMBOL",
                "FROM genome.GTI_ORTHOLOGS a, genome.GTI_GENES b",
                "WHERE a.RO_GENE_ID1=b.RO_GENE_ID AND a.TAX_ID1 ='9606'")
  ort <- querydbTmpTbl(comm,
                       "a.RO_GENE_ID2",
                       geneids, dbName(), binUser(), binPwd())
  colnames(ort) <- c("OrigGeneID", "OrigTaxID", "TaxID", "GeneID", "GeneSymbol")
  hGeneIDs <- as.integer(ort$GeneID)
  ort <- ort[order(hGeneIDs, decreasing=FALSE),] ## Smaller GeneIDs come first
  res <- matchColumn(geneids, ort, "OrigGeneID", multi=multiOrth)
  return(res)
}

annotateHumanOrthologsNoOrigTax <- function(...) {
  res <- annotateHumanOrthologs(...)
  return(res[,-2L])
}

#' Annotations of all features associated with the given chip type
#' 
#' The function returns annotations (see details below) of all or selected
#' features (probably probesets) associated with the given chip type.
#' 
#' The function is useful when someone wants to annotate all probesets on a
#' certain chip type, or some probesets of a known chip type. It connects to
#' the GTI Oracle database and performs the annotation.
#' 
#' It is also able to perform orthologue mapping by setting \code{orthologue}
#' to \code{TRUE}.
#' 
#' \code{annotateChip} is the synonymous function of \code{gtiChiptype}.
#' 
#' @aliases ORACLE.IN.NMAX gtiChipAnnotation annotateChip raceChipAnnotation
#' biosCurrentGeneSymbol
#' @param chip Character string, the chip type in interest. For a complete list
#' of supported chip names, use the \code{gtiChiptypes()} function.
#' @param ids Optional, character string vectors. IDs (probably probesets) to
#' be queried in the given chip
#' @param orthologue Logical, should genes be mapped to human orthologues?
#' Default is \code{FALSE}
#' @param multiOrth Logical, should multiple-mapping orthologues all returned.
#' Default is \code{FALSE}
#' @return A \code{data.frame} object.
#' 
#' If \code{orthologue=FALSE}, the data.frame contains following columns
#' \item{ProbeID}{Probeset ID} \item{GeneID}{Entrez Gene ID}
#' \item{GeneSymbol}{Entrez (HUGO) Gene Symbol} \item{GeneName}{Gene
#' descriptions} \item{Chip}{Chip type} \item{TaxID}{NCBI Taxonomy ID}
#' 
#' If \code{orthologue=TRUE}, the data.frame contains following columns:
#' \item{ProbeID}{Probeset ID} \item{GeneID}{Entrez Gene ID of human orthologue
#' gene (if existing)} \item{GeneSymbol}{Entrez (HUGO) Gene Symbol of human
#' orthologue gene (if existing)} \item{GeneName}{Gene descriptions of original
#' gene} \item{Chip}{Chip type} \item{TaxID}{Human taxonomy ID (9606)}
#' \item{OrigGeneID}{Entrez Gene ID of the gene to which queried probesets map}
#' \item{OrigGeneSymbol}{Gene Symbol of the gene to which queried probesets
#' map} \item{OrigTaxID}{Taxonomy ID of the gene to which queried probesetsmap}
#' @note From version 2.0-0 and on, \code{isSingleGeneID} has been removed.
#' 
#' From version 1.0-15 and on, \code{isSingleGeneID} is always set to
#' \code{TRUE} since only uniquely mapped probesets are stored. Provisionally
#' this column shall be removed in the future version.
#' @section Warning: If \code{multiOrth} is \code{TRUE}, the returning
#' data.frame may contain more rows than the number of input probesets.This is
#' not desired in situations where the output contains exactly the same number
#' of probesets as the input, e.g. annotating an \code{ExpressionSet} object.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{gtiChiptypes}} for a complete list of supported chip
#' types.
#' @examples
#' 
#' \dontrun{
#' mychip <- "HG-U133_PLUS_2"
#' mychipAnno <- gtiChipAnnotation(mychip)
#' dim(mychipAnno)
#' head(mychipAnno)
#' }
#' 
#' options(error=utils::recover)
#' gtiChipAnnotation(mychip, ids=c("1053_at","117_at", "121_at"))
#' gtiChipAnnotation(mychip, ids=c("1053_at","117_at", "121_at"),
#' orthologue=TRUE)
#' gtiChipAnnotation(mychip, ids=c("1053_at","117_at", "121_at"),
#' orthologue=TRUE, multiOrth=TRUE)
#' options(error=NULL)
#' 
#' @export gtiChipAnnotation
gtiChipAnnotation <- function(chip,ids, orthologue=FALSE, multiOrth=FALSE) {
  if(missing(chip)) stop("'chip' cannot be missing.Using 'gtiChiptypes()' to find supported chip types.")
  
  hasInColVal <- !missing(ids)
  woOrthCn <- c("ProbeID", "GeneID", "GeneSymbol", "GeneName", "Chip", "TaxID")
  wOrthCn <- c(woOrthCn, c("OrigTaxID", "OrigGeneID", "OrigGeneSymbol"))
  
  if(!hasInColVal) {
    state <- paste("SELECT a.PROBESET_ID, e.RO_GENE_ID, e.GENE_SYMBOL, e.DESCRIPTION, a.ARRAY_TYPE, e.TAX_ID ",
                   "FROM genome.chip_probeset_gene a, genome.GTI_GENES e ",
                   "where a.RO_GENE_ID(+) = e.RO_GENE_ID ",
                   "AND ARRAY_TYPE='",chip, "'", sep="")
    ann <- querydb(state, db=dbName(), user=binUser(), password=binPwd())
  } else {
    state <-  paste("SELECT a.PROBESET_ID, e.RO_GENE_ID, e.GENE_Symbol, e.DESCRIPTION, a.ARRAY_TYPE, e.TAX_ID ",
                    "FROM genome.chip_probeset_gene a, genome.GTI_GENES e ",
                    "WHERE a.RO_GENE_ID(+) = e.RO_GENE_ID ",
                    "AND ARRAY_TYPE='",chip, "'",sep="")
    ann <- querydbTmpTbl(state, inCol="a.PROBESET_ID",inValues=ids, db=dbName(), user=binUser(), password=binPwd())
  }

  if(!orthologue) {
    colnames(ann) <- woOrthCn
    res <- ann
  } else {
    colnames(ann) <- c("ProbeID", "OrigGeneID", "OrigGeneSymbol", "GeneName", "Chip", "OrigTaxID")
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, multiOrth=multiOrth)
    if(multiOrth) {
      res <- merge(ann, ort, by="OrigGeneID", all.x=TRUE)
    } else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", multi=FALSE)
      res <- cbind(ann, ort.re[,-1L])
    }
    res <- putColsFirst(res, wOrthCn)
  }
  
  if(hasInColVal) {
    res <- matchColumn(ids, res, "ProbeID", multi=orthologue && multiOrth)
  }
  rownames(res) <- id2rownames(res$ProbeID)
  return(res)
}

#' @export annotateChip
annotateChip <- gtiChipAnnotation

## Scenario 2: annotate probesets with unknown, but single chip type


#' Guess likely chiptype(s) from probesets
#' 
#' Given a vector of probesets, guess the most likely chiptype from which these
#' probesets come from. This function can be used in situations where probesets
#' are known to be all from a single, unknown chip type.
#' 
#' The function is used where a set of probesets have to be annotated, and they
#' are probably coming from the same chip type. If they may come from different
#' chips, use \code{\link{annotateAnyProbeset}} instead; if the chip type from
#' which they come is known, use \code{\link{gtiChipAnnotation}} which is
#' faster.
#' 
#' The function works as follows: first it samples a subset of input probesets
#' (given by the \code{sample} parameter, in order to speed up the function and
#' to save computational resources). Then it determines via SQL query which
#' chip types contain these probesets. All types containing at least one
#' probeset as in the input are returned, and they are sorted by the number of
#' matching probesets. If \code{maxOnly} is \code{TRUE}, only the most likely
#' chip type is returned.
#' 
#' @param ids Character vector, probeset IDs
#' @param maxOnly Logical, should only the chiptype to which most probesets
#' were mapped be returned? By default \code{FALSE}, namely all chiptypes are
#' returned which at least one probeset could be mapped.
#' @param sample Integer to sample probesets. See \code{details} below.
#' @return Character vector containing the chip types, with an attribute named
#' \code{count} recording number of matched probesets in the chiptypes.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{annotateProbesets}}, \code{\link{gtiChipAnnotation}},
#' \code{\link{annotateAnyProbeset}}.
#' @examples
#' 
#' options(error=utils::recover)
#' guessChiptype(ids=c("1399067_at","1378257_at", "1380581_at"))
#' guessChiptype(ids=c("1399067_at","1378257_at", "1380581_at"), maxOnly=TRUE)
#' options(error=NULL)
#' 
#' @export guessChiptype
guessChiptype <- function(ids, maxOnly=FALSE, sample=100) {
  ids <- unique(as.character(ids))
  if(!is.null(sample) && !is.na(sample) && is.numeric(sample)) {
    ids <- sample(ids, pmin(as.integer(sample), length(ids)), replace=FALSE)
  }
  state <- paste("SELECT a.PROBESET_ID, a.ARRAY_TYPE ",
                 "FROM genome.chip_probeset_gene a WHERE ")
  ann <- querydbSelectIn(state, inCol="a.PROBESET_ID", inValues=ids,
                         db=dbName(), user=binUser(), password=binPwd())
  res.raw <- sort(table(ann$ARRAY_TYPE), decreasing=TRUE)
  if(maxOnly) res.raw <- res.raw[1]
  res <- names(res.raw)
  attr(res, "count") <- unname(res.raw)
  return(res)
}

## Scenario 3: annotate probesets with mixed chip type
## TODO: add OrigGeneSymbol


#' Annotate any probeset without knowing chiptype
#' 
#' Annotate a set of probesets possibily of mixed chip types
#' 
#' The function annotates input probesets with essential information. These
#' probesets probably come from unknown and mixed chiptypes. In case they come
#' from a known, single chip type, consider using
#' \code{\link{gtiChipAnnotation}} which is faster; in case they are known to
#' come from a single chip but the type is unknown, consider using
#' \code{\link{guessChiptype}} first to guess the most likely chip type and
#' then use \code{\link{gtiChipAnnotation}}.
#' 
#' @param ids Character vector, probesets to be annotated
#' @param orthologue Logical, should genes be mapped to human orthologues.
#' @param multiOrth Logical, should multiple orthologues all be returned.
#' @return If \code{orthologue=FALSE}, the data.frame contains following
#' columns \item{ProbeID}{Probeset ID} \item{GeneID}{Entrez Gene ID}
#' \item{GeneSymbol}{Entrez (HUGO) Gene Symbol} \item{GeneName}{Gene
#' descriptions} \item{Chip}{Has NA as value} \item{TaxID}{NCBI Taxonomy ID}
#' 
#' If \code{orthologue=TRUE}, the data.frame contains following columns:
#' \item{ProbeID}{Probeset ID} \item{GeneID}{Entrez Gene ID of human orthologue
#' gene (if existing)} \item{GeneSymbol}{Entrez (HUGO) Gene Symbol of human
#' orthologue gene (if existing)} \item{GeneName}{Gene descriptions of original
#' gene} \item{Chip}{Has NA as value} \item{TaxID}{Human taxonomy ID (9606)}
#' \item{OrigGeneID}{Entrez Gene ID of the gene to which queried probesets map}
#' \item{OrigGeneSymbol}{Gene Symbol of the gene to which queried probesets
#' map} \item{OrigTaxID}{Taxonomy ID of the gene to which queried probesetsmap}
#' @section Warning: If \code{multiOrth} is \code{TRUE}, the returning
#' data.frame may contain more rows than the number of input probesets.This is
#' not desired in situations where the output contains exactly the same number
#' of probesets as the input, e.g. annotating an \code{ExpressionSet} object.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{gtiChipAnnotation}} and \code{\link{guessChiptype}}
#' @examples
#' 
#' options(error=utils::recover)
#' myProbesets <- gtiChipAnnotation(chip="HG-U133_PLUS_2")
#' freeProbesets <- annotateAnyProbeset(ids=myProbesets$ProbeID[1:100])
#' system.time(freeProbesetsSomeOrtho <- annotateAnyProbeset(ids=myProbesets$ProbeID[1:100],
#'                                                           orthologue=TRUE))
#' options(error=NULL)
#' 
#' @export annotateAnyProbeset
annotateAnyProbeset <- function(ids, orthologue=FALSE, multiOrth=FALSE) {
  comm <- paste("SELECT a.ANY_ID ,a.RO_GENE_ID ,c.GENE_SYMBOL, c.DESCRIPTION, 'NA' AS chip, c.TAX_ID ",
                " FROM genome.GTI_IDMAP a, ",
                " GTI_GENES c",
                " WHERE a.RO_GENE_ID=c.RO_GENE_ID", sep="")
  cnames <- c("ProbeID", "GeneID", "GeneSymbol", "GeneName", "Chip", "TaxID")
  conames <- c(cnames, c("OrigTaxID", "OrigGeneID", "OrigGeneSymbol"))
  ann <- querydbTmpTbl(comm,
                       "a.ANY_ID",
                       ids, dbName(), binUser(), binPwd())

  
  if(!orthologue) {
    colnames(ann) <- cnames
    res <- ann
  }
  if(orthologue) {
    colnames(ann) <- c("ProbeID", "OrigGeneID", "OrigGeneSymbol", "GeneName", "Chip", "OrigTaxID")
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, multiOrth=multiOrth)
    if(multiOrth) {
      res <- merge(ann, ort, by="OrigGeneID", all.x=TRUE)
    } else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", multi=FALSE)
      res <- cbind(ann, ort.re[,-1L])
    }
    res <- putColsFirst(res, conames)
  }
  res <- matchColumn(ids, res, "ProbeID", multi=orthologue && multiOrth)
  rownames(res) <- id2rownames(res$ProbeID)
  return(res)
}



#' Annotate probesets by their identifiers (ID) and chip type
#' 
#' Annotate probesets in chip types supported, as listed by the
#' \code{\link{gtiChiptypes}} function. If no chip type is given, the function
#' will automatically map probesets to genes via GTI.
#' 
#' The function is a wrapper of \code{gtiChiptype} (when \code{chip} is given)
#' and \code{annotateAnyProbeset} (when \code{chip} is missing).
#' 
#' \code{annotateProbeIDs} is a synonym of \code{annotateProbesets}.
#' 
#' The concept \sQuote{Probesets} here can either refer to a set of probes, or
#' a single probe, depending on manufacturing techniques used by the microarray
#' vendor.
#' 
#' @aliases annotateProbesets annotateProbeIDs
#' @param ids Probeset identifiers
#' @param chip Chip type identifier, for example \code{HG-U133_PLUS_2}. If
#' missing (or \code{NA}, \code{NULL}, \dQuote{any} and empty string),
#' automatic mapping is done via GTI.
#' @param orthologue Logical, should be human orthologues be returned?
#' @return A \code{data.frame} annotating probesets.
#' 
#' Rows are in the same order as of the given probesets. When the input ids are
#' unique, they are used as row names of the resulting \code{data.frame};
#' otherwise the row names is set to \code{NULL}, which when printed will be
#' shown as incremental integers.
#'
#' If \code{orthologue} is set to \code{TRUE}, orthologue information is
#' returned.
#' @note If the \code{chip} option is given, internally the function calls the
#' \code{gtiChipAnnotation} to first query all annotations about the given
#' chip. Subsequently they are matched against queried probesets, returning
#' their annotations only.
#' 
#' If no \code{chip} is given, the function calls \code{annotateAnyProbeset} to
#' annotate probesets.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso See \code{\link{gtiChiptypes}} for supported chip types.
#' 
#' See \code{\link{gtiChipAnnotation}} to get annotation for all probesets in a
#' chip, \code{\link{annotateAnyProbeset}} to get annotations for probesets
#' from unknown chiptype(s).
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' myprobes <- c("1000_at", "1004_at", "1002_f_at", "nonsense_at")
#' annotateProbesets(myprobes, chip="HG_U95A")
#' 
#' myprobes2 <- c( "Hs2.356162.1.S1_3p_at","g4885506_3p_at",
#' "Hs.210202.0.A1_3p_at", "Hs.166262.0.A1_3p_x_at",
#' "Hs.48499.1.S1_3p_a_at", "Hs.143268.0.A1_3p_x_at",
#' "Hs.86386.1.A1_3p_at", "g11497042_3p_a_at", "Hs.279556.1.S2_3p_at",
#' "Hs.83484.1.A1_3p_at")
#' annotateProbesets(myprobes2, chip="U133_X3P")
#' 
#' ## when using wrong chip it is likely that many (if not all)
#' ## features will not be annotated
#' annotateProbesets(myprobes2, chip="HGU_95A")
#' 
#' ## automatically mapping
#' annotateProbesets(c("1373268_at", "1002_f_at", "103416_at"))
#' annotateProbesets(c("1373268_at", "1002_f_at", "103416_at"), orthologue=TRUE)
#' 
#' options(error=NULL)
#' 
#' @export annotateProbesets
annotateProbesets <- function(ids, chip, orthologue=FALSE) {
  if(missing(chip) || notValid(chip)) {
    annotateAnyProbeset(ids, orthologue=orthologue)
  } else {
    gtiChipAnnotation(chip, ids=ids, orthologue=orthologue)
  }
}

#' @export annotateProbeIDs
annotateProbeIDs <- annotateProbesets



#' Annotate GeneIDs, GeneSymbols, or mRNA accession numbers.
#' 
#' Annotate GeneIDs, GeneSymbols or mRNA accesion numbers, optionally with
#' human orthologue mappings
#' 
#' Prior to version \code{1.2-0}, GeneIDs and GeneSymbols are mapped via chip
#' annotation files, which introduce unnecessary dependency and confusing results
#' when chip does not contain all genes.
#' 
#' From version \code{1.2-0}, \code{annotateGeneIDs} and
#' \code{annotateGeneSymbols} do not dependent on chip annotations, but on the
#' core tables of GTI, which can be more efficient. In addition, new
#' functionalities to map orthologues are added.
#' 
#' From version \code{2.0-0}, the speed of both functions are further
#' optimized.
#' 
#' @aliases annotateGeneIDs annotateGeneSymbols annotateRefSeqs annotatemRNAs
#' @param ids Character vector, GeneIDs, GeneSymbols or mRNAs of query. It can
#' contain \code{NA} or \code{NULL}
#' @param orthologue Logical, whether human orthologues should be returnd
#' @param multiOrth Logical, whether multiple mapped orthologues should be
#' returned or not. Only useful when \code{orthologue} is \code{TRUE}, By
#' default \code{FALSE}. Be cautious with the \code{TRUE} option: in this case
#' the returning \code{data.frame} may have different row numbers as the input
#' vector.
#' @return A \code{data.frame} object containing the annotations.
#' 
#' If \code{orthologue} is set to \code{FALSE}, the data frame contains
#' following columns: \code{GeneID, GeneSymbol, GeneName and TaxID}
#' 
#' If \code{orthologue} is \code{TRUE}, the data frame contains \code{GeneID,
#' GeneSymbol,TaxID, OrigTaxID, OrigGeneID, OrigGeneSymbol, OrigGeneName}. Note
#' that \code{GeneID, GeneSymbol, TaxID} contains the information of mapped
#' orthologues, while \code{OrigXXX} contains the information of queried genes.
#' 
#' If \code{multiOrth} is \code{TRUE} (only valid when \code{orthologue} is
#' \code{TRUE}), multiples orthologues mapped to the same gene are all
#' returned. This means that the result data frame may contain more rows than
#' the length of input vector. If set to \code{FALSE}, the resulting data frame
#' contains exact number of rows as the input vector length.
#' @note \code{annotatemRNAs} is an alias of \code{annotateRefSeqs}
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso See \code{\link{gtiChipAnnotation}} to get annotation for all
#' probesets in a chip.
#' 
#' See \code{\link{annotateProbesets}} to get annotation for probesets.
#' @examples
#' 
#' options(error=utils::recover)
#' 
#' ## normal use
#' annotateGeneIDs(ids=c(780, 5982, 3310, NA))
#' annotateGeneIDs(ids=c(780, 5982, 3310, NULL))
#' annotateGeneSymbols(ids=c("DDR1", "RFC2", "HSPA6",
#' "HSAP6"),organism="human", orthologue=FALSE)
#' 
#' ## note that GeneSymbols are not unique
#' annotateGeneSymbols(ids=c("DDR1", "RFC2", "HSPA6", "HSAP6"),organism="any", orthologue=FALSE)
#' 
#' ## orthologues
#' annotateGeneIDs(ids=c(1017, 12566, 362817, 26416, 81649),
#' orthologue=FALSE)
#' annotateGeneIDs(ids=c(1017, 12566, 362817, 26416, 81649),
#' orthologue=TRUE)
#' annotateGeneSymbols(ids=c("DDR1", "RFC2", "HSPA6",
#' "HSAP6"),organism="any", orthologue=TRUE)
#' annotateGeneIDs(ids=c(1017, 12566, 362817, 26416, 81649),
#' orthologue=TRUE, multiOrth=TRUE)
#' 
#' ## Following examples underlines the non-uniqueness of GeneSymbols
#' annotateGeneSymbols(ids=c("Cdk2", "CDK4", "Mapk14"),organism="rat",
#' orthologue=TRUE)
#' annotateGeneSymbols(ids=c("Cdk2", "CDK4", "Mapk14"),organism="mouse", orthologue=TRUE)
#' annotateGeneSymbols(ids=c("Cdk2", "CDK4", "Mapk14"),organism="any", orthologue=TRUE)
#' 
#' ## Annotate mRNA
#' annotatemRNAs(c("NM_007158", "NM_001007553"))
#' options(error=NULL)
#' 
#' @export annotateGeneSymbols
annotateGeneSymbols <- function(ids,
                                organism=c("human", "mouse", "rat", "any"),
                                orthologue=FALSE,
                                multiOrth=FALSE) {
  organism <- match.arg(organism)
  oid <- c("any"="", "human"="9606", "mouse"="10090", rat="10116")[organism]
  comm <- paste("SELECT c.RO_GENE_ID,c.GENE_SYMBOL, c.DESCRIPTION, c.TAX_ID ",
                " FROM GTI_GENES c ",
                ifelse(organism=="any", "", paste("WHERE c.TAX_ID='", oid, "' ",sep="")),
                sep="")
  cnames <- c("GeneID", "GeneSymbol", "GeneName", "TaxID")
  ann <- querydbTmpTbl(comm,
                       "c.GENE_SYMBOL",
                       ids, dbName(), binUser(), binPwd())

  if(!orthologue) {
    colnames(ann) <- cnames
    cn <- "GeneSymbol"
    res <- ann
  } else {
    colnames(ann) <- c("OrigGeneID", "OrigGeneSymbol", "OrigGeneName", "OrigTaxID")
    cn <- "OrigGeneSymbol"
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, multiOrth=multiOrth)
    if(multiOrth) {
      res <- merge(ann, ort, by="OrigGeneID", all.x=TRUE)
    } else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", multi=FALSE)
      res <- cbind(ann, ort.re[,-1L])
    }
    res <- putColsFirst(res, c("GeneID", "GeneSymbol", "TaxID",
                               "OrigTaxID", "OrigGeneID", "OrigGeneSymbol", "OrigGeneName"))
  }
                       
  res <- matchColumn(ids, res, cn, multi=multiOrth)
  rownames(res) <- id2rownames(res[,cn])
  return(res)
}

#' @export annotateRefSeqs
annotateRefSeqs <- function(ids, orthologue=FALSE, multiOrth=FALSE) {
  comm <- paste("SELECT a.item_id, c.RO_GENE_ID ,c.GENE_SYMBOL, c.DESCRIPTION, 'NA' AS chip, c.TAX_ID ",
                " FROM genome.GTI_GENE_ITEMS a, ",
                " GTI_GENES c",
                " WHERE a.item_type_id in (3,4) AND a.ro_gene_id=c.ro_gene_id", sep="")
  cnames <- c("mRNA", "GeneID", "GeneSymbol", "GeneName", "Chip", "TaxID")
  conames <- c(cnames, c("OrigTaxID", "OrigGeneID", "OrigGeneSymbol"))
  ann <- querydbTmpTbl(comm,
                       "lower(a.item_id)",
                       tolower(ids), dbName(), binUser(), binPwd())
  if(!orthologue) {
    colnames(ann) <- cnames
    res <- ann
  }
  if(orthologue) {
    colnames(ann) <- c("mRNA", "OrigGeneID", "OrigGeneSymbol", "GeneName", "Chip", "OrigTaxID")
    ort <- annotateHumanOrthologsNoOrigTax(ann$OrigGeneID, multiOrth=multiOrth)
    if(multiOrth) {
      res <- merge(ann, ort, by="OrigGeneID", all.x=TRUE)
    } else {
      ort.re <- matchColumn(ann$OrigGeneID, ort, "OrigGeneID", multi=FALSE)
      res <- cbind(ann, ort.re[,-1L])
    }
    res <- putColsFirst(res, conames)
  }
  res <- matchColumn(ids, res, "mRNA", multi=orthologue && multiOrth)
  rownames(res) <- id2rownames(res$mRNA)
  return(res)
}

#' @export annotatemRNAs
annotatemRNAs <- annotateRefSeqs
