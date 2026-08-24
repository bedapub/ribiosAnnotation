# Internal utilities for selecting annotation backend and resolving
# Bioconductor OrgDb packages by taxonomy ID.

normalizeAnnotationBackend <- function(backend = NULL) {
  if (is.null(backend) || is.na(backend) || backend == "") {
    backend <- getOption("ribiosAnnotation.backend", default = NA_character_)
  }
  if (is.null(backend) || is.na(backend) || backend == "") {
    backend <- Sys.getenv("RIBIOS_ANNOTATION_BACKEND", unset = "mongodb")
  }
  backend <- tolower(as.character(backend)[1])
  if (!(backend %in% c("mongodb", "bioconductor"))) {
    stop("backend must be one of: 'mongodb', 'bioconductor'", call. = FALSE)
  }
  backend
}

#' Set Global Annotation Backend
#'
#' Set the default annotation backend for the current R session.
#' This affects all annotation function calls where the \\code{backend}
#' argument is not explicitly provided.
#'
#' @param backend Character string, one of \\code{"mongodb"} or
#'   \\code{"bioconductor"}.
#' @return Invisibly returns the selected backend.
#' @export
setAnnotationBackend <- function(backend = c("mongodb", "bioconductor")) {
  backend <- normalizeAnnotationBackend(match.arg(backend))
  options(ribiosAnnotation.backend = backend)
  invisible(backend)
}

#' Get Effective Annotation Backend
#'
#' Returns the currently effective backend after resolving function argument,
#' option, and environment variable fallbacks.
#'
#' @return Character string, one of \\code{"mongodb"} or
#'   \\code{"bioconductor"}.
#' @export
getAnnotationBackend <- function() {
  normalizeAnnotationBackend(NULL)
}

biocOrgDbPackageByTaxId <- function(taxId) {
  taxId <- suppressWarnings(as.integer(taxId))
  pkg <- switch(as.character(taxId),
                "9606" = "org.Hs.eg.db",
                "10090" = "org.Mm.eg.db",
                "10116" = "org.Rn.eg.db",
                "9615" = "org.Cf.eg.db",
                "9823" = "org.Ss.eg.db",
                "9598" = "org.Pt.eg.db",
                "9541" = "org.Mmu.eg.db",
                NA_character_)
  pkg
}

biocSupportedOrgDbTaxIds <- function() {
  c(9606L, 10090L, 10116L, 9615L, 9823L, 9598L, 9541L)
}

biocRequireOrgDbByTaxId <- function(taxId) {
  pkg <- biocOrgDbPackageByTaxId(taxId)
  if (is.na(pkg)) {
    stop("Bioconductor backend supports taxId 9606 (human), 10090 (mouse), 10116 (rat), 9615 (dog), 9823 (pig), 9598 (chimpanzee), and 9541 (macaque via org.Mmu.eg.db).",
         call. = FALSE)
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required for backend='bioconductor'.", call. = FALSE)
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for backend='bioconductor'.",
         " Please install it with BiocManager::install('", pkg, "').",
         call. = FALSE)
  }
  getExportedValue(pkg, pkg)
}

biocEnsemblTaxId <- function(ensemblIDs) {
  x <- as.character(ensemblIDs)
  taxId <- rep(NA_integer_, length(x))
  taxId[grepl("^ENSG", x)] <- 9606L
  taxId[grepl("^ENSMUSG", x)] <- 10090L
  taxId[grepl("^ENSRNOG", x)] <- 10116L
  taxId[grepl("^ENSCAFG", x)] <- 9615L
  taxId[grepl("^ENSSSCG", x)] <- 9823L
  taxId[grepl("^ENSPTRG", x)] <- 9598L
  taxId[grepl("^ENSMMUG", x)] <- 9541L
  taxId
}

biocTaxIdToBabelSpecies <- function(taxId) {
  taxId <- suppressWarnings(as.integer(taxId))
  switch(as.character(taxId),
         "9606" = "Homo sapiens",
         "9598" = "Pan troglodytes",
         "9541" = "Macaca fascicularis",
         "10090" = "Mus musculus",
         "10116" = "Rattus norvegicus",
         "9986" = "Oryctolagus cuniculus",
         "36483" = "hamster",
         "9615" = "Canis lupus familiaris",
         "9823" = "Sus scrofa",
         "9031" = "Gallus gallus",
         "7955" = "Danio rerio",
         NA_character_)
}

biocBabelgeneOrthologs <- function(genes, species, multiOrth = FALSE) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    stop("Package 'babelgene' is required for non-human ortholog mapping in backend='bioconductor'.",
         call. = FALSE)
  }
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(data.frame())
  }
  out <- tryCatch(
    babelgene::orthologs(genes = genes,
                         species = species,
                         human = FALSE,
                         top = !multiOrth),
    error = function(e) data.frame()
  )
  out
}

biocSelectSafe <- function(orgdb, keys, keytype, columns) {
  keys <- unique(as.character(keys))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (length(keys) == 0) {
    return(data.frame())
  }
  validKeys <- tryCatch(AnnotationDbi::keys(orgdb, keytype = keytype),
                        error = function(e) character())
  keep <- intersect(keys, validKeys)
  if (length(keep) == 0) {
    return(data.frame())
  }
  AnnotationDbi::select(orgdb,
                        keys = keep,
                        keytype = keytype,
                        columns = columns)
}

biocAnnotateGeneIDs <- function(geneIds) {
  validIDs <- suppressWarnings(as.integer(geneIds))
  validIDs <- unique(validIDs[!is.na(validIDs)])
  if (length(validIDs) == 0) {
    return(data.frame())
  }

  resList <- list()
  for (taxId in biocSupportedOrgDbTaxIds()) {
    orgdb <- biocRequireOrgDbByTaxId(taxId)
    keys <- as.character(validIDs)
    dbRes <- biocSelectSafe(orgdb,
                keys = keys,
                keytype = "ENTREZID",
                columns = c("ENTREZID", "SYMBOL", "GENENAME"))
    if (nrow(dbRes) == 0) {
      next
    }
    dbRes <- dbRes[!is.na(dbRes$ENTREZID), , drop = FALSE]
    if (nrow(dbRes) == 0) {
      next
    }
    dbRes$TaxID <- taxId
    resList[[as.character(taxId)]] <- data.frame(
      GeneID = suppressWarnings(as.integer(dbRes$ENTREZID)),
      GeneSymbol = dbRes$SYMBOL,
      Description = dbRes$GENENAME,
      TaxID = dbRes$TaxID,
      Type = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (length(resList) == 0) {
    return(data.frame())
  }
  dplyr::bind_rows(resList) %>% unique()
}

biocAnnotateGeneSymbols <- function(symbols, taxId) {
  orgdb <- biocRequireOrgDbByTaxId(taxId)
  keys <- unique(as.character(symbols[!is.na(symbols)]))
  if (length(keys) == 0) {
    return(data.frame())
  }
  dbRes <- biocSelectSafe(orgdb,
                          keys = keys,
                          keytype = "SYMBOL",
                          columns = c("SYMBOL", "ENTREZID", "GENENAME"))
  if (nrow(dbRes) == 0) {
    return(data.frame())
  }
  dbRes <- dbRes[!is.na(dbRes$SYMBOL), , drop = FALSE]
  if (nrow(dbRes) == 0) {
    return(data.frame())
  }
  data.frame(
    GeneID = suppressWarnings(as.integer(dbRes$ENTREZID)),
    GeneSymbol = dbRes$SYMBOL,
    Description = dbRes$GENENAME,
    TaxID = as.integer(taxId),
    Type = NA_character_,
    stringsAsFactors = FALSE
  ) %>% unique()
}

biocAnnotateEnsemblGeneIDs <- function(ids) {
  EnsemblID <- GeneID <- GeneSymbol <- Description <- TaxID <- NULL

  ids <- as.character(ids)
  uvids <- removeEnsemblVersion(ids)
  taxVec <- biocEnsemblTaxId(uvids)

  resList <- list()
  for (taxId in unique(taxVec[!is.na(taxVec)])) {
    idx <- which(taxVec == taxId)
    keys <- unique(uvids[idx])
    orgdb <- biocRequireOrgDbByTaxId(taxId)
    dbRes <- biocSelectSafe(orgdb,
                keys = keys,
                keytype = "ENSEMBL",
                columns = c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"))
    if (nrow(dbRes) == 0) {
      next
    }
    dbRes <- dbRes[!is.na(dbRes$ENSEMBL), , drop = FALSE]
    if (nrow(dbRes) == 0) {
      next
    }
    dbRes$TaxID <- taxId
    resList[[as.character(taxId)]] <- data.frame(
      UVID = dbRes$ENSEMBL,
      GeneID = suppressWarnings(as.integer(dbRes$ENTREZID)),
      GeneSymbol = dbRes$SYMBOL,
      Description = dbRes$GENENAME,
      TaxID = dbRes$TaxID,
      stringsAsFactors = FALSE
    )
  }
  if (length(resList) == 0) {
    return(data.frame(EnsemblID = ids,
                      GeneID = NA,
                      GeneSymbol = NA,
                      Description = NA,
                      TaxID = NA,
                      stringsAsFactors = FALSE))
  }

  input <- data.frame(EnsemblID = ids, UVID = uvids, stringsAsFactors = FALSE)
  res <- dplyr::bind_rows(resList) %>%
    unique() %>%
    dplyr::left_join(input, by = "UVID") %>%
    dplyr::select(EnsemblID, GeneID, GeneSymbol, Description, TaxID)
  res <- sortAnnotationByQuery(res, ids, "EnsemblID", multi = FALSE)
  res
}

biocAnnotateUniProt <- function(accessions) {
  acc <- unique(as.character(accessions[!is.na(accessions)]))
  if (length(acc) == 0) {
    return(data.frame())
  }
  resList <- list()
  for (taxId in biocSupportedOrgDbTaxIds()) {
    orgdb <- biocRequireOrgDbByTaxId(taxId)
    dbRes <- biocSelectSafe(orgdb,
                keys = acc,
                keytype = "UNIPROT",
                columns = c("UNIPROT", "ENTREZID"))
    if (nrow(dbRes) == 0) {
      next
    }
    dbRes <- dbRes[!is.na(dbRes$UNIPROT), , drop = FALSE]
    if (nrow(dbRes) == 0) {
      next
    }
    resList[[as.character(taxId)]] <- data.frame(
      Accession = dbRes$UNIPROT,
      EntryName = NA_character_,
      GeneID = suppressWarnings(as.integer(dbRes$ENTREZID)),
      TaxID = taxId,
      EnsemblGeneID = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  if (length(resList) == 0) {
    return(data.frame())
  }
  dplyr::bind_rows(resList) %>% unique()
}

biocAnnotateHumanOrthologs <- function(geneids, multiOrth = FALSE) {
  geneAnno <- annotateGeneIDsWithoutHumanOrtholog(geneids, backend = "bioconductor")
  app <- biocAppendHumanOrthologs(geneAnno, multiOrth = multiOrth)
  sortAnnotationByQuery(app[, c("GeneID", "TaxID", "HumanGeneID"), drop = FALSE],
                        geneids,
                        id_column = "GeneID",
                        multi = multiOrth)
}

biocAppendHumanOrthologs <- function(anno, multiOrth = FALSE) {
  Description <- GeneID <- GeneSymbol <- TaxID <- Type <- NULL
  .rowid <- HumanGeneID.mapped <- HumanGeneSymbol.mapped <- HumanGeneSymbolAnno <- NULL

  res <- anno
  res$.rowid <- seq_len(nrow(res))

  if (!"GeneSymbol" %in% colnames(res)) {
    res$GeneSymbol <- NA_character_
  }

  # Self-map known human genes
  isHuman <- !is.na(res$TaxID) & as.character(res$TaxID) == "9606"
  res$HumanGeneID <- NA
  res$HumanGeneSymbol <- NA_character_
  res$HumanDescription <- NA_character_
  res$HumanType <- NA_character_
  res$HumanGeneID[isHuman] <- res$GeneID[isHuman]
  res$HumanGeneSymbol[isHuman] <- as.character(res$GeneSymbol[isHuman])

  # Map non-human genes to human orthologs via babelgene.
  nonHuman <- res[!isHuman & (!is.na(res$GeneID) | !is.na(res$GeneSymbol)),
                  c(".rowid", "GeneID", "GeneSymbol", "TaxID"), drop = FALSE]
  orthoParts <- list()
  if (nrow(nonHuman) > 0) {
    for (tax in unique(nonHuman$TaxID)) {
      species <- biocTaxIdToBabelSpecies(tax)
      if (is.na(species)) {
        warning("No babelgene species mapping configured for taxId ", tax,
                "; orthologs set to NA for these rows.")
        next
      }

      sub <- nonHuman[as.character(nonHuman$TaxID) == as.character(tax), , drop = FALSE]
      sub$key <- ifelse(!is.na(sub$GeneID), as.character(sub$GeneID), as.character(sub$GeneSymbol))
      sub$source <- ifelse(!is.na(sub$GeneID), "entrez", "symbol")

      for (src in c("entrez", "symbol")) {
        srcSub <- sub[sub$source == src & !is.na(sub$key), c(".rowid", "key"), drop = FALSE]
        if (nrow(srcSub) == 0) {
          next
        }
        bg <- biocBabelgeneOrthologs(srcSub$key, species = species, multiOrth = multiOrth)
        if (nrow(bg) == 0) {
          next
        }
        if (!(src %in% colnames(bg))) {
          next
        }
        map <- data.frame(key = as.character(bg[[src]]),
                          HumanGeneID = suppressWarnings(as.integer(bg$human_entrez)),
                          HumanGeneSymbol = as.character(bg$human_symbol),
                          stringsAsFactors = FALSE)
        map <- map[!is.na(map$key), , drop = FALSE]
        if (nrow(map) == 0) {
          next
        }
        orthoParts[[paste0(tax, "_", src)]] <- dplyr::left_join(srcSub, map, by = "key")
      }
    }
  }

  if (length(orthoParts) > 0) {
    orthoDf <- dplyr::bind_rows(orthoParts)
    res <- dplyr::left_join(res,
                            orthoDf[, c(".rowid", "HumanGeneID", "HumanGeneSymbol"), drop = FALSE],
                            by = ".rowid",
                            suffix = c("", ".mapped"))
    mappedID <- suppressWarnings(as.integer(res$HumanGeneID.mapped))
    hasMappedID <- !is.na(mappedID)
    res$HumanGeneID[hasMappedID] <- mappedID[hasMappedID]
    mappedSym <- as.character(res$HumanGeneSymbol.mapped)
    hasMappedSym <- !is.na(mappedSym) & mappedSym != ""
    res$HumanGeneSymbol[hasMappedSym] <- mappedSym[hasMappedSym]
    res <- res %>% dplyr::select(-HumanGeneID.mapped, -HumanGeneSymbol.mapped)
  }

  # Fill human annotation columns from human GeneIDs.
  humanIDs <- unique(suppressWarnings(as.integer(res$HumanGeneID)))
  humanIDs <- humanIDs[!is.na(humanIDs)]
  if (length(humanIDs) > 0) {
    humanAnno <- annotateGeneIDsWithoutHumanOrtholog(humanIDs,
                                                     backend = "bioconductor") %>%
      dplyr::select(HumanGeneID = GeneID,
                    HumanGeneSymbolAnno = GeneSymbol,
                    HumanDescription = Description,
                    HumanType = Type)
    res <- dplyr::left_join(res, humanAnno, by = "HumanGeneID")
    needSym <- is.na(res$HumanGeneSymbol) | res$HumanGeneSymbol == ""
    res$HumanGeneSymbol[needSym] <- res$HumanGeneSymbolAnno[needSym]
    res <- res %>% dplyr::select(-HumanGeneSymbolAnno)
  }

  res <- res %>% dplyr::select(-.rowid)
  unique(res)
}

biocAnnotateAnyIDsFeatureMap <- function(ids) {
  ids <- as.character(ids)
  res <- data.frame(Input = ids,
                    GeneID = NA,
                    TaxID = NA,
                    IDType = NA,
                    stringsAsFactors = FALSE)

  for (i in seq_along(ids)) {
    id <- ids[i]
    if (is.na(id) || id == "") {
      next
    }
    if (likeGeneID(id)) {
      res$GeneID[i] <- suppressWarnings(as.integer(id))
      res$IDType[i] <- "GeneID"
      ganno <- annotateGeneIDsWithoutHumanOrtholog(id, backend = "bioconductor")
      res$TaxID[i] <- ganno$TaxID[1]
    } else if (likeEnsemblGeneID(id)) {
      res$IDType[i] <- "EnsemblGeneID"
      eanno <- annotateEnsemblGeneIDsWithoutHumanOrtholog(id, backend = "bioconductor")
      res$GeneID[i] <- eanno$GeneID[1]
      res$TaxID[i] <- eanno$TaxID[1]
    } else if (likeGeneSymbol(id)) {
      res$IDType[i] <- "GeneSymbol"
      sanno <- annotateGeneSymbolsWithoutHumanOrtholog(id,
                                                       taxId = 9606,
                                                       backend = "bioconductor")
      res$GeneID[i] <- sanno$GeneID[1]
      res$TaxID[i] <- sanno$TaxID[1]
    } else if (likeUniProt(id)) {
      res$IDType[i] <- "UniProt"
      uanno <- biocAnnotateUniProt(id)
      if (nrow(uanno) > 0) {
        res$GeneID[i] <- uanno$GeneID[1]
        res$TaxID[i] <- uanno$TaxID[1]
      }
    }
  }

  res
}