library(ribiosAnnotation)
library(testthat)

test_that("bioconductor backend annotates GeneIDs with stable shape", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  ids <- c("1", "2", NA, "999999999")
  res <- annotateGeneIDsWithoutHumanOrtholog(ids, backend = "bioconductor")

  expect_identical(res$GeneID, ids)
  expect_true(all(c("GeneID", "GeneSymbol", "Description", "TaxID", "Type") %in% colnames(res)))
  expect_identical(nrow(res), length(ids))
})

test_that("bioconductor backend annotates GeneSymbols", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  ids <- c("AKT1", "NOT_A_GENE")
  res <- annotateGeneSymbolsWithoutHumanOrtholog(ids, taxId = 9606, backend = "bioconductor")

  expect_identical(res$GeneSymbol, ids)
  expect_identical(nrow(res), length(ids))
  expect_true(all(c("GeneID", "GeneSymbol", "Description", "TaxID", "Type") %in% colnames(res)))
})

test_that("bioconductor backend annotates Ensembl GeneIDs", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  ids <- c("ENSG00000142208", NA)
  res <- annotateEnsemblGeneIDsWithoutHumanOrtholog(ids, backend = "bioconductor")

  expect_identical(res$EnsemblID, ids)
  expect_identical(nrow(res), length(ids))
  expect_true(all(c("EnsemblID", "GeneID", "GeneSymbol", "Description", "TaxID", "Type") %in% colnames(res)))
})

test_that("bioconductor backend works in annotateAnyIDs", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  ids <- c("1", "AKT1", "ENSG00000142208", NA)
  res <- annotateAnyIDs(ids, backend = "bioconductor")
  
  expect_identical(res$Input, ids)
  expect_identical(nrow(res), length(ids))
  expect_true(all(c("Input", "GeneID", "TaxID", "IDType") %in% colnames(res)))
})

test_that("both backends expose compatible core columns for GeneID annotation", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  mongo_ok <- TRUE
  mongo_res <- tryCatch({
    annotateGeneIDsWithoutHumanOrtholog(c(1, 2), backend = "mongodb")
  }, error = function(e) {
    mongo_ok <<- FALSE
    NULL
  })
  skip_if_not(mongo_ok)

  bioc_res <- annotateGeneIDsWithoutHumanOrtholog(c(1, 2), backend = "bioconductor")

  expect_identical(colnames(mongo_res), colnames(bioc_res))
  expect_identical(mongo_res$GeneID, bioc_res$GeneID)
})

test_that("bioconductor backend maps mouse symbols to human orthologs via babelgene", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("babelgene")

  ids <- c("Akt1", "Erbb2", "Tlr7")
  res <- annotateGeneSymbols(ids,
                             taxId = 10090,
                             orthologue = TRUE,
                             backend = "bioconductor")

  expect_identical(res$GeneSymbol, ids)
  expect_true(all(!is.na(res$HumanGeneID)))
  expect_true(all(!is.na(res$HumanGeneSymbol)))
})

test_that("bioconductor backend maps dog and pig symbols to human orthologs", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("babelgene")

  dog_ids <- c("AKT1", "ERBB2")
  pig_ids <- c("AKT1", "ERBB2")

  dog_res <- annotateGeneSymbols(dog_ids,
                                 taxId = 9615,
                                 orthologue = TRUE,
                                 backend = "bioconductor")
  pig_res <- annotateGeneSymbols(pig_ids,
                                 taxId = 9823,
                                 orthologue = TRUE,
                                 backend = "bioconductor")

  expect_identical(dog_res$GeneSymbol, dog_ids)
  expect_identical(pig_res$GeneSymbol, pig_ids)
  expect_true(all(!is.na(dog_res$HumanGeneSymbol)))
  expect_true(all(!is.na(pig_res$HumanGeneSymbol)))
})

test_that("global backend option applies when backend argument is omitted", {
  skip_if_not_installed("AnnotationDbi")
  skip_if_not_installed("org.Hs.eg.db")

  old_opt <- getOption("ribiosAnnotation.backend", default = NULL)
  on.exit(options(ribiosAnnotation.backend = old_opt), add = TRUE)

  setAnnotationBackend("bioconductor")
  expect_identical(getAnnotationBackend(), "bioconductor")

  res <- annotateGeneIDsWithoutHumanOrtholog(c("1", "2"))
  expect_identical(res$GeneID, c("1", "2"))
})
