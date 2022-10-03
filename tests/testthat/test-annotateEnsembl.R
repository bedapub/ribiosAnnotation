library(ribiosAnnotation)
library(testthat)

ids <- c("ENSG00000197535.1", "ENSG00000105221.2", "ENSG00000112062", NA, NULL)
idsWov <- removeEnsemblVersion(ids)

test_that("removeEnsemblVersion works", {
  expect_identical(idsWov,
                   c("ENSG00000197535", "ENSG00000105221", "ENSG00000112062", NA, NULL))
})

test_that("annotateEnsembl works without version number", {
  ids <- c("ENSG00000197535", "ENSG00000105221", NA,
               "ENSG00000112062", "ENST00000229795", "ENSP00000229795")
  skip('updating annotation backend from Oracle to MongoDB')
  res <- annotateEnsembl(ids)
  expect_identical(res$EnsemblID, ids)
  expect_identical(res$GeneID, c("4644", "208", NA, "1432", "1432", "1432"))
})

test_that("annotateEnsembl works with version number", {
  ids <- c("ENSMUSG00000031710.1", 
           "ENSMUSG00000031710",
           "ENSMUSG00000041272.1")
  skip('updating annotation backend from Oracle to MongoDB')
  res <- annotateEnsembl(ids)
  expect_identical(res$EnsemblID, ids)
  expect_identical(res$GeneID, c("22227", "22227", "252838"))
})
