library(ribiosAnnotation)
library(testthat)

test_that("guessAndAnnotate works with gene symbols", {
  gs <- c("AKT1", "AKT2", "MAPK14")
  gsRes <- guessAndAnnotate(gs)
  expect_identical(gs, gsRes$GeneSymbol)
  expect_identical(as.character(c(207, 208, 1432)), gsRes$GeneID)
})

test_that("guessAndAnnotate works with Entrez Gene IDs", {
  gid <- as.character(c(1, 2, 14, 149))
  gidRes <- guessAndAnnotate(gid)
  expect_identical(gid, gidRes$GeneID)
})

test_that("guessAndAnnotate works with EnsemblIDs", {
  eid <- c("ENST00000613858.4", "ENST00000553916.5", "ENST00000399229.6")
  eidRes <- guessAndAnnotate(eid)
  expect_identical(eidRes$EnsemblID, eid)
  expect_identical(eidRes$GeneID, as.character(rep(4644, 3)))
})

test_that("guessAndAnnotate works with UniProtIDs", {
  uid <- c("O60583", "P05997", "Q7Z624")
  uidRes <- guessAndAnnotate(uid)
  expect_identical(uid, uidRes$Input)
  expect_identical(uidRes$GeneID, as.character(c(905, 1290, 79823)))
})

test_that("guessAndAnnotate returns an empty data.frame in case of unknown identifiers", {
  ukid <- c("CM000677.2", "AB003434.2")
  ukidRes <- guessAndAnnotate(ukid)
  expect_identical(ukidRes, data.frame(row.names=ukid))
})