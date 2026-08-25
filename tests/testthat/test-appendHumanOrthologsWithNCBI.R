library(ribiosAnnotation)
library(testthat)

test_that("appendHumanOrthologsWithNCBI works", {
  skip_if_no_mongodb_backend()

  anno <- data.frame(GeneID=c(780, 1506, 114483548, 102129055, NA),
                    TaxID=c(9606, 9606, 10116, 9541, NA))
  annoHoApp <- appendHumanOrthologsWithNCBI(anno)

  testthat::expect_equal(annoHoApp$GeneID, anno$GeneID)
  testthat::expect_equal(annoHoApp$TaxID, anno$TaxID)

  testthat::expect_true(all(c("HumanGeneID", "HumanGeneSymbol") %in% colnames(annoHoApp)))
  testthat::expect_equal(annoHoApp$HumanGeneID[1:2], c(780, 1506))
  testthat::expect_equal(annoHoApp$HumanGeneSymbol[1:2], c("DDR1", "CTRL"))
  testthat::expect_true(is.na(annoHoApp$HumanGeneID[3]))
  testthat::expect_true(is.na(annoHoApp$HumanGeneSymbol[3]))
  testthat::expect_true(is.na(annoHoApp$HumanGeneID[5]))
  testthat::expect_true(is.na(annoHoApp$HumanGeneSymbol[5]))

  if ("HumanType" %in% colnames(annoHoApp)) {
    testthat::expect_equal(annoHoApp$HumanType[1:2], c("protein-coding", "protein-coding"))
  }
})
