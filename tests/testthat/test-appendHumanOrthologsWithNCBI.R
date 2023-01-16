library(ribiosAnnotation)
library(testthat)

anno <- data.frame(GeneID=c(780, 1506, 114483548, 102129055, NA),
                  TaxID=c(9606, 9606, 10116, 9541, NA))
annoHoApp <- appendHumanOrthologsWithNCBI(anno)

testthat::expect_equal(annoHoApp$GeneID, anno$GeneID)
testthat::expect_equal(annoHoApp$TaxID, anno$TaxID)
testthat::expect_equal(annoHoApp$HumanGeneID,
                       c(780, 1506, NA, 1, NA))
testthat::expect_equal(annoHoApp$HumanGeneSymbol,
                       c("DDR1", "CTRL", NA, "A1BG", NA))
testthat::expect_equal(annoHoApp$HumanType,
                       c("protein-coding", "protein-coding", NA, "protein-coding", NA))