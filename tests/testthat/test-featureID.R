library(ribiosAnnotation)
library(testthat)

## test isValidFeatureID
test_that("isValidFeatureID and validFeatureIDs work with characters",
          {
            myIDs <- c("ACTB", NA, "MYC", "-", "", "C21XORF9")
            isMyIDvalid <- isValidFeatureID(myIDs)
            myValidIDs <- validFeatureIDs(myIDs)

            expect_identical(isMyIDvalid, 
                             c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
            expect_identical(myValidIDs,
                             c("ACTB", "MYC", "C21XORF9"))
          })


test_that("isValidFeatureID and validFeatureIDs work with factors",
          {
            myIDs <- factor(c("ACTB", NA, "MYC", "-", "", "C21XORF9"))
            isMyIDvalid <- isValidFeatureID(myIDs)
            myValidIDs <- validFeatureIDs(myIDs)
            expect_identical(isMyIDvalid,
                             c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE))
            expect_identical(myValidIDs,
                             factor(c("ACTB", "MYC", "C21XORF9")))
          })


test_that("isValidFeatureID and validFeatureIDs work with integers",
          {
            myNumberIDs <- c(3, 9, 12, NA)
            isMyNumberIDvalid <- isValidFeatureID(myNumberIDs)
            myValidNumberIDs <- validFeatureIDs(myNumberIDs)
            expect_identical(isMyNumberIDvalid,
                             c(TRUE, TRUE, TRUE, FALSE))
            expect_identical(myValidNumberIDs,
                             c(3, 9, 12))
          })

test_that("likeGeneID and similar functions work",
          {
             geneids <- c("1234", "345")
             invalidids <- c("-", "", NA)
             genesymbols <- c("LOX", "NKX-1", "CXorf21", "MT-ATP6")
             nonHumanGeneSymbols <- c("snail")
             uniprots <- c("A2BC19", "P12345", "A0A023GPI8")
             refseqs <- c("NM_000259", "NM_000259.3")
             ensembls <- c("ENSG00000197535", "ENST00000399231.7",
                           "ENSP00000418960.2")
             
             feats <- c(geneids, invalidids, genesymbols,
                        nonHumanGeneSymbols, uniprots, refseqs, ensembls)
             
             expect_identical(likeGeneID(feats),
                              feats %in% geneids)
             expect_identical(likeGeneSymbol(feats),
                              feats %in% c(genesymbols, nonHumanGeneSymbols))
             expect_identical(likeRefSeq(feats),
                              feats %in% refseqs) 
             expect_identical(likeEnsembl(feats),
                              feats %in% ensembls)
             expect_identical(likeUniProt(feats),
                              feats %in% uniprots)
             expect_identical(likeHumanGeneSymbol(feats),
                              feats %in% genesymbols)
})
