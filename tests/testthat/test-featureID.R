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

myNumberIDs <- c(3,9,12,NA)
isMyNumberIDvalid <- isValidFeatureID(myNumberIDs)
myValidNumberIDs <- validFeatureIDs(myNumberIDs)

test_that("isValidFeatureID and validFeatureIDs work with integers",
          {
            expect_identical(isMyNumberIDvalid, 
                             c(TRUE, TRUE, TRUE, FALSE))
            expect_identical(myValidNumberIDs,
                             c(3,9,12))
          })
