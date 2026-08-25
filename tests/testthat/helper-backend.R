skip_if_no_mongodb_backend <- function() {
  testthat::skip_if_not(
    ribiosAnnotation::mongodbConnectionAvailable(),
    "MongoDB backend not requested or unavailable"
  )
}
