library(testthat)
library(ribiosAnnotation)

options(
	ribiosAnnotation.backend = Sys.getenv(
		"RIBIOS_ANNOTATION_BACKEND",
		unset = "bioconductor"
	)
)

test_check("ribiosAnnotation")
