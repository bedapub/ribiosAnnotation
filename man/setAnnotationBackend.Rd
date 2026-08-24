% Generated manually
\name{setAnnotationBackend}
\alias{setAnnotationBackend}
\title{Set Global Annotation Backend}
\usage{
setAnnotationBackend(backend = c("mongodb", "bioconductor"))
}
\arguments{
\item{backend}{Character string, one of \code{"mongodb"} or \code{"bioconductor"}.}
}
\value{
\code{setAnnotationBackend} invisibly returns the selected backend.
}
\description{
Set the default annotation backend for the current R session.
When annotation functions are called without an explicit \code{backend}
argument, this setting is used.
}
\examples{
setAnnotationBackend("bioconductor")
setAnnotationBackend("mongodb")
}
