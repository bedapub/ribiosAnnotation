% Generated manually
\name{setAnnotationBackend}
\alias{setAnnotationBackend}
\alias{getAnnotationBackend}
\title{Set or Get Global Annotation Backend}
\usage{
setAnnotationBackend(backend = c("mongodb", "bioconductor"))

getAnnotationBackend()
}
\arguments{
\item{backend}{Character string, one of \code{"mongodb"} or \code{"bioconductor"}.}
}
\value{
\code{setAnnotationBackend} invisibly returns the selected backend.

\code{getAnnotationBackend} returns the effective backend as a character string.
}
\description{
Set or get the default annotation backend for the current R session.
When annotation functions are called without an explicit \code{backend}
argument, this setting is used.
}
\examples{
setAnnotationBackend("bioconductor")
getAnnotationBackend()
setAnnotationBackend("mongodb")
}
