## format IN syntax

#' Formatting a vector for SQL SELECT query with IN syntax
#' 
#' Prepare a vector for SQL SELECT query with the IN syntax
#' 
#' 
#' @param x A vector to be queried with the IN syntax
#' @return A character string to be used after IN. See examples.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @examples
#' 
#' myvec <- c("HH", "HM", "TH")
#' formatIn(myvec)
#' mysel <- "SELECT * FROM table WHERE city IN"
#' paste(mysel,formatIn(myvec))
#' 
#' 
#' @export formatIn
formatIn <- function(x) paste("(",paste("'", x, "'", sep="", collapse=","),")", sep="")



