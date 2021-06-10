RIBIOS_TMP_TBL <- "RIBIOS_ID_TMP"
RIBIOS_JDBC_TMP_TBL <- "RIBIOS_JDBC_ID_TMP"

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
tmpTbl <- function(forceJDBC=FALSE) ifelse(hasOracle() & !forceJDBC, RIBIOS_TMP_TBL, RIBIOS_JDBC_TMP_TBL)



#' Query database instances with prepared SQL commands
#' 
#' The function takes one SQL command (probably a SELECT query) and returns the
#' results. It takes care of database connection, sending queries, fetching
#' results, and disconnection.
#' 
#' The \emph{sqlComm} parameter must not end with comma (;).
#' 
#' The function stops at SQL error by printing error messages.
#' 
#' @param sqlComm A SQL query, a SELECT command.
#' @param db Characters, which Oracle database should be connected?
#' @param user Characters, user name for login
#' @param password Characters, password for login
#' @param forceJDBC Logical, forcing the user of JDBC interface to fetch data
#' from Oracle servers. By default it is set to \code{FALSE}, since the JDBC
#' interface is slower than native ROracle implementation. This option is
#' largely for debugging purposes.
#' @return A data frame queried by the user
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso Internally the function calls following comamnds:
#' \code{dbSendQuery}, \code{dbClearResult} and \code{dbDisconnect}.
#' @references \url{http://bioinfo.bas.roche.com:8080/bioracle/}
#' @examples
#' 
#' options(error=utils::recover)
#' querydb("SELECT * FROM genome_sequence WHERE DB='HUMANN'", db=dbName(),
#' user="genome", password="genome")
#' options(error=NULL)
#' 
#' @export querydb
querydb <- function(sqlComm, db=dbName(), user=biaroUser(), password=biaroPwd(), forceJDBC=FALSE) {
  isORA <- hasOracle() & !forceJDBC
  con <- ribiosCon(db=db, user=user, password=password, forceJDBC=forceJDBC)
  rs <- dbSendQuery(con, sqlComm)
  if(isORA) {
    while (!dbHasCompleted(rs))
      ann <- fetch(rs, n = -1)
  } else {
    ann <- fetch(rs, n=-1)
  }
  dbClearResult(rs)
  dbDisconnect(con)
  ann
}


## select in: large IN queries


#' Query database by selecting with SQL IN syntax
#' 
#' The function builds SELECT SQL commands with IN syntax. By default the
#' Oracle SQL server supports IN syntax up to 1000 items. The function handles
#' such cases by splitting the query. Therefore it is capable to let users
#' perform SELECT operation with IN syntax without worrying about the number of
#' items.
#' 
#' The function builds SELECT SQL command with IN syntax, and splits the
#' database query if needed. This is necessary when there are over 1000 IN
#' items.
#' 
#' @param sqlComm Character string. It must be a \sQuote{SELECT} clause, and
#' should be complete without the \dQuote{IN (sQuote...)} part. See example
#' below.
#' @param inCol Character string, the input items are in which column?
#' @param inValues Character vector, input items
#' @param db Database
#' @param user User name
#' @param password User password
#' @param forceJDBC Logical, forcing the user of JDBC interface to fetch data
#' from Oracle servers. By default it is set to \code{FALSE}, since the JDBC
#' interface is slower than native ROracle implementation. This option is
#' largely for debugging purposes.
#' @return As \code{\link{querydb}}, a \code{data.frame} of the query results.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{querydb}} to query databases.
#' @examples
#' 
#' options(error=utils::recover)
#' querydbSelectIn(sqlComm="SELECT * FROM genome_sequence WHERE DB='HUMANN' AND ",
#'                 inCol="SEQ", inValues=c("CHR1", "CHR5", "CHR16", "CHRY"),
#'                 db=dbName(), user="genome", password="genome")
#' 
#' options(error=NULL)
#' 
#' @export querydbSelectIn
querydbSelectIn <- function(sqlComm, inCol, inValues,
                            db=dbName(), user=biaroUser(), password=biaroPwd(),
                            forceJDBC=FALSE) {
  isORA <- hasOracle() & !forceJDBC
  inValues <- unique(as.character(inValues))
  if (length(inValues) <= oracleInNmax()) {
    state <- paste(sqlComm, inCol, "IN", formatIn(inValues), 
                   collapse = " ")
    querydb(state, db = db, user = user, password = password, forceJDBC=forceJDBC)
  } else {
    con <- ribiosCon(db = db, user = user, password = password, forceJDBC=forceJDBC)
    nob <- ceiling(length(inValues)/oracleInNmax())
    res <- vector("list", nob)
    for (i in 1:nob) {
      ind <- seq((i - 1) * oracleInNmax() + 1L, i * oracleInNmax())
      state <- paste(sqlComm, inCol, "IN", formatIn(inValues[ind]))
      rs <- dbSendQuery(con, state)
      if(isORA) {
        while (!dbHasCompleted(rs))
          ann <- fetch(rs, n = -1)
      } else {
        ann <- fetch(rs, n=-1)
      }
      dbClearResult(rs)
      res[[i]] <- ann
    }
    dbDisconnect(con)
    do.call(rbind, res)
  }
}


## The temporary table is created one a new Oracle() instance (e.g. a new R session) is made
## The table content is private to the session
## ON OCMMIT DELETE ROWS: the rows are cleared after dbCommit() is run
fillOneColTmpTbl <- function(con,  values) {
  values <- as.character(values)
  values[is.na(values)] <- "NA"
  values[values==""] <- "NA"
  if(any(nchar(values)>100)) {
    warning("Identifiers longer than 100 characters are truncated!")
  }
  values <- substr(values, 1, 100)
  isORA <- inherits(con, "OraConnection")
  if(isORA) {
    if (!dbExistsTable(con, RIBIOS_TMP_TBL)) {
      state <- paste("CREATE GLOBAL TEMPORARY TABLE", RIBIOS_TMP_TBL, 
                     "(ID VARCHAR2(100) NOT NULl PRIMARY KEY) ON COMMIT DELETE ROWS")
      rs <- dbSendQuery(con, state)
    }
    inputDf <- data.frame(ID = values)
    state2 <- paste("insert into", RIBIOS_TMP_TBL, "(ID) values (:1)")
    rs <- dbSendQuery(con, state2, data = inputDf)
    return(dbHasCompleted(rs))
  } else {
    if(!dbExistsTable(con, RIBIOS_JDBC_TMP_TBL)) {
      state <- paste("CREATE GLOBAL TEMPORARY TABLE", RIBIOS_JDBC_TMP_TBL, 
                     "(ID VARCHAR2(100) NOT NULl PRIMARY KEY) ON COMMIT PRESERVE ROWS")
      rs <- RJDBC::dbSendUpdate(con, state)
    }
    state2 <- paste("insert into",RIBIOS_JDBC_TMP_TBL, " (ID) values (?)")
    ## TODO: SLOW: batch insert is desired
    for(i in seq(along=values))
      rs <- RJDBC::dbSendUpdate(con, state2, values[i])
    return(TRUE)
  }
}

## querydbTmpTbl shows principles of using temporary table. The SQL building is not finished: currently it only supports WHERE-free syntax


#' Query database with single-column temporary table
#' 
#' Low-level function to query database with filling single-column temporary
#' table with user inptu values.
#' 
#' This function is intended to be used by developer or advanced users who wish
#' to access the database directly. End-users should not use this function
#' unless 100\% clear what he or she is doing.
#' 
#' This function uses temporary table to perform query with user input. This
#' can be alternatively done with \code{\link{querydbSelectIn}}. Howevever
#' \code{\link{querydbSelectIn}} has the limitation of 1000 elements per SQL
#' query, and the current implementation (splitting input values in
#' 1000-element blocks) can be less efficient compared to \code{queryTmpTbl}.
#' 
#' @param sqlComm Character string of SQL command. The command can be written
#' normally as if no temporary table is used.
#' @param inCol Which column stated in the SQL command should be joined with
#' the temporary table values?
#' @param inValues Character vector, values to be filled into the temporary
#' table
#' @param db Database
#' @param user User name
#' @param password User password
#' @param forceJDBC Logical, forcing the user of JDBC interface to fetch data
#' from Oracle servers. By default it is set to \code{FALSE}, since the JDBC
#' interface is slower than native ROracle implementation. This option is
#' largely for debugging purposes.
#' @return A \code{data.frame} depending on the query.
#' @author Jitao David Zhang <jitao_david.zhang@@roche.com>
#' @seealso \code{\link{querydbSelectIn}}
#' @examples
#' 
#' options(error=utils::recover)
#' hcIn <- querydbTmpTbl("SELECT * FROM genome_sequence WHERE DB='HUMANN' ",
#'                       inCol="SEQ", inValues=c("CHR1", "CHR5", "CHRX"),
#'                       db=dbName(), user="genome", password="genome")
#' options(error=NULL)
#' 
#' @export querydbTmpTbl
querydbTmpTbl <- function(sqlComm, inCol, inValues,
                          db=dbName(), user=biaroUser(), password=biaroPwd(), 
                          forceJDBC=FALSE) {
  isORA <- hasOracle() & !forceJDBC
  con <- ribiosCon(db=db, user=user, password=password, forceJDBC=forceJDBC)
  inValues <- setdiff(unique(as.character(inValues)), "")
  
  TMP_TBL <- tmpTbl(forceJDBC=forceJDBC)
  fillOneColTmpTbl(con = con, values = inValues)
  hasFrom <- grepl("from", sqlComm, ignore.case = TRUE)
  hasWhere <- grepl("where", sqlComm, ignore.case = TRUE)
  if (!hasFrom) 
    stop("Cannot find 'from' in the SQL command line\n")
  if (hasWhere) {
    state <- gsub("WHERE", paste(",", TMP_TBL, " t WHERE t.ID=", 
                                 inCol, " AND ", sep = ""), sqlComm, ignore.case = TRUE)
  } else {
    sqlComm <- gsub("FROM", paste("FROM ", TMP_TBL, 
                                  " t,", sep = ""), sqlComm, ignore.case = TRUE)
    state <- paste(sqlComm, " WHERE ", inCol, "=t.ID", sep = "")
  }
  rs <- dbSendQuery(con, state)
  if(isORA) {
    while (!dbHasCompleted(rs))
      ann <- fetch(rs, n = -1)
  } else {
    ann <- fetch(rs, n=-1)
  }
  dbClearResult(rs)
  dbDisconnect(con)
  ann
}


