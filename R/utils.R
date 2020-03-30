#' @importFrom DBI dbConnect dbExistsTable dbSendQuery dbHasCompleted  fetch dbClearResult dbDisconnect
#' @importFrom ribiosUtils matchColumn putColsFirst

##--------------------##
## crendentials
##--------------------##
#' @importFrom jsonlite read_json
parseSecrets <- function(
    path=file.path("/pstore/apps/bioinfo",
                   "ribios/secrets",
                   "ribiosAnnotation-secrets.json")) {
  if (!file.exists(path)) {
    warning("Can't find secret file: '", path, "'.")
    path <- system.file("secrets/secrets-template.json", 
                        package="ribiosAnnotation")
  
  return(jsonlite::read_json(path))
}

loadSecrets <- function(path=file.path("/pstore/apps/bioinfo",
                   "ribios/secrets",
                   "ribiosAnnotation-secrets.json")) {
    credentials <- parseSecrets(path)
    assign(ORACLE.BIA.USER, credentials$bia$username, pos=sys.frame())
    assign(ORACLE.BIA.PWD, credentials$bia$password, pos=sys.frame())
    assign(ORACLE.BIA2.USER, credentials$bia2$username, pos=sys.frame())
    assign(ORACLE.BIA2.PWD, credentials$bia2$password, pos=sys.frame())
    assign(ORACLE.BIARO.USER, credentials$biread$username, pos=sys.frame())
    assign(ORACLE.BIARO.PWD, credentials$biread$password, pos=sys.frame())
    assign(ORACLE.BIN.USER, credentials$bin$username, pos=sys.frame())
    assign(ORACLE.BIN.PWD, credentials$bin$password, pos=sys.frame())
}

##--------------------##
## constants
##--------------------##
## ORACLE.LIB <- ":/opt/oracle/client/10/run_1/lib"
## maximum vector length in the IN syntax
#' @export ORACLE.IN.NMAX
ORACLE.IN.NMAX <- 1000L

## function to test whether Oracle is available
#' @export hasOracle
hasOracle <- function() {
  return(requireNamespace("ROracle", quietly=TRUE))
}

## onload / onAttach
## .onLoad <- function(libname, pkgname) {}
.onAttach <- function(libname, pkgname) {
  if(hasOracle()) {
    if(!"package:ROracle" %in% search())
      attachNamespace("ROracle")
    assign("ORA", ROracle::Oracle(), pos=sys.frame())
  } 
  loadSecrets()
}

## automatically establish a connection, depending on whether Oracle client is installed
#' @export ribiosCon
## TODO: 1g memorz at least is hard-coded. How to infer it?
ribiosCon <- function(db="bia", user=ORACLE.BIARO.USER, password=ORACLE.BIARO.PWD, forceJDBC=FALSE) {
  if(hasOracle() & !forceJDBC) {
    con <- dbConnect(ORA, user = user, password = password, db = db)
  } else {
    options(java.parameters = "-Xmx1g" ) ## increase the heap size before the RJDBC package is loaded
    suppressWarnings(suppressMessages(hasJDBC <- requireNamespace("RJDBC")))
    if(!hasJDBC)
      stop("No JDBC package installed: please run 'install.packages('RJDBC')' first and then load ribiosAnnotation again.")
    drv <- RJDBC::JDBC("oracle.jdbc.OracleDriver",
                       system.file("drivers", "ojdbc7.jar", package="ribiosAnnotation"))
    port <- switch(EXPR=db, bia=15210, bin=15001)
    str <- paste("jdbc:oracle:thin:", user, "/", password, "@orse04p-scan.kau.roche.com:", port, "/", db, ".kau.roche.com", sep="")
    con <- dbConnect(drv,str)
  }
  return(con)
}

## shortcuts for common connections
newconBIA <- function() ribiosCon(db="bia", user=ORACLE.BIA.USER, password=ORACLE.BIA.PWD)
newconBIA2 <- function() ribiosCon(db="bia", user=ORACLE.BIA2.USER, password=ORACLE.BIA2.PWD)
newconBIARO <- function() ribiosCon(db="bia", user=ORACLE.BIARO.USER, password=ORACLE.BIARO.PWD)
newconBIN <- function() ribiosCon(db="bia", user=ORACLE.BIN.USER, password=ORACLE.BIN.PWD)

##----------------------------------------##
## deprecated functions
##----------------------------------------##
biosCurrentGeneSymbol <- function(...) {
  .Deprecated("gtiChipAnnotation",
              package="ribiosAnnotation")
  gtiChipAnnotation(...)
}

raceChipAnnotation <- function(...) {
  .Deprecated("gtiChipAnnotation",
              package="ribiosAnnotation")
  gtiChipAnnotation(...)
}
