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
  }
  
  return(jsonlite::read_json(path))
}

loadSecrets <- function(path=file.path("/pstore/apps/bioinfo",
                   "ribios/secrets",
                   "ribiosAnnotation-secrets.json")) {
    credentials <- parseSecrets(path)
    options("ribiosAnnotation"=list(credentials=credentials))
    return(invisible(credentials))
}

#' Quick access of credentials
#' @export
biaUser <- function() options()$ribiosAnnotation$credentials$bi$username

#' @describeIn biaUser BIA password
#' @export
biaPwd <- function() options()$ribiosAnnotation$credentials$bi$password

#' @describeIn biaUser BIA2 username
#' @export
bia2User <- function() options()$ribiosAnnotation$credentials$bi2$username

#' @describeIn biaUser BIA2 password
#' @export
bia2Pwd <- function() options()$ribiosAnnotation$credentials$bi2$password

#' @describeIn biaUser BIA READONLY username
#' @export
biaroUser <- function() options()$ribiosAnnotation$credentials$biaro$username

#' @describeIn biaUser BIA READONLY password
#' @export
biaroPwd <- function() options()$ribiosAnnotation$credentials$biaro$password

#' @describeIn biaUser BIN username
#' @export 
binUser <- function() options()$ribiosAnnotation$credentials$bin$username

#' @describeIn biaUser BIN password
#' @epxport
binPwd <- function() options()$ribiosAnnotation$credentials$bin$password

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
ribiosCon <- function(db="bia", user=biaroUser(), password=ORACLE.BIARO.PWD, forceJDBC=FALSE) {
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
newconBIA <- function() ribiosCon(db="bia", user=biaUser(), password=biaPwd())
newconBIA2 <- function() ribiosCon(db="bia", user=bia2User(), password=bia2Pwd())
newconBIARO <- function() ribiosCon(db="bia", user=biaroUser(), password=biaroPwd())
newconBIN <- function() ribiosCon(db="bia", user=binUser(), password=binPwd())

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
