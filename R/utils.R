#' @rdname loadSecrets
#' @export
ribiosAnnotationSecretFile <- file.path("~",
			       ".credentials", "ribiosAnnotation-secrets.json")

#' Load secrets from a file and set them in options
#'
#' ribiosAnnotation needs to access databases to fetch annotations, the process
#' of which requires credentials for these databases. The package looks for a
#' file \sQuote{\code{~/.credentials/ribiosAnnotation-secrets.json}} which
#' contains the credentials. If this file is not found, no queries can be made.
#'
#' A template of the credential file is provided in
#' \sQuote{\code{secrets/secrets-template.json}} file of the package.
#'
#' @param path Path to the secret file
#'
#' @return The current options of \code{ribiosAnnotation}
#' The function writes the \code{credentials} field of the options
#' After running this function, database names and passwords can be accessed
#' @importFrom rjson fromJSON
#' @export
loadSecrets <- function(path=ribiosAnnotationSecretFile) {
  if (!file.exists(path)) {
    ## the message below must not be sent with warning, otherwise
    ## packages that depdend on ribiosAnnotation will not load properly.
    message("Can't find secret file: '", path, "'. Queries will not work.")
    path <- system.file("secrets/secrets-template.json", 
                        package="ribiosAnnotation")
  }
  secrets <- rjson::fromJSON(file=path)
  opts <- options("ribiosAnnotation")[[1]]
  opts$credentials <- secrets$credentials
  opts$connection <- secrets$connection
  options("ribiosAnnotation"=opts)
  return(invisible(opts))
}


## onload / onAttach
.onAttach <- function(libname, pkgname) {
  loadSecrets()
}

##----------------------------------------##
## MongoDB functions
##----------------------------------------##

#' Get secrets for MongoDB connections
#' @param file The secret JSON file.
#' @param instance String, which must be found under the \code{mongodb} section
#'  of the JSON file
#' @return A list of the following items:
#' \itemize{
#'   \item{\code{hostname}}{Hostname of the MongoDB}
#'   \item{\code{port}}{Port of the MongoDB}
#'   \item{\code{dbname}}{Database of the MongoDB}
#'   \item{\code{username}}{User name}
#'   \item{\code{password}}{Password}
#' }
#' @examples 
#' loadMongodbSecrets(instance="bioinfo_read")
#' \dontrun{
#'     loadMongodbSecrets(instance="decoy")
#' }
#' @export
loadMongodbSecrets <- function(file=ribiosAnnotationSecretFile,
                              instance="bioinfo_read") {
  secrets <-  rjson::fromJSON(file=file)
  dbSecrets <- secrets$mongodb[[instance]]
  if(is.null(dbSecrets)) {
    warning("Secrets for the instance '", instance, 
            "' was not found. Make sure that the secret file '", file, 
            "' contains a field 'mongodb' with the instance as a list.\n",
            "A decoy is returned. Subsequent queries will not work")
    dbSecrets <- list(hostname=NULL, port=NULL, dbname=NULL,
                      username=NULL, password=NULL)
  }
  return(dbSecrets)
}

#' Connect to a MongoDB instance
#' @param instance Character string, the MongoDB instance to connect to
#' @param collection Character string, the collection to be used
#' @return A pointer to a collection on the server, as returned by
#'  \code{\link[mongolite]{mongo}}.
#' @examples 
#' giCon <- connectMongoDB(instance="bioinfo_read",
#'                         collection="ncbi_gene_info")
#' @seealso \code{\link{loadMongodbSecrets}}
#' @export
connectMongoDB <- function(instance="bioinfo_read",
                           collection="ncbi_gene_info") {
  bioinfoReadSecrets <- loadMongodbSecrets(instance=instance)
  bioinfoReadURL <- sprintf("mongodb://%s:%s@%s:%s/%s?authSource=%s", 
                            bioinfoReadSecrets$username,
                            bioinfoReadSecrets$password,
                            bioinfoReadSecrets$hostname, 
                            bioinfoReadSecrets$port,
                            bioinfoReadSecrets$dbname,
                            bioinfoReadSecrets$dbname)
  giCon <- mongolite::mongo(collection=collection,
                            db=bioinfoReadSecrets$dbname,
                            url=bioinfoReadURL)
  return(giCon)
}

#' Construct a JSON string to indicate returned fields from a MongoDB query
#' @param fields A vector of character strings that should be included
#' @param include_id Logical, whether \code{_id} should be returned. Default 
#'   is \code{FALSE}
#' @return A JSON string that represents the fields to be returned
#' @examples 
#' returnFieldsJson(c("name", "birthday"))
#' returnFieldsJson(c("name", "birthday"), include_id=TRUE)
#' @importFrom rjson toJSON
#' @export
returnFieldsJson <- function(fields, include_id=FALSE) {
  logvec <- rep(TRUE, length(fields))
  names(logvec) <- fields
  if(!include_id) {
    logvec <- c(logvec, "_id"=FALSE)
  }
  res <- rjson::toJSON(logvec)
  return(res)
}


id2rownames <- function(ids) {
  if(identical(anyDuplicated(ids), 0L) & !any(is.na(ids))) {
    return(ids)
  } else {
    return(NULL)
  }
}
notValid <- function(x)  is.null(x) || is.na(x) || tolower(x)=="any" || x==""

