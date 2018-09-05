#' @title Show all supported DBs
#'
#' @description Returns a table containing all supported DBs with DB name and ID. The ID can
#' used for example as targetDbs parameter in the \code{query} and \code{convertFromTo} function.
#'
#' @return A \preformatted{data.frame} object with the columns DatabaseDisplayName and DatabaseID.
#'
#' @examples
#' \dontrun{
#'   getValidDatabases();
#' }
#'
#' @export
#' @seealso \code{\link{query}} and \code{\link{convertFromTo}}
getValidDatabases = function() {

              x = .postNeo4jRequest(paste("CALL rcsi.params.getValidDatabases()"));
              return(x);
          }

#' @title Available attributes
#'
#' @description Returns a table containing all supported attributes that can serve as
#' 'params' input parameter in the \code{\link{query}} function.
#'
#' @return A \preformatted{data.frame} object with the columns `DatabaseDisplayName` and `DatabaseID`.
#'
#' @examples
#' \dontrun{
#'   getValidGngAttributes();
#' }
#' @seealso \code{\link{query}}
getValidGngAttributes = function() {
              x = .postNeo4jRequest("CALL rcsi.params.getValidAttributes()");
              return(x);
          }

#' @title Available attributes
#'
#' @description Returns a table containing all supported attributes that can serve as
#' 'params' input parameter in the \code{\link{query}} function.
#'
#' @return A \preformatted{data.frame} object with the columns DatabaseDisplayName and DatabaseID.
#'
#' @examples
#' \dontrun{
#'   getValidMirnaMetadataValues();
#' }
#' @seealso \code{\link{query}}
#' @export
getValidMirnaMetadataValues = function() {
              x = .postNeo4jRequest("CALL rcsi.params.getValidMirnaMetadataValues()");
              return(x);
          }


#' @title Get highest supported miRBase release version information
#'
#' @examples
#' \dontrun{
#'   getCurrentMirbaseVersion();
#' }
#' @export
getCurrentMirbaseVersion = function() {
              x = .postNeo4jRequest("MATCH (db :MirnaDB {id: 'miRBase_mature_name'}) RETURN db.release AS release;")$release[1];
              # x1 = strsplit(x, ",")[[1]][1];
              # m <- regexpr("\\d+\\.?\\d+", x1, perl=TRUE)
              # releaseVersion = regmatches(strsplit(x, ",")[[1]][1], m);
              return(x);
}

#' @title Return the Ensembl DB version of GeneNameGenie.
#'
#' @examples
#' \dontrun{
#'   getEnsemblVersion();
#' }
#' @export
getEnsemblVersion = function() {
    x = .postNeo4jRequest("MATCH (db :EnsemblDB {id: 'ArrayExpress'}) RETURN db.release AS release;")$release[1];
    return(x);
}

