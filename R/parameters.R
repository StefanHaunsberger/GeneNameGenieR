#' @title Show all supported DBs
#'
#' @description Returns a table containing all supported DBs with DB name and ID. The ID can
#' used for example as targetDbs parameter in the \code{query} and \code{convertFromTo} function.
#'
#' @param gng An instance of the GeneNameGenieR object (e.g. \code{gng = GeneNameGenieR()})
#'
#' @return A \preformatted{data.frame} object with the columns DatabaseDisplayName and DatabaseID.
#'
#' @examples
#' \dontrun{
#'   gng = GeneNameGenieR();
#'   getValidDatabases(gng);
#' }
#'
#' @export
#' @seealso \code{\link{query}} and \code{\link{convertFromTo}}
setGeneric("getValidDatabases", function(gng) standardGeneric("getValidDatabases"))
# @export
setMethod("getValidDatabases",
          signature(gng = "GeneNameGenieR"),
          function(gng) {

              x = .postNeo4jRequest(gng, paste("CALL rcsi.params.getValidDatabases() YIELD value",
                                          "RETURN value.DatabaseDisplayName AS DatabaseDisplayName, value.DatabaseId AS DatabaseId",
                                          "ORDER BY DatabaseDisplayName"));
              return(x);
          })

#' @title Available attributes
#'
#' @description Returns a table containing all supported attributes that can serve as
#' 'params' input parameter in the \code{\link{query}} function.
#'
#' @param gng An instance of the \code{GeneNameGenieR} object (e.g. \preformatted{gng = new("GeneNameGenieR")})
#'
#' @return A \preformatted{data.frame} object with the columns DatabaseDisplayName and DatabaseID.
#'
#' @examples
#' \dontrun{
#'   gng = GeneNameGenieR();
#'   getValidGngAttributes(gng);
#' }
#' @seealso \code{\link{query}}
# @export
setGeneric("getValidGngAttributes", function(gng) standardGeneric("getValidGngAttributes"))

setMethod("getValidGngAttributes",
          signature(gng = "GeneNameGenieR"),
          function(gng) {
              x = .postNeo4jRequest(gng, "CALL rcsi.params.getValidAttributes() YIELD value RETURN value AS Parameter");
              return(x);
          })

#' @title Available attributes
#'
#' @description Returns a table containing all supported attributes that can serve as
#' 'params' input parameter in the \code{\link{query}} function.
#'
#' @param gng An instance of the \code{GeneNameGenieR} object (e.g. \preformatted{gng = new("GeneNameGenieR")})
#'
#' @return A \preformatted{data.frame} object with the columns DatabaseDisplayName and DatabaseID.
#'
#' @examples
#' \dontrun{
#'   gng = GeneNameGenieR();
#'   getValidMirnaMetadataValues(gng);
#' }
#' @seealso \code{\link{query}}
#' @export
setGeneric("getValidMirnaMetadataValues", function(gng) standardGeneric("getValidMirnaMetadataValues"))
# @export
setMethod("getValidMirnaMetadataValues",
          signature(gng = "GeneNameGenieR"),
          function(gng) {
              x = .postNeo4jRequest(gng, "CALL rcsi.params.getValidMirnaMetadataValues() YIELD value RETURN value AS Parameter");
              return(x);
          })


#' @title Get highest supported miRBase release version information
#'
#' @param gng An instance of the \code{GeneNameGenieR} object (e.g. \preformatted{gng = GeneNameGenieR()})
#'
#' @examples
#' \dontrun{
#'   getCurrentMirbaseVersion(GeneNameGenieR());
#' }
#' @export
setGeneric("getCurrentMirbaseVersion", function(gng) standardGeneric("getCurrentMirbaseVersion"))
# @export
setMethod("getCurrentMirbaseVersion",
          signature(gng = "GeneNameGenieR"),
          function(gng) {
              x = .postNeo4jRequest(gng, "MATCH (db :MirnaDB {id: 'miRBase_mature_name'}) RETURN db.release AS release;")$release[1];
              # x1 = strsplit(x, ",")[[1]][1];
              # m <- regexpr("\\d+\\.?\\d+", x1, perl=TRUE)
              # releaseVersion = regmatches(strsplit(x, ",")[[1]][1], m);
              return(x);
          })
