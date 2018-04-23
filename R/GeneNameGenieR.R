#' @title GeneNameGenieR: Lightning fast graph based gene name converter
#' @author Stefan J. Haunsberger
#' @name GeneNameGenieR
#' @description
#'
#' @docType package

library(dplyr)
library(magrittr)

DEFAULT_URL = "http://localhost:7474/db/data/";
URL = DEFAULT_URL;
setOldClass("graph")

#' @title Instantiate from MiRNANameConverter class
#'
#' @description This function returns back an instance of a
#' \emph{MiRNANAmeConverter} object.
#'
#' @slot url Database connection
#' @slot graph Valid/Supported miRBase versions
#' @slot toAddDbCols Number of different organisms supported
#' @author Stefan J. Haunsberger
#'
#' For bug reports, please
#' \href{https://github.com/StefanHaunsberger/GeneNameGenieR/issues}{create a new GitHub issue}.
GeneNameGenieR = setClass("GeneNameGenieR",
    slots = list(
        baseUrl = "character",
        cypherEndpoint = "character",
        verbose = "logical",
        toAddDbCols = "character"
    ),
    prototype = list(
        url = DEFAULT_URL,
        cypherEndpoint = "cypher",
        verbose = FALSE,
        toAddDbCols = c(ArrayExpress = "EnsemblGeneId",
                        Ens_Hs_transcript = "EnsemblTranscriptId",
                        Ens_Hs_translation = "EnsemblProteinId")
    )
);
# #' @exportClass GeneNameGenieR

#' @param url : Graph database URL, such as "http://localhost:7474/db/data/"
#' @export
GeneNameGenieR = function(url = NA_character_) {
    if (is.na(url)) {
        url = DEFAULT_URL;
    }

    g = new("GeneNameGenieR", url = url);
    g@url = url;
    return(g);
}

##### initialization ####
# Initialization function (constructor)
# setMethod(
#     f = "initialize",
#     signature("GeneNameGenieR"),
#     definition = function(.Object, url) {
#
#         args = as.list(match.call());
#         url = "";
#         print(args)
#         if (is.na(args["url"])) {
#             url = DEFAULT_URL;
#         } else {
#             url = args["url"];
#         }
#         .Object@graph = RNeo4j::startGraph(url);
#
#         return(.Object);
#     }
# )

# ##### GeneNameGenieR - constructor ####
# #' @title GeneNameGenieR constructor
# #'
# #' @description This function returns an instance of a GeneNameGenieR class.
# #'
# #' @return an object of class 'GeneNameGenieR'
# #' @examples
# #' \dontrun{
# #' gng = GeneNameGenieR() # Instance of class 'GeneNameGenieR'
# #' }
# #' @details
# #' This function initializes an object of the class GeneNameGenieR It is a
# #' wrapper for \code{new()}.
# #'
# #' @param url : Graph database URL, such as "http://localhost:7474/db/data/"
# #'
# #' @seealso \code{\link{new}}
# #' @author Stefan J. Haunsberger
# #' @export GeneNameGenieR
# setMethod(
#     f = "GeneNameGenieR",
#     signature(),
#     definition = function(url = NA_character_) {
#         return(new('GeneNameGenieR', url));
#     }
# )
source("R/parameters.R")
source("R/utils.R")
# source("R/mirnameconversion.R")

#' @title Convert input id to its corresponding official gene symbol
#'
#' @description Compute the official gene symbol for a single or multiple `queryId`(s).
#' If `sourceDb` is missing, the database will automatically try to detect the queryId's
#'   source database. Per default `chromosomal` is set to `true`.
#'
#' @param gng A GeneNameGenieR object (e.g. \code{GeneNameGenieR()})
#' @param queryId A character vector of molecular identifiers that needs to be converted
#' @param targetDb A character vector of target databases (default equals to 'EntrezGene' and 'GeneSymbolDb' for
#' the Entrez gene ID and official gene symbol respectively)
#' @param sourceDb String containing the origin database of the input ID(s) (optional).
#' @param chromosomal Boolean parameter (default TRUE). Whether the search and output shall include features located on
#' chromosomes other than the standard 23 {1,...21, X, Y}. Ensembl also includes features that are located on so called
#' alternative human chromosomes, such as CHR_HSCHR5_6_CTG1. In the case where \preformatted{chromosomal} is \preformatted{FALSE}, not even input
#' identifiers that are associated with this location will be found.
#'
#' @importFrom magrittr %>%
#'
#' @export
setGeneric(
    "getOfficialGeneSymbol",
    signature = c("gng"),
    function(gng, queryId, sourceDb = NA_character_, chromosomal = TRUE) {
        standardGeneric("getOfficialGeneSymbol")
    })

#' @export
setMethod("getOfficialGeneSymbol",
          signature(gng = "GeneNameGenieR"),
  function(gng, queryId, sourceDb, chromosomal) {

  q = paste0("CALL rcsi.convert.table.getOfficialGeneSymbol(",
      ifelse(length(queryId) == 1, "[{ids}], ", "{ids}, "),
      ifelse(is.na(sourceDb), "null, ", "{sourceDb}, "),
      "{chromosomal}) YIELD value ",
      "RETURN ",
      "   value.InputId AS InputId, ",
      "   value.InputSourceDb AS InputSourceDb, ",
      "   value.OfficialGeneSymbol AS OfficialGeneSymbol");


  x = .postNeo4jRequest(gng, q,
                     ids = queryId,
                     sourceDb = sourceDb,
                     chromosomal = chromosomal);

  .postCheckInput(x);

  return(x);

})

#' @title Convert molecular identifier to target identifiers
#'
#' @description Converts molecular input identifier to identifier from target database. Thereby, the source database of the
#' input identifier does not have to be specified as the database detects the corresponding database automatically.
#'
#' @param gng A GeneNameGenieR object (e.g. \code{GeneNameGenieR()})
#' @param queryId A character vector of molecular identifiers that needs to be converted
#' @param targetDb A character vecotr of target databases (default equals to 'EntrezGene' and 'GeneSymbolDb' for
#' the Entrez gene ID and official gene symbol respectively)
#' @param sourceDb String containing the origin database of the input ID(s) (optional).
#' @param longFormat Boolean parameter (default TRUE). Indicates if returned dataframe is in the form of a longFormat
#' @param chromosomal Boolean parameter (default TRUE). Whether the search and output shall include features located on
#' chromosomes other than the standard 23 {1,...21, X, Y}. Ensembl also includes features that are located on so called
#' alternative human chromosomes, such as CHR_HSCHR5_6_CTG1. In the case where \preformatted{chromosomal} is \preformatted{FALSE}, not even input
#' identifiers that are associated with this location will be found.
#'
#' @return A data.frame object with the columns InputId, InputSourceDb and columns containing TargetDb and TargetId values.
#'
#' @export
#' @seealso \code{\link{query}}
setGeneric(
   "convertFromTo",
   signature = c("gng"),
   function(gng, queryId, targetDb = c("GeneSymbolDB"), sourceDb = NA, longFormat = TRUE, chromosomal = TRUE) {
      standardGeneric("convertFromTo")
})

#' @export
setMethod("convertFromTo",
   signature(gng = "GeneNameGenieR"),
    function(gng, queryId, targetDb, sourceDb, longFormat, chromosomal) {

    tDb = targetDb;
    if (!longFormat) {
        tDb = c(tDb, names(gng@toAddDbCols));
    }

    q = paste0("CALL rcsi.convert.table.convertIds(",
               ifelse(length(queryId) == 1, "[{ids}], ", "{ids}, "),
               ifelse(length(tDb) == 1, "[{dbs}], ", "{dbs}, "),
               # "{dbs}, ",
               ifelse(is.na(sourceDb), "null, ", "{sourceDb}, "),
               "{chromosomal}) YIELD value ",
               "RETURN value.InputId AS InputId, ",
               "value.InputSourceDb AS InputSourceDb, ",
               # "value.EnsemblGeneId AS EnsemblGeneId, ",
               "value.TargetDb AS TargetDb, ",
               "value.TargetId AS TargetId;");

    x = .postNeo4jRequest(gng, q,
               ids = queryId,
               dbs = tDb,
               sourceDb = sourceDb,
               chromosomal = chromosomal);

    if (!longFormat) {
        x = .unstackDf(x);
        # remove artifially added columns
        colDiff = setdiff(tDb, targetDb);
        x = x[,!(names(x) %in% gng@toAddDbCols[colDiff])];
    }
    x = dplyr::distinct(x);
    .postCheckInput(x);

    return(x);

})





