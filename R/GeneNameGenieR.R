#' @title GeneNameGenieR: Lightning fast graph based gene name converter
#' @author Stefan J. Haunsberger
#' @name GeneNameGenieR
#' @description
#'
#' @importFrom magrittr %>%
#' @import httr
#' @importFrom jsonlite fromJSON validate
#'
#' @docType package

library(dplyr)
library(magrittr)
library(jsonlite)
library(httr)

DEFAULT_PORT = 7474;
DEFAULT_PATH = "db/data/";
DEFAULT_HOST = "http://localhost";
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
        url = "character",
        host = "character",
        path = "character",
        port = "numeric",
        cypherEndpoint = "character",
        verbose = "logical",
        toAddDbCols = "character"
    ),
    prototype = list(
        url = paste(paste(DEFAULT_HOST, DEFAULT_PORT, sep = ":"), DEFAULT_PATH, sep = "/"),
        host = "localhost",
        path = "db/data/",
        port = 7474,
        cypherEndpoint = "cypher",
        verbose = FALSE,
        toAddDbCols = c(ArrayExpress = "Ensembl.Human.Gene",
                        Ens_Hs_transcript = "Ensembl.Human.Transcript",
                        Ens_Hs_translation = "Ensembl.Human.Translation")
    )
);

validityGeneNameGenieObj = function(object) {
    errors = character(0);
    for (slotName in c("host", "path", "port")) {
        if (length(slot(object, slotName)) != 1) {
            errors = c(errors, sprintf("Length of slot %s must be 1", slotName));
        }
    }

    if (!endsWith(object@path, "/")) {
        errors = c(errors, sprintf("Slot 'path' must end on '/' (slash)"));
    }

    if (length(errors) > 0) {
        errors;
    } else {
        TRUE
    }
}

setValidity("GeneNameGenieR", validityGeneNameGenieObj)

# #' @exportClass GeneNameGenieR

#' @param host : Graph database URL, such as "localhost"
#' @param port : Port where the DB is exposed, such as 7474
#' @param path : Location of the data (defaults to `/db/data/`)
#' @export
GeneNameGenieR = function(host = DEFAULT_HOST, port = DEFAULT_PORT, path = DEFAULT_PATH) {
    # if (is.na(url)) {
    #     url = DEFAULT_HOST;
    # }

    g = new("GeneNameGenieR", host = host, port = port, path = path);
    # g@baseUrl = url;
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
#' @examples
#'
#' \dontrun{
#'   gng = GeneNameGenieR();
#'   getOfficialGeneSymbol(gng, "Bcl-2")
#'     InputId     InputSourceDb OfficialGeneSymbol
#'   1   Bcl-2 Gene Symbol Alias               BCL2
#' }
#'
#' @export
setGeneric(
    "getOfficialGeneSymbol",
    signature = c("gng"),
    function(gng, queryId, sourceDb = NA_character_, chromosomal = TRUE) {
        standardGeneric("getOfficialGeneSymbol")
    })

# @export
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
#' @examples
#' \dontrun{
#'   gng = GeneNameGenieR();
#'   convertFromTo(gng, c("BCL2", "AMPK"), "EntrezGene");
#'   InputId        InputSourceDb  TargetDb TargetId
#' 1    BCL2 Official Gene Symbol NCBI gene      596
#' 2    AMPK    Gene Symbol Alias NCBI gene     5563
#' }
#'
#' @export
#' @seealso \code{\link{query}}
setGeneric(
   "convertFromTo",
   signature = c("gng"),
   function(gng, queryId, targetDb = c("GeneSymbolDB"), sourceDb = NA, longFormat = TRUE, chromosomal = TRUE) {
      standardGeneric("convertFromTo")
})

# @export
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
    if (nrow(x) > 0) {
        x = dplyr::distinct(x);
        .postCheckInput(x);
    }

    return(x);

})





