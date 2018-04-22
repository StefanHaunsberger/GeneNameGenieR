#' @title GeneNameGenieR: Lightning fast graph based gene name converter
#' @author Stefan J. Haunsberger
#' @name GeneNameGenieR
#' @description
#'
#' @docType package

library(RNeo4j)
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
        url = "character",
        verbose = "logical",
        graph = "graph",
        toAddDbCols = "character"
    ),
    prototype = list(
        url = DEFAULT_URL,
        verbose = FALSE,
        graph = RNeo4j::startGraph(URL),
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
    g@graph = RNeo4j::startGraph(url);
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

  # x = RNeo4j::cypher(gng@graph, q,
  #            ids = queryId,
  #            sourceDb = sourceDb,
  #            chromosomal = chromosomal);

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

    # x = RNeo4j::cypher(gng@graph, q,
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


#' @title Translate identifiers and receive attributes for Ensembl genes, transcripts and proteins
#'
#' @description Similar to the \code{\link{convertFromTo}} function, the \code{\link{query}} function converts any matchable
#' input identifier to the desired target database identifier. Moreover, with the \preformatted{params} parameter
#' \code{\link{query}} returns attributes for Ensembl genes, transcripts and proteins, such as location and biotype.
#'
#' @param gng An instance of the GeneNameGenieR object (e.g. \code{gng = GeneNameGenieR()})
#' @param queryId A \preformatted{character} vector of molecular identifiers that needs to be converted
#' @param targetDbs A \preformatted{character} vector of target databases (default equals to 'EntrezGene' and 'GeneSymbolDb' for
#' the Entrez gene ID and official gene symbol respectively)
#' @param params String or vector containing retrievable attributes, such as 'engs:region' for the Ensembl gene region (\code{\link{getValidAttributes}}).
#' @param sourceDb String containing the origin database of the input ID(s) (optional).
#'
#' @return A \preformatted{data.frame} object with the columns InputId, InputSourceDb, TargetDb and TargetId.
#'
#' @examples
#'
#' \dontrun{
#' query(
#'    GeneNameGenieR(),
#'    c('AMPK'))
#' }
#' \dontrun{
#' query(
#'    GeneNameGenieR(),
#'    c('AMPK', 'Bcl-2', '596', 'NM_000657'),
#'    targetDb = c('EntrezGene', 'ArrayExpress', 'RefSeq_mRNA'),
#'    params = c('ensg:region', 'ensg:regionStart', 'ensg:regionEnd'),
#'    longFormat = FALSE)
#'
# A tibble: 1 x 2
#'      InputId     n
#'        <chr> <int>
#'    1     596     2
#'      InputId     InputSourceDb  GeneEnd GeneStart   EnsemblGeneId GeneRegion NCBI.gene RefSeq.mRNA
#' 1        596          WikiGene 63320128  63123346 ENSG00000171791         18       596   NM_000657
#' 2        596          WikiGene 63320128  63123346 ENSG00000171791         18       596   NM_000633
#' 3        596          WikiGene 63320128  63123346 ENSG00000171791         18       596        <NA>
#' 4        596         NCBI gene 63320128  63123346 ENSG00000171791         18       596   NM_000657
#' 5        596         NCBI gene 63320128  63123346 ENSG00000171791         18       596   NM_000633
#' 6        596         NCBI gene 63320128  63123346 ENSG00000171791         18       596        <NA>
#' 7       AMPK Gene Symbol Alias 56715335  56645322 ENSG00000162409          1      5563        <NA>
#' 8       AMPK Gene Symbol Alias 56715335  56645322 ENSG00000162409          1      5563   NM_006252
#' 9      Bcl-2 Gene Symbol Alias 63320128  63123346 ENSG00000171791         18       596   NM_000657
#' 10     Bcl-2 Gene Symbol Alias 63320128  63123346 ENSG00000171791         18       596   NM_000633
#' 11     Bcl-2 Gene Symbol Alias 63320128  63123346 ENSG00000171791         18       596        <NA>
#' 12 NM_000657       RefSeq mRNA 63320128  63123346 ENSG00000171791         18       596   NM_000657
#' Warning message:
#' In postCheckInput(x) :
#' Some input identifiers match to more than one input source database
#' }
#' @export
#' @seealso \code{\link{convertFromTo}}
setGeneric(
   "query",
   signature = c("gng"),
   function(gng,
            queryId,
            targetDb = c("GeneSymbolDB"),
            params = NA,
            sourceDb = NA,
            longFormat = TRUE,
            chromosomal = TRUE) {
      standardGeneric("query")
})

#' @export
setMethod("query",
          signature(gng = "GeneNameGenieR"),
          function(gng, queryId, targetDb, params, sourceDb, longFormat, chromosomal) {

    tDb = targetDb;

    if (!longFormat) {
        tDb = c(tDb, names(gng@toAddDbCols));
    }

    q = paste0("CALL rcsi.convert.json.convertIdsExtended(",
             ifelse(length(queryId) == 1, "[{ids}], ", "{ids}, "),
             ifelse(length(tDb) == 1, "[{dbs}], ", "{dbs}, "),
             # "{dbs}, ",
             ifelse(length(params) == 1, ifelse(is.na(params), "null, ", "[{params}], "), "{params}, "),
             ifelse(is.na(sourceDb), "null, ", "{sourceDb}, "), "{chromosomal});");

    res = RNeo4j::cypherToList(gng@graph, q,
                       ids = queryId,
                       dbs = tDb,
                       params = params,
                       sourceDb = sourceDb,
                       chromosomal = chromosomal);

    #
    x = .parseGngList(res);

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






