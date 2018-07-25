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

pkg.env = new.env();

DEFAULT_PORT = 7474;
DEFAULT_PATH = "db/data/";
DEFAULT_HOST = "http://localhost";
pkg.env$port = DEFAULT_PORT;
pkg.env$path = DEFAULT_PATH;
pkg.env$host = DEFAULT_HOST;
pkg.env$url = "http://localhost:7474/db/data/";
pkg.env$cypherEndpoint = "cypher";
pkg.env$toAddDbCols = c(ArrayExpress = "Ensembl.Human.Gene",
                Ens_Hs_transcript = "Ensembl.Human.Transcript",
                Ens_Hs_translation = "Ensembl.Human.Translation");

setOldClass("graph")

#' @title Set GeneNameGenie-Neo4j connection
#'
#' @description The connection defaults to 'localhost' if not specified differently.
#'
#' @param host A string value of the host, such as 'localhost' or an IP address
#' @param port An integer value for the port, such as 7474 (default)
#' @param path String value containing the path to the database files (defaults to `db/data/`). The value must end on a `/` (slash).
#'
#' @export setNeo4jConnection
setNeo4jConnection = function(host = DEFAULT_HOST, port = DEFAULT_PORT, path = DEFAULT_PATH) {
    errors = character(0);
    for (varName in c("host", "path", "port")) {
        if (length(varName) != 1) {
            errors = c(errors, sprintf("Length of variable %s must be 1", varName));
        }
    }

    # Attach trailing '/' in case it is missing
    if (!endsWith(path, "/")) {
        path = paste0(path, "/");
    }

    if (length(errors) > 0) {
        errors;
    } else {
    pkg.env$host = host;
    pkg.env$port = port;
    pkg.env$path = path;
    pkg.env$url = paste(paste(host, port, sep = ":"), path, sep = "/");

    }
}

#' @title Show GeneNameGenie-Neo4j connection details
#'
#' @export showGeneNameGenieConfig
showGeneNameGenieConfig = function() {
    cat(paste0(
        "GeneNameGenie-Neo4j configuration:\n",
        "----------------------------------\n",
        "host: ", pkg.env$host, "\n",
        "port: ", pkg.env$port, "\n",
        "path: ", pkg.env$path, "\n",
        "url: ", pkg.env$url, "\n",
        "cypher-endpoint: ", pkg.env$cypherEndpoint, "\n"
    ));
}

source("R/parameters.R")
source("R/utils.R")
# source("R/mirnameconversion.R")

#' @title Convert input id to its corresponding official gene symbol
#'
#' @description Compute the official gene symbol for a single or multiple `queryId`(s).
#' If `sourceDb` is missing, the database will automatically try to detect the queryId's
#'   source database. Per default `chromosomal` is set to `true`.
#'
#' @param queryId A character vector of molecular identifiers that needs to be converted
#' @param targetDb A character vector of target databases (default equals to 'EntrezGene' and 'GeneSymbolDb' for
#' the Entrez gene ID and official gene symbol respectively)
#' @param sourceDb String containing the origin database of the input ID(s) (optional).
#' @param chromosomal Boolean parameter (default TRUE). Whether the search and output shall include features located on
#' chromosomes other than the standard 23 {1,...21, X, Y}. Ensembl also includes features that are located on so called
#' alternative human chromosomes, such as CHR_HSCHR5_6_CTG1. In the case where `chromosomal = FALSE`, not even input
#' identifiers that are associated with this location will be found.
#'
#' @examples
#'
#' \dontrun{
#'   getOfficialGeneSymbol("Bcl-2")
#'     InputId     InputSourceDb OfficialGeneSymbol
#'   1   Bcl-2 Gene Symbol Alias               BCL2
#' }
#'
#' @export
getOfficialGeneSymbol = function(queryId, sourceDb = NA_character_, chromosomal = TRUE) {

  q = paste0("CALL rcsi.convert.table.getOfficialGeneSymbol(",
      ifelse(length(queryId) == 1, "[{ids}], ", "{ids}, "),
      ifelse(is.na(sourceDb), "null, ", "{sourceDb}, "),
      "{chromosomal}) YIELD value ",
      "RETURN ",
      "   value.InputId AS InputId, ",
      "   value.InputSourceDb AS InputSourceDb, ",
      "   value.OfficialGeneSymbol AS OfficialGeneSymbol");

  x = .postNeo4jRequest(q,
                     ids = queryId,
                     sourceDb = sourceDb,
                     chromosomal = chromosomal);

  if (nrow(x) > 0) {
      .postCheckInput(x, queryId);
  }

  return(x);

}

#' @title Convert molecular identifier to target identifiers
#'
#' @description Converts molecular input identifier to identifier from target database. Thereby, the source database of the
#' input identifier does not have to be specified as the database detects the corresponding database automatically.
#'
#' @param queryId A character vector of molecular identifiers that needs to be converted
#' @param targetDb A character vecotr of target databases (default equals to 'EntrezGene' and 'GeneSymbolDb' for
#' the Entrez gene ID and official gene symbol respectively)
#' @param sourceDb String containing the origin database of the input ID(s) (optional).
#' @param longFormat Boolean parameter (default TRUE). Indicates if returned dataframe is in the form of a longFormat
#' @param chromosomal Boolean parameter (default TRUE). Whether the search and output shall include features located on
#' chromosomes other than the standard 23 {1,...21, X, Y}. Ensembl also includes features that are located on so called
#' alternative human chromosomes, such as CHR_HSCHR5_6_CTG1. In the case where `chromosomal = FALSE`, not even input
#' identifiers that are associated with this location will be found.
#'
#' @return A data.frame object with the columns InputId, InputSourceDb and columns containing TargetDb and TargetId values.
#'
#' @examples
#' \dontrun{
#'   convertFromTo(c("BCL2", "AMPK"), "EntrezGene");
#'   InputId        InputSourceDb  TargetDb TargetId
#' 1    BCL2 Official Gene Symbol NCBI gene      596
#' 2    AMPK    Gene Symbol Alias NCBI gene     5563
#' }
#'
#' @export
#' @seealso \code{\link{query}}
convertFromTo = function(queryId, targetDb = c("GeneSymbolDB"), sourceDb = NA, longFormat = TRUE, chromosomal = TRUE) {

    tDb = targetDb;
    if (!longFormat) {
        tDb = c(tDb, names(pkg.env$toAddDbCols));
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

    x = .postNeo4jRequest(q,
               ids = queryId,
               dbs = tDb,
               sourceDb = sourceDb,
               chromosomal = chromosomal);

    if (!longFormat) {
        x = .unstackDf(x);
        # remove artifially added columns
        colDiff = setdiff(tDb, targetDb);
        x = x[,!(names(x) %in% pkg.env$toAddDbCols[colDiff])];
    }
    if (nrow(x) > 0) {
        x = dplyr::distinct(x);
        .postCheckInput(x, queryId);
    }

    return(x);

}
