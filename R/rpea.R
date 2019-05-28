

#' @title Perform Reactome Pathway Enrichment Analysis (RPEA)
#'
#' @description Run the Cypher stored procedure `rcsi.enrichment.reactomePathwayEnrichment`.
#'
#' @param queryIds A character vector of molecular identifiers.
#' @param sourceDb String containing the origin database of the input ID(s) (optional).
#' @param pValueCutoff Floating point value (default 0.05) for p-value threshold.
#' @param minNAnnotatedGenes Integer value (default 10) denoting the minimum number of annotated genes a Reactome
#' pathway needs to be included in the analysis.
#' @param maxNAnnotatedGenes Integer value (default 500) denoting the maximum number of annotated genes a Reactome
#' pathway needs to be included in the analysis.
#' @param pAdjustmentMethod String (default 'BF' - Bonferroni) containing the p-value correction for multiple testing
#' (other options: 'BH' - Benjamini and Hochberg)
#'
#'
#' @return A data.frame object with the columns InputId, InputSourceDb and columns containing TargetDb and TargetId values.
#'
#' @examples
#' \dontrun{
#'   GeneNameGenieR::rpea(c('PRKAA1', 'PRKAB2', 'PRKAG2', 'STK11'))
#' }
#'
#'
#' @export
rpea = function(queryIds, sourceDb = 'GeneSymbolDB', pValueCutoff = 0.05, minNAnnotatedGenes = 10, maxNAnnotatedGenes = 500,
                pAdjustmentMethod = 'BF', qValueCutoff = 0.2) {

  q = paste0("CALL rcsi.enrichment.reactomePathwayEnrichment(",
             ifelse(length(queryIds) == 1, "[{queryIds}], ", "{queryIds}, "),
             "{sourceDb}, {pValueCutoff}, {minNAnnotatedGenes}, {maxNAnnotatedGenes}, {pAdjustmentMethod}, {qValueCutoff}",
             ") YIELD reactomeID, description, bgRatio, geneRatio, geneIds, pValue, pValueAdjusted, numberOfPopulationSuccesses, populationSize, sampleSize, numberOfSampleSuccesses ",
             "RETURN reactomeID, description, bgRatio, geneRatio, geneIds, pValue, pValueAdjusted, numberOfPopulationSuccesses, populationSize, sampleSize, numberOfSampleSuccesses")

  x = .postNeo4jRequest(q,
                        src = "rpea",
                        queryIds = queryIds,
                        sourceDb = sourceDb,
                        pValueCutoff = pValueCutoff,
                        minNAnnotatedGenes = minNAnnotatedGenes,
                        maxNAnnotatedGenes = maxNAnnotatedGenes,
                        pAdjustmentMethod = pAdjustmentMethod,
                        qValueCutoff = qValueCutoff);

  df = data.frame(
      reactomeID = character(),
      description = character(),
      bgRatio = character(),
      geneRatio = character(),
      geneIds = vector(mode="character"),
      pValue = double(),
      pValueAdjusted = double(),
      numberOfPopulationSuccesses = integer(),
      populationSize = integer(),
      sampleSize = integer(),
      numberOfSampleSuccesses = integer(),
      stringsAsFactors = FALSE);

  if (length(x$data) > 0) {
      for (i in 1:length(x$data)) {
          entry = x$data[[i]];
          for (j in 1:length(x$columns)) {
              if (x$columns[j] == "geneIds") {
                  df[i, x$columns[j]] = paste0(entry[[j]], collapse = ",");
              } else {
                  df[i, x$columns[j]] = entry[[j]];
              }
          }
      }
  }

  return(df);

}
