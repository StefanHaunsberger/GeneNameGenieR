
Metadata_Map = c(confidence = "value.Conficence AS Confidence",
                 type = "value.Type AS Type",
                 sequence = "value.Sequence AS Sequence",
                 comments = "value.Comments AS Comments",
                 previousIds = "value.PreviousIds AS PreviousIds",
                 url= "value.URL AS Url",
                 chromosome = "value.Chromosome AS Chromosome",
                 regionStart = "value.RegionStart AS RegionStart",
                 regionEnd = "value.RegionEnd AS RegionEnd",
                 strand = "value.Strand AS Strand",
                 communityAnnotation = "value.CommunityAnnotation AS CommunityAnnotation",
                 nExperiments = "value.nExperiments AS nExperiments",
                 reads = "value.Reads AS Reads",
                 evidenceType = "value.EvidenceType AS EvidenceType");

#' @title Get current miRNA name for mature input miRNAs
#'
#' @description Return the current miRNA name for any mature input identifier `queryId`(s).
#' Thereby, one can pass in either MIMAT IDs or miRNAs or a mix of both. If desired,
#' the `queryType` can be specified passing through either `"mature"` or `"precursor"`.
#'
#' @param queryId A string vector of miRNA identifiers which needs to be converted
#' @param species Any of the 271 supported, three letter KEGG species, such as 'hsa' for _homo sapiens_ (species codes: https://www.genome.jp/kegg/catalog/org_list.html)
#' @param metadata A string vector with parameter values, such as `'sequence'` or `'type'`.
#'
#' @seealso \code{\link{getValidMirnaMetadataValues}}
#'
#' @examples
#'
#' \dontrun{
#'   convertToCurrentMirbaseVersion('hsa-miR-29a');
#'         InputId    Accession   CurrentMirna CurrentVersion
#'   1 hsa-miR-29a MIMAT0000086 hsa-miR-29a-3p             22
#' }
#'\dontrun{
#'   queryIds = c('MIMAT0000447', 'hsa-mir-134', 'hsa-miR-134', 'MI0000146', 'hsa-miR-29a')
#'   convertToCurrentMirbaseVersion(queryIds, 'hsa', 'sequence');
#'        InputId    Accession   CurrentMirna CurrentVersion                                                                  Sequence
#' 1 MIMAT0000447 MIMAT0000447 hsa-miR-134-5p             22                                                    UGUGACUGGUUGACCAGAGGGG
#' 2  hsa-miR-134 MIMAT0000447 hsa-miR-134-5p             22                                                    UGUGACUGGUUGACCAGAGGGG
#' 3  hsa-miR-29a MIMAT0000086 hsa-miR-29a-3p             22                                                    UAGCACCAUCUGAAAUCGGUUA
#' 4  hsa-mir-134    MI0000474    hsa-mir-134             22 CAGGGUGUGUGACUGGUUGACCAGAGGGGCAUGCACUGUGUUCACCCUGUGGGCCACCUAGUCACCAACCCUC
#' }
#' \dontrun{
#'   convertToCurrentMirbaseVersion(queryIds, metadata = 'sequence');
#'        InputId    Accession   CurrentMirna CurrentVersion                                                                  Sequence
#' 1    MI0000146    MI0000146    mmu-mir-99a             22         CAUAAACCCGUAGAUCCGAUCUUGUGGUGAAGUGGACCGCGCAAGCUCGUUUCUAUGGGUCUGUG
#' 2 MIMAT0000447 MIMAT0000447 hsa-miR-134-5p             22                                                    UGUGACUGGUUGACCAGAGGGG
#' 3  hsa-miR-134 MIMAT0000447 hsa-miR-134-5p             22                                                    UGUGACUGGUUGACCAGAGGGG
#' 4  hsa-miR-29a MIMAT0000086 hsa-miR-29a-3p             22                                                    UAGCACCAUCUGAAAUCGGUUA
#' 5  hsa-mir-134    MI0000474    hsa-mir-134             22 CAGGGUGUGUGACUGGUUGACCAGAGGGGCAUGCACUGUGUUCACCCUGUGGGCCACCUAGUCACCAACCCUC
#' }
#' \dontrun{
#'   convertToCurrentMirbaseVersion('hsa-mir-29a', metadata = c('sequence', 'type', 'previousIds'));
#'       InputId Accession CurrentMirna CurrentVersion                                                         Sequence      Type PreviousIds
#' 1 hsa-mir-29a MI0000087  hsa-mir-29a             22 AUGACUGAUUUCUUUUGGUGUUCAGAGUCAAUAUAAUUUUCUAGCACCAUCUGAAAUCGGUUAU antisense  hsa-mir-29
#' }
#'\dontrun{
#'   queryIds = c('hsa-mir-134', 'hsa-mir-29a')
#'   convertToCurrentMirbaseVersion(queryIds, metadata = c('type', 'previousIds'));
#'       InputId Accession CurrentMirna CurrentVersion      Type PreviousIds
#' 1 hsa-mir-134 MI0000474  hsa-mir-134             22 antisense
#' 2 hsa-mir-29a MI0000087  hsa-mir-29a             22 antisense  hsa-mir-29
#' }
#'
#' @export
convertToCurrentMirbaseVersion = function(queryId, species = "", metadata = NA_character_) {

      if (anyNA(metadata)) {
          metadata = "";
      }

      q = paste0("CALL rcsi.mirna.convertToCurrentMirbaseVersion(",
                 ifelse(length(queryId) == 1, "[{queryIds}], ", "{queryIds}, "),
                 "{species}",
                 ifelse(any(metadata == ""), "", ifelse(length(metadata) == 1, ", [{metadata}]", ", {metadata}")),
                 ") YIELD value ",
                 "RETURN DISTINCT value.InputId AS InputId, ",
                 "value.Accession AS Accession, ",
                 "value.CurrentMirna AS CurrentMirna, ",
                 "value.CurrentVersion AS CurrentVersion"
      );

      if (any(metadata != "")) {
          q = paste0(q, ", ", paste0(sapply(metadata, function(m) {
              if (is.na(Metadata_Map[m])) {
                  stop(paste0("requested metadata value '", m, "' not recognised"))
              } else {
                  Metadata_Map[m]
              }}), collapse = ", "));
      }

      x = .postNeo4jRequest(q,
                 queryIds = queryId,
                 species = species,
                 metadata = metadata);

      x = dplyr::distinct(x);
      if (nrow(x) > 0) {
        .postCheckMirnaTranslation2(x, queryId);
      }

      return(x);

}

#' @title Translate mature miRNA names to different versions
#'
#' @description Taking single or multiple `queryId`s return the name(s) for a single or
#' multiple `targetVersion`(s). Optionally the `species` can be provided and if
#' `sequence` == TRUE the sequence for each version respecively returned.
#'
#' @param queryId A character vector of miRNA identifiers which needs to be converted
#' @param targetVersion A numeric vector with the target miRBase version(s), such as `c(20, 21)` (default: 22).
#' @param species A string containing the three letter species, such as `'hsa'` or `'mmu'`.
#' @param sequence A boolean depicting if sequence for respective miRBase versions is returned or not (default: FALSE).
#'
#' @examples
#'
#' \dontrun{
#'   convertMatureMirnasToVersions('hsa-miR-29a');
#'         InputId MatureAccession miRBaseVersion    TargetMirna
#'   1 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p
#' }
#' \dontrun{
#'   convertMatureMirnasToVersions('hsa-miR-29a', c(17, 21, 22))
#'        InputId MatureAccession miRBaseVersion    TargetMirna
#'  1 hsa-miR-29a    MIMAT0000086             17    hsa-miR-29a
#'  2 hsa-miR-29a    MIMAT0000086             21 hsa-miR-29a-3p
#'  3 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p
#' }
#' \dontrun{
#'   convertMatureMirnasToVersions('hsa-miR-29a', c(17, 21, 22), sequence = TRUE)
#'       InputId MatureAccession miRBaseVersion    TargetMirna         TargetSequence
#' 1 hsa-miR-29a    MIMAT0000086             17    hsa-miR-29a UAGCACCAUCUGAAAUCGGUUA
#' 2 hsa-miR-29a    MIMAT0000086             21 hsa-miR-29a-3p UAGCACCAUCUGAAAUCGGUUA
#' 3 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p UAGCACCAUCUGAAAUCGGUUA
#' }
#'
#' \dontrun{
#'   convertMatureMirnasToVersions(c('hsa-miR-29a', 'hsa-miR-378'), c(17, 21, 22))
#'           InputId     n
#'             <chr> <int>
#'     1 hsa-miR-378     2
#' # --------
#'       InputId MatureAccession miRBaseVersion     TargetMirna
#' 1 hsa-miR-29a    MIMAT0000086             17     hsa-miR-29a
#' 2 hsa-miR-29a    MIMAT0000086             21  hsa-miR-29a-3p
#' 3 hsa-miR-29a    MIMAT0000086             22  hsa-miR-29a-3p
#' 4 hsa-miR-378    MIMAT0000732             17     hsa-miR-378
#' 5 hsa-miR-378    MIMAT0000731             17    hsa-miR-378*
#' 6 hsa-miR-378    MIMAT0000732             21 hsa-miR-378a-3p
#' 7 hsa-miR-378    MIMAT0000731             21 hsa-miR-378a-5p
#' 8 hsa-miR-378    MIMAT0000732             22 hsa-miR-378a-3p
#' 9 hsa-miR-378    MIMAT0000731             22 hsa-miR-378a-5p
#' Warning message:
#'     In postCheckMirnaTranslation(x) :
#'     Some input identifiers match to more than one MIMAT accession!
#' }
#'
#' @export
convertMatureMirnasToVersions = function(queryId, targetVersion = NA_integer_, species = "", sequence = FALSE) {

              if (anyNA(targetVersion)) {
                  targetVersion = 0.0;
              }

              q = paste0("CALL rcsi.mirna.convertMatureMirnas(",
                         ifelse(length(queryId) == 1, "[{queryIds}], ", "{queryIds}, "),
                         ifelse(length(targetVersion) == 1, "[{versions}], ", "{versions}, "),
                         "{species}, ",
                         "{sequence}",
                         ") YIELD value ",
                         "RETURN ",
                         "value.inputId AS InputId, ",
                         "value.matureAccession AS MatureAccession, ",
                         "value.targetMirna AS TargetMirna, ",
                         "value.miRBaseVersion AS miRBaseVersion"
              );

              if (sequence) {
                  q = paste0(q, ", value.sequence AS Sequence");
              }

              x = .postNeo4jRequest(q,
                         queryIds = queryId,
                         versions = targetVersion,
                         species = species,
                         sequence = sequence);

              x = dplyr::distinct(x);

              if (nrow(x) > 0) {
                .postCheckMirnaTranslation(x, queryId);
              }

              return(x);

          }

