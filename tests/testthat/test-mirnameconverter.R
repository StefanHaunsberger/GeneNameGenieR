
require(GeneNameGenieR);

testthat::context("mirnameconversion");

gng = GeneNameGenieR();

testthat::test_that("convertToCurrentMirbaseVersion of 'hsa-miR-29a' and get some metadata", {

    inputIds = 'hsa-miR-29a';
    target = data.frame(
        InputId = inputIds,
        Accession = "MIMAT0000086",
        CurrentMirna = 'hsa-miR-29a-3p',
        CurrentVersion = "22",
        nExperiments = "77",
        EvidenceType = "experimental",
        stringsAsFactors = FALSE);

    testthat::expect_equal(convertToCurrentMirbaseVersion(gng, inputIds, metadata = c('nExperiments', 'evidenceType')), target);
})

testthat::test_that("convertMatureMirnasToVersions(gng, 'hsa-miR-29a', c(17, 21, 22))", {

    inputIds = 'hsa-miR-29a';
    versions = c(17, 21, 22)
    target = data.frame(
        InputId = inputIds,
        MatureAccession = rep("MIMAT0000086", 3),
        miRBaseVersion = c("17", "21", "22"),
        TargetMirna = c("hsa-miR-29a", "hsa-miR-29a-3p", "hsa-miR-29a-3p"),
        stringsAsFactors = FALSE);

    testthat::expect_equal(convertMatureMirnasToVersions(gng, inputIds, c(17, 21, 22)), target);
})
