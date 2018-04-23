
require(GeneNameGenieR);

testthat::context("mirnameconversion");

gng = GeneNameGenieR();

testthat::test_that("", {

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
