
library(GeneNameGenieR);

testthat::context("getOfficialGeneSymbol");

testthat::test_that("get EntrezGene ID for 'BCL2' and 'AMPK'", {

    inputIds = c("BCL2", "AMPK");
    inputSourceDbs = c("Official Gene Symbol", "Gene Symbol Alias");
    targetDbs = c("NCBI gene", "NCBI gene");
    targetIds = c("596", "5563");
    target = data.frame(InputId = inputIds, InputSourceDb = inputSourceDbs, TargetDb = targetDbs, TargetId = targetIds, stringsAsFactors = FALSE);

    testthat::expect_equal(convertFromTo(inputIds, "EntrezGene"), target);
})
