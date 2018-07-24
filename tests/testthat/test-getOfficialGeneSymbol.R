
require(GeneNameGenieR);

testthat::context("getOfficialGeneSymbol");

testthat::test_that("With 'Bcl-2'", {
    alias = "Bcl-2";

    target = data.frame(InputId = alias, InputSourceDb = "Gene Symbol Alias", OfficialGeneSymbol = "BCL2", stringsAsFactors = FALSE);

    testthat::expect_equal(getOfficialGeneSymbol(alias), target);
})
