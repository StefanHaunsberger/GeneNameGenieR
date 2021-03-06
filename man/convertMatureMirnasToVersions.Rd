% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mirnameconversion.R
\name{convertMatureMirnasToVersions}
\alias{convertMatureMirnasToVersions}
\title{Translate mature miRNA names to different versions}
\usage{
convertMatureMirnasToVersions(queryId, targetVersion = NA_integer_,
  species = "", sequence = FALSE)
}
\arguments{
\item{queryId}{A character vector of miRNA identifiers which needs to be converted}

\item{targetVersion}{A numeric vector with the target miRBase version(s), such as `c(20, 21)` (default: 22).}

\item{species}{A string containing the three letter species, such as `'hsa'` or `'mmu'`.}

\item{sequence}{A boolean depicting if sequence for respective miRBase versions is returned or not (default: FALSE).}
}
\description{
Taking single or multiple `queryId`s return the name(s) for a single or
multiple `targetVersion`(s). Optionally the `species` can be provided and if
`sequence` == TRUE the sequence for each version respecively returned.
}
\examples{

\dontrun{
  convertMatureMirnasToVersions('hsa-miR-29a');
        InputId MatureAccession miRBaseVersion    TargetMirna
  1 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p
}
\dontrun{
  convertMatureMirnasToVersions('hsa-miR-29a', c(17, 21, 22))
       InputId MatureAccession miRBaseVersion    TargetMirna
 1 hsa-miR-29a    MIMAT0000086             17    hsa-miR-29a
 2 hsa-miR-29a    MIMAT0000086             21 hsa-miR-29a-3p
 3 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p
}
\dontrun{
  convertMatureMirnasToVersions('hsa-miR-29a', c(17, 21, 22), sequence = TRUE)
      InputId MatureAccession miRBaseVersion    TargetMirna         TargetSequence
1 hsa-miR-29a    MIMAT0000086             17    hsa-miR-29a UAGCACCAUCUGAAAUCGGUUA
2 hsa-miR-29a    MIMAT0000086             21 hsa-miR-29a-3p UAGCACCAUCUGAAAUCGGUUA
3 hsa-miR-29a    MIMAT0000086             22 hsa-miR-29a-3p UAGCACCAUCUGAAAUCGGUUA
}

\dontrun{
  convertMatureMirnasToVersions(c('hsa-miR-29a', 'hsa-miR-378'), c(17, 21, 22))
          InputId     n
            <chr> <int>
    1 hsa-miR-378     2
# --------
      InputId MatureAccession miRBaseVersion     TargetMirna
1 hsa-miR-29a    MIMAT0000086             17     hsa-miR-29a
2 hsa-miR-29a    MIMAT0000086             21  hsa-miR-29a-3p
3 hsa-miR-29a    MIMAT0000086             22  hsa-miR-29a-3p
4 hsa-miR-378    MIMAT0000732             17     hsa-miR-378
5 hsa-miR-378    MIMAT0000731             17    hsa-miR-378*
6 hsa-miR-378    MIMAT0000732             21 hsa-miR-378a-3p
7 hsa-miR-378    MIMAT0000731             21 hsa-miR-378a-5p
8 hsa-miR-378    MIMAT0000732             22 hsa-miR-378a-3p
9 hsa-miR-378    MIMAT0000731             22 hsa-miR-378a-5p
Warning message:
    In postCheckMirnaTranslation(x) :
    Some input identifiers match to more than one MIMAT accession!
}

}
