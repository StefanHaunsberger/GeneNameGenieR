


.geneCols = list(
  aliases = "GeneSymbolAlias",
  symbol = "OfficialGeneSymbol",
  region = "GeneRegion",
  regionStart = "GeneStart",
  regionEnd = "GeneEnd",
  id = "EnsemblGeneId",
  version = "EnsemblGeneVersion",
  biotype = "EnsemblGeneBiotype",
  name = "GeneName"
)

.transcriptCols = list(
  region = "TranscriptRegion",
  regionStart = "TranscriptStart",
  regionEnd = "TranscriptEnd",
  id = "EnsemblTranscriptId",
  version = "EnsemblTranscriptVersion",
  biotype = "EnsemblTranscriptBiotype"
);
.proteinCols = list(
  id = "EnsemblProteinId",
  version = "EnsemblProteinVersion"
);

.elTypeListFun = function(type, n) {call(as.character(type), n)};
.elTypeList = list(
  aliases = "character",
  symbol = "character",
  region = "character",
  regionStart = "integer",
  regionEnd = "integer",
  id = "character",
  version = "integer",
  biotype = "character",
  name = "character"
);

.processParameter = function(parameterName, parameterValue) {
    if (length(parameterValue) == 1 && is.na(parameterValue)) {
        parameterValue = "";
    }
    if (typeof(parameterValue) == "character") {
        return(paste0("\"", parameterName, "\" : ",
                    ifelse(length(parameterValue) == 1,
                            paste0("\"", parameterValue, "\""),
                            paste0("[\"", paste0(parameterValue, collapse = "\", \""), "\"]")
                           )
                    )
               );
    } else {
        if (typeof(parameterValue) == "logical") {
            parameterValue = tolower(parameterValue);
        }
        return(paste0("\"", parameterName, "\" : ",
                      ifelse(length(parameterValue) == 1,
                             parameterValue,
                             paste0("[", paste0(parameterValue, collapse = ", "), "]")
                            )
                      )
               );

    }
           # return(paste0("\"", parameterName, "\" : \"", ifelse(length(parameterValue) == 1, parameterValue, paste0("[\"", paste0(parameterValue, collapse = "\", \""), "\"]"))));
}

.createParameterBody = function(params) {
    paramBody = "";
    if (length(params) > 0) {
        paramNames = names(params);
        paramBody = ", \n\"params\" : {";
        paramBody = paste0(paramBody, .processParameter(paramNames[1], params[[1]]));
        if (length(paramNames) > 1) {
            for (i in 2:(length(paramNames))) {
                paramBody = paste0(paramBody, ", \n", .processParameter(paramNames[i], params[[i]]));
            }
        }
        paramBody = paste0(paramBody, "}");
    }
    return(paramBody);
}

.postNeo4jRequest = function(query, ...) {

    args = as.list(match.call());

    x = data.frame();

    if (!is.null(args$query)) {
        dots = list(...);

        paramBody = .createParameterBody(dots);

        body = paste0('{ \"query\": \"', query, "\"", paramBody, "\n}");

        # req = httr::POST("http://localhost:7474/db/data/transaction/commit",
        req = httr::POST(paste0(pkg.env$url, pkg.env$cypherEndpoint, sep = ""),
                         body = body);

        if (req$status_code != 200) {
            stop("Unexpected status code returned by requst:\n",
                 paste(capture.output(print(req)), collapse="\n"))
        }

        httr::stop_for_status(req);
        json = httr::content(req, "text");

        if (!jsonlite::validate(json)) {
            stop("Malformatted JSON object was returned!");
        }

        # parse JSON object
        obj = jsonlite::fromJSON(json);

        x = as.data.frame(x = obj$data, stringsAsFactors = FALSE);
        if (nrow(x) == 0) {
            return(x)
        } else {
            colnames(x) = obj$columns;
        }
    } else {
        print("no query provided");
    }

    return(x);

}

.unstackDf = function(x, keyCols = c()) {

    if (length(keyCols) == 0) {
        keyColsBool = !(names(x) %in% c("TargetDb", "TargetId"));
        keyCols = names(x)[keyColsBool];
        keyColsElTypes = as.character(sapply(x[keyColsBool], function(i) class(i)));
    }
    colNames = keyCols;

    # Add all DB names as column names
    colNames = c(colNames, sapply(unique(x$TargetDb), "make.names"));
    # keyColsElTypes = c(keyColsElTypes, rep("character", (length(colNames) - length(keyColsElTypes))));
    colElTypes = c(keyColsElTypes, rep("character", (length(colNames) - length(keyColsElTypes))));

    nRows = nrow(unique(x[keyCols]));

    xL = list()
    for (idx in 1:length(colElTypes)) {
        xL[[colNames[idx]]] = eval(.elTypeListFun(colElTypes[idx], nRows))
    }
    xDf = as.data.frame(xL, stringsAsFactors = FALSE);
    xDf[,] = NA;

    x = dplyr::arrange_(x, c(keyCols, "TargetDb"));

    idx = 0;
    key = "";
    keyDb = "";
    for (rowIdx in 1:nrow(x)) {
        row = x[rowIdx,];
        keyCurr = paste0(row[keyCols], collapse = "");
        # new entry
        if (key != keyCurr) {
            idx = idx + 1;
            key = keyCurr;
            # fill values from each key column into new list
            xDf[idx,keyCols] = row[keyCols];
        }
        name = make.names(row[1,"TargetDb"]);
        if (is.na(xDf[idx,name])) {
            xDf[idx,name] = row[1,"TargetId"];
        } else {
            xDf[idx,name] = paste(xDf[idx,name], row[1,"TargetId"], sep = ",");
        }
    }

    return(xDf);

}

.parseGngList = function(gng) {

  #########################
  # Count number of rows and create columns
  #  In the case of the columns only look at the first entry
  #  First four columns are inputDB, inputID, targetID and targetDB
  nCols = 4;
  colNames = c("InputId", "InputSourceDb", "TargetId", "TargetDb");
  colNamesUnstack = c("InputId", "InputSourceDb");
  colElTypes = c(rep("character", 4));
  # Gene column names and types
  #   For sapply <<- would enable to manipulate variable from outer scope
  # sapply(names(gng[[1]]$value[["gene"]]), function(name) {
  #     if (name %in% names(elTypeList)) {
  #         colNames = c(colNames, geneCols[[name]]);
  #         colNamesUnstack = c(colNamesUnstack, geneCols[[name]]);
  #         colElTypes = c(colElTypes, elTypeList[[name]]);
  #     }
  # })
  for (name in names(gng[[1]]$value$gene)) {
    if (name %in% names(.elTypeList)) {
      colNames = c(colNames, .geneCols[[name]]);
      colNamesUnstack = c(colNamesUnstack, .geneCols[[name]]);
      colElTypes = c(colElTypes, .elTypeList[[name]]);
    }
  }
  # Transcript
  for (name in names(gng[[1]]$value$gene$transcript[[1]])) {
    if (name %in% names(.elTypeList)) {
      colNames = c(colNames, .transcriptCols[[name]]);
      colNamesUnstack = c(colNamesUnstack, .transcriptCols[[name]]);
      colElTypes = c(colElTypes, .elTypeList[[name]]);
    }
  }
  # Peptide
  for (name in names(gng[[1]]$value$gene$transcript[[1]]$peptide)) {
    if (name %in% names(.elTypeList)) {
      colNames = c(colNames, .proteinCols[[name]]);
      colNamesUnstack = c(colNamesUnstack, .proteinCols[[name]]);
      colElTypes = c(colElTypes, .elTypeList[[name]]);
    }
  }

  # compute number of rows
  nRows = 0;
  for (entry in gng) {
    # For each gene entry
    # nRows += 1;
    # Number of external gene identifiers
    gene = entry$value$gene;
    nRows = nRows + (length(gene$targetIds) * length(gene$transcript));
    for (transcript in gene$transcript) {
      nRows = nRows + length(transcript$targetIds);
      nRows = nRows + length(transcript$peptide$targetIds);
    }
  }
  nRows

  # Initialise list which later is converted to dataframe
  xL = list()
  for (idx in 1:length(colElTypes)) {
    xL[[colNames[idx]]] = eval(.elTypeListFun(colElTypes[idx], nRows))
  }

  geneColsNames = names(geneCols);
  transcriptColsNames = names(transcriptCols);
  proteinColsNames = names(proteinCols);

  # Fill list
  lowG = 1;
  highG = 1;
  lowT = 1;
  highT = 1;
  for (entry in gng) {
    # For each gene entry
    lowT = lowG;
    highT = lowT;
    highG = lowG;
    # Number of external gene identifiers
    gene = entry$value$gene;
    externalIds = character();
    externalDbs = character();
    if (length(gene$targetIds) > 0) {
      externalDbs = c(externalDbs, sapply(gene$targetIds, '[[', "dbName"));
      externalIds = c(externalIds, sapply(gene$targetIds, '[[', "id"));
    }
    externalGeneIds = externalIds;
    externalGeneDbs = externalDbs;
    for (transcript in gene$transcript) {
      # add external transcript information
      if (length(transcript$targetIds) > 0) {
        externalDbs = c(externalDbs, sapply(transcript$targetIds, '[[', "dbName"));
        externalIds = c(externalIds, sapply(transcript$targetIds, '[[', "id"));
      }
      # add external protein information
      if (length(transcript$peptide$targetIds) > 0) {
        externalDbs = c(externalDbs, sapply(transcript$peptide$targetIds, '[[', "dbName"));
        externalIds = c(externalIds, sapply(transcript$peptide$targetIds, '[[', "id"));
      }

      highT = highT + (length(externalIds) - 1);
      # Add attribute info, such as transcript ID, version,...
      for (transcriptKey in names(transcript)) {
        if (transcriptKey %in% transcriptColsNames) {
          xL[[transcriptCols[[transcriptKey]]]][lowT:highT] = transcript[[transcriptKey]];
        }
      }
      for (peptideKey in names(transcript$peptide)) {
        if (peptideKey %in% proteinColsNames) {
          xL[[proteinCols[[peptideKey]]]][lowT:highT] = transcript$peptide[[peptideKey]];
        }
      }
      # add external IDs and DBs
      xL$TargetId[lowT:highT] = externalIds;
      xL$TargetDb[lowT:highT] = externalDbs;

      # reset external IDs/DBs with gene IDs/DBs (they are needed for next transcript, if there is one)
      externalIds = externalGeneIds;
      externalDbs = externalGeneDbs;

      lowT = highT + 1;
      highT = lowT;
    }
    highG = highT - 1;
    for (name in names(gene)) {
      if ((name %in% geneColsNames) && length(gene[[name]]) > 0) {
        val = gene[[name]];
        if (name == "aliases") {
          val = paste0(val, collapse = ",");
        }
        xL[[geneCols[[name]]]][lowG:highG] = val;
      }
    }
    xL$InputId[lowG:highG] = entry$value$inputId;
    xL$InputSourceDb[lowG:highG] = entry$value$inputDb;
    lowG = highG + 1;
  }

  colElTypesUnstack = colElTypes[sapply(colNames, function(name) name %in% colNamesUnstack)];

  return(as.data.frame(xL, stringsAsFactors = FALSE));

}

.postCheckInput = function(x) {
    # Find all input values that resolve to multiple input sources
    # e.g. '596' is an EntrezGene ID as well as a WikiGene ID
    x2 = x %>%
        dplyr::select(InputId, InputSourceDb) %>%
        unique() %>%
        dplyr::group_by(InputId) %>%
        dplyr::summarise(nEl = n()) %>%
        dplyr::filter(nEl > 1);

    if (nrow(x2) != 0) {
        warning("Some input identifier(s) match to more than one input source database");
        print(x2);
    }

}

.postCheckMirnaTranslation = function(x, input) {
    x2 = x %>%
        dplyr::select(InputId, MatureAccession) %>%
        unique() %>%
        dplyr::group_by(InputId) %>%
        dplyr::summarise(nEl = n()) %>%
        dplyr::filter(nEl > 1);

    if (nrow(x2) != 0) {
        warning("Input identifier(s) match to more than one MIMAT accession!");
        print(x2);
    }

    dif = setdiff(input, x$InputId);
    if (length(dif) != 0) {
        warning("Input identifier(s) not found in DB.");
        print(dif);
    }
}

.postCheckMirnaTranslation2 = function(x, input) {
    x2 = x %>%
        dplyr::select(InputId, Accession) %>%
        unique() %>%
        dplyr::group_by(InputId) %>%
        dplyr::summarise(nEl = n()) %>%
        dplyr::filter(nEl > 1);

    if (nrow(x2) != 0) {
        warning("Some input identifiers match to more than one MIMAT accession!");
        print(x2);
    }

    dif = setdiff(input, x$InputId);
    if (length(dif) != 0) {
        warning("Input identifier(s) not found in DB.");
        print(dif);
    }

}
