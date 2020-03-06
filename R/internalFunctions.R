# Specifies the internally used functions and data

.validSpecies <- function() {
    c("mouse", "rat", "human")
}

.validRegions <- function() {
    c("HPC", "vHPC", "dHPC", "DG", "vDG", "dDG", "PFC", "NACC", "STR", "CA", "AMY", "VTA")
}

.validSex <- function() {
    c("female", "male")
}

.validPlatforms <- function() {
    c("MicroArray", "RNASeq")
}

.revDT <- function(theData, comparisons) {
    colsToRev <- c("logFC", "t", "cohensD", "hedgesG")

    # first split by comparisonsID
    theData <- split(theData, by="comparisonID")
    theData <- lapply(theData, function(x) {
        if (!(x[["comparisonID"]][1] %in% comparisons)) {
            return(x)
        } else {
            for (tr in colsToRev) {
                x[[tr]] <- -x[[tr]]
            }
            groups <- unlist(strsplit(x[["comparison"]][1], " - "))
            x[["comparison"]] <- paste(groups[2], "-", groups[1])
            temp <- x[["n.group1"]]
            x[["n.group1"]] <- x[["n.group2"]]
            x[["n.group2"]] <- temp
            x[["comparisonID"]] <- paste0("-", x[["comparisonID"]])

            return(x)
        }
    })

    return(rbindlist(theData))
}

.createIds <- function(o) {

    # id's are set on the number of species present. If there is only one species, set it to the ensembl id
    species <- getSpecies(o)
    names(species) <- species

    ens <- getGeneEns(o, species)

    # If we have only one species, use the ensembl ids as id
    if (length(ens) == 1) {
        ids <- ens[[1]]
        names(ids) <- ids

        # now set this id in all the experiments.
        o@experiments <- lapply(o@experiments, function(x) {
            geneId(x) <- data.table(ensembl = ids, id = ids)
            return(x)
        })

        # get the symbols
        cols <- paste0(species, c("ID", "GeneName"))
        ids <- homologs[[species]][, ..cols]
        ids <- unique(ids)
        ens <- na.omit(data.table(ensembl = ens[[1]], key = "ensembl"))

        ids <- merge(ens, ids, by.y = "mouseID", by.x = "ensembl", all.x=TRUE)

        ids$id <- ids$ensembl
        names(ids)[1:2] <- paste(species, c("ensembl", "symbol"), sep = ".")

        o$ids <- ids
        return(o)
    }

    # if we have homologs for more than the species we have, remove them
    if (length(species) < length(.validSpecies())) {
        colsToKeep <- as.character(sapply(species, paste0, c("ID", "GeneName", "Identical")))

        homologs <- lapply(homologs, function(x) {
            ctk <- colsToKeep[colsToKeep %in% names(x)]
            x <- x[, ..ctk]
            unique(x)
        })
    }

    # if there are more species, the id will be a combination of the ensembl id's of the species.
    # convert the ensemble from a character vecotr to a data.table with a key
    ens <- Map(function(x, nx) {
        d <- na.omit(data.table(x, key = "x"))
        names(d) <- paste0(nx, "ID")
        return(d)
    }, ens, names(ens))

    # join with homologs... data.table makes it FAST! Also, if an id doesn't exist in the other species, set it to NA
    ensJoined <- Map(function(x, nx) {
        nx2 <- paste0(nx, "ID")

        # In the homologs, only keep all the ensemble ids we have
        temp <- homologs[[nx]][x]

        otherSpecies <- species[!(species %in% nx)]

        # if the other species ID doesn't exist in the database, set it to NA
        for (os in otherSpecies) {
            osID <- paste0(os, "ID")

            # since we are using "NA" and not NA, set the other ones which are NA to "NA"
            temp[[osID]][is.na(temp[[osID]])] <- "NA"

            # now set the ones for which there are no homologs to "NA"
            toNA <- temp[[osID]] %in% ens[[os]][[osID]]
            temp[[osID]][!toNA] <- "NA"

            # if any are still missing, as in "", set them to "NA"
            temp[[osID]][temp[[osID]] == ""] <- "NA"
        }

        # now, remove those rows for which two conditions are fulfilled:
        # 1: An otherSpecies gene exists for the gene
        # 2: then, rows to remove = those where otherSpecies gene = "NA"
        # 3: If an otherSpecies gene doesn't exist, keep the gene with the highest otherSpeciesIdentical percentage
        temp <- split(temp, temp[[nx2]])
        temp <- lapply(temp, function(y) {
            if (nrow(y) == 1)
                return(y)  # if only one row, return it

            osID <- paste0(otherSpecies, "ID")
            identicalCol <- paste0(otherSpecies, "Identical")
            isNA <- y[, ..osID] == "NA"

            # exclude if both are NA
            isNA <- apply(isNA, 1, all)

            # Next condition: all "NA"... keep the one which is most identical
            if (sum(isNA) == nrow(y))
                y[which.max(rowSums(y[, ..identicalCol]))]
            else
                y[!isNA]
        })

        rbindlist(temp)

    }, ens, names(ens))

    # now, create the ids
    ensJoined <- Map(function(x, nx) {
        species <- species[match(species, .validSpecies())]
        allCols <- paste0(species, "ID")
        x$id <- do.call(paste, c(x[, ..allCols], sep="_"))
        setkey(x, "id")

        # cols to keep
        ctk <- c("id", paste0(nx, c("ID", "GeneName")))
        x <- x[, ..ctk]
        names(x) <- c("id", paste(nx, c("ensembl", "symbol"), sep="."))
        x
    }, ensJoined, names(ensJoined))

    # merge them
    ids <- Reduce(function(l1, l2) merge(l1, l2, by = "id", all=TRUE), ensJoined)

    # now add it to the experiments
    o@experiments <- lapply(o@experiments, function(x) {
        nx <- getSpecies(x)

        tCols <- c("id", paste(nx, "ensembl", sep="."))
        tIDs <- ids[, ..tCols]
        names(tIDs) <- c("id", "ensembl")

        geneId(x) <- tIDs

        return(x)
    })

    # make the compound symbol... there has to be a better way to do this, but my brain isn't working right now.
    cCols <- names(ids)[grepl("symbol", names(ids))]
    compound.symbol <- do.call(paste, c(ids[, ..cCols], sep="-"))

    compound.symbol <- strsplit(compound.symbol, "-")
    compound.symbol <- sapply(compound.symbol, function(cs) {
        if (any(cs %in% "NA")) {
            cs <- cs[!(cs %in% "NA")]
        }

        cs <- unique(cs)
        if (length(cs) == 1)
            return(cs)
        else
            return(paste(cs, collapse="-"))
    })

    ids$compound.symbol <- compound.symbol
    ids$compound.symbol <- ifelse(ids$compound.symbol != "", ids$compound.symbol, ids$id)

    o@ids <- ids

    return(o)

}
