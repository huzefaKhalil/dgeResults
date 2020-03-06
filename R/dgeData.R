#' @include dgeExperiment.R
#' @importFrom DBI dbGetQuery
#' @importClassesFrom RSQLite SQLiteConnection
#' @importFrom RSQLite SQLite
NULL

# this file specifies the dgeData class which is just a collection of dgeExperiments with associated methods.

#' Class `dgeData`
#'
#' Has only one slot, a list with all the experiments that form this DGE data set.
#'
#' @slot experiments list of `dgeExperiments`
#' @slot ids A data.table of unique DGE ids present in this data set along with the species ensembl ids and symbols
#'
#' @export
#'
#' @examples
setClass(
    "dgeData",

    slots = list(
        experiments = "list",
        ids = "data.table"
    )
)

validDGEData <- function(object) {

    hClass <- sapply(object@experiments, class)
    if (!all(hClass == "dgeExperiment"))
        return("All objects are not of class dgeExperiment.")
    else
        TRUE
}

setValidity("dgeData", validDGEData)

# other generics

#' Print the experiments names in the data set.
#'
#' @param dgeData
#'
#' @export
#'
setMethod("show", signature = "dgeData", function(object) {
    cat("Experiment list:\n")
    for (exp in object@experiments) {
        cat(getName(exp), "\n")
    }
})

#' Get Species in this data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector of species in the data set. Currently limited to "mouse", "rat" and "human"
#' @export
#'
#' @examples getSpecies(dgeDataSet)
setMethod("getSpecies", signature = "dgeData", function(o) unique(sapply(o@experiments, getSpecies)))

#' Get the brain regions present in the DGE data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector of brain regions present in this DGE data set.
#' @export
#'
#' @examples
setMethod("getRegion", signature = "dgeData", function(o) unique(unlist(sapply(o@experiments, getRegion))))

#' Get the sexes present in this DGE data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector of the sexes which are included in this data set. Can include "both" if some comparisons are across sexes.
#' @export
#'
#' @examples
setMethod("getSex", signature = "dgeData", function(o) unique(unlist(sapply(o@experiments, getSex))))

#' Get the treatments present in this data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector of the treatments in this data set.
#' @export
#'
#' @examples
setMethod("getTreatment", signature = "dgeData", function (o) unique(unlist(sapply(o@experiments, getTreatment))))

#' Get the treatments present in this data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector of the timepoint in this data set.
#' @export
#'
#' @examples
setMethod("getTimepoint", signature = "dgeData", function (o) unique(unlist(sapply(o@experiments, getTimepoint))))

#' Get the platforms present in this data set.
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector of platforms, either Microarray or RNASeq or both.
#' @export
#'
#' @examples
setMethod("getPlatform", signature = "dgeData", function(o) unique(unlist(sapply(o@experiments, getPlatform))))

#' Get the models present in this data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector with the model names.
#' @export
#'
#' @examples
setMethod("getModel", signature = "dgeData", function(o) unique(unlist(sapply(o@experiments, getModel))))

#' Get the experiment names present in this data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A character vector with the experiment names.
#' @export
#'
#' @examples
setMethod("getName", signature = "dgeData", function(o) unique(unlist(sapply(o@experiments, getName))))

#' Get a vector of groups
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A string vector of groups which are a part of the data set.
#' @export
#'
#' @examples
setMethod("getGroups", signature = "dgeData", function(o) unique(unlist(sapply(o@experiments, getGroups))))

#' Get the internal comparison IDs of the comparisons in this data base
#'
#' @param o An object of class \code{dgeData} from which to get the IDs
#'
#' @return A character vector of IDs
#' @export
#'
#' @examples
setMethod("getIds", signature = "dgeData", function(o) unlist(lapply(o@experiments, getIds)))

#' Print the comparisons in this data object.
#'
#' @param o An object of class \code{dgeData} for which to print the comparisons
#' @param name Logical, indicating if the name of the experiment should be printed. Default \code{FALSE}.
#' @param region Logical, indicating if the brain region of the comparison should be printed. Default \code{FALSE}.
#' @param treatment Logical, indicating if the treatment of the comparison should be printed. Default \code{FALSE}.
#' @param timepoint Logical, indicating if the timepoint of the comparison should be printed. Default \code{FALSE}.
#'
#' @return A vector of strings with the comparisons. The vector is named by the ids of the comparisons
#' @export
#'
#' @examples
setMethod("printComparison", signature = "dgeData", function(o,
                                                             name = FALSE,
                                                             region = FALSE,
                                                             treatment = FALSE,
                                                             timepoint = FALSE)
    unlist(
        lapply(
            o@experiments,
            printComparison,
            name = name,
            region = region,
            treatment = treatment,
            timepoint = timepoint
        )
    ))

#' Get the ensemble id's for a particular species in a DGE data set
#'
#' @param o An object of class \code{dgeData} from which to get id's
#' @param species The species for which to get id's. Can be a vector.
#'
#' @return A list with one element for each species giving the unique Ensembl id's present in each experiment.
#' @export
#'
#' @examples
#' getGeneEnsembl(dataSet, "mouse")
#' getGeneEnsembl(dataSet, c("Mouse", "Rat"))
setMethod("getGeneEns", signature = "dgeData", function(o, species) {

    species <- tolower(species)
    if (!all(species %in% .validSpecies())) {
        stop("One or more species are not supported.")
    }

    ens <- lapply(species, function(sp) {
        o1 <- getExperimentsBySpecies(o, sp)

        if (is.null(o1))
            NULL
        else
            unique(unlist(lapply(o1@experiments, getGeneEns), use.names = FALSE))
    })

    names(ens) <- species

    return(ens)
})

#' Get the gene symbols for a particular species in a DGE data set
#'
#' @param o An object of class \code{dgeData} from which to get the symbols
#' @param species The species for which to get the symbols. Can be a vector.
#'
#' @return A list with one element for each species giving the unique gene symbols present in each experiment.
#' @export
#'
#' @examples
#' getGeneEnsembl(dataSet, "mouse")
#' getGeneEnsembl(dataSet, c("Mouse", "Rat"))
setMethod("getGeneSymbol", signature = "dgeData", function(o, species) {

    species <- tolower(species)
    if (!all(species %in% .validSpecies())) {
        stop("One or more species are not supported.")
    }

    ens <- lapply(species, function(sp) {
        o1 <- getExperimentsBySpecies(o, sp)

        if (is.null(o1))
            NULL
        else
            unique(unlist(lapply(o1@experiments, getGeneSymbol), use.names = FALSE))
    })

    names(ens) <- species

    return(ens)
})

#' Get ids for a data set
#'
#' @param o An object of class \code{dgeData}
#'
#' @return a character vector of internal DGE ids
#' @export
#'
#' @examples
setMethod("geneId", signature = "dgeData", function(o) o@ids$id)

# generic functions we will use
getExperiments <- function(o) stop("Method undefined for an object of this class.")
getExperimentsBySpecies <- function(o, species) stop("Method undefined for an object of this class.")
getColumnData <- function(o, column, comparisonID, conn) stop("Method undefined for an object of this class.")
symbolToId <- function(o, ...) stop("Method undefined for an object of this class.")

setGeneric("getExperiments", function(o, species) { standardGeneric("getExperiments") })
setGeneric("getExperimentsBySpecies", function(o, species) { standardGeneric("getExperimentsBySpecies") })
setGeneric("getColumnData", function(o, column, comparisonID, conn) { standardGeneric("getColumnData") })
setGeneric("symbolToId", function(o, ...) { standardGeneric("symbolToId") })

#' Get all experiments in a data set.
#'
#' @param o An object of class \code{dgeData}
#'
#' @return A list of objects of class \code{dgeExperiment}
#' @export
#'
#' @examples
setMethod("getExperiments", signature = "dgeData", function(o) o@experiments )

#' Gets the experiments of a particular species in an DGE data set
#'
#' @param o An object of class \code{dgeData}
#' @param species ONE species for which to get the experiments
#'
#' @return \code{NULL} if there are no experiments for that species. Else an object of class \code{dgeData} with experiments from \code{species}
#' @export
#'
#' @examples
setMethod("getExperimentsBySpecies", signature = "dgeData", function(o, species) {

    species <- tolower(species)
    if (length(species) != 1 || !(species %in% .validSpecies()))
        stop("Species not valid or more than one specified.")

    experiments <- lapply(o@experiments, function(e) {
        if (getSpecies(e) == species)
            e
        else
            NULL
    })

    # remove NULL values
    experiments <- experiments[lengths(experiments) != 0]

    if (is.null(experiments))
        NULL
    else
        new("dgeData", experiments = experiments, ids = o@ids)
})

#' Get column data for particular comparisons
#'
#' @param o An object of class \code{dgeData}
#' @param column This can be any column name in the \code{data.table} returned by \code{flattenDge}. These include
#'               \enumerate{
#'                   \item logFC
#'                   \item se
#'                   \item pvalue
#'                   \item fdr
#'                   \item cohensD
#'                   \item varD
#'                   \item hedgesG
#'                   \item varG
#'               }
#' @param comparisonID The comparison IDs for which the data is needed
#' @param conn The connection to the SQLite DB if the data is stored there
#'
#' @return A data.table with columns as comparisonIDs
#' @export
#'
#' @examples
setMethod("getColumnData", signature = "dgeData", function(o, column, comparisonID, conn = NULL) {

    validCols <- c("logFC", "se", "pvalue", "fdr", "cohensD", "varD", "hedgesG", "varG")

    if (length(column) != 1)
        stop("The column supplied to method 'getColumnData' isn't of length 1.")
    if (!(column %in% validCols))
        stop("In method 'getColumnData', the supplied column isn't valid. Must be one of ", paste(validCols, collapse=", "), " but is ", column)

    # get the data using flattenDge
    if (!is.null(comparisonID) || length(comparisonID) != 0) # if any of these are true, assume that the data set has only required ids
        o <- getComparisonById(o, comparisonID)

    fd <- flattenDge(o, conn = conn)

    # make it wide
    pMat <- dcast(fd, id ~ comparisonID, value.var = column)
    pMat$symbol <- o@ids[pMat]$compound.symbol

    return(pMat)
})

#' Converts human gene symbols to the internal IDs used by the package.
#'
#' @param o An object of class \code{dgeData} from which to extract the IDs
#' @param geneSymbols The symbols for which we need the IDs
#' @param species The species which the symbols correpond to. If unspecified, the function returns matches across all species.
#'
#' @return A character vector with the gene ids which can be used to subset the data
#' @export
#'
#' @examples
setMethod("symbolToId", signature = "dgeData", function(o, geneSymbols, species = NULL) {

    if (is.null(species))
        species <- getSpecies(o)

    species <- tolower(species)
    if (!all(species %in% .validSpecies()))
        stop("In the function 'symbolToId', the species specified were invalid.")

    if (!all(species %in% getSpecies(o))) {
        warning("Species specified to function 'symbolToId' are not present in this data set.")
        return(NULL)
    }

    sCols <- paste(species, "symbol", sep=".")

    ids <- rbindlist(lapply(sCols, function(s) o@ids[tolower(o@ids[[s]]) %in% tolower(geneSymbols)] ))

    # remove dups
    return(unique(ids))

})

#' Subset by gene id
#'
#' This function returns a subsetted dgeData object containing all the experiments but only the selected genes for which an id is given.
#' i can be a vector of id's or a logical vector with the same number of rows as the data set.
#'
#' @param x An object of class \code{dgeData}
#'
#' @return and object of class `dgeData`
#' @export
#'
#' @examples `dge[c(TRUE, FALSE, FALSE, TRUE)]`
setMethod("[", signature = "dgeData", function(x, i) {
    x@experiments <- lapply(x@experiments, function(exp) {
        exp[i, ]
    })
    return(x)
})

#' Reverses a set of comparisons in the data
#'
#' @param o An object of class \code{dgeData}
#' @param comparisons A character vector with the list of comparisons to be reversed
#'
#' @return
#' @export
#'
#' @examples
setMethod("reverseComparison", signature = "dgeData", function(o, comparisons) {
    o@experiments <- lapply(o@experiments, reverseComparison, comparisons = comparisons)
    return(o)
})

#' Flatten a data set to a data table
#'
#' @param o An object of class \code{dgeData}
#' @param geneId Information for only these id's is returned
#' @param conn A connection to an SQLite DB. Has to be of class \code{SQLiteConnection}.
#'
#' @return A data.table conatining all the information in the experiment.
#' @export
#'
#' @examples
setMethod("flattenDge", signature = "dgeData", function(o, geneId = NULL, conn = NULL) {
    if (is.null(conn))
        return(rbindlist(lapply(o@experiments, flattenDge, geneId = geneId)))
    else {

        # now we have to get it from a DB...

        if (!is(conn, "SQLiteConnection"))
            stop("In function 'flattenDge', a connection object is given but it is not of the class 'SQLiteConnection'.")

        # so, connection exists and the object 'o' is assumed to be a lite object
        comparisons <- names(printComparison(o))

        # check which ones are reversed
        toRev <- which(grepl("^-", comparisons))

        # first, get the comparisons.
        selectedComparisons <- gsub("^-", "", comparisons)
        cPlaceholder <- paste0("$comparison", seq_along(selectedComparisons))

        # check if genes are specified... if not, then we select only by comparison ID
        if (is.null(geneId) || length(geneId) == 0) {
            sqlQuery <- paste0("SELECT *
                           FROM comparisonData
                               WHERE comparisonID IN (", paste(cPlaceholder, collapse = ","), ")")
            cParams <- as.list(setNames(selectedComparisons, gsub("\\$", "", cPlaceholder)))
            paramList <- cParams
        } else {
            selectedGenes <- geneId   # Not checking here to ensure that these are valid ids... Assuming they will be
            gPlaceholder <- paste0("$gene", seq_along(selectedGenes))
            sqlQuery <- paste0("SELECT *
                           FROM comparisonData
                               WHERE comparisonID IN (", paste(cPlaceholder, collapse = ","), ")
                               AND id IN (", paste(gPlaceholder, collapse = ","), ")")

            cParams <- as.list(setNames(selectedComparisons, gsub("\\$", "", cPlaceholder)))
            gParams <- as.list(setNames(selectedGenes, gsub("\\$", "", gPlaceholder)))
            paramList <- c(cParams, gParams)
        }

        outDT <- data.table(dbGetQuery(conn, sqlQuery, param = paramList))
        setkey(outDT, id)

        # now, reverse the comparisons in the data
        if (length(toRev) == 0) {
            return(outDT)
        } else {
            return(.revDT(outDT, selectedComparisons[toRev]))
        }
    }
})

#' Gets the comparisons within this data set by id.
#'
#' @param o An object of class \code{dgeData} which we want to subset
#' @param id A character vector with the ids we need
#'
#' @return An \code{dgeData} object with the needed ids. If none found, NULL is returned.
#' @export
#'
#' @examples
setMethod("getComparisonById", signature = "dgeData", function(o, id) {
    exp <- lapply(o@experiments, getComparisonById, id)

    # remove NULL values
    exp <- exp[lengths(exp) != 0]

    o@experiments <- exp

    if (length(exp) == 0)
        return(NULL)
    else
        return(o)
})

#' Get one or many comparisons from a `dgeData` object.
#'
#' The comparisons matching the given parameters are returned from the given data set. If groups are specified, they need to be a
#' character vector of length two. If need be, the comparison will be reversed. One of groups, treatment or region
#' needs to be specified. If both treatment and region are specified, the comparison will be matched using both parameters.
#'
#' @param o An object of class \code{dgeData} from which the comparison is to be extracted
#' @param groups The groups for which the comparison is needed. If specified, the other two parameters are ignored
#' @param treatment Optional parameter specifying the treatment needed for the comparison. Can be used in conjunction with region
#' @param region The brain region of the desired comparison. Can be specified in conjunction with treatment
#' @param sex The sex of the desired comparison
#' @param species The desired species
#' @param model Select only specified models
#' @param name Select only experiments with specified name
#'
#' @return If no comparison is matched, `NULL` is returned, else an object of class `dgeData` with the specified comparisons.
#' @export
#'
#' @examples
setMethod("getComparison", signature = "dgeData",
          function(o, groups=NULL, treatment=NULL, timepoint = NULL, region=NULL, sex=NULL, species = NULL, model=NULL, name=NULL)  {
              #TODO Add select by model and experiment name as well.

              # we either need both groups or treatment or region or sex)nd
              stopifnot((!is.null(groups) && length(groups) == 2) ||
                            !is.null(treatment) ||
                            !is.null(timepoint) ||
                            (!is.null(region) && region %in% .validRegions()) ||
                            (!is.null(sex) && tolower(sex) %in% .validSex()) ||
                            (!is.null(species) && tolower(species) %in% .validSpecies()) ||
                            !is.null(model) ||
                            !is.null(name)
              )

              ex <- o@experiments

              if (!is.null(species)) {
                  t <- sapply(ex, getSpecies)
                  t <- t %in% tolower(species)
                  if (!any(t))
                      return(NULL)

                  ex <- ex[t]
              }

              if (!is.null(model)) {
                  t <- sapply(ex, getModel)
                  t <- t %in% model
                  if (!any(t))
                      return(NULL)

                  ex <- ex[t]
              }

              if (!is.null(name)) {
                  t <- sapply(ex, getName)
                  t <- t %in% name
                  if (!any(t))
                      return(NULL)

                  ex <- ex[t]
              }

              # if groups, treatment, region and sex are unspecified,
              # we don't neet to go through the individual experiments
              if (!is.null(groups) || !is.null(treatment) || !is.null(region) ||
                  !is.null(sex) || !is.null(timepoint)) {
                  # get a list of experiments back with the correct comparisons
                  ex <- lapply(ex, getComparison,
                               groups = groups,
                               treatment = treatment,
                               timepoint = timepoint,
                               region = region,
                               sex = sex)
              }

              experiments <- as.list(unlist(ex, use.names = FALSE))  #Faster way to remove NULL values

              # If no comparisons are found, return NULL
              if (length(experiments) == 0)
                  return(NULL)

              # get the unique ids
              #ids <- unique(unlist(lapply(experiments, geneId), use.names = FALSE))
              #ids <- o@ids[ids]

              # return the new dgeData object.
              #TODO use the constructor and create new IDs here. If we remove one species, the ids will change.
              #but, that might not be the best idea. Let's keep the id's stable once the entire data set has been created.
              o@experiments <- experiments
              return(o)
          })

#' Constructor for the `dgeData` class. Creates the dataset ids when called.
#'
#' @param experiments  The list of experiments. Must be of class `dgeExperiment`.
#'
#' @return A new `dgeData` object.
#' @export
#'
#' @examples
dgeData <- function(experiments) {
    # Just a note! Since this calls the .createIds method, NEVER USE IT INTERNALLY!
    # to initialize a new `dgeData` internally, always call the `new` method.
    hd <- new("dgeData",
              experiments = experiments,
              ids = data.table())
    .createIds(hd)
}
