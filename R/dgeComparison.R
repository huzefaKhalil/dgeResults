#' @include internalFunctions.R
#' @import data.table
#' @importFrom digest digest
#' @importFrom methods is
#' @importFrom methods new
#' @importFrom stats na.omit
#' @importFrom stats qnorm
#' @importFrom stats setNames
#' @importFrom utils read.csv
NULL

# Defines the formal class dgeComparison and associated functions

#' The \code{dgeComparison} class
#'
#' @description The basic class holds the DE results for a particular comparison along with meta-data.
#' The meta-data is stored in specific slots and getter methods can be used to access them.
#' The DE results assume output from limma-voom and calculate effect sizes and their variances
#' when the object is created. Methods are present to subset this object as well. Subsetting the object
#' means subsetting the DE results only; all the meta-data stays the same.
#'
#' @slot deRes A data.frame which contains the DE results. Has the following columns:
#'                 \enumerate{
#'                   \item id: internal gene id. Usually a concatenation of mouse and rat Ensembl
#'                   \item ensembl: Ensemble id
#'                   \item symbol: gene symbol
#'                   \item logFC: log fold-change
#'                   \item se: standard error
#'                   \item t: t-statistic
#'                   \item pvalue: associated pvalue
#'                   \item df.residual: residual degrees of freedom
#'                   \item df.total: total degrees of freedom
#'                   \item cohensD: effect size statistic Cohen's D
#'                   \item varD: variance associated with Cohen's D
#'                   \item hedgesG: effect size statistic Hedge's G
#'                   \item varG: variance associated with Hedge's G
#'                 }
#' @slot name The internal name of the experiment that this comparison is from
#' @slot model The model of the experiment. Possible values include: GRov, FSL, HR-LR, CSDS, Cort, Stress
#' @slot treatment If any treatments have been used in this particular comparison.
#' @slot timepoint The timepoint for an experiment with various timepoints.
#' @slot species The species include mouse, rat, human
#' @slot region The brain region. possible values include HPC, vHPC, dHPC, DG, vDG, dDG, PFC, NACC, STR
#' @slot sex Female or Male
#' @slot platform MicroArray or RNASeq
#' @slot group1 The first group name. The comparison ALWAYS givew the results for Group1 - Group2
#' @slot group2 The second group name
#' @slot n.group1 Group 1 N
#' @slot n.group2 GRoup 2 N
#' @slot n.effective Internally calculated in order to calculate effect sizes
#' @slot owner The system lab name of the owner. Only useful if it is a private data set
#' @slot public Logical flag. If true, then the data set can be viewed by all with access to the system.
#' @slot comparisonID The internal comparison ID
#'
#' @export
#'
setClass(
    "dgeComparison",
    slots = list(
        # First the per-gene data points
        deRes = "data.table",
        # this table has the following columns:
        # id: the internal DGE gene ID, made up of the concatenation of rat and mouse ensembl
        # ensembl: the ensemble id
        # symbol: gene symbol
        # logFC, se, t, pvalue, df.residual, df.total: DE results
        # cohensD, varD, hedgesG, varG: effect size measures calculated when the object is created

        # Now, the meta-data for the experiment
        name = "character",
        # The names by which the experiment is known in the DGE database
        model = "character",
        # GRov, HR-LR, FSL, CSDS, Cort, Env Enrichment, BDNF, etc.
        treatment = "character",
        # Any treatments... BDNF, WT, standard or enriched housing, drug treatments...
        timepoint = "character",
        # timepoint, if any
        species = "character",
        # one of Mouse or Rat (add Human)?
        region = "character",
        # enforce the naming convention of HPC, vHPC, dHPC, DG, vDG, dDG, PFC, NACC, STR... which others?
        sex = "character",
        # one of male / female
        platform = "character",
        # RNASeq or microarray
        group1 = "character",
        # The two groups. The comparison is always group1 - group2
        group2 = "character",
        n.group1 = "numeric",
        # The n's of the two groups
        n.group2 = "numeric",
        n.effective = "numeric",  # calculated internally. Used for effect size measures.

        owner = "character",
        public = "logical",

        comparisonID = "character"  # A unique id identifying this comparison
    )
)

# let's check the validity of the object
validDGEComparison <- function(object) {
    errors <- character()

    if (is.null(object@deRes)) {
        errors <-
            c(errors, "The differential expression results cannot be NULL.")
    }

    if (!any(grepl(tolower(object@species), .validSpecies()))) {
        errors <-
            c(errors, paste(
                "The species is invalid. Must be one of",
                paste(.validSpecies(), collapse = ", ")
            ))
    }

    if (!any(grepl(object@region, .validRegions()))) {
        errors <-
            c(errors,
              paste(
                  "The brain region is invalid. Must be one of",
                  paste(.validRegions(), collapse = ", ")
              ))
    }

    if (!any(grepl(tolower(object@sex), .validSex()))) {
        errors <-
            c(errors, paste("The sex is invalid. Must be one of", paste(sex, collapse = ", ")))
    }

    if (!any(grepl(object@platform, .validPlatforms()))) {
        errors <-
            c(errors,
              paste(
                  "The platform is invalid. Must be one of",
                  paste(.validPlatforms(), collapse = ", ")
              ))
    }

    # the others shouldn't be null
    if (any(
        is.null(object@group1),
        is.null(object@group2),
        is.null(object@n.group1),
        is.null(object@n.group2)
    )) {
        errors <- c(errors, "One of the groups is NULL.")
    }

    if (length(errors) == 0)
        TRUE
    else
        errors
}

setValidity("dgeComparison", validDGEComparison)

# getter methods for this class. Once set, we should need to have setter methods, unless its for id.
# First, create the generic functions we will need
getModel <- function(o) stop("Method undefined for an object of this class.")
getName <- function(o) stop("Method undefined for an object of this class.")
getTreatment <- function(o) stop("Method undefined for an object of this class.")
getTimepoint <- function(o) stop("Method undefined for an object of this class.")
getSpecies <- function(o) stop("Method undefined for an object of this class.")
getRegion <- function(o) stop("Method undefined for an object of this class.")
getPlatform <- function(o) stop("Method undefined for an object of this class.")
getSex <- function(o) stop("Method undefined for an object of this class.")
getGroups <- function(o) stop("Method undefined for an object of this class.")
getNs <- function(o) stop("Method undefined for an object of this class.")
getDERes <- function(o, groups = NULL) stop("Method undefined for an object of this class.")
getGeneEns <- function(o, ...) stop("Method undefined for an object of this class.")
getGeneId <- function(o) stop("Method undefined for an object of this class.")
getGeneSymbol <- function(o, ...) stop("Method undefined for an object of this class.")
getOwner <- function(o) stop("Method undefined for an object of this class.")
getPublic <- function(o) stop("Method undefined for an object of this class.")
geneId <- function(o) stop("Method undefined for an object of this class.")
`geneId<-` <- function(o, value) stop("Method undefined for an object of this class.")
reverseComparison <- function(o, ...) stop("Method undefined for an object of this class.")
flattenDge <- function(o, ...) stop("Method undefined for an object of this class.")
printComparison <- function(o, name=FALSE, region=FALSE, treatment=FALSE, timepoint=FALSE) stop("Method undefined for an object of this class.")

# now we set them as generic
setGeneric("getModel", function(o) { standardGeneric("getModel") })
setGeneric("getName", function(o) { standardGeneric("getName") })
setGeneric("getTreatment", function(o) { standardGeneric("getTreatment") })
setGeneric("getTimepoint", function(o) { standardGeneric("getTimepoint") })
setGeneric("getSpecies", function(o) { standardGeneric("getSpecies") })
setGeneric("getRegion", function(o) { standardGeneric("getRegion") })
setGeneric("getPlatform", function(o) { standardGeneric("getPlatform") })
setGeneric("getSex", function(o) { standardGeneric("getSex") })
setGeneric("getGroups", function(o) { standardGeneric("getGroups") })
setGeneric("getNs", function(o) { standardGeneric("getNs") })
setGeneric("getDERes", function(o, ...) { standardGeneric("getDERes") })
setGeneric("getGeneEns", function(o, ...) { standardGeneric("getGeneEns") })
setGeneric("getGeneId", function(o) { standardGeneric("getGeneId") })
setGeneric("getGeneSymbol", function(o, ...) { standardGeneric("getGeneSymbol") })
setGeneric("getOwner", function(o) { standardGeneric("getOwner") })
setGeneric("getPublic", function(o) { standardGeneric("getPublic") })
setGeneric("geneId", function(o) { standardGeneric("geneId") })
setGeneric("geneId<-", function(o, value) { standardGeneric("geneId<-") })
setGeneric("reverseComparison", function(o, ...) { standardGeneric("reverseComparison") })
setGeneric("flattenDge", function(o, ...) { standardGeneric("flattenDge") })
setGeneric("printComparison", function(o, name=FALSE, region=FALSE, treatment=FALSE, timepoint=FALSE) { standardGeneric("printComparison") })

# first the show method to print it on screen
#' Print an DGE Comparison
#'
#' Prints the basic meta-data from an dgeComparison
#'
#' @param dgeComparison
#'
#' @export
#'
setMethod("show", signature = "dgeComparison",
          function(object) {
              cat("Experiment Name:", object@name, "\n")
              cat("Model:", object@model, "\n")
              cat("Species:", object@species, "\n")
              cat("Brain Region:", object@region, "\n")
              cat("Treatment:", object@treatment, "\n")
              cat("Timepoint:", object@timepoint, "\n")
              cat("Sex:", object@sex, "\n")

              if (grepl("^-", object@comparisonID)) {
                  cat("Comparison:", object@group2, "-", object@group1, "\n")
                  cat("Ns: ", object@n.group2, " (", object@group2, ") ",
                      object@n.group1, " (", object@group1, ")\n", sep = "")
              } else {
                  cat("Comparison:", object@group1, "-", object@group2, "\n")
                  cat("Ns: ", object@n.group1, " (", object@group1, ") ",
                      object@n.group2, " (", object@group2, ")\n", sep = "")
              }

              cat("Number of Genes:", sum(!is.na(object@deRes$logFC)), "\n")
          })

# now the getters
#' Get the model name
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return the model name
#' @export
#'
setMethod("getModel", signature = "dgeComparison", function(o) o@model)

#' Get the experiment name
#'
#' Usually of the form LabName-experimentNumber. Is the internal DGE identifier for this particular experiment.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return the experiment name
#' @export
#'
setMethod("getName", signature = "dgeComparison", function(o) o@name)

#' Get the treatment for this comparison. If no treatment, "None" is returned.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return Character with the treatment name
#' @export
#'
setMethod("getTreatment", signature = "dgeComparison", function(o) o@treatment)

#' Get the timepoint for this comparison. If no timepoint, "None" is returned.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return Character with the timepoint name
#' @export
#'
setMethod("getTimepoint", signature = "dgeComparison", function(o) o@timepoint)

#' Get the species.
#'
#' One of "mouse", "rat", or "human". The last is not yet implemented.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return Character with the species name
#' @export
#'
setMethod("getSpecies", signature = "dgeComparison", function(o) o@species)

#' Gets the brain region.
#'
#' One of "HPC", "vHPC", "dHPC", "DG", "vDG", "dDG", "PFC", "NACC", "STR" or "VTA".
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return Character with brain region.
#' @export
#'
setMethod("getRegion", signature = "dgeComparison", function(o) o@region)

#' Gets the experiment platform.
#'
#' Either MicroArray or RNASeq. Maybe extend to others like ChIPSeq?
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return Character with the platform name.
#' @export
#'
setMethod("getPlatform", signature = "dgeComparison", function(o) o@platform)

#' Sex of the animals used in this comparison.
#'
#' "male" or "female" or "both"
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return Character vector with the sex.
#' @export
#'
setMethod("getSex", signature = "dgeComparison", function(o) o@sex)

#' Gets the groups for this comparison.
#'
#' Returns the names of the groups in the comparison. The comparison is always Group1 - Group2
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A names character vector of length two with the group names.
#' @export
#'
setMethod("getGroups", signature = "dgeComparison", function(o) {
    setNames(c(o@group1, o@group2), c(o@group1, o@group2))
})

#' Gets the number of animals per group.
#'
#' Is used to calculate the effect size. Not sure how useful it will be for users, but its there just in case.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A named integer vector of length two with the group sizes.
#' @export
#'
setMethod("getNs", signature = "dgeComparison", function(o) {
    setNames(c(o@n.group1, o@n.group2), c(o@group1, o@group2))
})

#' Returns the lab where this experiment was conducted. If it is private,
#' only lab members will have access to this data.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A character string with the owner of the data. Returns \code{NULL} if not set.
#' @export
#'
setMethod("getOwner", signature = "dgeComparison", function(o) o@owner)

#' Returns a boolean value which is TRUE if the data is universally accessible.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A logical value, \code{TRUE} if it is public, \code{FALSE} otherwise.
#' @export
#'
setMethod("getPublic", signature = "dgeComparison", function(o) o@public)

#' Get a vector of ensembl id's for the genes in this comparison.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A vector of ensemble ids
#' @export
#'
setMethod("getGeneEns", signature = "dgeComparison", function(o) o@deRes$ensembl)

#' Get a vector of gene symbols for this comparison.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A vector of gene symbols
#' @export
#'
setMethod("getGeneSymbol", signature = "dgeComparison", function(o) o@deRes$symbol)


#' Gets the vector of DGE id's for the genes in this comparison.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A vector of DGE ids
#' @export
#'
setMethod("getGeneId", signature = "dgeComparison", function(o) o@deRes$id)

#' Gets the vector of DGE id's for the genes in this comparison.
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A vector of DGE ids
#' @export
#'
setMethod("geneId", signature = "dgeComparison", function(o) o@deRes$id)

#' Sets the DGE gene ids for a `dgeComparison`
#'
#' @param o An object of class \code{dgeComparison} for which the values need to be set
#' @param value A data.table with ids and ensembl.
#'
#' @return An object of class `dgeComparison` with the ids set
#' @export
#'
setMethod("geneId<-", signature = "dgeComparison", function(o, value) {

    # we expect a data frame as a value
    if (!("data.table" %in% class(value)))
        stop("The value must be a data table.")

    if (!all(c("id", "ensembl") %in% names(value)))
        stop("The value must have columns ensembl and id.")

    # add rows... basically a right join... This will keep all ids for which no gene data exists.
    # TODO: a better way would be to remove the NA rows and keep only the data for which gene data exists.
    # Then, we would have to add that data in everytime a call is made to get data from this.
    # This would be much more memory efficient and wouldn't require any API changes.

    # starting to do the above TODO... if ensembl is NA, remove
    value <- value[!(value$ensembl == "NA"), ] ### Doesn't work!!

    temp <- o@deRes[value, on = c("ensembl")]
    temp$id <- temp$i.id

    temp <- temp[, -c("i.id")]

    o@deRes <- temp
    setkey(o@deRes, id)

    return(o)

})

#' Gets the DE Results
#'
#' This is a data frame containing the per gene DE results for this particular comparison.
#'
#' @param o An object of class \code{dgeComparison} from which to return the DE result
#' @param groups An optional groups vector. If the groups don't correspond to the comparison groups, an error is thrown.
#' However, if the groups are present and correct and the order of the groups is reversed, the DE results returned are reversed also.
#'
#' @return A data frame containing the DE results and effect sizes for all genes.
#' @export
#'
setMethod("getDERes", signature = "dgeComparison", function(o, groups = NULL) {
    if (is.null(groups)) return(o@deRes)

    if (length(groups) != 2) stop("If specified, there can be only two groups.")

    g <- getGroups(o)
    names(g) <- names(groups) <- NULL

    if (identical(groups, g)) {
        return(o@deRes)
    } else if (identical(groups, rev(g))) {
        deRes <- o@deRes
        deRes$logFC <- deRes$logFC *  -1
        deRes$cohensD <- deRes$cohensD * -1
        deRes$hedgesG <- deRes$hedgesG * -1
        return(deRes)
    } else {
        stop("Unknown groups.")
    }
})

#' Convenience method to reverse the direction of a comparison
#'
#' @param o An object of class \code{dgeComparison}
#' @param comparisons Ignored for this class
#'
#' @return An object of class dgeComparison
#' @export
#'
setMethod("reverseComparison", signature = "dgeComparison", function(o) {

    # groups <- getGroups(o)
    # Ns <- getNs(o)
    #
    # groups <- rev(groups)
    # Ns <- rev(Ns)
    #
    # o@deRes <- getDERes(o, groups)
    # o@group1 <- groups[1]
    # o@group2 <- groups[2]
    # o@n.group1 <- Ns[1]
    # o@n.group2 <- Ns[2]

    # just set id to negative here
    if (grepl("^-", o@comparisonID)) {
        o@comparisonID <- gsub("^-", "", o@comparisonID)
    } else {
        o@comparisonID <- paste0("-", o@comparisonID)
    }
    return(o)
})

#' Returns a data.table with all relevant information
#'
#' @param o An object of class \code{dgeComparison} from which to get the data
#' @param geneId An optional vector of geneId for which to get the information
#'
#' @return A data.table with all relevant information
#' @export
#'
setMethod("flattenDge", signature = "dgeComparison", function(o, geneId = NULL) {

    if (is.null(geneId)) {
        deRes <- o@deRes
    } else {
        if (!haskey(o@deRes)) stop("The id's haven't yet bee created!")
        deRes <- o@deRes[geneId, ]
    }

    if (nrow(deRes) == 0) {
        return(data.table)
    }

    if (grepl("^-", o@comparisonID)) {
        # change the sign of the numbers
        toRev <- c("logFC", "t", "cohensD", "hedgesG")
        for (t in toRev)
            deRes[[t]] <- -deRes[[t]]

        out <- data.table(deRes,
                          name = o@name,
                          model = o@model,
                          platform = o@platform,
                          treatment = o@treatment,
                          timepoint = o@timepoint,
                          region = o@region,
                          sex = o@sex,
                          species = o@species,
                          comparison = paste(o@group2, "-", o@group1),
                          n.group1 = o@n.group2,
                          n.group2 = o@n.group1,
                          n.effective = o@n.effective,
                          owner = o@owner,
                          public = o@public,
                          comparisonID = o@comparisonID)
    } else {
        out <- data.table(deRes,
                          name = o@name,
                          model = o@model,
                          platform = o@platform,
                          treatment = o@treatment,
                          timepoint = o@timepoint,
                          region = o@region,
                          sex = o@sex,
                          species = o@species,
                          comparison = paste(o@group1, "-", o@group2),
                          n.group1 = o@n.group1,
                          n.group2 = o@n.group2,
                          n.effective = o@n.effective,
                          owner = o@owner,
                          public = o@public,
                          comparisonID = o@comparisonID)
    }

    # remove any rows for which the estimate doesn't exist
    return(out[!is.na(out$logFC)])
})

#' Print the comparison, namely Group1 - Group2
#'
#' @param o An object of class \code{dgeComparison}
#' @param name Logical, indicating if the name of the experiment should be printed. Default \code{FALSE}.
#' @param region Logical, indicating if the brain region of the comparison should be printed. Default \code{FALSE}.
#' @param treatment Logical, indicating if the treatment of the comparison should be printed. Default \code{FALSE}.
#' @param timepoint Logical, indicating if the timepoint of the comparison should be printed. Default \code{FALSE}.
#'
#' @return String with the specified comparison
#' @export
#'
setMethod("printComparison", signature = "dgeComparison", function(o,
                                                                   name = FALSE,
                                                                   region = FALSE,
                                                                   treatment = FALSE,
                                                                   timepoint = FALSE) {
    addedVars <- c(name, region, treatment, timepoint)
    if (!is.logical(addedVars))
        stop("In method 'printComparison', some non-logical variables were passed as arguments.")

    if (grepl("^-", o@comparisonID)) {
        outName <- paste(o@group2, "-", o@group1)
    } else {
        outName <- paste(o@group1, "-", o@group2)
    }

    if (any(addedVars)) {

        names(addedVars) <- c(o@name, o@region, o@treatment, o@timepoint)
        toAdd <- names(addedVars)[addedVars]

        # None should also be printed...
        #toAdd <- toAdd[toAdd != "None"]

        if (length(toAdd) > 0)
            outName <- paste0(outName, " (", paste(toAdd, collapse = ", "), ")")
    }

    return(setNames(outName, o@comparisonID))
})

#' Returns the number of genes in the DE results
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return Integer with the number of genes.
#' @export
#'
setMethod("nrow", signature = "dgeComparison", function(x) nrow(x@deRes))

#' Overloading names to get the name of the comparison
#'
#' @param o An object of class \code{dgeComparison}
#'
#' @return A string with the name
#' @export
#'
setMethod("names", signature = "dgeComparison", function(x) getName(x))

#' Subset the DE results of this comparison.
#'
#' It doesn't effect the meta-data in any way. It mirrors base R's subset functions for a data.frame.
#'
#' @param x An object of class \code{dgeComparison}
#'
#' @return An dgeComparison object with the specified DE results.
#' @export
#'
setMethod("[", signature = "dgeComparison", function(x, i) {

    # This only works if the key has been set or if i is logical
    if (!(is.logical(i) || haskey(x@deRes)))
        stop("The 'id' for the data set has not been set OR the arument is not logical.")

    x@deRes <- x@deRes[i, nomatch=NULL]
    return(x)
})

#' Returns a column of the DE results data frame.
#'
#' Perhaps think about making it convenient to also return metadata this way??
#'
#' @param x An object of class \code{dgeComparison}
#'
#' @return A column of the DE results data frame.
#' @export
#'
setMethod("$", signature = "dgeComparison", function(x, name) {
    x@deRes[[name, exact = TRUE]]
})

#' Constructor function for the dgeComparison class.
#'
#' The function takes in the DE results or the file name containing the DE results, along with
#' the meta-data and checks it for any errors. The effect sizes are also calculated at this stage and
#' the resultant dgeComparison object is returned. The effect sizes calculated are Cohen's D and Hedge's G
#' along with their variances.
#'
#' @param deData A file name or a data frame containing the DE results. if a file name, it must be a CSV file.
#'               The following columns are expected in the file or data frame:
#'               \itemize{
#'                 \item ensembl: The Ensembl id for the gene
#'                 \item symbol: The gene symbol
#'                 \item logFC: the lof fold change
#'                 \item se: the standard error for the log fold change
#'                 \item t: the t-statistic for the gene
#'                 \item pvalue: the associated pvalue
#'                 \item df.residual: the residual degrees of freedom
#'                 \item df.total: the total degrees of freedom
#'               }
#' @param name The name of the experiment.
#' @param model The model name.
#' @param treatment The treatment, if any. If absent, treatment will be specified as "None".
#' @param species The species. One of "mouse", "rat" or "human"
#' @param region The brain region
#' @param sex The sex. "female", "male" or "both"
#' @param platform The experimental platform. MIcroArray or RNASeq.
#' @param groups The two groups. A character vector of length two is expected.
#' @param Ns The number of samples in each group. Integer vector of length two is expected.
#'
#' @return A dgeComparison object.
#' @export
#'
dgeComparison <- function(deData, name, model, treatment = "None",
                          timepoint = "None", species, region, sex,
                          platform, groups, Ns, owner = "", public = TRUE) {

    # treatment is the only optional argument. All others must be specified and be strings and vectors. make sure they are
    stopifnot(class(deData) %in% c("data.frame", "character", "data.table"),
              typeof(name) == "character",
              typeof(model) == "character",
              typeof(species) == "character",
              typeof(region) == "character",
              typeof(sex) == "character",
              typeof(platform) == "character",
              typeof(treatment) == "character",
              typeof(timepoint) == "character",
              length(groups) == 2,
              typeof(groups) == "character",
              length(Ns) == 2,
              class(Ns) == "numeric")

    # assume the DE Data is correct...
    if (is(deData, "character")) {
        stopifnot(file.exists(deData))

        deData <- read.csv(deData)
    }

    # Calculate the effect size measures
    n.effective <- (Ns[1] * Ns[2]) / (Ns[1] + Ns[2])
    nFactor <- (gamma(deData$df.total / 2) /
                    (sqrt(deData$df.total / 2) *
                         (gamma((deData$df.total - 1) / 2))))

    cohensD <- deData$t / sqrt(n.effective)
    varD <- (deData$df.total / ((deData$df.total - 2) * n.effective) *
                 (1 + (n.effective * (cohensD ^ 2)))) -
        ((cohensD ^ 2) / (nFactor ^ 2))

    hedgesG <- cohensD * nFactor
    varG <- varD * (nFactor ^ 2)

    comparisonID <- digest::digest(deData)

    # I asusme a csv file or data.frame with the following fields
    # LogFC, SE, df.residual. t, F, p.value. F.p.value, df.total, Symbol, Ensembl
    new("dgeComparison",
        deRes = data.table(
            id = rep(NA, nrow(deData)),
            ensembl = as.character(deData$ensembl),
            symbol = as.character(deData$symbol),
            logFC = deData$logFC,
            se = deData$SE,
            t = deData$t,
            pvalue = deData$pvalue,
            fdr = deData$fdr,
            df.residual = deData$df.residual,
            df.total = deData$df.total,
            cohensD = cohensD,
            varD = varD,
            hedgesG = hedgesG,
            varG = varG,
            stringsAsFactors = FALSE
        ),

        name = name,
        model = model,
        treatment = treatment,
        timepoint = timepoint,
        species = tolower(species),
        region = region,
        sex = tolower(sex),
        platform = platform,
        group1 = groups[1],
        group2 = groups[2],
        n.group1 = Ns[1],
        n.group2 = Ns[2],
        n.effective = n.effective,

        owner = owner,
        public = public,

        comparisonID = comparisonID
    )
}



