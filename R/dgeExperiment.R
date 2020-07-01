#' @include dgeComparison.R
NULL

# Define the class dgeExperiment. This will contain all the comparisons associated with a
# given experiment as well as an optional description of the experiment.

#' The \code{dgeExperiment} class
#'
#' @description The dgeExperiment class holds a list of dgeComparisons. The comparisons belong to the same
#' experiment, so have the same name. They should also have the same species and the same model name
#' and platform. What can change is the treatments, the brain regions and the sex. This class should list ALL
#' The treatments, regions and sexes involved in this experiment. In fact, this is a requirement when creating
#' this class from a list of comparisons.
#'
#' @slot comparisons A list of objects of the class dgeComparison which comprise this experiment. Needs to be at lest 1.
#' @slot name The name of the experiment. Usually lab name with an hyphenated number.
#' @slot model The model name.
#' @slot species The species.
#' @slot platform The platform.
#' @slot treatments A character vector with the list of treatments for each comparison. If there aren't any then "None"
#' @slot timepoints A character vector with the timepoints for each comparison. If there aren't any, then "None"
#' @slot regions A character vector which specifies the brain regions in this particular experiment
#' @slot sexes The sex: male or female or both
#' @slot groups The different groups in the experiment
#' @slot totalN The total number of animals in this experiment
#' @slot description Text describing the experiment. This is optional but highly recommended.
#' @slot owner The system username of the owner. Only useful if it is a private data set
#' @slot public Logical flag. If true, then the data set can be viewed by all with access to the system.
#'
#' @export
#'
setClass(
    "dgeExperiment",

    slots = list(
        comparisons = "list",      # This is the list of comparisons including the DE results from the experiment

        name = "character",        # The experiment attributes. can be copied from the comparisons themselves
        model = "character",
        species = "character",
        platform = "character",

        treatments = "character",  # Vectors of all treatmenst, regions and sexes involved in the experiment
        timepoints = "character",
        regions = "character",
        sexes = "character",
        groups = "character",
        totalN = "numeric",

        description = "character",  # An optional description of the experiment
        owner = "character",  # The username of the owner, if any. Only important if private
        public = "logical"    # If the dataset is public or not.
    )
)

validDGEExperiment <- function(object) {
    errors <- c()

    # all comparisons must be of the class dgeComparison
    hClass <- sapply(object@comparisons, class)
    if (!all(hClass == "dgeComparison")) {
        errors <- c(errors, "All comparisons must be of class \'dgeComparison\'.")
    }

    # Check to ensure that the name, model, species and platform matches that of all the comparisons
    hName <- sapply(object@comparisons, getName)
    if (!all(hName == object@name)) {
        errors <- c(errors, "The comparisons have different names.")
    }

    hSpecies <- sapply(object@comparisons, getSpecies)
    if (!all(hSpecies == object@species)) {
        errors <- c(errors, "The comparisons have different species.")
    }

    hPlatform <- sapply(object@comparisons, getPlatform)
    if (!all(hPlatform == object@platform)) {
        errors <- c(errors, "The comparisons have different platforms.")
    }

    totalN <- sum(sapply(object@comparisons, getNs))
    if (!identical(totalN, object@totalN)) {
        errors <- c(errors, "The total N doen't correspond with the N in the comparisons.")
    }

    # next, check the groups, regions, treatments, sexes... but I don't think we need to check them because
    # the constructor only uses the comparisons to construct an experiment and so the variables will
    # only be derived from the comparisons in the first place.

    if (length(errors) == 0)
        TRUE
    else
        errors
}

setValidity("dgeExperiment", validDGEExperiment)

# now the generic methods which work on objects of this class.

#' Show method.
#'
#' Prints the experiment
#'
#' @param dgeExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("show", signature = "dgeExperiment",
          function(object) {
              cat("Experiment Name:", object@name, "\n")
              cat("Model:", object@model, "\n")
              cat("Species:", object@species, "\n")
              cat("Brain Regions:", object@regions, "\n")
              cat("Treatments:", object@treatments, "\n")
              cat("Timepoints:", object@timepoints, "\n")
              cat("Sex:", object@sexes, "\n")
              cat("Number of comparisons:", length(object@comparisons), "\n")
              cat("Groups:", object@groups, "\n")
              if (!is.null(object@owner)){
                  cat("Owner:", object@owner, "\n")
              }
              cat("Public:", ifelse(object@public, "Yes", "No"), "\n")
          })

#' Returns the experiment model
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character with the experiment model
#' @export
#'
#' @examples
#' getModel(dgeExperiment)
setMethod("getModel", signature = "dgeExperiment", function(o) o@model)

#' Returns the experiment name
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return The experiment name
#' @export
#'
#' @examples
#' getName(dgeExperiment)
setMethod("getName", signature = "dgeExperiment", function(o) o@name)

#' The species.
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return The species is returned
#' @export
#'
#' @examples
#' getSpecies(dgeExperiment)
setMethod("getSpecies", signature = "dgeExperiment", function(o) o@species)

#' The platform: microarray or rna-seq
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return String specifying the platform
#' @export
#'
#' @examples
#' getPlatform(dgeExperiment)
setMethod("getPlatform", signature = "dgeExperiment", function(o) o@platform)

#' Get a vector of groups
#'
#' A possible complication. Should we ensure that each and every combination of groups has a comparison present?
#' Don't think we need to make this a requirement, but just make it highly recommended that a two-way comparison
#' with all group combinations is present.
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A string vector of groups which are a part of the experiment.
#' @export
#'
#' @examples
setMethod("getGroups", signature = "dgeExperiment", function(o) o@groups)

#' Get the treatments in the experiment
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector of the treatments. "None" is returned if there are none.
#' @export
#'
#' @examples
#' getTreatments(dgeExperiment)
setMethod("getTreatment", signature = "dgeExperiment", function(o) o@treatments)

#' Get the timepoints in the experiment
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector of the timepoints. "None" is returned if there are none.
#' @export
#'
#' @examples
#' getTreatments(dgeExperiment)
setMethod("getTimepoint", signature = "dgeExperiment", function(o) o@timepoints)

#' Get the brain regions of the experiment
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector of the brain regions.
#' @export
#'
#' @examples
#' getRegions(dgeExperiment)
setMethod("getRegion", signature = "dgeExperiment", function(o) o@regions)

#' Get the sexes in the experiment.
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector of the sexes of which we have comparisons. Can include "both"
#' @export
#'
#' @examples getSex(dgeExperiment)
setMethod("getSex", signature = "dgeExperiment", function(o) o@sexes)

#' Get the gene ensembl id's for this experiment
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector with ensembl id's
#' @export
#'
#' @examples getGeneEns(dgeExp)
setMethod("getGeneEns", signature = "dgeExperiment", function(o) o@comparisons[[1]]@deRes$ensembl)

#' Get the gene symbols for this experiment
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector with gene symbols
#' @export
#'
#' @examples getGeneSymbol(dgeExp)
setMethod("getGeneSymbol", signature = "dgeExperiment", function(o) o@comparisons[[1]]@deRes$symbol)

#' Get the gene DGE id's for this experiment
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector with DGE id's
#' @export
#'
#' @examples getGeneId(dgeExp)
setMethod("getGeneId", signature = "dgeExperiment", function(o) o@comparisons[[1]]@deRes$id)

#' Get the gene DGE id's for this experiment
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A character vector with DGE id's
#' @export
#'
#' @examples geneId(dgeExp)
setMethod("geneId", signature = "dgeExperiment", function(o) o@comparisons[[1]]$id)

#' Set the gene DGE id's for this experiment
#'
#' @param o An object of class \code{dgeExperiment}
#' @param value A data.table with ids and ensembl values.
#'
#' @return
#' @export
#'
#' @examples geneId(dgeExp) <- dgeIDs
setMethod("geneId<-", signature = "dgeExperiment", function(o, value) {
    o@comparisons <- lapply(o@comparisons, function(x) {
        geneId(x) <- value
        return(x)
    })
    return(o)
})

#' Flatten an experiment to a data table
#'
#' @param o An object of class \code{dgeExperiment}
#' @param geneId Information for only these id's is returned
#'
#' @return A data.table conatining all the information in the experiment.
#' @export
#'
#' @examples
setMethod("flattenDge", signature = "dgeExperiment", function(o, geneId = NULL) {
    rbindlist(lapply(o@comparisons, flattenDge, geneId = geneId))
})

#' Reverse comparisons in this experiment
#'
#' @param o An object of class \code{dgeExperiment}
#' @param comparisons A character vector of comparisonIDs which are to be reversed
#'
#' @return The \code{dgeExperiment} object with the IDs reversed
#' @export
#'
#' @examples
setMethod("reverseComparison", signature = "dgeExperiment", function(o, comparisons) {
    # check if any comparison is to be reversed here
    o@comparisons <- lapply(o@comparisons, function(x) {

        if (x@comparisonID %in% comparisons) {
            return(reverseComparison(x))
        } else {
            return(x)
        }
    })

    return(o)
})

#' Print the comparisons in this experiment object.
#'
#' @param o An object of class \code{dgeExperiment} for which to print the comparisons
#' @param name Logical, indicating if the name of the experiment should be printed. Default \code{FALSE}.
#' @param region Logical, indicating if the brain region of the comparison should be printed. Default \code{FALSE}.
#' @param treatment Logical, indicating if the treatment of the comparison should be printed. Default \code{FALSE}.
#' @param timepoint Logical, indicating if the timepoint of the comparison should be printed. Default \code{FALSE}.
#'
#' @return A vector of strings with the comparisons
#' @export
#'
#' @examples
setMethod("printComparison", signature = "dgeExperiment", function(o,
                                                                   name = FALSE,
                                                                   region = FALSE,
                                                                   treatment = FALSE,
                                                                   timepoint = FALSE)
    sapply(
        o@comparisons,
        printComparison,
        name = name,
        region = region,
        treatment = treatment,
        timepoint = timepoint
    ))

# Now define the generics which haven't been defined before
getAllComparisons <- function(o) stop("Method undefined for an object of this class.")
getComparison <- function(o, groups=NULL, model=NULL, treatment=NULL, timepoint=NULL, region=NULL, sex=NULL, ...) stop("Method undefined for an object of this class.")
getComparisonById <- function(o, id) stop("Method undefined for an object of this class.")
getTotalN <- function(o) stop("Method undefined for an object of this class.")
getIds <- function(o) stop("Method undefined for an object of this class.")
getOwner <- function(o) stop("Method undefined for an object of this class.")
getPublic <- function(o) stop("Method undefined for an object of this class.")
removeComparison <- function(o, id) stop("Method undefined for an object of this class.")

setGeneric("getAllComparisons", function(o) { standardGeneric("getAllComparisons") })
setGeneric("getComparison", function(o, groups=NULL, model=NULL, treatment=NULL, timepoint=NULL, region=NULL, sex=NULL, ...) { standardGeneric("getComparison") })
setGeneric("getComparisonById", function(o, id) { standardGeneric("getComparisonById") })
setGeneric("getTotalN", function(o) { standardGeneric("getTotalN") })
setGeneric("getIds", function(o) { standardGeneric("getIds") })
setGeneric("getOwner", function(o) { standardGeneric("getOwner") })
setGeneric("getPublic", function(o) { standardGeneric("getPublic") })
setGeneric("removeComparison", function(o, id) { standardGeneric("removeComparison") })

#' Returns the lab where this experiment was conducted. If it is private,
#' only lab members will have access to this data.
#'
#' @param dgeExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("getOwner", signature = "dgeExperiment", function(o) o@owner)

#' Title
#'
#' @param dgeExperiment
#'
#' @return
#' @export
#'
#' @examples
setMethod("getPublic", signature = "dgeExperiment", function(o) o@public)

#' Get all the internal comparison IDs in this experiment
#'
#' @param o An object of class \code{dgeExperiment} from which we want the IDs
#'
#' @return A character vector of IDs
#' @export
#'
#' @examples
setMethod("getIds", signature = "dgeExperiment", function(o) sapply(o@comparisons, function(x) x@comparisonID))

#' Gets the comparisons within this experiment by id.
#'
#' @param o An object of class \code{dgeExperiment} which we want to subset
#' @param id A character vector with the ids we need
#'
#' @return An \code{dgeExperiment} object with the needed ids. If none found, NULL is returned.
#' @export
#'
#' @examples
setMethod("getComparisonById", signature = "dgeExperiment", function(o, id) {
    sc <- lapply(o@comparisons, function(x) {
        if (x@comparisonID %in% id) {
            return(x)
        } else if (reverseComparison(x)@comparisonID %in% id) {
            return(reverseComparison(x))
        } else {
            return(NULL)
        }
    })

    sc <- sc[lengths(sc) != 0]

    o@comparisons <- sc

    if (length(sc) == 0)
        return(NULL)
    else
        return(o)
})

#' Remove comparisons from an experiment
#'
#' @param o The expermient from which to remove the comparisons
#' @param id A character vector of the ids of the comparison to remove.
#'
#' @return The experiment with the comparison removed
#' @export
#'
#' @examples
setMethod("removeComparison", signature = "dgeExperiment", function(o, id) {
    sc <- lapply(o@comparisons, function(x) if (x@comparisonID %in% id) return(NULL) else return(x))

    sc <- sc[lengths(sc) != 0]

    o@comparisons <- sc

    if (length(sc) == 0)
        return(NULL)
    else
        return(o)
})

#' Get all comparisons
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return A list of objects of class `dgeComparison`.
#' @export
#'
#' @examples
setMethod("getAllComparisons", signature = "dgeExperiment", function(o) o@comparisons)

#' Get total N
#'
#' @param o An object of class \code{dgeExperiment}
#'
#' @return The total number of animals in the experiment.
#' @export
#'
#' @examples
#' getTotalN(dgeExperiment)
setMethod("getTotalN", signature = "dgeExperiment", function(o) o@totalN)

#' Get one or many comparisons from an experiment object
#'
#' The comparisons matching the given parameters are returned from the experiment. If groups are specified, they need to be a
#' character vector of length two. If need be, the comparison will be reversed. One of groups, treatment or region
#' needs to be specified. If both treatment and region are specified, the comparison will be matched using both parameters.
#' In general, this performs an AND, not an OR. Comparisons matching ALL parameters will be returned.
#'
#' @param o An object of class \code{dgeExperiment} from which the comparison is to be extracted
#' @param groups The groups for which the comparison is needed. If specified, the other two parameters are ignored
#' @param model Select only specified models.
#' @param treatment Optional parameter specifying the treatments needed for the comparison. Can be used in conjunction with others
#' @param timepoint Optional parameter specifying the timepoints needed for the comparison. Can be used in conjunction with others
#' @param region The brain region of the desired comparison. Can be specified in conjunction with treatment
#' @param sex The sex of the desired comparison
#'
#' @return If no comparison is matched, `NULL` is returned, else an object of class `dgeExperiment` with the specified comparisons.
#' @export
#'
#' @examples
#' #getComparison(fslExperiment, groups = c("FSL", "FRL"))
#' #getComparison(fslExperiment, region = "HPC")
setMethod("getComparison", signature = "dgeExperiment",
          function(o, groups=NULL, model=NULL, treatment=NULL, timepoint=NULL, region=NULL, sex=NULL) {

              # we either need both groups or treatment or region or sex)nd
              stopifnot((!is.null(groups) && length(groups) == 2) ||
                            !is.null(treatment) ||
                            !is.null(timepoint) ||
                            !is.null(model) ||
                            (!is.null(region) && region %in% .validRegions()) ||
                            (!is.null(sex) && tolower(sex) %in% .validSex())
              )

              sc <- o@comparisons

              # let's start with sex
              if (!is.null(sex)) {
                  t <- sapply(sc, function(x) any(getSex(x) %in% tolower(sex)))
                  if (!any(t))
                      return(NULL)

                  sc <- sc[t]
              }

              # now moving on to region
              if (!is.null(region)) {
                  t <- sapply(sc, function(x) any(getRegion(x) %in% region))
                  if (!any(t))
                      return(NULL)

                  sc <- sc[t]
              }

              # now treatment
              if (!is.null(treatment)) {
                  t <- sapply(sc, function(x) any(getTreatment(x) %in% treatment))
                  if (!any(t))
                      return(NULL)

                  sc <- sc[t]
              }

              # now timepoint
              if (!is.null(timepoint)) {
                  t <- sapply(sc, function(x) any(getTimepoint(x) %in% timepoint))
                  if (!any(t))
                      return(NULL)

                  sc <- sc[t]
              }

              # now model
              if (!is.null(model)) {
                  t <- sapply(sc, function(x) any(getModel(x) %in% model))
                  if (!any(t))
                      return(NULL)

                  sc <- sc[t]
              }


              # now look for the groups
              if (!is.null(groups)) {
                  g <- lapply(sc, getGroups)
                  t <- sapply(g, function(x) all(x %in% groups))
                  if (!any(t))
                      return(NULL)

                  sc <- sc[t]
                  # check if we need to reverse it
                  sc <- lapply(sc, function(x) {
                      if (all(getGroups(x) == groups))
                          x
                      else
                          reverseComparison(x)
                  })
              }

              if (length(sc) == 0)
                  return(NULL)
              else {
                  o@comparisons <- sc
                  return(o)
              }
          })

#' Subset operator for DGE Experiment
#'
#' Rather simplistic. Don't specify `j`... it won't do anything. Just `i` as a logical vector will work. It would be
#' best to use `getGeneEns` to get the Ensembl id's and subset using those.
#'
#' @param x An object of class \code{dgeExperiment}
#'
#' @return A subsetted object of class `dgeExperiment`.
#' @export
#'
#' @examples
setMethod("[", signature = "dgeExperiment", function(x, i) {

    x@comparisons <- lapply(x@comparisons, function(a) a[i, ] )

    return(x)
})

# make the constructor
#' Constructor for the `dgeExperiment` class.
#'
#' Takes a list of `dgeComparisons` belonging to the same experiment and creates an object of `dgeExperiment` from them.
#'
#' @param comparisons A `list` of comparisons belonging to the same experiment
#' @param description An optional description of the experiment.
#'
#' @return An object of class`dgeExperiment`
#' @export
#'
#' @examples
dgeExperiment <- function(comparisons, owner = NULL, public = TRUE, description = "") {

    # all comparisons must be of the class dgeComparison
    stopifnot(all(sapply(comparisons, class) == "dgeComparison"))

    n <- unique(sapply(comparisons, getName))
    stopifnot(length(n) == 1)

    s <- unique(sapply(comparisons, getSpecies))
    stopifnot(length(s) == 1)

    p <- unique(sapply(comparisons, getPlatform))
    stopifnot(length(p) == 1)

    # make sure that all the ensembl ids are the same. if not, add the missing ones and fill with NA values
    # also, make sure that all of they are ordered by ensembl ids. Might make sense to switch this to data.table
    # and use ensemble as the key at some point
    ens <- lapply(comparisons, function(x) x$ensembl)

    # get the unique elements only
    ens <- unique(unlist(ens, use.names = FALSE))

    comparisons <- lapply(comparisons, function(comp) {
        toAdd <- ens[which(!(ens %in% comp$ensembl))]

        if (length(toAdd) > 0) {
            # add the ens which are not present in the comparison
            toAdd <- data.frame(ensembl = toAdd, stringsAsFactors = FALSE)
            comp@deRes[(nrow(comp) + 1):(nrow(comp) + nrow(toAdd)), names(toAdd)] <- toAdd
        }

        # now sort it
        comp@deRes <- comp@deRes[order(ens), ]

        return(comp)
    })

    new("dgeExperiment",
        comparisons = comparisons,
        name = n,
        model = unique(sapply(comparisons, getModel)),
        species = s,
        platform = p,
        treatments = unique(sapply(comparisons, getTreatment)),
        timepoints = unique(sapply(comparisons, getTimepoint)),
        regions = unique(sapply(comparisons, getRegion)),
        sexes = unique(sapply(comparisons, getSex)),
        groups = unique(c(sapply(comparisons, getGroups))),
        totalN = sum(sapply(comparisons, getNs)),
        owner = owner,
        public = public,
        description = description
    )

}


