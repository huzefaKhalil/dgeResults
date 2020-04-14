#' @include dgeExperiment.R
#' @include internalPlottingFunctions.R
#' @import ggplot2
#' @importFrom ggpubr as_ggplot
#' @importFrom gridExtra arrangeGrob
#' @importFrom heatmaply heatmaply
#' @importFrom grid textGrob
#' @importFrom grid gpar
NULL

#' Title
#'
#' @param hData Needs to be a \code{data.table} object.
#'              The data can contain information for a single gene or multiple genes.
#'              However many genes there are, those many plots will be returned.
#'              The function expects the columns to be the same as those returned by the
#'              \code{flattenDge} method.
#' @param includeModel logical. Includes the model if \code{TRUE}.
#' @param includeName logical. Includes the experiment name if \code{TRUE}.
#' @param includeRegion logical. Includes the brain region if \code{TRUE}.
#' @param includeEstimate logical. Includes the estimate if \code{TRUE}.
#' @param includeTreatment logical. Includes the treatment if \code{TRUE}.
#' @param includeTimepoint logical. Includes the ttimepoint if \code{TRUE}.
#' @param includeSex logical. Includes the sex if \code{TRUE}.
#' @param includePval logical. Includes the pvalue if \code{TRUE}.
#' @param comparisonNames A character vector. If specified, the graph will show these as the comparison names.
#' @param orderBy One of \code{asis}, \code{estimate}, \code{experiment}, \code{pval}.
#' @param estimate Which estimate should be plotted. Options are \code{hedgesG} (the default), \code{cohensD}, \code{logFC}.
#' @param splitBy Generate plots either by "gene" or by "comparison".
#' @param ids The data.table of mappings between ids and symbols for genes
#'
#' @return A list of ggplot2 objects which can be plotted.
#' @export
#'
#' @examples
smoothForest <- function(hData,
                         includeModel = FALSE,
                         includeName = FALSE,
                         includeRegion = FALSE,
                         includeEstimate = FALSE,
                         includeTreatment = FALSE,
                         includeTimepoint = FALSE,
                         includeSex = FALSE,
                         includePval = FALSE,
                         comparisonNames = NULL,
                         orderBy = "estimate",
                         estimate = "hedgesG",
                         splitBy = "gene",
                         ids) {

    # function to draw the smooth forest plot.
    # The logic is as follows:
    # 1) flatten any data or experiment object
    # 2) Split it by individual genes; do the following for each gene
    # 3) Sort the table
    # 4) Create a text data table
    # 5) Create a numeric data table
    # 6) Plot both
    # 7) Combine and return the plots

    if (!is.data.table(hData))
        stop("Input to the 'smoothForest' method must be a 'data.table'.")

    if (!(splitBy %in% c("gene", "comparison")))
        splitBy <- "gene"

    hData$compound.symbol <- ids[hData]$compound.symbol

    # step 2: split by unique genes or comparison
    if (splitBy == "gene")
        hData <- split(hData, by = "compound.symbol")
    else
        hData <- split(hData, by = "comparisonID")

    # next, loop over each genes
    outPlots <- Map(function(x, nx) {

        x <- .order(x, orderBy, estimate)

        # make the data table
        metaData <- .getMetaData(x,
                                 includeModel = includeModel,
                                 includeName = includeName,
                                 includeRegion = includeRegion,
                                 includeEstimate = includeEstimate,
                                 includeTreatment = includeTreatment,
                                 includeTimepoint = includeTimepoint,
                                 includeSex = includeSex,
                                 includePval = includePval,
                                 estimate = estimate,
                                 comparisonNames = comparisonNames,
                                 splitBy = splitBy)

        metaPlot <- .plotMetaData(metaData)
        forestPlot <- .plotForest(x, estimate)

        if (splitBy == "gene") {
            title <- ids[id %in% x$id[1]]$compound.symbol
            fSize <- 15
            layout <- matrix(c(1, 1, 1, 1, 1, 2, 2, 2), nrow=1)
        } else {
            title <- paste(x$comparison[1], x$name[1], x$region[1], sep = ", ")
            fSize <- 10
            layout <- matrix(c(1, 2), nrow=1)
        }

        # layout <- matrix(c(rep(1, times = 1 + ncol(metaData)), 2, 2, 2, 2, 2), nrow = 1)
        # as_ggplot(gridExtra::arrangeGrob(metaPlot, forestPlot, layout_matrix = layout))
        ggpubr::as_ggplot(gridExtra::arrangeGrob(metaPlot, forestPlot, layout_matrix = layout,
                                                 top = grid::textGrob(label = title,
                                                                      gp = grid::gpar(fontsize=fSize, fontface=2)),
                                                 padding = unit(1, "line")))
    }, hData, names(hData))

    return(outPlots)
}


#' Plots an interactive heatmap suitable for webpages
#'
#' @param hData The \code{data.table} from which to plot the heatmap.
#'              All the data points in this object will be plotted, so any subsetting must happen before.
#' @param estimate The estimate which should be plotted. Must be one of \code{hedgesG} (default),
#'                 \code{cohensD} or \code{logFC}.
#' @param ids The ids which will be printed. Can be a character vector or data.table
#'            with the ids as well as the compound symbol
#' @param Rowv Should the rows be reordered? Defualt is \code{TRUE}.
#' @param Colv Should the columns be reordered? Defualt is \code{TRUE}.
#'
#' @return An object of class heatmapr
#' @export
#'
#' @examples
heatmap.3 <- function(hData, estimate = "hedgesG", ids = NULL, Rowv = TRUE, Colv = TRUE) {

    if (!is.data.table(hData))
        stop("Input to the 'heatmap.3' method must be an 'data.table' object.")

    if (!(estimate %in% c("hedgesG", "cohensD", "logFC")))
        estimate <- "hedgesG"

    # If there are too many NA values, the function will throw an error.
    # It might make senese to remove comparisons / gene combinations with a lot of NA vals

    .heatmap.3(hData, estimate, ids)
}
