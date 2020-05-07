#' @include dgeExperiment.R
#' @include internalPlottingFunctions.R
#' @import ggplot2
#' @importFrom ggpubr as_ggplot
#' @importFrom gridExtra arrangeGrob
#' @importFrom heatmaply heatmaply
#' @importFrom grid textGrob
#' @importFrom grid gpar
#' @importFrom metafor rma
NULL

#' Title
#'
#' @param hData Needs to be a \code{data.table} object.
#'              The data can contain information for a single gene or multiple genes.
#'              However many genes there are, those many plots will be returned.
#'              The function expects the columns to be the same as those returned by the
#'              \code{flattenDge} method.
#' @param metaStatistic logical. A random-effect meta-analysis statistic is included if \code{TRUE}.
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
#' @param fontSize The font size when printing. Defaults to 3.
#' @param ids The data.table of mappings between ids and symbols for genes
#'
#' @return A list of ggplot2 objects which can be plotted.
#' @export
#'
#' @examples
smoothForest <- function(hData,
                         metaStatistic = FALSE,
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
                         fontSize = 3,
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
    if (!requireNamespace("parallel", quietly = TRUE)) {
        outPlots <- mapply(function(x, nx) {
            .forestPlot(x, nx,
                        metaStatistic,
                        includeModel,
                        includeName,
                        includeRegion,
                        includeEstimate,
                        includeTreatment,
                        includeTimepoint,
                        includeSex,
                        includePval,
                        comparisonNames,
                        orderBy,
                        estimate,
                        splitBy,
                        fontSize,
                        ids)
        }, hData, names(hData), SIMPLIFY = FALSE)
    } else {
        outPlots <- parallel::mcmapply(function(x, nx) {
            .forestPlot(x, nx,
                        metaStatistic,
                        includeModel,
                        includeName,
                        includeRegion,
                        includeEstimate,
                        includeTreatment,
                        includeTimepoint,
                        includeSex,
                        includePval,
                        comparisonNames,
                        orderBy,
                        estimate,
                        splitBy,
                        fontSize,
                        ids)
        }, hData, names(hData), SIMPLIFY = FALSE,
        mc.silent = TRUE, mc.cores = getOption("mc.cores", 4L))
    }

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
