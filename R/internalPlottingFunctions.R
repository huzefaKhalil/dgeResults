
# plot the heatmap
.heatmap.3 <- function(hd, estimate, ids) {

    # get only unique rows... This shouldn't happen often, but
    hd <- hd[!duplicated(hd[, -c("id")])]
    # convert to wide
    pMat <- dcast(hd, id ~ comparisonID, value.var = estimate)

    # remove rows with a lot of NA values... Arbitrarily choosing 70% as the cut-off
    # so, the gene must be present in 70% of the comparisons to be plotted
    # nNA <- rowSums(is.na(pMat[, -1]))
    # nNA <- (nNA / (ncol(pMat) - 1)) < 0.3
    #
    # pMat <- pMat[nNA, ]

    genes <- pMat$id
    pMat <- as.matrix(pMat[, -1])

    if (is.null(ids))
        row.names(pMat) <- unique(hd[id %in% genes]$symbol)
    else if (is.character(ids))
        row.names(pMat) <- ids
    else if (is.data.frame(ids))
        row.names(pMat) <- ids[id %in% genes]$compound.symbol

    # make the hover text
    cNames <- unique(hd[, c("comparisonID", "name", "model", "region", "comparison", "treatment", "sex", "timepoint")])
    cNames <- cNames[match(colnames(pMat), cNames$comparisonID)]

    hoverTextCols <- paste0("Experiment: ", cNames$name, "\n",
                            "Model: ", cNames$model, "\n",
                            "Brain region: ", cNames$region, "\n",
                            "Comparison: ", cNames$comparison, "\n",
                            "Sex: ", cNames$sex, "\n",
                            "Treatment: ", cNames$treatment, "\n",
                            "Timepoint: ", cNames$timepoint, "\n")
    hoverTextRows <- paste0("Gene Name: ", row.names(pMat), "\n")

    hoverTextMat <- matrix(paste0(rep(hoverTextRows, each=length(hoverTextCols)), hoverTextCols),
                           nrow = length(hoverTextRows),
                           ncol = length(hoverTextCols), byrow = TRUE)

    # add the gene estimate
    hoverTextMat <- matrix(paste0(hoverTextMat, "Estimate: ", round(pMat, digits = 4)), ncol = ncol(pMat), nrow = nrow(pMat))

    #colnames(pMat) <- cNames$name

    tryCatch(
        heatmaply::heatmapr(pMat,
                            custom_hovertext = hoverTextMat),
        error = function(e) {
            heatmaply::heatmapr(pMat,
                                custom_hovertext = hoverTextMat,
                                Rowv = FALSE,
                                Colv = FALSE,
                                show_dendrogram = FALSE)
        }
    )
}

# putting this here so that we can make use of parallel and clusters to deal with plotting
.forestPlot <- function(x, nx,
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
                        ids) {
    x <- .order(x, orderBy, estimate)

    ms <- NULL
    # do the meta analysis
    if (metaStatistic) {
        if (estimate == "cohensD") {
            estM <- "cohensD"
            varM <- "varD"
        } else {
            estM <- "hedgesG"
            varM <- "varG"
        }

        tryCatch({
            meta <- metafor::rma(x[[estM]], x[[varM]])
            ms <- data.table("statistic" = meta$b[1, ],
                             "se" = meta$se,
                             "pval" = meta$pval)
        }, finally = {})
    }

    # make the data table
    metaData <- .getMetaData(x,
                             ms,
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

    metaPlot <- .plotMetaData(metaData, fontSize)
    forestPlot <- .plotForest(x, ms, estimate)

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
}

# returns a forest plot only
.plotForest <- function(fd, ms = NULL, estimate) {

    # basic idea by dspaks @ https://gist.github.com/dsparks/818997..
    # That one was way too convoluted... This simplifies things quite a bit
    # set it up for a smooth coef plot

    # first get the estimate we need, along with the se
    fd <- fd[!is.na(fd$logFC), ]    # remove NA rows
    fd$estimate <- switch(estimate,
                          hedgesG = fd$hedgesG,
                          cohensD = fd$cohensD,
                          logFC = fd$logFC,
                          fd$hedgesG)
    fd$stdErr <- switch(estimate,
                        hedgesG = sqrt(fd$varG),
                        cohensD = sqrt(fd$varD),
                        logFC = fd$se,
                        sqrt(fd$varG))

    alphas <- seq(1, 99, 2) / 100
    multiplier <- qnorm(1 - alphas/2)
    emphasis <- 1 - seq(0, 1, length = length(multiplier) +1)[-1]
    transparency <- 1 / (length(multiplier) / 4)

    # get the data we need to plot.. What should we plot? Hedger's g
    fd <- fd[, c("estimate", "stdErr")]

    if (!is.null(ms)) {
        fd <- rbind(data.table("estimate" = c(ms$statistic, NA),
                               "stdErr" = c(ms$se, NA)),
                    fd)
    }

    fd$id <- 1:nrow(fd)
    nr <- nrow(fd)

    # now, we have to expand the data to include the multiplier
    fd <- fd[rep(seq(1, nrow(fd)), length(multiplier))]

    # and now add this to the data
    multiplier <- rep(multiplier, each = nr)
    emphasis <- rep(emphasis, each = nr)

    fd <- cbind(fd, multiplier, emphasis)
    fd <- na.omit(fd)

    xLim <- range(fd$estimate - fd$multiplier*fd$stdErr, fd$estimate + fd$multiplier*fd$stdErr)
    yLim <- c(0, nr+1.5)

    ggplot(data = fd, aes(x = estimate, y = id, xmin = estimate - multiplier * stdErr, xmax = estimate + multiplier * stdErr)) +
        geom_vline(xintercept = 0, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12)) +
        geom_errorbarh(data = fd, aes(size = 1/emphasis), alpha = I(transparency),
                       colour = I(gray(0)), lwd = I(7/12), height = 0) +
        scale_y_continuous(name = "", breaks = c(1:nr), labels = NULL) +
        geom_point(aes(x = estimate, y = id), colour = I(gray(0))) + theme_bw() +
        coord_cartesian(xlim = xLim, ylim = yLim, expand = F) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.border = element_blank(),
              axis.text.x = element_text(colour="black"),
              axis.text.y = element_blank(),
              axis.ticks.x = element_line(colour="black"),
              axis.ticks.y = element_line(colour="black"),
              axis.line.x = element_line(colour="black"),
              axis.line.y = element_line(colour = "black"),
              plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 0, unit = "pt")) +
        labs(x = "", y = "")
}

.plotMetaData <- function(md, fontSize) {

    # this function plots the text summary and returns it as a ggplot2 plot object..
    # kudos to metaviz.. I'm gonna borrow much of their code here.
    # TODO: Don't forget to cite them

    # The coulmns should all be characters!

    # get the size of each column
    tmd <- rbind(matrix(colnames(md), nrow = 1), md, use.names = FALSE)
    colArea <- cumsum(c(1, apply(tmd, 2, function(y) max(round(max(nchar(y))/100, 2), 0.03))))
    names(colArea) <- NULL
    # Add a little to columns after the first
    colArea <- c(colArea[1:2], colArea[-c(1,2)] + cumsum(rep(0.02, length(colArea) - 2)))

    xVals <- colArea[1:ncol(md)]
    xLim <- range(colArea)
    yLim <- c(0, nrow(md) + 1.5)

    # if FDR is present, color red
    if ("FDR" %in% colnames(md)) {
        clr <- as.numeric(md$FDR) < 0.05
    } else {
        clr <- rep(FALSE, nrow(md))
    }

    # create the data to plot
    df <- data.frame(y = rep(1:nrow(md), ncol(md)),  # Get the number of rows as facotrs to plot separately
                     x = rep(xVals, each=nrow(md)),
                     value = unlist(lapply(md, as.vector), use.names = FALSE),  # the text to print
                     stringsAsFactors = FALSE)
    df$clr <- rep(clr, ncol(md))

    dfTitle <- data.frame(y = rep(nrow(md) + 1, ncol(md)),
                          x = xVals,
                          value = names(md))

    text_size <- fontSize
    ggplot(df, aes(x = x, y = y)) +
        geom_text(aes(label = value, color = clr), size = text_size, hjust = 0, vjust = 0.5)  +
        scale_color_manual(values = c("black", "red")) +
        geom_text(data = dfTitle, aes(x = x, y = y, label = value), size = text_size, hjust = 0, vjust = -0.5) +
        coord_cartesian(xlim = xLim, ylim = yLim, expand = F) +
        geom_hline(yintercept = nrow(md) + 0.5) +
        theme_bw() +
        theme(text = element_text(size = 1/0.352777778*text_size),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.border = element_blank(),
              axis.text.x = element_text(colour="white"),
              axis.text.y = element_blank(),
              axis.ticks.x = element_line(colour="white"),
              axis.ticks.y = element_blank(),
              axis.line.x = element_line(colour="white"),
              axis.line.y = element_blank(),
              plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 5.5, unit = "pt")) +
        labs(x = "", y = "")
}


.getMetaData <- function(x, ms = NULL, includeModel=F, includeName=F,
                         includeRegion=F, includeEstimate=F,
                         includeTreatment=F, includeTimepoint=F,
                         includeSex=F,
                         includePval=F, estimate="hedgesG",
                         splitBy="gene", comparisonNames = NULL) {

    # keep only the rows for which we have a logFC
    x <- x[!is.na(x$logFC), ]

    if (splitBy == "gene") {

        outDF <- data.table("Comparison" = x$comparison)

        if (!is.null(comparisonNames)) {
            # check if it has the same length as the data
            if (length(comparisonNames) != nrow(x))
                stop("In internal plotting function, the comparison names are not equal to the number of rows of data.table")

            outDF$Comparison <- comparisonNames
        }

        if (isTRUE(includeModel))
            outDF <- data.table(outDF, "Model" = x$model)

        if (isTRUE(includeName))
            outDF <- data.table(outDF, "Name" = x$name)

        if (isTRUE(includeRegion))
            outDF <- data.table(outDF, "Region" = x$region)

        if (isTRUE(includeTreatment))
            outDF <- data.table(outDF, "Treatment" = x$treatment)

        if (isTRUE(includeTimepoint))
            outDF <- data.table(outDF, "Timepoint" = x$timepoint)

        if (isTRUE(includeSex))
            outDF <- data.table(outDF, "Sex" = x$sex)

        if (isTRUE(includeEstimate))
            outDF <- data.table(outDF,
                                "Estimate" = format(switch(estimate,
                                                           hedgesG = x$hedgesG,
                                                           cohensD = x$cohensD,
                                                           logFC = x$logFC,
                                                           x$hedgesG), digits = 2))

        if (isTRUE(includePval))
            outDF <- data.table(outDF, "P-value" = format(x$pvalue, digits = 3))

        if (!is.null(ms)) {
            pl <- outDF[(nrow(outDF) + 1):(nrow(outDF) + 2),]
            pl[1:2, ] <- ""
            pl[1, "Comparison"] <- "Summary"
            pl[1, "Estimate"] <- format(ms$statistic, digits = 2)
            pl[1, "P-value"] <- format(ms$pval, digits = 3)

            outDF <- rbind(pl, outDF)
        }

    } else {
        # just show the gene name, logFC, p-value and fdr
        outDF <- data.table("Gene Symbol" = x$compound.symbol)
        outDF <- data.table(outDF, "LogFC" = format(x$logFC, digits = 2))

        # make them zero if very low
        x$pvalue <- ifelse(x$pvalue <= 0.0001, 0, x$pvalue)
        x$fdr <- ifelse(x$fdr <= 0.0001, 0, x$fdr)

        outDF <- data.table(outDF, "P-value" = format(x$pvalue, digits = 3,
                                                      scientific = FALSE))
        outDF <- data.table(outDF, "FDR" = format(x$fdr, digits = 3,
                                                  scientific = FALSE))
    }
    return(outDF)
}

# change the order of the data.table
.order <- function(x, orderBy, estimate = "hedgesG") {

    if (!(orderBy %in% c("asis", "estimate", "experiment", "pval")))
        stop("Error in the internal function '.order' in file 'internalPlottingFunctions.R'. The 'orderBy' value is not valid.")

    if (orderBy == "asis") {
        return(x[nrow(x):1])
    }

    if (orderBy == "estimate") {
        return(x[order(x[[estimate]])])
    }

    if (orderBy == "experiment") {
        return(x[order(-name, x[[estimate]])])
    }

    if (orderBy == "pval") {
        return(x[order(-pvalue)])
    }
}

