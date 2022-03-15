#' Plot ligand-receptor gene expression in reduced dimensions
#'
#' Plot ligand and receptor gene expression on low-dimensional projections
#' stored in a \code{SingleCellExperiment} object.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param dimred A string or integer scalar indicating the reduced dimension 
#' result in \code{reducedDims(object)} to plot.
#' @param lr_pair A character vector of length 2 containing the ligand and 
#' receptor gene symbol.
#' @param lr_desc A character vector of length 2 containing short description 
#' to change legend title. (default: \code{c("Ligand","Receptor")})
#' @param lr_color A character vector of length 2 containing colour aesthetics. 
#' (default: \code{c("blue","red")})
#' @param lr_sep A character string to define how the 2 genes terms are 
#' separated. (default: "-")
#' @param oneplot Logical scalar indicating whether to overlay expressions 
#' in a single plot or generate 2 side-by-side plots. (default: TRUE)
#' @param by_exprs_values A string or integer scalar specifying which assay 
#' to obtain expression values from, for use in point aesthetics. 
#' (default: "logcounts")
#' @param point_size A numeric scalar specifying the size of the points. 
#' (default: 2)
#' @param point_alpha A numeric scalar (between 0 and 1) specifying the 
#' transparency. (default: 0.4)
#' @param point_shape An integer scalar (between 0 and 25) specifying the 
#' shape aesthetics. (default: 16)
#' @param text_by A string specifying the column metadata field with which to 
#' add text labels on the plot.
#' @param text_size A numeric scalar specifying the size of added text. 
#' (default: 8)
#' @param text_color A string specifying the colour of the added text. 
#' (default: "black")
#' @param theme_size A numeric scalar specifying the base font size. 
#' (default: 14)
#'
#' @return A ggplot object
#'
#' @details
#' The function is based on \code{plotReducedDim()} from the 
#' \code{\link[scater]{scater}} package. It uses the \code{new_scale_colour()}
#' from the \code{\link[scater]{ggnewscale}} package to add an additional \emph{layer} 
#' where a second \code{geom} will use another colour scale to show the gene 
#' expression intensity.
#'
#' @author I-Hsuan Lin
#'
#' @name plotReducedDimLR
#'
#' @seealso \code{\link[scater::plotReducedDimLR()]{scater::plotReducedDimLR}}
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom scater retrieveCellInfo
#' @importFrom stats quantile
#' @importFrom scales squish
#' @importFrom ggnewscale new_scale_colour
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom ggplot2 guide_colorbar
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 calc_element
#' @importFrom cowplot theme_cowplot
#' @importFrom cowplot draw_label
#' @importFrom cowplot ggdraw
#' @importFrom cowplot plot_grid
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' sce <- runPCA(sce, ncomponents = 5)
#' plotReducedDimLR(sce, "PCA", c("Gene_0001","Gene_1111"))
plotReducedDimLR <- function(object, dimred, lr_pair, lr_desc = c("Ligand","Receptor"),
                             lr_color = c("blue","red"), lr_sep = "-", oneplot = TRUE,
                             by_exprs_values = "logcounts", point_size = 2,
                             point_alpha = 0.4, point_shape = 16, text_by = NULL,
                             text_size = 8, text_color = "black", theme_size = 14) {
    red_dim <- as.data.frame(reducedDim(object, dimred))

    if(length(lr_pair) == 1L) {
        lr_pair <- unlist(strsplit(lr_pair, "[-]+"))
    }

    if(length(lr_pair) != 2L) {
        stop("Wrong LR pair format.")
    }

    if(length(lr_color) != 2L) {
        stop("Wrong LR pair color format.")
    }

    if(length(lr_desc) != 2L) {
        stop("Wrong LR pair description format.")
    }

    if(length(point_shape) != 2L) {
        if(length(point_shape) == 1L) {
            point_shape <- c(point_shape, point_shape)
        } else {
            stop("Wrong LR pair shape format.")
        }
    }

    gene1 <- retrieveCellInfo(object, lr_pair[1], exprs_values = by_exprs_values)
    gene2 <- retrieveCellInfo(object, lr_pair[2], exprs_values = by_exprs_values)
    label1 <- paste0(lr_desc[1], "\n", gene1$name)
    label2 <- paste0(lr_desc[2], "\n", gene2$name)
    mid1 <- mean(gene1$value)
    mid2 <- mean(gene2$value)
    max_val <- max(quantile(gene1$value, 0.95), quantile(gene2$value, 0.95))
    max_val <- ifelse(max_val > 0, max_val, min(c(max(gene1$value), max(gene2$value))))
    max_val <- ifelse(max_val > 0, max_val, c(max(gene1$value), max(gene2$value)))
    limits <- c(0, max_val)
    title <- paste(lr_pair, collapse = lr_sep)

    dat <- cbind(red_dim, data.frame(L = gene1$value, R = gene2$value))

    # Set up plot
    p <- ggplot(mapping = aes_string(x = "V1", y = "V2"))

    if(oneplot) {
        p <- p + geom_point(data = subset(dat, L > mid1), aes_string(color = "L"),
                            shape = point_shape[1], size = point_size, alpha = point_alpha) +
        scale_colour_gradientn(label1, colours = c("white", lr_color[1]), limits = limits,
                               oob = squish, guide = guide_colorbar(order = 1)) +
        new_scale_colour() +
        geom_point(data = subset(dat, R > mid2), aes_string(color = "R"),
                   shape = point_shape[2], size = point_size, alpha = point_alpha) +
        scale_colour_gradientn(label2, colours = c("white", lr_color[2]), limits = limits,
                               oob = squish, guide = guide_colorbar(order = 2)) +
        theme_cowplot(theme_size) + labs(title = title, x = paste(dimred, "1"), y = paste(dimred, "2"))

        if (!is.null(text_by)) {
            p <- p + add_label(object, dimred, text_by = text_by, text_size = text_size,
                               text_color = text_color)
        }

        p
    } else {
        p1 <- p + geom_point(data = subset(dat, L > mid1), aes_string(color = "L"),
                             shape = point_shape[1], size = point_size, alpha = point_alpha) +
        scale_colour_gradientn(label1, colours = c("white", lr_color[1]), limits = limits,
                               oob = squish) +
        theme_cowplot(theme_size) + labs(x = paste(dimred, "1"), y = paste(dimred, "2"))

        p2 <- p + geom_point(data = subset(dat, R > mid2), aes_string(color = "R"),
                             shape = point_shape[2], size = point_size, alpha = point_alpha) +
        scale_colour_gradientn(label2, colours = c("white", lr_color[2]), limits = limits,
                               oob = squish) +
        theme_cowplot(theme_size) + labs(x = paste(dimred, "1"), y = paste(dimred, "2"))

        if (!is.null(text_by)) {
            p1 <- p1 + add_label(object, dimred, text_by = text_by, text_size = text_size,
                                 text_color = text_color)
            p2 <- p2 + add_label(object, dimred, text_by = text_by, text_size = text_size,
                                 text_color = text_color)
        }

        title_theme <- calc_element("plot.title", theme_cowplot())

        title <- ggdraw() + draw_label(title, x = 0.05,
                                       hjust = title_theme$hjust,
                                       vjust= title_theme$vjust,
                                       fontface = title_theme$face,
                                       color = title_theme$colour,
                                       size = title_theme$size * 1.3,
                                       lineheight = title_theme$lineheight,
                                       angle = title_theme$angle)

        plot_grid(title, plot_grid(p1, p2), ncol = 1, rel_heights = c(0.1, 1))
    }
}

#' Add labels to reduced dimension plots
#'
#' Add labels to reduced dimension plots.
#'
#' @param object A \code{SingleCellExperiment} object.
#' @param dimred A string or integer scalar indicating the reduced dimension 
#' result in \code{reducedDims(object)} to plot.
#' @param text_by A string specifying the column metadata field with which to
#' add text labels on the plot.
#' @param text_size A numeric scalar specifying the size of added text.
#' (default: 8)
#' @param text_color A string specifying the colour of the added text.
#' (default: "black")
#'
#' @return A geom (geometric object)
#'
#' @details
#' The function is designed to use with \code{\link[scater::plotReducedDim]{scater::plotReducedDim()}}. 
#' It overrides the repel in the \code{plotReducedDim()} with force = 0, hence placing 
#' the labels centrally.
#'
#' The repel away from center bug in \code{plotReducedDim()} will be fixed in
#' \code{\link[scater]{scater}} package v1.23.5.
#'
#' @author I-Hsuan Lin
#'
#' @name add_label
#'
#' @seealso \code{\link[scater::plotReducedDimLR]{scater::plotReducedDimLR()}}
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom scater retrieveCellInfo
#' @importFrom stats median
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 aes_string
#' @examples
#' library(scater)
#' sce <- mockSCE()
#' sce <- logNormCounts(sce)
#' sce <- runPCA(sce, ncomponents = 5)
#' plotReducedDim(sce, "PCA", colour_by = "Treatment") + 
#'      add_label(sce, "PCA", text_by = "Treatment")
add_label <- function(object, dimred, text_by = "label", text_size = 8,
                      text_color = "black") {
    text_out <- retrieveCellInfo(object, text_by, search = "colData")
    text_out$val <- as.factor(text_out$val)
    red_dim <- as.matrix(reducedDim(object, dimred))
    df_to_plot <- data.frame(red_dim[, seq_len(2), drop = FALSE])
    colnames(df_to_plot)[seq_len(2)] <- c("X", "Y")
    by_text_x <- vapply(split(df_to_plot$X, text_out$val), median, FUN.VALUE = 0)
    by_text_y <- vapply(split(df_to_plot$Y, text_out$val), median, FUN.VALUE = 0)

    geom_text_repel(data = data.frame(x = by_text_x, y = by_text_y,
                                      label = names(by_text_x)),
                    mapping = aes_string(x = "x", y = "y", label = "label"),
                    inherit.aes = FALSE, size = text_size, colour = text_color,
                    force = 0, max.overlaps = Inf)
}
