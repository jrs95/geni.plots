#' @title Regional plot
#'
#' @description `fig_region` creates a regional plot, i.e. a scatter graph of
#'   genomic markers associations (e.g. log10(p-values)) with a gene bar
#'   underneath.
#'
#' @name fig_region
#'
#' @import dplyr
#'
#' @import tidyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @import ggrepel
#'
#' @import ggnewscale
#'
#' @import patchwork
#'
#' @param data a `data.frame` containing the association statistics for each
#'   marker within a genomic region with the following columns:
#'   \itemize{
#'     \item{\code{marker}} {
#'       the genomic marker identifier (e.g. rsID)
#'     }
#'     \item{\code{chr}} {
#'       the chromosome for each genomic marker
#'     }
#'     \item{\code{pos}} {
#'       the genomic position for each genomic marker
#'     }
#'     and one of the following:
#'     \item{\code{pvalue}} {
#'       the association p-value for each genomic marker
#'     }
#'     \item{\code{z}} {
#'       the association z-statistic for each genomic marker
#'     }
#'     \item{\code{prob}} {
#'       the association probability for each genomic marker
#'     }
#'   }
#'
#' @param corr a `numeric` `matrix` of correlation statistics between
#'   the markers (default: `NULL`)
#'
#' @param corr_top a `numeric` `vector` of correlation statistics between the
#'   top marker and the rest of the markers (default: `NULL`)
#'
#' @param top_marker a `character` value depicting the marker to plot the
#'   correlation statistics of the rest of the markers against (default: `NULL`)
#'
#' @param r2 a `logical` value indicating whether the set of correlation
#'   statistics entered in `corr` or `corr_top` are squared (default: `FALSE`)
#'
#' @param build a `numeric` value indicating the genome build used to determine
#'   genomic position (default: `38` representing human assembly GRCh38)
#'
#' @param prob a `logical` value indicating whether probability statistics
#'   should be plotted instead of -log10(p-values) (default: `FALSE`)
#'
#' @param interactive a `logical` value indicating whether the plot
#'   should be interactive (default: `FALSE`)
#'
#' @param thresh a `numeric` `vector` providing the p-value thresholds to be
#'   plotted (default: `NULL`)
#'
#' @param thresh_colour a `character` `vector` indicating the colours of the
#'   lines indicating the p-value thresholds (default: `"grey50"`)
#'
#' @param x_min a `numeric` value depicting the minimum plotted x-axis value
#'   representing the start of the genomic region (default: `NULL`)
#'
#' @param x_max a `numeric` value depicting the maximum plotted x-axis value
#'   representing the end of the genomic region (default: `NULL`)
#'
#' @param y_title a `character` string defining the title of the y-axis
#'   (default: `NULL`)
#'
#' @param point_size a `numeric` value indicating the size of each point
#'   (default: `3`)
#'
#' @param alpha a `numeric` value adjusting the opacity of colours
#'   representing the correlation statistics (default: `1`)
#'
#' @param genebar a `logical` value indicating whether bars representing
#'   the genes should be included in the plot (default: `TRUE`)
#'
#' @param genebar_ntracks an `integer` value indicating the number of
#'   tracks to be included in the gene bar (default: `NULL`)
#'
#' @param genebar_label_pos a `numeric` value indicating the relative
#'   position of gene labels with respect to each gene bar (default: `3.6`)
#'
#' @param genebar_label_size a `numeric` value defining the size of each
#'   gene label (default: `4.25`)
#'
#' @param genebar_line_size a `numeric` value defining the line size of
#'   each gene bar (default: `0.8`)
#'
#' @param label_size a `numeric` value indicating the size of each label
#'   (default: `3.5`)
#'
#' @param highlights a `character` `vector` defining a set of markers to
#'   highlight in the plot (default: `NULL`)
#'
#' @param highlights_cat a `character` `vector` defining the category for
#'   each highlighted marker (default: `NULL`)
#'
#' @param highlights_label a `logical` value indicating whether highlighted
#'   points should be labelled (default: `TRUE`)
#'
#' @param highlights_shape a value defining the shape for highlighted points
#'   (default: `22`)
#'
#' @param highlights_nolabel_shape a value defining the shape for points
#'   which are not highlighted (default: `21`)
#'
#' @param highlights_sort a `logical` value indicating whether to sort highlight
#'   group label levels (default: `TRUE`)
#'
#' @param highlights_colours a `character` `vector` specifying colours for
#'   highlighted points (default: `NULL`)
#'
#' @param highlights_title a `character` string providing a title for the legend
#'   corresponding to the highlighted points (default: `"Group"`)
#'
#' @param title a `character` string providing a title for the plot
#'   (default: `NULL`)
#'
#' @param title_size a `numeric` value indicating the size of the title text
#'   for the plot (default: `NULL`)
#'
#' @param title_center a `logical` value indicating whether the plot title
#'   should be centered (default: `FALSE`)
#'
#' @param axis_text_size a `numeric` value indicating the size of the axis
#'   text for the plot (default: `14`)
#'
#' @param axis_title_size a `numeric` value indicating the size of the axis
#'   title text for the plot (default: `16`)
#'
#' @param legend a `logical`  value indicating whether a legend corresponding
#'   to the displayed groups should be included (default: `TRUE`)
#'
#' @param legend_text_size a `numeric` value indicating the size of the legend
#'   text (default: `12`)
#'
#' @param legend_title_size a `numeric` value indicating the size of the legend
#'   title (default: `12`)
#'
#' @param point_padding a `numeric` value indicating the relative distance of
#'   labels from plotted points (default: `0`
#'
#' @param nudge_x a `numeric` value indicating the degree to which label
#'   placement on the x-axis should be adjusted (default: `0`)
#'
#' @param nudge_y a `numeric` value indicating the degree to which label
#'   placement on the y-axis should be adjusted (default: `0`)
#'
#' @param nudge_y_top a `numeric` value indicating the degree to which the top
#'   marker should be adjusted on the y-axis by a proportion of the y-axis limit
#'   (default: `0.06`)
#'
#' @param ylim_prob a `numeric` value defining the upper y-axis limit for
#'   probability plots (default: `1`)
#'
#' @param assoc_plot_size a `numeric` value determining the size of the
#'   association plot (default: `NULL`)
#'
#' @param genebar_plot_size a `numeric` value determining the size of the gene
#'   bar plot (default: `NULL`)
#'
#' @param legend_plot_dist a `numeric` value defining the distance and size of
#'   the legend from the bottom of the regional plot (default = `NULL`)
#'
#' @param plot_width a `numeric` value indicating the width of the plot
#'   (default: `9`)
#'
#' @param plot_height a `numeric` value indicating the height of the plot
#'   (default: `7`)
#'
#' @param girafe a `logical` value indicating whether an interactive plot
#'   should be turned into an interactive graphic using
#'   [girafe()][ggiraph::girafe()] (default = `TRUE`)
#'
#' @return `fig_region` returns a regional plot visualising associations
#'   of markers within a genomic region.
#'
#' @examples
#' fig_region(
#'   data = geni.plots::geni_test_region$assoc,
#'   corr = geni.plots::geni_test_region$corr,
#'   build = 37,
#'   axis_text_size = 11,
#'   axis_title_size = 12,
#'   genebar_label_size = 3.5,
#'   legend_text_size = 10,
#'   legend_title_size = 10
#' )
#'
#' # Notes:
#' #   (i) corr has to have the same markers as assoc in the same order,
#' #   (ii) by default fig_region assumes corr contains correlation
#'          statistics that have not been squared
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
fig_region <- function(data, corr = NULL, corr_top = NULL, top_marker = NULL,
  r2 = FALSE, build = 38, prob = FALSE, interactive = FALSE, thresh = NULL,
  thresh_colour = "grey50", x_min = NULL, x_max = NULL, y_title = NULL,
  point_size = 3, alpha = 1, genebar = TRUE, genebar_ntracks = NULL,
  genebar_label_pos = 3.6, genebar_label_size = 4, genebar_line_size = 0.8,
  label_size = 3.5, highlights = NULL, highlights_cat = NULL,
  highlights_label = TRUE, highlights_shape = 22, highlights_nolabel_shape = 21,
  highlights_sort = TRUE, highlights_colours = NULL, highlights_title = "Group",
  title = NULL, title_size = 16, title_center = FALSE, axis_text_size = 14,
  axis_title_size = 16, legend = TRUE, legend_text_size = 12,
  legend_title_size = 12, point_padding = 0, nudge_x = 0, nudge_y = 0,
  nudge_y_top = 0.06, ylim_prob = 1, assoc_plot_size = NULL,
  genebar_plot_size = NULL, legend_plot_dist = NULL, plot_width = 9,
  plot_height = 7, girafe = TRUE) {

  # Errors
  if (any(is.na(data)))
    stop("there are missing values in the dataset")
  if (!is.logical(prob) || length(prob) > 1)
    stop("the type of plot has to be either log10p or prob")
  if (prob == FALSE) {
    if (!all(c("marker", "chr", "pos") %in% names(data)))
      stop("dataset needs to include marker, chr and pos columns")
    if (!any(c("pvalue", "z") %in% names(data)))
      stop("dataset needs to include either a pvalue or z column")
  } else {
    if (!all(c("marker", "chr", "pos", "prob") %in% names(data)))
      stop("dataset needs to include marker, chr, pos and prob columns")
  }
  if (!is.null(corr)) {
    if (ncol(corr) != nrow(data) || nrow(corr) != nrow(data))
      stop(
        "corr has to have the same dimensions as the number of rows in ",
        "the markers dataset"
      )
    if (any(rownames(corr) != data$marker))
      stop("corr has to have the same markers in the same order as the dataset")
  }
  if (length(unique(data$chr)) > 1)
    stop(
      "there should only be markers from one chromosome in the markers dataset"
    )
  if (!(data$chr[1] %in% 1:22))
    stop("the plotting tool is only for autosomal chromosomes")
  if (class(data$pos) != "integer")
    stop("the pos variable has to be an integer")
  if ("pvalue" %in% names(data)) {
    if (class(data$pvalue) != "numeric")
      stop("the pvalue variable has to be an numeric")
  }
  if ("z" %in% names(data)) {
    if (class(data$z) != "numeric")
      stop("the z variable has to be an numeric")
  }
  if ("prob" %in% names(data)) {
    if (class(data$prob) != "numeric")
      stop("the prob variable has to be an numeric")
  }
  if (is.null(corr) && !is.null(corr_top) && is.null(top_marker))
    stop("top_marker must be defined if corr_top is provided")
  if (is.null(corr) && !is.null(corr_top)) {
    if (length(corr_top) != nrow(data))
      stop(
        "corr_top has to have the same length as the number of rows in the ",
        "markers dataset"
      )
  }
  if (!is.null(top_marker) && length(which(top_marker == data$marker)) == 0)
    stop("top_marker is not contained in the markers dataset")
  if (!is.null(top_marker) && length(which(top_marker == data$marker)) > 1)
    stop("top_marker maps to multiple markers in the markers dataset")
  if (!(build %in% c(37, 38)))
    stop("genome build can only be 37 or 38")
  if (!is.null(highlights) && any(is.na(highlights)))
    stop("highlights cannot contain missing values")
  if (
    !is.null(highlights_cat) && any(is.na(highlights_cat)) &&
      length(highlights) != length(highlights_cat)
  )
    stop(
      "highlights_cat cannot contain missing values and has to be the same ",
      "length as highlights"
    )
  if (
    !is.null(highlights_label) && any(is.na(highlights_label)) &&
      !is.logical(highlights_label) &&
      (
        length(highlights) != length(highlights_label) ||
          length(highlights_label) == 1
      )
  )
    stop(
      "highlights_label cannot contain missing values, has to be a logical ",
      "vector and has either be the same length as highlights or of length 1"
    )
  if (!is.null(x_min) && !is.numeric(x_min))
    stop("x_min has to be an integer")
  if (!is.null(x_max) && !is.numeric(x_max))
    stop("x_max has to be an integer")

  # Dataset

  ## Markers
  data <- data %>%
    mutate(
      marker = as.character(marker),
      chr = as.character(chr),
      pos = as.integer(pos)
    )

  ## Compute statistic
  if (prob == FALSE) {
    if ("pvalue" %in% names(data)) {
      data <- data %>%
        mutate(
          stats = -log10(as.numeric(pvalue)),
          stats = if_else(stats > 300, 300, stats)
        )
    } else {
      data <- data %>%
        mutate(
          stats = -(log(2) + pnorm(-abs(as.numeric(z)), log.p = TRUE)) /
            log(10),
          stats = if_else(stats > 1000, 1000, stats)
        )
    }
  } else {
    data <- data %>%
      mutate(stats = as.numeric(prob))
  }

  ## Hover text
  if (interactive == TRUE) {
    if (!("text" %in% names(data))) {
      if (prob == FALSE) {
        if (any(grepl("^rs", data$marker))) {
          if ("pvalue" %in% names(data)) {
            data <- data %>%
              mutate(
                text = paste0(
                  "SNP: ", marker,
                  "<br>p-value: ", signif(as.numeric(pvalue), 3)
                )
              )
          } else {
            data <- data %>%
              mutate(
                text = paste0(
                  "SNP: ", marker,
                  "<br>p-value: ",
                  signif(2 * pnorm(-abs(as.numeric(z))), 3)
                )
              )
          }
        } else {
          if ("pvalue" %in% names(data)) {
            data <- data %>%
              mutate(
                text = paste0(
                  "Marker: ", marker,
                  "<br>p-value: ", signif(as.numeric(pvalue), 3)
                )
              )
          } else {
            data <- data %>%
              mutate(
                text = paste0(
                  "Marker: ", marker, "<br>p-value: ",
                  signif(2 * pnorm(-abs(as.numeric(z))), 3)
                )
              )
          }
        }
      } else {
        if (any(grepl("^rs", data$marker))) {
          data <- data %>%
            mutate(
              text = paste0(
                "SNP: ", marker, "<br>Probability: ", round(stats, 4)
              )
            )
        } else {
          data <- data %>%
            mutate(
              text = paste0(
                "Marker: ", marker, "<br>Probability: ", round(stats, 4)
              )
            )
        }
      }
    } else {
      data <- data %>%
        mutate(text = as.character(text))
    }
    data <- data %>%
      select(marker, chr, pos, text, stats)
  } else {
    data <- data %>%
      select(marker, chr, pos, stats)
  }

  # Plot

  ## Genomic location
  chr <- as.character(first(pull(data, chr)))
  if (!is.null(x_min)) {
    x_min <- as.integer(x_min)
  } else {
    x_min <- min(as.integer(data$pos))
  }
  if (!is.null(x_max)) {
    x_max <- as.integer(x_max)
  } else {
    x_max <- max(as.integer(data$pos))
  }
  if ((x_max - x_min) > 5000000)
    stop("the plotting tool can plot a maximum of 5MB")
  if (!is.null(thresh) && prob == TRUE)
    thresh <- NULL

  ## Max and min
  x_min <- x_min - 0.02 * (x_max - x_min)
  x_max <- x_max + 0.02 * (x_max - x_min)

  ## Gene bar
  if (genebar == TRUE) {
    genes_df <- fig_gene_bar_data(
      chr, x_min, x_max, build, genebar_ntracks, interactive,
      genebar_label_pos
    )
    gene_bar <- fig_gene_bar_plot(
      genes_df$genes_df, chr, x_min, x_max, genes_df$ntracks, interactive,
      genebar_label_size, genebar_line_size, axis_text_size, axis_title_size
    )
  }

  ## Association plot
  stack <- FALSE
  if (genebar == TRUE) {
    x_labels <- FALSE
  } else {
    x_labels <- TRUE
  }
  assoc_df <- fig_region_data(
    data, corr, corr_top, top_marker, r2, interactive, highlights,
    highlights_cat, highlights_label, highlights_sort
  )
  assoc <- fig_region_plot(
    assoc_df, chr, x_min, x_max, build, prob, interactive, stack, thresh,
    thresh_colour, x_labels, y_title, point_size, alpha, label_size,
    highlights_shape, highlights_nolabel_shape, highlights_colours,
    highlights_title, title, title_size, title_center, axis_text_size,
    axis_title_size, legend_text_size, legend_title_size, point_padding,
    nudge_x, nudge_y, nudge_y_top, ylim_prob
  )

  ## Combined plot
  if (genebar == TRUE) {
    if (legend == TRUE) {
      if (!is.null(legend_plot_dist)) {
        fig <- assoc / gene_bar / (ggplot() + theme_void()) &
          theme(
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.position = "bottom", legend.box = "vertical",
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      } else {
        fig <- assoc / gene_bar &
          theme(
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.position = "bottom", legend.box = "vertical",
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      }
    } else {
      fig <- assoc / gene_bar &
        theme(
          legend.position = "none",
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size)
        )
    }
  } else {
    if (legend == TRUE) {
      if (!is.null(legend_plot_dist)) {
        fig <- assoc / (ggplot() + theme_void()) &
          theme(
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.position = "bottom", legend.box = "vertical",
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      } else {
        fig <- assoc &
          theme(
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.position = "bottom", legend.box = "vertical",
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      }
    } else {
      fig <- assoc &
        theme(
          legend.position = "none",
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size)
        )
    }
  }

  ## Title
  if (!is.null(title)) {
    title_hjust <- 0
    if (title_center == TRUE)
      title_hjust <- 0.5
    fig <- fig &
      theme(plot.title = element_text(size = title_size, hjust = title_hjust))
  }

  ## Plot size
  if (is.null(assoc_plot_size))
    assoc_plot_size <- 50
  if (genebar == TRUE && is.null(genebar_plot_size)) {
    if (genes_df$ntracks == 2) {
      genebar_plot_size <- 14
    } else if (genes_df$ntracks > 2 && genes_df$ntracks <= 8) {
      genebar_plot_size <- 24
    } else if (genes_df$ntracks > 8) {
      genebar_plot_size <- 40
    }
  }

  ## Layout
  if (genebar == TRUE) {
    if (legend == TRUE && !is.null(legend_plot_dist)) {
      fig <- fig +
        plot_layout(
          heights = c(assoc_plot_size, genebar_plot_size, legend_plot_dist),
          guides = "collect"
        )
    } else {
      fig <- fig +
        plot_layout(
          heights = c(assoc_plot_size, genebar_plot_size),
          guides = "collect"
        )
    }
  } else {
    if (legend == TRUE && !is.null(legend_plot_dist)) {
      fig <- fig +
        plot_layout(
          heights = c(assoc_plot_size, legend_plot_dist),
          guides = "collect"
        )
    } else {
      fig <- fig +
        plot_layout(heights = c(assoc_plot_size))
    }
  }

  ## Interactive
  if (interactive == TRUE && girafe == TRUE)
    fig <- girafe(
      print(fig),
      width_svg = plot_width, height_svg = plot_height,
      options = list(
        opts_tooltip(
          css = paste0(
            "background-color:",
            "#EEEEEE;color:black;",
            "padding:10px;",
            "border-radius:5px;"
          )
        ),
        opts_sizing(width = 0.85)
      )
    )

  # Return
  return(fig)

}

#' @title Stacked regional plot
#'
#' @description `fig_region_stack` creates a stacked regional association plot.
#'
#' @name fig_region_stack
#'
#' @import dplyr
#'
#' @import tidyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @import ggrepel
#'
#' @import ggnewscale
#'
#' @import patchwork
#'
#' @param data a `data.frame` containing the association statistics for each
#'   marker within a genomic region with the following columns:
#'   \itemize{
#'     \item{\code{marker}} {
#'       the genomic marker identifier (e.g. rsID)
#'     }
#'     \item{\code{chr}} {
#'       the chromosome for each genomic marker
#'     }
#'     \item{\code{pos}} {
#'       the genomic position for each genomic marker
#'     }
#'     and one of the following:
#'     \item{\code{pvalue_1}, \code{pvalue_2}, \code{pvalue_3}, etc.} {
#'       the association p-values for each genomic marker and trait
#'     }
#'     \item{\code{z_1}, \code{z_2}, \code{z_3}, etc.} {
#'       the association z-statistics for each genomic marker and trait
#'     }
#'     \item{\code{prob_1}, \code{prob_2}, \code{prob_3}, etc.} {
#'       the association probabilities for each genomic marker and trait
#'     }
#'   }
#'
#' @param traits a `character` `vector` of trait names
#'
#' @param corr a `numeric` `matrix` of correlation statistics between
#'   the markers (default: `NULL`)
#'
#' @param corr_top a `numeric` `vector` of correlation statistics between the
#'   top marker and the rest of the markers (default: `NULL`)
#'
#' @param top_marker a `character` value depicting the marker to plot the
#'   correlation statistics of the rest of the markers against (default: `NULL`)
#'
#' @param r2 a `logical` value indicating whether the set of correlation
#'   statistics entered in `corr` or `corr_top` are squared (default: `FALSE`)
#'
#' @param build a `numeric` value indicating the genome build used to determine
#'   genomic position (default: `38` representing human assembly GRCh38)
#'
#' @param prob a `logical` value indicating whether probability statistics
#'   should be plotted instead of -log10(p-values) (default: `FALSE`)
#'
#' @param interactive a `logical` value indicating whether the plot
#'   should be interactive (default: `FALSE`)
#'
#' @param thresh a `numeric` `vector` providing the p-value thresholds to be
#'   plotted (default: `NULL`)
#'
#' @param thresh_colour a `character` `vector` indicating the colours of the
#'   lines indicating the p-value thresholds (default: `"grey50"`)
#'
#' @param x_min a `numeric` value depicting the minimum plotted x-axis value
#'   representing the start of the genomic region (default: `NULL`)
#'
#' @param x_max a `numeric` value depicting the maximum plotted x-axis value
#'   representing the end of the genomic region (default: `NULL`)
#'
#' @param y_title a `character` string defining the title of the y-axis
#'   (default: `NULL`)
#'
#' @param point_size a `numeric` value indicating the size of each point
#'   (default: `3`)
#'
#' @param alpha a `numeric` value adjusting the opacity of colours
#'   representing the correlation statistics (default: `1`)
#'
#' @param genebar a `logical` value indicating whether bars representing
#'   the genes should be included in the plot (default: `TRUE`)
#'
#' @param genebar_ntracks an `integer` value indicating the number of
#'   tracks to be included in the gene bar (default: `NULL`)
#'
#' @param genebar_label_pos a `numeric` value indicating the relative
#'   position of gene labels with respect to each gene bar (default: `3.6`)
#'
#' @param genebar_label_size a `numeric` value defining the size of each
#'   gene label (default: `4.25`)
#'
#' @param genebar_line_size a `numeric` value defining the line size of
#'   each gene bar (default: `0.8`)
#'
#' @param label_size a `numeric` value indicating the size of each label
#'   (default: `3.5`)
#'
#' @param highlights a `character` `vector` defining a set of markers to
#'   highlight in the plot (default: `NULL`)
#'
#' @param highlights_cat a `character` `vector` defining the category for
#'   each highlighted marker (default: `NULL`)
#'
#' @param highlights_label a `logical` value indicating whether highlighted
#'   points should be labelled (default: `TRUE`)
#'
#' @param highlights_shape a value defining the shape for highlighted points
#'   (default: `22`)
#'
#' @param highlights_nolabel_shape a value defining the shape for points
#'   which are not highlighted (default: `21`)
#'
#' @param highlights_sort a `logical` value indicating whether to sort highlight
#'   group label levels (default: `TRUE`)
#'
#' @param highlights_colours a `character` `vector` specifying colours for
#'   highlighted points (default: `NULL`)
#'
#' @param highlights_title a `character` string providing a title for the legend
#'   corresponding to the highlighted points (default: `"Group"`)
#'
#' @param title_size a `numeric` value indicating the size of the title text
#'   for the plot (default: `NULL`)
#'
#' @param title_center a `logical` value indicating whether the plot title
#'   should be centered (default: `FALSE`)
#'
#' @param axis_text_size a `numeric` value indicating the size of the axis
#'   text for the plot (default: `14`)
#'
#' @param axis_title_size a `numeric` value indicating the size of the axis
#'   title text for the plot (default: `16`)
#'
#' @param legend a `logical`  value indicating whether a legend corresponding
#'   to the displayed groups should be included (default: `TRUE`)
#'
#' @param legend_text_size a `numeric` value indicating the size of the legend
#'   text (default: `NULL`)
#'
#' @param legend_title_size a `numeric` value indicating the size of the legend
#'   title (default: `NULL`)
#'
#' @param point_padding a `numeric` value indicating the relative distance of
#'   labels from plotted points (default: `0`
#'
#' @param nudge_x a `numeric` value indicating the degree to which label
#'   placement on the x-axis should be adjusted (default: `0`)
#'
#' @param nudge_y a `numeric` value indicating the degree to which label
#'   placement on the y-axis should be adjusted (default: `0`)
#'
#' @param nudge_y_top a `numeric` value indicating the degree to which the top
#'   marker should be adjusted on the y-axis by a proportion of the y-axis limit
#'   (default: `0.06`)
#'
#' @param ylim_prob a `numeric` value defining the upper y-axis limit for
#'   probability plots (default: `1`)
#'
#' @param assoc_plot_size a `numeric` value determining the size of the
#'   association plot (default: `NULL`)
#'
#' @param genebar_plot_size a `numeric` value determining the size of the gene
#'   bar plot (default: `NULL`)
#'
#' @param legend_plot_dist a `numeric` value defining the distance and size of
#'   the legend from the bottom of the regional plot (default = `NULL`)
#'
#' @param plot_width a `numeric` value indicating the width of the plot
#'   (default: `NULL`)
#'
#' @param plot_height a `numeric` value indicating the height of the plot
#'   (default: `NULL`)
#'
#' @param girafe a `logical` value indicating whether an interactive plot
#'   should be turned into an interactive graphic using
#'   [girafe()][ggiraph::girafe()] (default = `TRUE`)
#'
#' @return `fig_region_stack` returns a stacked regional plot.
#'
#' @examples
#' fig <- fig_region_stack(
#'   data = geni.plots::geni_test_stack_region$assoc,
#'   traits = c("Interleukin-6 levels", "Interleukin-6 receptor levels"),
#'   corr = geni.plots::geni_test_stack_region$corr,
#'   build = 37,
#'   highlights = "rs11265611",
#'   title_center = TRUE
#' )
#'
#' # Notes:
#' #   (i) corr has to have the same markers as assoc in the same order,
#' #   (ii) by default fig_region_stack assumes corr contains correlation
#'          statistics that have not been squared
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
fig_region_stack <- function(data, traits, corr = NULL, corr_top = NULL,
  top_marker = NULL, r2 = FALSE, build = 38, prob = FALSE, interactive = FALSE,
  thresh = NULL, thresh_colour = "grey50", x_min = NULL, x_max = NULL,
  y_title = NULL, point_size = 3, alpha = 1, genebar = TRUE,
  genebar_ntracks = NULL, genebar_label_pos = 3.6, genebar_label_size = 4,
  genebar_line_size = 0.8, label_size = 3.5, highlights = NULL,
  highlights_cat = NULL, highlights_label = TRUE, highlights_shape = 22,
  highlights_nolabel_shape = 21, highlights_sort = TRUE,
  highlights_colours = NULL, highlights_title = "Group", title_size = 16,
  title_center = FALSE, axis_text_size = 14, axis_title_size = 16,
  legend = TRUE, legend_text_size = 12, legend_title_size = 12,
  point_padding = 0, nudge_x = 0, nudge_y = 0, nudge_y_top = 0.06,
  ylim_prob = 1, assoc_plot_size = NULL, genebar_plot_size = NULL,
  legend_plot_dist = NULL, plot_width = NULL, plot_height = NULL,
  girafe = TRUE) {

  # Errors
  if (!is.logical(prob) || length(prob) > 1)
    stop("the type of plot has to be either log10p or prob")
  if (prob == FALSE) {
    if (!all(c("marker", "chr", "pos") %in% names(data)))
      stop("dataset needs to include marker, chr and pos columns")
    if (
      !any(grepl("z_", names(data))) &&
        !any(grepl("pvalue_", names(data)))
    )
      stop("dataset needs to include pvalue or z columns")
    if (any(grepl("z_", names(data))) && any(grepl("pvalue_", names(data))))
      stop("dataset cannot include both pvalue and z columns")
    if (any(grepl("z_", names(data)))) {
      if (!all(paste0("z_", seq_along(traits)) %in% names(data)))
        stop("the number of z columns is not equal to the number of traits")
    }
    if (any(grepl("pvalue_", names(data)))) {
      if (!all(paste0("pvalue_", seq_along(traits)) %in% names(data)))
        stop(
          "the number of pvalue columns is not equal to the number of traits"
        )
    }
  } else {
    if (!all(c("marker", "chr", "pos") %in% names(data)))
      stop("dataset needs to include marker, chr and pos columns")
    if (!any(grepl("prob_", names(data))))
      stop("dataset needs to include prob columns")
    if (!any(grepl("prob_", names(data)))) {
      if (!all(paste0("prob_", seq_along(traits)) %in% names(data)))
        stop("the number of prob columns is not equal to the number of traits")
    }
  }
  if (!is.null(corr)) {
    if (ncol(corr) != nrow(data) || nrow(corr) != nrow(data))
      stop(
        "corr has to have the same dimensions as the number of rows in the ",
        "markers dataset"
      )
    if (any(rownames(corr) != data$marker))
      stop("corr has to have the same markers in the same order as the dataset")
  }
  if (length(unique(data$chr)) > 1)
    stop(
      "there should only be markers from one chromosome in the markers dataset"
    )
  if (!(data$chr[1] %in% 1:22))
    stop("the plotting tool is only for autosomal chromosomes")
  if (any(is.na(select(data, marker, chr, pos))))
    stop("there are missing values in the dataset")
  if (class(data$pos) != "integer")
    stop("the pos variable has to be an integer")
  if (is.null(corr) && !is.null(corr_top) && is.null(top_marker))
    stop("top_marker must be defined if corr_top is provided")
  if (is.null(corr) && !is.null(corr_top)) {
    if (length(corr_top) != nrow(data))
      stop(
        "corr_top has to have the same length as the number of rows in the ",
        "markers dataset"
      )
  }
  if (!is.null(top_marker) && length(which(top_marker == data$marker)) == 0)
    stop("top_marker is not contained in the markers dataset")
  if (!is.null(top_marker) && length(which(top_marker == data$marker)) > 1)
    stop("top_marker maps to multiple markers in the markers dataset")
  if (!(build %in% c(37, 38)))
    stop("genome build can only be 37 or 38")
  if (!is.null(highlights) && any(is.na(highlights)))
    stop("highlights cannot contain missing values")
  if (
    !is.null(highlights_cat) && any(is.na(highlights_cat)) &&
      length(highlights) != length(highlights_cat)
  )
    stop(
      "highlights_cat cannot contain missing values and has to be the same ",
      "length as highlights"
    )
  if (
    !is.null(highlights_label) && any(is.na(highlights_label)) &&
      !is.logical(highlights_label) &&
      (
        length(highlights) != length(highlights_label) ||
          length(highlights_label) == 1
      )
  )
    stop(
      "highlights_label cannot contain missing values, has to be a ",
      "logical vector and has either be the same length as highlights ",
      "or of length 1"
    )
  if (!is.null(x_min) && !is.numeric(x_min))
    stop("x_min has to be an integer")
  if (!is.null(x_max) && !is.numeric(x_max))
    stop("x_max has to be an integer")

  # Markers
  data <- data %>%
    mutate(
      marker = as.character(marker),
      chr = as.character(chr),
      pos = as.integer(pos)
    )

  # Plot

  ## Genomic location
  chr <- as.character(first(pull(data, chr)))
  if (!is.null(x_min)) {
    x_min <- as.integer(x_min)
  } else {
    x_min <- min(as.integer(data$pos))
  }
  if (!is.null(x_max)) {
    x_max <- as.integer(x_max)
  } else {
    x_max <- max(as.integer(data$pos))
  }
  if ((x_max - x_min) > 5000000)
    stop("the plotting tool can plot a maximum of 5MB")
  if (!is.null(thresh) && prob == TRUE)
    thresh <- NULL

  ## Max and min
  x_min <- x_min - 0.02 * (x_max - x_min)
  x_max <- x_max + 0.02 * (x_max - x_min)

  ## Genes
  if (genebar == TRUE) {
    genes_df <- fig_gene_bar_data(
      chr, x_min, x_max, build, genebar_ntracks,
      interactive, genebar_label_pos
    )
    gene_bar <- fig_gene_bar_plot(
      genes_df$genes_df, chr, x_min, x_max, genes_df$ntracks, interactive,
      genebar_label_size, genebar_line_size, axis_text_size, axis_title_size
    )
  }

  ## Association plot
  for (i in seq_along(traits)) {

    ### Dataset
    df <- data

    ### Compute statistic
    if (prob == FALSE) {
      if (paste0("pvalue_", i) %in% names(data)) {
        df <- df %>%
          mutate(
            stats = -log10(as.numeric(df[[paste0("pvalue_", i)]])),
            stats = if_else(stats > 300, 300, stats, stats)
          )
      } else {
        df <- df %>%
          mutate(
            stats = -(
              log(2) +
                pnorm(-abs(as.numeric(df[[paste0("z_", i)]])), log.p = TRUE)
            ) / log(10),
            stats = if_else(stats > 300, 300, stats, stats)
          )
      }
    } else {
      df <- df %>%
        mutate(stats = as.numeric(df[[paste0("prob_", i)]]))
    }

    ### Hover text
    if (interactive == TRUE) {
      if (!(paste0("text_", i) %in% names(data))) {
        if (prob == FALSE) {
          if (any(grepl("^rs", data$marker))) {
            df <- df %>%
              mutate(
                text = paste0(
                  "SNP: ", marker, "<br>p-value: ", signif(10^(-1 * stats), 3)
                )
              )
          } else {
            df <- df %>%
              mutate(
                text = paste0(
                  "Marker: ", marker,
                  "<br>p-value: ", signif(10^(-1 * stats), 3)
                )
              )
          }
        } else {
          if (any(grepl("^rs", data$marker))) {
            df <- df %>%
              mutate(
                text = paste0(
                  "SNP: ", marker,
                  "<br>Probability: ", round(stats, 4)
                )
              )
          } else {
            df <- df %>%
              mutate(
                text = paste0(
                  "Marker: ", marker,
                  "<br>Probability: ", round(stats, 4)
                )
              )
          }
        }
      } else {
        df <- df %>%
          mutate(text = as.character(df[[paste0("text_", i)]]))
      }
      df <- df %>%
        select(marker, chr, pos, text, stats)
    } else {
      df <- df %>%
        select(marker, chr, pos, stats)
    }

    ### Association plot
    stack <- TRUE
    if (genebar == TRUE || i < length(traits)) {
      x_labels <- FALSE
    } else {
      x_labels <- TRUE
    }
    assoc_df <- fig_region_data(
      df, corr, corr_top, top_marker, r2, interactive, highlights,
      highlights_cat, highlights_label, highlights_sort
    )
    assoc <- fig_region_plot(
      assoc_df, chr, x_min, x_max, build, prob, interactive, stack, thresh,
      thresh_colour, x_labels, y_title, point_size, alpha, label_size,
      highlights_shape, highlights_nolabel_shape, highlights_colours,
      highlights_title, bquote(italic(.(traits[i]))), title_size,
      title_center, axis_text_size, axis_title_size, legend_text_size,
      legend_title_size, point_padding, nudge_x, nudge_y, nudge_y_top,
      ylim_prob
    )
    if (i == 1) {
      fig <- assoc
    } else {
      fig <- fig / assoc
    }

  }

  ## Combined plot
  if (genebar == TRUE) {
    if (legend == TRUE) {
      if (!is.null(legend_plot_dist)) {
        fig <- fig / gene_bar / (ggplot() + theme_void()) &
          theme(
            legend.position = "bottom", legend.box = "vertical",
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      } else {
        fig <- fig / gene_bar &
          theme(
            legend.position = "bottom", legend.box = "vertical",
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      }
    } else {
      fig <- fig / gene_bar &
        theme(
          legend.position = "none",
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size)
        )
    }
  } else {
    if (legend == TRUE) {
      if (!is.null(legend_plot_dist)) {
        fig <- fig / (ggplot() + theme_void()) &
          theme(
            legend.position = "bottom", legend.box = "vertical",
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      } else {
        fig <- fig &
          theme(
            legend.position = "bottom", legend.box = "vertical",
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.title = element_text(size = legend_title_size),
            legend.text = element_text(size = legend_text_size)
          )
      }
    } else {
      fig <- fig &
        theme(
          legend.position = "none",
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size)
        )
    }
  }

  ## Title
  title_hjust <- 0
  if (title_center == TRUE)
    title_hjust <- 0.5
  fig <- fig &
    theme(plot.title = element_text(size = title_size, hjust = title_hjust))

  ## Plot size
  if (is.null(assoc_plot_size))
    assoc_plot_size <- 50
  if (genebar == TRUE && is.null(genebar_plot_size)) {
    if (genes_df$ntracks == 2) {
      genebar_plot_size <- 8 * length(traits)
    } else if (genes_df$ntracks > 2 && genes_df$ntracks <= 8) {
      genebar_plot_size <- 16 * length(traits)
    } else if (genes_df$ntracks > 8) {
      genebar_plot_size <- 20 * length(traits)
    }
  }

  ## Layout
  if (genebar == TRUE) {
    if (legend == TRUE && !is.null(legend_plot_dist)) {
      fig <- fig +
        plot_layout(
          heights = c(
            rep(assoc_plot_size, length(traits)),
            genebar_plot_size,
            legend_plot_dist
          ),
          guides = "collect"
        )
    } else {
      fig <- fig +
        plot_layout(
          heights = c(rep(assoc_plot_size, length(traits)), genebar_plot_size),
          guides = "collect"
        )
    }
  } else {
    if (legend == TRUE && !is.null(legend_plot_dist)) {
      fig <- fig +
        plot_layout(
          heights = c(rep(assoc_plot_size, length(traits)), legend_plot_dist),
          guides = "collect"
        )
    } else {
      fig <- fig +
        plot_layout(heights = c(rep(assoc_plot_size, length(traits))))
    }
  }

  ## Interactive
  if (interactive == TRUE && girafe == TRUE) {
    if (is.null(plot_width)) {
      plot_width <- 9
    }
    if (is.null(plot_height)) {
      plot_height <- 3 + 4 * length(traits)
    }
    fig <- girafe(
      print(fig),
      width_svg = plot_width, height_svg = plot_height,
      options = list(
        opts_tooltip(
          css = paste0(
            "background-color:#EEEEEE;",
            "color:black;",
            "padding:10px;",
            "border-radius:5px;"
          )
        ),
        opts_sizing(width = 0.85)
      )
    )
  }

  # Return
  return(fig)

}

#' @title Figure regional plot data
#'
#' @description `fig_region_data` creates the `data.frame` used to plot the
#'   regional plot.
#'
#' @name fig_region_data
#'
#' @import dplyr
#'
#' @import tidyr
#'
#' @param data `data.frame` containing the columns:
#'   `marker` (marker name),
#'   `chr` (chromosome)
#'   `pos` (position),
#'   `stats` (either -log10(p-values) or probabilities)
#'
#' @param traits trait names
#'
#' @param corr correlation matrix between markers; default: `NULL`
#'
#' @param corr_top correlation statistics between the top marker and the rest of
#'   the markers; default: `NULL`
#'
#' @param top_marker the marker to plot the correlation statistics in relation
#'   to against the rest of the markers; default: `NULL`
#'
#' @param r2 correlation statistics are squared; default: `FALSE`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param highlights markers to highlight; default: `NULL`
#'
#' @param highlights_cat highlight categories; default: `NULL`
#'
#' @param highlights_label label the highlighted points; default: `TRUE`
#'
#' @param highlights_sort sort highlight group label levels; default: `TRUE`
#'
#' @return `fig_region_data` returns the `data.frame` used to plot the regional
#'   plot.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_region_data <- function(data, corr = NULL, corr_top = NULL,
  top_marker = NULL, r2 = FALSE, interactive = FALSE, highlights = NULL,
  highlights_cat = NULL, highlights_label = TRUE, highlights_sort = TRUE) {

  # Dataset
  data <- data %>%
    as_tibble() %>%
    mutate(
      x = as.integer(data$pos),
      y = as.numeric(data$stats),
      marker = as.character(data$marker)
    )
  if (interactive == TRUE) {
    data <- data %>%
      mutate(text = as.character(text)) %>%
      select(x, y, marker, text)
  } else {
    data <- data %>%
      select(x, y, marker)
  }

  # Highlight points
  nrow_df <- nrow(data)
  if (!is.null(highlights)) {
    data <- data %>%
      mutate(highlight = if_else(marker %in% !!highlights, TRUE, FALSE, FALSE))
    if (!is.null(highlights_cat)) {
      highlights_df <- tibble(
        marker = highlights,
        highlight_cat = as.character(highlights_cat)
      )
      if (highlights_sort == TRUE) {
        highlights_df <- highlights_df %>%
          mutate(
            highlight_cat = factor(
              highlight_cat,
              levels = sort(unique(highlight_cat))
            )
          )
      } else {
        highlights_df <- highlights_df %>%
          mutate(
            highlight_cat = factor(
              highlight_cat,
              levels = unique(highlight_cat)
            )
          )
      }
      highlights_df <- highlights_df %>%
        distinct(marker, .keep_all = TRUE)
      data <- data %>%
        left_join(highlights_df, by = "marker")
    }
    if (length(highlights_label) > 1) {
      highlights_df <- tibble(
        marker = highlights,
        highlight_label = highlights_label
      )
      highlights_df <- highlights_df %>%
        distinct(marker, .keep_all = TRUE)
      data <- data %>%
        left_join(highlights_df, by = "marker") %>%
        replace_na(list(highlight_label = FALSE))
    } else {
      if (highlights_label == TRUE) {
        data <- data %>%
          mutate(
            highlight_label = if_else(
              marker %in% !!highlights,
              TRUE,
              FALSE,
              FALSE
            )
          )
      } else {
        data <- data %>%
          mutate(highlight_label = FALSE)
      }
    }
  } else {
    data <- data %>%
      mutate(highlight = FALSE, highlight_label = FALSE)
  }
  if (nrow(data) != nrow_df)
    stop("number of rows in data different after adding highlight points")

  # Keep non-missing data
  keep <- data %>%
    select(which(!(names(data) %in% "highlight_cat"))) %>%
    complete.cases(.)
  if (!is.null(corr)) {
    corr <- corr[which(keep), which(keep), drop = FALSE]
  }
  if (!is.null(corr_top)) {
    corr_top <- corr_top[keep]
  }
  data <- data %>%
    filter(!!keep)

  # Top marker
  if (!is.null(corr) || !is.null(corr_top)) {
    if (length(top_marker) != 0) {
      top_marker <- which(data$marker == top_marker)
      if (length(top_marker) > 1) {
        top_marker <- sample(top_marker, 1)
        if (is.null(corr) && !is.null(corr_top))
          warning("top_marker maps to multiple markers")
      }
      if (length(top_marker) == 0) {
        top_marker <- which.max(pull(data, y))
      }
      data <- data %>%
        mutate(top = FALSE)
      data$top[top_marker] <- TRUE
    } else {
      top_marker <- which.max(pull(data, y))
      data <- data %>%
        mutate(top = FALSE)
      data$top[top_marker] <- TRUE
    }
  } else {
    data <- data %>%
      mutate(top = FALSE)
  }

  # Correlation statistics
  if (!is.null(corr) || !is.null(corr_top)) {
    if (!is.null(corr)) {
      if (r2 == TRUE) {
        data$r2 <- corr[, top_marker]
      } else {
        data$r2 <- corr[, top_marker]^2
      }
    } else {
      if (r2 == TRUE) {
        data$r2 <- corr_top
      } else {
        data$r2 <- corr_top^2
      }
    }
    if (interactive == TRUE) {
      data <- data %>%
        mutate(text = paste0(text, "<br>r2: ", round(r2, 4)))
    }
    data <- data %>%
      mutate(
        r2_cat = "miss",
        r2_cat = if_else(r2 >= 0 & r2 < 0.2, "0.0-0.2", r2_cat, r2_cat),
        r2_cat = if_else(r2 >= 0.2 & r2 < 0.4, "0.2-0.4", r2_cat, r2_cat),
        r2_cat = if_else(r2 >= 0.4 & r2 < 0.6, "0.4-0.6", r2_cat, r2_cat),
        r2_cat = if_else(r2 >= 0.6 & r2 < 0.8, "0.6-0.8", r2_cat, r2_cat),
        r2_cat = if_else(r2 >= 0.8 & r2 <= 1, "0.8-1.0", r2_cat, r2_cat),
        r2_cat = factor(
          r2_cat,
          levels = c(
            "miss", "0.0-0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1.0"
          )
        )
      )
  }

  # Return
  return(data)

}

#' @title Figure regional plot
#'
#' @description `fig_region_plot` plots a regional plot, i.e. a scatter graph of
#'   genetic associations (e.g. log10 p-values).
#'
#' @name fig_region_plot
#'
#' @import dplyr
#'
#' @import tidyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @import ggrepel
#'
#' @import ggnewscale
#'
#' @param df data.frame containing the plotting variables:
#'   `x` (x-axis),
#'   `y` (y-axis),
#'   `marker` (marker name),
#'   `text` (hover text),
#'   `top` (top marker indicator),
#'   `highlight` (highlight indicator),
#'   `highlight_label` (highlight_label indicator),
#'   `r2` (r-squared),
#'   `r2_cat` (r-squared category)
#'
#' @param chr chromosome
#'
#' @param x_min start of region
#'
#' @param x_max end of region
#'
#' @param build genome build; default: `38`
#'
#' @param prob plot probability statistics; default: `FALSE`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param stack the right-hand axis title for stacked plot; default: `FALSE`
#'
#' @param thresh p-value threshold; default: `NULL`
#'
#' @param thresh_colour p-value threshold colour; default: `"grey50"`
#'
#' @param y_title y-axis left-hand title; default: `NULL`
#'
#' @param point_size point size; default: `3`
#'
#' @param alpha opacity of correlation statistics colours; default: `1`
#'
#' @param label_size label text size; default: `3.5`
#'
#' @param highlights_shape shape for highlighted points; default: `22`
#'
#' @param highlights_nolabel_shape shape for highlights not labelled;
#'   default: `21`
#'
#' @param highlights_colours highlight point colours; default: `NULL`
#'
#' @param highlights_title highlight points legend title; default: `"Group"`
#'
#' @param title figure title; default: `NULL`
#'
#' @param title_size title text size; default: `NULL`
#'
#' @param title_center center title; default: `FALSE`
#'
#' @param axis_text_size axis text size; default: `14`
#'
#' @param axis_title_size axis title size; default: `16`
#'
#' @param legend display legend for groups; default: `TRUE`
#'
#' @param legend_text_size legend text size; default: `NULL`
#'
#' @param legend_title_size legend title size; default: `NULL`
#'
#' @param point_padding point padding for labels; default: `0`
#'
#' @param nudge_x nudge x-axis of labels; default: `0`
#'
#' @param nudge_y nudge y-axis of labels; default: `0`
#'
#' @param nudge_y_top nudge y-axis of top marker label by proportion of
#'   y-axis limit; default: `0.06`
#' @param ylim_prob upper y-axis limit for probability plots; default: `1`
#'
#' @return `fig_region_plot` returns a regional plot.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_region_plot <- function(df, chr, x_min, x_max, build = 38, prob = FALSE,
  interactive = FALSE, stack = FALSE, thresh = NULL, thresh_colour = "grey50",
  x_labels = FALSE, y_title = NULL, point_size = 3, alpha = 1, label_size = 3.5,
  highlights_shape = 22, highlights_nolabel_shape = 21,
  highlights_colours = NULL, highlights_title = "Group", title = NULL,
  title_size = 16, title_center = FALSE, axis_text_size = 14,
  axis_title_size = 16, legend_text_size = 12, legend_title_size = 12,
  point_padding = 0, nudge_x = 0, nudge_y = 0, nudge_y_top = 0.06,
  ylim_prob = 1) {

  # Plot

  ## Base plot
  fig <- ggplot(
    data = filter(df, top == FALSE & highlight == FALSE),
    mapping = aes(x = x, y = y)
  )

  ## Y-axis scale
  if (prob == TRUE) {

    ### Probability scale
    if (!is.null(ylim_prob)) {
      ylim <- ylim_prob
    } else {
      ylim <- max(max(df$y) + 0.1 * max(df$y), 1)
    }
    if (is.null(y_title))
      y_title <- "Probability"

  } else {

    ### mlog10(p) scale
    ylim <- max(pull(df, y), na.rm = TRUE) +
      0.1 * max(pull(df, y), na.rm = TRUE)
    if (!is.null(thresh)) {
      ylim <- max(ylim, (-log10(thresh) + -log10(thresh) * 0.1))
    }
    if (is.null(y_title))
      y_title <- expression("-log"["10"] * paste("(", italic("p"), ")"))

  }

  ## Recombination rate
  recomb_df <- fig_recombination_rate_data(chr, x_min, x_max, build)
  fig <- fig +
    geom_line(
      mapping = aes(x = x, y = y2 / (100 / ylim)),
      data = recomb_df,
      colour = "steelblue1"
    )

  ## Threshold
  if (!is.null(thresh) && prob == FALSE)
    fig <- fig + geom_hline(
      yintercept = -log10(thresh),
      color = thresh_colour, linetype = "dashed"
    )

  ## Point colouring
  if ("r2_cat" %in% names(df)) {

    ### Correlation colouring scale
    if (interactive == TRUE) {
      fig <- fig +
        geom_point_interactive(
          mapping = aes(fill = r2_cat, tooltip = text),
          size = point_size, stroke = point_size / 8, shape = 21, alpha = alpha
        )
    } else {
      fig <- fig +
        geom_point(
          mapping = aes(fill = r2_cat),
          size = point_size, stroke = point_size / 8, shape = 21, alpha = alpha
        )
    }
    fig <- fig +
      scale_fill_manual(
        values = c(
          "#E8E8E8", "#66FFFF", "#66FF66", "#FFCC00", "#FF9933",
          "#CC3300"
        ),
        guide = guide_legend(
          order = 1, nrow = 1, override.aes = list(alpha = 1)
        ),
        drop = FALSE
      )
    fig <- fig +
      labs(fill = "r2")

  } else {

    ### Light grey scale
    if (interactive == TRUE) {
      fig <- fig +
        geom_point_interactive(
          mapping = aes(tooltip = text),
          fill = "#E8E8E8", size = point_size,
          stroke = point_size / 8, shape = 21
        )
    } else {
      fig <- fig +
        geom_point(
          fill = "#E8E8E8", size = point_size,
          stroke = point_size / 8, shape = 21
        )
    }

  }

  ## Highlight points
  if (any(df$highlight == TRUE & df$top == FALSE)) {

    ### Plot highlights
    if (!("highlight_cat" %in% names(df))) {

      #### Set highlights colour
      if (!is.null(highlights_colours)) {
        highlights_colour <- highlights_colours[1]
      } else {
        highlights_colour <- "blue3"
      }

      #### Plot highlights
      if (interactive == TRUE) {
        fig <- fig +
          geom_point_interactive(
            mapping = aes(tooltip = text),
            data = filter(df, highlight == TRUE & top == FALSE),
            fill = highlights_colour, size = point_size,
            stroke = point_size / 8, shape = highlights_shape
          )
      } else {
        fig <- fig +
          geom_point(
            data = filter(df, highlight == TRUE & top == FALSE),
            fill = highlights_colour, size = point_size,
            stroke = point_size / 8, shape = highlights_shape
          )
      }

    } else {

      #### Set new fill scale
      if ("r2_cat" %in% names(df))
        fig <- fig +
          new_scale_fill()

      #### Create highlights data frame
      df_highlight <- df %>%
        filter(highlight == TRUE & !is.na(highlight_cat) & top == FALSE)

      #### Plot highlights
      if (interactive == TRUE) {
        fig <- fig +
          geom_point_interactive(
            mapping = aes(fill = highlight_cat, tooltip = text),
            data = df_highlight,
            size = point_size, stroke = point_size / 8, shape = highlights_shape
          )
      } else {
        fig <- fig +
          geom_point(
            mapping = aes(fill = highlight_cat),
            data = df_highlight,
            size = point_size, stroke = point_size / 8, shape = highlights_shape
          )
      }

      #### Highlight colour scheme & legend point type
      if (!is.null(highlights_colours)) {
        fig <- fig +
          scale_fill_manual(
            values = highlights_colours, drop = FALSE,
            guide = guide_legend(
              order = 2,
              override.aes = list(shape = highlights_shape)
            )
          )
      }

      #### Highlight legend label
      fig <- fig +
        labs(fill = highlights_title)

    }

  }

  ## Top marker
  if (any(df$top)) {

    ### Plot top marker with purple diamond
    if (interactive == TRUE) {
      fig <- fig +
        geom_point_interactive(
          mapping = aes(tooltip = text),
          data = filter(df, top == TRUE),
          fill = "purple", size = point_size,
          stroke = point_size / 8, shape = 23
        )
    } else {
      fig <- fig +
        geom_point(
          data = filter(df, top == TRUE),
          fill = "purple", size = point_size,
          stroke = point_size / 8, shape = 23
        )
    }

  }

  ## Labels
  if (any(df$top | df$highlight_label)) {

    ### Label top marker and the highlighted points specified
    if (any(df$highlight_label) == FALSE) {
      if (max(df$y[df$top == TRUE]) / max(df$y) > 0.99) {
        fig <- fig +
          geom_text_repel(
            mapping = aes(label = marker),
            data = filter(df, top == TRUE),
            size = label_size, segment.color = NA, direction = "y",
            point.padding = point_padding, nudge_y = nudge_y_top * ylim
          )
      } else {
        fig <- fig +
          geom_label_repel(
            mapping = aes(label = marker),
            data = filter(df, top == TRUE),
            size = label_size, segment.size = 0.5, segment.color = "grey",
            point.padding = point_padding, nudge_x = nudge_x, nudge_y = nudge_y
          )
      }
    } else {
      fig <- fig +
        geom_label_repel(
          mapping = aes(label = marker),
          data = filter(df, top == TRUE | highlight_label == TRUE),
          size = label_size, segment.size = 0.5, segment.color = "grey",
          point.padding = point_padding, nudge_x = nudge_x, nudge_y = nudge_y
        )
    }

  }

  ## Theme
  fig <- fig +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      plot.margin = unit(c(0.5, 0, 0, 0), "cm")
    )

  ## Axes limits
  if (x_labels == TRUE) {
    fig <- fig +
      scale_x_continuous(limits = c(x_min, x_max))
  } else {
    fig <- fig +
      scale_x_continuous(limits = c(x_min, x_max), breaks = NULL)
  }
  if (stack == TRUE) {
    fig <- fig +
      scale_y_continuous(
        limits = c(0, ylim),
        sec.axis = sec_axis(
          ~. * (100 / ylim),
          name = "Recomb. rate",
          breaks = c(0, 25, 50, 75, 100)
        )
      )
  } else {
    fig <- fig +
      scale_y_continuous(
        limits = c(0, ylim),
        sec.axis = sec_axis(
          ~. * (100 / ylim),
          name = "Recombination rate (cM/Mb)",
          breaks = c(0, 25, 50, 75, 100)
        )
      )
  }

  ## Axes labels
  if (x_labels == TRUE) {
    fig <- fig +
      xlab(paste0("Position on chromosome ", chr))
  } else {
    fig <- fig +
      xlab(NULL)
  }
  fig <- fig +
    ylab(y_title) +
    theme(
      axis.title.y = element_text(vjust = 2.25, size = axis_title_size),
      axis.title.y.right = element_text(vjust = 1.5, size = axis_title_size),
      axis.text = element_text(size = axis_text_size)
    )

  ## Title
  if (!is.null(title)) {
    title_hjust <- 0
    if (title_center == TRUE)
      title_hjust <- 0.5
    fig <- fig +
      ggtitle(title) +
      theme(plot.title = element_text(size = title_size, hjust = title_hjust))
  }

  ## Legend
  fig <- fig +
    theme(
      legend.text = element_text(size = legend_text_size),
      legend.title = element_text(size = legend_title_size),
      legend.background = element_rect(colour = "black"),
      panel.background = element_rect(fill = NA),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  # Return
  return(fig)

}

#' @title Figure recombination rate data
#'
#' @description `fig_recombination_rate_data` creates the `data.frame` used to
#'   plot the recombination rate. The build 37 (hg19) recombination rate data
#'   is from IMPUTE 2 and the build 38 (hg38) recombination rate data is from
#'   LocusZoom. Both are based on HapMap Phase II.
#'
#' @name fig_recombination_rate_data
#'
#' @import dplyr
#'
#' @param chrom chromosome
#'
#' @param x_min start of region
#'
#' @param x_max end of region
#'
#' @param build genome build; default: `38`
#'
#' @return `fig_recombination_rate_data` returns the `data.frame` used to
#'   plot the recombination rate.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_recombination_rate_data <- function(chrom, x_min, x_max, build = 38) {

  # Recombination data
  if (build == 37) {
    recomb <- geni.plots::recomb_b37 %>%
      filter(chr == chrom & pos >= x_min & pos <= x_max)
  } else if (build == 38) {
    recomb <- geni.plots::recomb_b38 %>%
      filter(chr == chrom & pos >= x_min & pos <= x_max)
  } else {
    stop("genome build has to be either 37 or 38")
  }
  recomb <- recomb %>%
    mutate(
      x = as.integer(pos),
      y2 = as.numeric(recomb_rate)
    ) %>%
    select(x, y2)

  # Return
  return(recomb)

}

#' @title Figure gene bar data
#'
#' @description `fig_gene_bar_data` creates the `data.frame` used to plot
#'   the gene bar.
#'
#' @name fig_gene_bar_data
#'
#' @import dplyr
#'
#' @param chrom chromosome
#'
#' @param x_min start of region
#'
#' @param x_max end of region
#'
#' @param build genome build; default: `38`
#'
#' @param ntracks number of tracks in gene bar; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param label_pos gene label position in relation to the gene lines;
#'   default: `3.6`
#'
#' @return `fig_gene_bar_data` returns the `data.frame` used to plot the
#'   gene bar.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_gene_bar_data <- function(chrom, x_min, x_max, build = 38, ntracks = NULL,
  interactive = FALSE, label_pos = 3.6) {

  # Gene information
  if (build == 37) {
    genes <- geni.plots::genes_pos_b37 %>%
      filter(chr == chrom & !(end < x_min) & !(start > x_max))
  } else if (build == 38) {
    genes <- geni.plots::genes_pos_b38 %>%
      filter(chr == chrom & !(end < x_min) & !(start > x_max))
  } else {
    stop("genome build has to be either 37 or 38")
  }
  ngenes <- nrow(genes)

  # Process gene information
  if (ngenes > 0) {

    ## Start and end of window
    genes <- genes %>%
      mutate(
        start = as.numeric(start),
        start = if_else(start < x_min, x_min, start, start),
        end = as.numeric(end),
        end = if_else(end > x_max, x_max, end, end)
      ) %>%
      arrange(start)

    ## Extend small genes
    small_gene <- (genes$end - genes$start) < (x_max - x_min) / 190
    genes <- genes %>%
      mutate(
        mid_point = start + (end - start) / 2,
        start = if_else(
          !!small_gene,
          mid_point - (x_max - x_min) / 380, start, start
        ),
        end = if_else(!!small_gene, mid_point + (x_max - x_min) / 380, end, end)
      )

    ## Variables
    genes <- genes %>%
      mutate(id = row_number())
    if (interactive == TRUE) {
      genes <- genes %>%
        mutate(
          text = paste0("Gene: ", gene, "<br>Ensembl ID: ", ensembl),
          onclick = sprintf(
            "window.open(\"%s%s%s\")",
            "https://www.ensembl.org/Multi/Search/Results?q=",
            ensembl, ";site=ensembl_all"
          )
        ) %>%
        select(id, gene, ensembl, text, onclick, chr, mid_point, start, end)
    } else {
      genes <- genes %>%
        select(id, gene, ensembl, chr, mid_point, start, end)
    }

    ## Dataset
    genes_df <- genes %>%
      pivot_longer(
        cols = c("start", "end"),
        names_to = "start_end",
        values_to = "x"
      )
    genes_df$overlap <- c(
      FALSE,
      (
        genes_df$x[2:length(genes_df$x)] -
          genes_df$x[1:(length(genes_df$x) - 1)]
      ) <= 0
    )
    if (interactive == TRUE) {
      genes_df <- genes_df %>%
        select(id, gene, ensembl, text, onclick, overlap, chr, mid_point, x)
    } else {
      genes_df <- genes_df %>%
        select(id, gene, ensembl, overlap, chr, mid_point, x)
    }

    ## Tracks
    if (is.null(ntracks)) {
      if (ngenes > 0 && ngenes <= 8) {
        ntracks <- 2
        if (ngenes > 2 && any(pull(genes_df, overlap)))
          ntracks <- 4
      } else if (ngenes > 8 && ngenes <= 16) {
        ntracks <- 4
        if (any(pull(genes_df, overlap)))
          ntracks <- 6
      } else if (ngenes > 16 && ngenes <= 24) {
        ntracks <- 6
        if (any(pull(genes_df, overlap)))
          ntracks <- 8
      } else if (ngenes > 24 && ngenes <= 32) {
        ntracks <- 8
        if (any(pull(genes_df, overlap)))
          ntracks <- 10
      } else if (ngenes > 32) {
        ntracks <- 10
        if (any(pull(genes_df, overlap)))
          ntracks <- 12
      }
    }
    genes_df <- genes_df %>%
      select(-overlap)

    ## Y-axis
    genes_df$y <- (
      8 * ntracks - 8 *
        rep(
          rep(1:ntracks, each = 2),
          times = ceiling(nrow(genes_df) / (2 * ntracks))
        )[seq_len(nrow(genes_df))]
    )
    genes_df <- genes_df %>%
      mutate(ylab = y + label_pos)

  } else {
    genes_df <- tibble(x = c(x_min, x_max), y = c(16, 8))
    ntracks <- 2
  }

  # Output
  out <- list(genes_df = genes_df, ntracks = ntracks)

  # Return
  return(out)

}

#' @title Figure gene bar plot
#'
#' @description `fig_gene_bar_plot` plots the gene bar.
#'
#' @name fig_gene_bar_plot
#'
#' @import dplyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @param df genes `data.frame`
#'
#' @param chr chromosome
#'
#' @param x_min start of region
#'
#' @param x_max end of region
#'
#' @param ntracks number of tracks in gene bar; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param label_size gene label size; default: `4.25`
#'
#' @param line_size gene line size; default: `0.8`
#'
#' @param label_size label text size; default: `3.5`
#'
#' @param axis_text_size axis text size; default: `14`
#'
#' @param axis_title_size axis title size; default: `16`
#'
#' @return `fig_gene_bar_plot` returns the gene bar.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_gene_bar_plot <- function(df, chr, x_min, x_max, ntracks,
  interactive = FALSE, label_size = 4, line_size = 0.8,
  axis_text_size = 14, axis_title_size = 16) {

  # Plot

  ## Base plot
  fig <- ggplot(
    data = df,
    mapping = aes(x = x, y = y)
  )

  ## Gene bar
  if (nrow(df) > 0 && ("id" %in% names(df))) {
    if (interactive == TRUE) {
      fig <- fig +
        geom_line_interactive(
          mapping = aes(group = id, tooltip = text, onclick = onclick),
          colour = "blue4", size = line_size
        )
    } else {
      fig <- fig +
        geom_line(
          mapping = aes(group = id),
          colour = "blue4", size = line_size
        )
    }
  }

  ## Axes
  fig <- fig +
    xlab(paste0("Position on chromosome ", chr)) +
    scale_x_continuous(limits = c(x_min, x_max)) +
    scale_y_continuous(limits = c(-5, (8 * ntracks + 1)))

  ## Gene labels
  if (nrow(df) > 0 && ("id" %in% names(df))) {
    fig <- fig +
      geom_text(
        mapping = aes(x = mid_point, y = ylab, label = gene),
        data = df %>%
          distinct(id, .keep_all = TRUE) %>%
          select(id, gene, mid_point, ylab),
        color = "black", size = (label_size - 0.5 * log((ntracks - 1))),
        fontface = 3
      )
  }

  ## Theme
  fig <- fig +
    theme_bw() +
    theme(
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size),
      axis.title.x = element_text(vjust = -0.5),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )

  # Return
  return(fig)

}

#' @title Regional plot save
#'
#' @description `fig_region_save` saves the regional plots to file.
#'
#' @name fig_region_save
#'
#' @import dplyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @import ggrepel
#'
#' @import patchwork
#'
#' @param fig a regional plot object
#'
#' @param file a `character` string of the file path to save the plot
#'
#' @param ntraits an `integer` value of the number of traits plotted
#'   (default: `1`)
#'
#' @param interactive a `logical` value indicating whether the plot
#'   is interactive (default: `FALSE`)
#'
#' @param width a `numeric` value indicating the width of the plot
#'   in inches (default: `NULL`)
#'
#' @param height a `numeric` value indicating the height of the plot
#'   in inches (default: `NULL`)
#'
#' @param dpi the resolution of the plot (default: `500`)
#'
#' @examples
#' fig <- fig_region(
#'   data = geni.plots::geni_test_region$assoc,
#'   corr = geni.plots::geni_test_region$corr,
#'   build = 37
#' )
#' fig_region_save(fig, "test.png")
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
fig_region_save <- function(fig, file, ntraits = 1, interactive = FALSE,
  width = NULL, height = NULL, dpi = 500) {

  # Plot dimensions (in inches)
  if (is.null(width)) {
    width <- 9
  }
  if (is.null(height)) {
    if (ntraits == 1) {
      height <- 7
    } else {
      height <- 3 + 4 * ntraits
    }
  }

  # Save
  if (interactive == TRUE) {
    htmlwidgets::saveWidget(fig, file)
  } else {
    ggsave(
      file, plot = fig,
      width = width, height = height,
      units = "in", limitsize = FALSE, dpi = dpi
    )
  }

}
