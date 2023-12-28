#' @title Manhattan plot
#'
#' @description `fig_manhattan` creates a Manhattan plot for genomic markers
#'   from across the genome, e.g. results from genome-wide association studies.
#'
#' @details This plotting function plots a Manhattan plot for genomic markers
#'   from across the genome. The default is to truncate these results to
#'   p-value cut-off of `1e-30`.
#'
#' @name fig_manhattan
#'
#' @import dplyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @import ggrepel
#'
#' @param data a `data.frame` containing the chromosome-position information and
#'   association statistics for each genomic marker (e.g., genetic variants)
#'   with the following columns:
#'   \itemize{
#'     \item{\code{chr}} {
#'       the chromosome for each genomic marker
#'     }
#'     \item{\code{pos}} {
#'       the genomic position for each genomic marker
#'     }
#'     \item{\code{pvalue}} {
#'       the association p-value for each genomic marker
#'     }
#'     \item{\code{highlight}} {
#'       the optional highlight points variable, where
#'       `0` = point not highlighted,
#'       `1` = first highlight colour,
#'       `2` = second highlight colour, etc.
#'     }
#'     \item{\code{highlight_shape}} {
#'       the optional highlight point shape variable, where
#'       `0` = standard circle,
#'       `1` = standard circle with border,
#'       `2` = standard rectangle,
#'       `3` = standard rectangle with border,
#'       `4` = standard diamond,
#'       `5` = standard diamond with border
#'     }
#'     \item{\code{label}} {
#'       the optional point labelling variable (e.g. gene name),
#'       if `label = ""` for a point then no
#'       label is presented for that point
#'     }
#'     \item{\code{text}} {
#'       the optional hover text variable for interactive plots to
#'       display further information,
#'       if `text = ""` for a point then no hover text is
#'       presented for that point
#'     }
#'   }
#'
#' @param colours a `vector` of colours used to differentiate chromosomes,
#'   these can either be a separate colour for each chromosome or
#'   a pair of colours (default: `c("#7395D3", "#1D345D")`)
#'
#' @param rank_pos a `logical` value whether genomic markers should be
#'   plotted by their rank position (default: `FALSE`)
#'
#' @param thin_thresh a `numeric` value representing the minimum p-value
#'   threshold for genomic markers to be displayed (default: `NULL`)
#'
#' @param block_thresh a `numeric` value for representing a p-value threshold,
#'   above which genomic markers are represented using blocks of colour
#'   (default: `NULL`)
#'
#' @param interactive a `logical` value indicating whether the plot
#'   should be interactive (default: `FALSE`)
#'
#' @param interactive_n a `numeric` value indicating the number of top
#'   associated points (max = 100,000) to present before using blocks of colour
#'   to minimise file size (default: `NULL`)
#'
#' @param thresh a `numeric` `vector` providing p-value thresholds to be plotted
#'   (default: `c(1e-5, 5e-8)`)
#'
#' @param thresh_size a `numeric` value indicating the width of the lines
#'   indicating the p-value thresholds (default: `0.5`)
#'
#' @param thresh_colours a `character` `vector` indicating the colours of the
#'   lines indicating the p-value thresholds (default: `c("grey50", "red")`)
#'
#' @param trunc a `numeric` value representing the maximum p-value for which
#'   results are displayed (default: `1e-30`)
#'
#' @param highlight_colours a `vector` indicating the colours for the
#'   highlighted genomic markers (default: `NULL`)
#'
#' @param point_size a `numeric` value indicating the size of each point
#'   (default: `2`)
#'
#' @param chr_dist a `numeric` value indicating the gap between different
#'   chromosomes (default: `10000000`)
#'
#' @param x_labels a `logical` value whether the x-axis should be labelled
#'   (default: `TRUE`)
#'
#' @param label_top a `logical` value whether the top associated points should
#'   be labelled (default: `TRUE`)
#'
#' @param label_thresh a `numeric` value providing a p-value threshold for
#'   labelling points (default: `1e-5`)
#'
#' @param label_size a `numeric` value indicating the size of each label
#'   (default: `3`)
#'
#' @param label_ylim a `numeric` value indicating maximum y-axis value at which
#'   labels can be displayed (default: `-log10(1e-5)`)
#'
#' @param label_nudge_y a `numeric` value indicating the degree to which label
#'   placement on the y-axis should be adjusted (default: `0`)`
#'
#' @param label_box a `logical` value indicating whether labels should be
#'   surrounded by a box (default: `FALSE`)
#'
#' @param title a `character`` string providing a title for the plot
#'   (default: `NULL`)
#'
#' @param title_size a `numeric` value indicating the size of the title text
#'   for the plot (default: `NULL`)
#'
#' @param title_center a `numeric` value indicating whether the plot title
#'   should be centered (default: `FALSE`)
#'
#' @param axis_text_size a `numeric` value indicating the size of the axis
#'   text for the plot (default: `NULL`)
#'
#' @param axis_title_size a `numeric` value indicating the size of the axis
#'   title text for the plot (default: `NULL`)
#'
#' @param plot_width a `numeric` value indicating the width of the plot
#'   (default: `16`)
#'
#' @param plot_height a `numeric` value indicating the height of the plot
#'   (default: `8`)
#'
#' @param girafe a `logical` value indicating whether an interactive plot
#'   should be turned into an interactive graphic using
#'   [girafe()][ggiraph::girafe()] (default = `TRUE`)
#'
#' @return `fig_manhattan` returns a Manhattan plot for genomic markers
#'   from across the genome, e.g. results from genome-wide association studies.
#'
#' @examples
#' fig_manhattan(
#'   data = geni.plots::geni_test_manhattan,
#'   block_thresh = 1e-4,
#'   label_box = TRUE
#' )
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
fig_manhattan <- function(data, colours = c("#7395D3", "#1D345D"),
  rank_pos = FALSE, thin_thresh = NULL, block_thresh = NULL,
  interactive = FALSE, interactive_n = NULL, thresh = c(1e-5, 5e-8),
  thresh_size = 0.5, thresh_colours = c("grey50", "red"), trunc = 1e-30,
  highlight_colours = NULL, point_size = 2, chr_dist = 10000000,
  x_labels = TRUE, label_top = TRUE, label_thresh = 1e-5, label_size = 3,
  label_ylim = -log10(1e-5), label_nudge_y = 0, label_box = FALSE,
  title = NULL, title_size = NULL, title_center = FALSE,
  axis_text_size = NULL, axis_title_size = NULL, plot_width = 16,
  plot_height = 8, girafe = TRUE) {

  # Errors
  if (any(is.na(data)))
    stop("data cannot contain missing values")
  if (!all(c("chr", "pos", "pvalue") %in% names(data)))
    stop(
      "the variables chr, pos & pvalue are required to be in the data ",
      "data.frame"
    )
  if (!all(as.character(data$chr) %in% as.character(c(1:22, "X", "Y"))))
    stop("chromosomes have to either be autosomal or sex")
  if (!is.integer(data$pos))
    stop("position has to be an integer vector")
  if ("highlight" %in% names(data)) {
    if (is.logical(data$highlight))
      stop("highlight has to be a logical vector")
    if (!is.numeric(data$highlight) && !is.integer(data$highlight))
      stop("highlight has to be a either a numeric or integer vector")
    if (sum(data$highlight != 0) > 0.1 * nrow(data))
      stop("more than 10% of the data is highlighted")
  }
  if (!is.null(colours)) {
    if (any(is.na(colours)))
      stop("colours cannot contain missing values")
    if (any(duplicated(colours)))
      stop("colours cannot contain duplicates")
    if (length(colours) != 2 && length(colours) != length(unique(data$chr)))
      stop(
        "colours has to be the same length as the number of chromosomes or ",
        "of length 2"
      )
  } else {
    colours <- c("lightblue3", "dodgerblue4")
  }
  if (!is.null(thresh)) {
    if (length(thresh) != length(thresh_colours) && length(thresh_colours) != 1)
      stop(
        "thresh_colours needs to either be the length of the thresh or ",
        "of length 1"
      )
  }
  if (!is.null(highlight_colours)) {
    if ((length(highlight_colours) + 1) != length(unique(data$highlight)))
      stop(
        "highlight_colours has to be the same length as the non-zero ",
        "highlight entries"
      )
  }
  if ("highlight_shape" %in% names(data) && !all(data$highlight_shape %in% 0:5))
    stop("highlight_shape can only contain the values 0 to 5")
  if (interactive == TRUE && !("text" %in% names(data)))
    stop(
      "the text variable must be provided in the data data.frame if the plot ",
      "is specified to be interactive"
    )
  if (label_top == TRUE && !("label" %in% names(data)))
    stop(
      "the label variable must be provided in the data data.frame if the plot ",
      "is specified to have labels"
    )
  if (
    !is.null(block_thresh) &&
      block_thresh > 0.5 * max(as.numeric(data$pvalue), na.rm = TRUE)
  )
    stop("block_thresh has to be less than 0.5 * max(pvalue)")

  # Data
  data <- fig_manhattan_data(
    data, rank_pos, thin_thresh, interactive, trunc, chr_dist,
    label_top, label_thresh
  )

  # Plot
  fig <- fig_manhattan_plot(
    data$df, data$df_plot, colours, rank_pos, block_thresh,
    interactive, interactive_n, thresh, thresh_size, thresh_colours,
    highlight_colours, point_size, x_labels, label_top, label_size,
    label_ylim, label_nudge_y, label_box, title, title_size,
    title_center, axis_text_size, axis_title_size
  )
  if (interactive == TRUE && girafe == TRUE)
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

  # Return
  return(fig)

}

#' @title Figure Manhattan data
#'
#' @description `fig_manhattan_data` creates the `data.frames` used to plot the
#'   Manhattan plot.
#'
#' @name fig_manhattan_data
#'
#' @import dplyr
#'
#' @param data data.frame containing the columns:
#'   `chr` (chromosome),
#'   `pos` (position),
#'   `pvalue`,
#'   `highlight` (highlight points in increasing priority, an integer variable
#'     where `0` = point not highlighted, `1` = first highlight colour,
#'     `2` = second highlight colour, etc.),
#'   `highlight_shape` (highlight point shape, an integer variable where
#'     `0` = standard circle, `1` = standard circle with border,
#'     `2` = standard rectangle, `3` = standard rectangle with border,
#'     `4` = standard diamond, `5` = standard diamond with border),
#'   `label` (point labelling variable, if `label = ""` for a point then no
#'     label is presented for that point),
#'   `text` (hover text variable, if `text = ""` for a point then no hover
#'     is presented for that point)
#'
#' @param rank_pos rank the positions; default: `FALSE`
#'
#' @param thin_thresh thin the data to only include p-values below the
#'   threshold; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param trunc truncate the association results to this values;
#'   default: `1e-30`
#'
#' @param chr_dist distance between chromosomes; default: `10000000`
#'
#' @param label_top label the top associated points; default: `TRUE`
#'
#' @param label_thresh label all points with p-value less than this threshold;
#'   default: `1e-5`
#'
#' @return `fig_manhattan_data` returns a list with the `data.frames` used
#'   to plot the Manhattan plot.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_manhattan_data <- function(data, rank_pos = FALSE, thin_thresh = NULL,
  interactive = FALSE, trunc = 1e-30, chr_dist = 10000000, label_top = TRUE,
  label_thresh = 1e-5) {

  # Data
  data <- data %>%
    as_tibble()

  ## Thin threshold
  data <- data %>%
    filter(!is.na(as.numeric(pvalue)))
  if (!is.null(thin_thresh))
    data <- data %>%
      filter(as.numeric(pvalue) < !!thin_thresh)

  ## Chromosome count
  chr_count <- data %>%
    count(chr) %>%
    filter(n >= 100)
  data <- data %>%
    filter(chr %in% !!pull(chr_count, chr))

  ## Variables

  ### mlog10p
  data <- data %>%
    mutate(y = round(-log10(as.numeric(pvalue)), 4)) %>%
    select(-pvalue)
  if (!is.null(trunc))
    data <- data %>%
      mutate(
        y = if_else(
          y > -log10(!!trunc),
          round(-log10(!!trunc), 4),
          y,
          y
        )
      )

  ### Chromosome
  data <- data %>%
    mutate(chr = as.character(chr))
  chromosomes <- as.character(c(1:22, "X", "Y"))
  chromosomes <- chromosomes[chromosomes %in% unique(pull(data, chr))]
  chromosomes_int <- chromosomes
  chromosomes_int[chromosomes_int == "X"] <- 23
  chromosomes_int[chromosomes_int == "Y"] <- 24
  chromosomes_int <- as.integer(chromosomes_int)
  data <- data %>%
    mutate(
      chr_cat = factor(chr, levels = !!chromosomes, labels = !!chromosomes_int),
      chr_int = as.integer(chr_cat)
    )

  ### Position
  data <- data %>%
    mutate(x = NA)
  max_position <- chr_dist
  if (max_position > 1000000 && rank_pos == TRUE) {
    max_position <- 50000
  }
  for (i in sort(unique(pull(data, chr_int)))) {
    if (rank_pos == TRUE) {
      data <- data %>%
        mutate(
          x = if_else(
            chr_int == !!i,
            rank(pos, na.last = "keep") + !!max_position,
            x,
            x
          )
        )
      max_position <- max_position +
        max(rank(data$pos[data$chr_int == i], na.last = "keep")) -
        min(rank(data$pos[data$chr_int == i], na.last = "keep")) +
        chr_dist
    } else {
      data$x[data$chr_int == i] <- data$pos[data$chr_int == i] -
        min(data$pos[data$chr_int == i]) +
        max_position
      max_position <- max_position +
        max(data$pos[data$chr_int == i]) -
        min(data$pos[data$chr_int == i]) +
        chr_dist
    }
  }

  ### Highlights
  if (!("highlight" %in% names(data)))
    data <- data %>%
      mutate(highlight = 0)
  if (!("highlight_shape" %in% names(data)))
    data <- data %>%
      mutate(highlight_shape = 0)

  ### Labels & hover text
  if (label_top == TRUE) {
    data <- data %>%
      mutate(label = as.character(label))
    if (
      !is.null(label_thresh) && length(label_thresh) &&
        is.numeric(label_thresh)
    ) {
      data <- data %>%
        mutate(label = if_else(y < -log10(!!label_thresh), "", label, ""))
    } else {
      data <- data %>%
        mutate(label = if_else(y < -log10(1e-5), "", label, ""))
    }
  }
  if (interactive == TRUE) {
    data <- data %>%
      mutate(text = as.character(text))
    if (
      !is.null(label_thresh) && length(label_thresh) &&
        is.numeric(label_thresh)
    ) {
      data <- data %>%
        mutate(text = if_else(y < -log10(!!label_thresh), "", text, ""))
    } else {
      data <- data %>%
        mutate(text = if_else(y < -log10(1e-5), "", text, ""))
    }
  }

  # Datasets

  ## Re-order columns
  data <- data %>%
    relocate(x, y)

  ## Order by chromosome-position and point highlighting
  data <- data %>%
    arrange(chr_int, pos) %>%
    arrange(highlight, highlight_shape)

  ## Chromosome mid-point
  data_plot <- data %>%
    group_by(chr_int) %>%
    summarize(
      chr = first(chr),
      mid_point = min(x) + (max(x) - min(x)) / 2,
      .groups = "drop"
    )

  # Output
  out <- list(df = data, df_plot = data_plot)

  # Return
  return(out)

}

#' @title Figure Manhattan plot
#'
#' @description `fig_manhattan_plot` plots a Manhattan plot for genome-wide
#'   association studies.
#'
#' @name fig_manhattan_plot
#'
#' @import dplyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @import ggrepel
#'
#' @param df `data.frame` containing the plotting variables:
#'   `x` (x-axis),
#'   `y` (y-axis),
#'   `chr` (chromosome),
#'   `chr_cat` (chromosome category),
#'   `chr_int` (chromosome integer),
#'   `pos` (position),
#'   `highlight` (highlight points),
#'   `label` (point labels, if `label = ""` for a point then no label is
#'     presented for that point),
#'   `text` (hover text, if `text = ""` for a point then no hover is
#'     presented for that point)
#'
#' @param df_plot `data.frame` containing the x-axis labelling:
#'   `chr_int` (chromosome integer),
#'   `chr` (chromosome),
#'   `mid_point` (x-axis mid-point)
#'
#' @param colours colours of the groups; default: `c("#7395D3", "#1D345D")`
#'
#' @param rank_pos rank the positions; default: `FALSE`
#'
#' @param block_thresh thin the data to present blocks of colour for
#'   p-values above the threshold for static plots; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param interactive_n number of top associated points (max = 100,000) to
#'   present prior to presenting blocks of colour, required to prevent huge
#'   file size; default: `NULL`
#'
#' @param thresh p-value thresholds; default: `c(1e-5, 5e-8)`
#'
#' @param thresh_size threshold size; default: `0.5`
#'
#' @param thresh_colours threshold colours; default: `c("grey50", "red")`
#'
#' @param highlight_colours highlight colours; `NULL`
#'
#' @param point_size point size; default: `2`
#'
#' @param x_labels label the x-axis; default: `TRUE`
#'
#' @param label_top label the top associated points; default: `TRUE`
#'
#' @param label_size label text size; default: `3`
#'
#' @param label_ylim y-axis label limit maximum; default: `-log10(1e-5)`
#'
#' @param label_nudge_y nudge y-axis of labels; default: `0`
#'
#' @param label_box put the label in a box; default: `FALSE`
#'
#' @param title figure title; default: `NULL`
#'
#' @param title_size title text size; default: `NULL`
#'
#' @param title_center center title; default: `FALSE`
#'
#' @param axis_text_size axis text size; default: `NULL`
#'
#' @param axis_title_size axis title size; default: `NULL`
#'
#' @return `fig_manhattan_plot` returns a Manhattan plot for genome-wide
#'   association studies.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_manhattan_plot <- function(df, df_plot, colours = c("#7395D3", "#1D345D"),
  rank_pos = FALSE, block_thresh = NULL, interactive = FALSE,
  interactive_n = NULL, thresh = c(1e-5, 5e-8), thresh_size = 0.5,
  thresh_colours = c("grey50", "red"), highlight_colours = NULL, point_size = 2,
  x_labels = TRUE, label_top = TRUE, label_size = 3, label_ylim = -log10(1e-5),
  label_nudge_y = 0, label_box = FALSE, title = NULL, title_size = NULL,
  title_center = FALSE, axis_text_size = NULL, axis_title_size = NULL) {

  # Plot

  ## Base plot
  fig <- ggplot(
    data = filter(df, highlight == 0),
    mapping = aes(x = x, y = y)
  )

  ## Points
  if (interactive == TRUE) {
    if (!is.null(interactive_n)) {
      df_cut_off <- df %>%
        filter(highlight == 0) %>%
        arrange(desc(y))
      if (interactive_n > 100000) {
        interactive_n <- 100000
      }
      cut_off <- df_cut_off$y[interactive_n]
      if (rank_pos == TRUE) {
        fig <- fig +
          geom_ribbon(
            mapping = aes(
              x = x, ymin = 0, ymax = cut_off,
              colour = chr_cat, fill = chr_cat
            ),
            data = filter(df, y >= !!cut_off & highlight ==  0),
            size = point_size,
            outline.type = "full"
          )
      } else {
        df_cut_off <- df %>%
          filter(y <= cut_off & highlight == 0)
        df_cut_off_window <- tibble()
        for (i in unique(pull(df_cut_off, chr_int))) {
          df_cut_off_chr <- df_cut_off %>%
            filter(chr_int == !!i)
          df_cut_off_chr$window <- cut(
            df_cut_off_chr$x,
            seq(
              min(df_cut_off_chr$x),
              max(df_cut_off_chr$x),
              100000
            ),
            include.lowest = TRUE,
            labels = FALSE
          )
          df_cut_off_chr_window <- tibble(
            window = seq_along(
              zoo::rollmean(
                seq(
                  min(df_cut_off_chr$x),
                  max(df_cut_off_chr$x),
                  100000
                ),
                2
              )
            ),
            mid_point = zoo::rollmean(
              seq(
                min(df_cut_off_chr$x),
                max(df_cut_off_chr$x),
                100000
              ),
              2
            )
          )
          df_cut_off_chr_window <- df_cut_off_chr %>%
            group_by(window) %>%
            summarize(
              y_cut_off = max(y),
              .groups = "drop"
            ) %>%
            left_join(
              x = df_cut_off_chr_window,
              y = .,
              by = "window"
            ) %>%
            mutate(
              y_cut_off = if_else(is.na(y_cut_off), 0, y_cut_off),
              chr_cat = first(pull(df_cut_off_chr, chr_cat))
            )
          df_cut_off_window <- df_cut_off_window %>%
            bind_rows(df_cut_off_chr_window)
        }
        fig <- fig +
          geom_ribbon(
            mapping = aes(
              x = mid_point, ymin = 0, ymax = y_cut_off,
              colour = chr_cat, fill = chr_cat
            ),
            data = df_cut_off_window,
            size = point_size,
            outline.type = "full",
            inherit.aes = FALSE
          )
      }
      fig <- fig +
        geom_point_interactive(
          mapping = aes(colour = chr_cat, fill = chr_cat, tooltip = text),
          data = filter(df, y >= !!cut_off & highlight == 0),
          shape = 21, size = point_size, stroke = 0
        )
      warning(
        paste0(
          "only the top ",
          interactive_n,
          " variants are plotted interactively in addition to highlighted
          points, the rest of the plot consists of blocks of colour"
        )
      )
    } else {
      if (nrow(df) > 100000)
        stop(
          "a maximum of 100,000 points can be plotted interactively in ",
          "addition to highlighted points, specify interactive_n to plot ",
          "a block of colour below the top 100,000 points or a thin_thresh"
        )
      fig <- fig +
        geom_point_interactive(
          mapping = aes(colour = chr_cat, fill = chr_cat, tooltip = text),
          shape = 21, size = point_size, stroke = 0
        )
    }
  } else {
    if (!is.null(block_thresh)) {
      df_cut_off <- df %>%
        filter(highlight == 0)
      if (rank_pos == TRUE) {
        fig <- fig +
          geom_ribbon(
            mapping = aes(
              x = x, ymin = 0, ymax = -log10(block_thresh),
              colour = chr_cat
            ),
            data = filter(df, y >= -log10(!!block_thresh) & highlight == 0),
            size = point_size, outline.type = "full"
          )
      } else {
        df_cut_off <- df %>%
          filter(y <= -log10(!!block_thresh) & highlight == 0)
        df_cut_off_window <- tibble()
        for (i in unique(pull(df_cut_off, chr_int))) {
          df_cut_off_chr <- df_cut_off %>%
            filter(chr_int == !!i)
          df_cut_off_chr$window <- cut(
            df_cut_off_chr$x,
            seq(
              min(df_cut_off_chr$x),
              max(df_cut_off_chr$x),
              100000
            ),
            include.lowest = TRUE,
            labels = FALSE
          )
          df_cut_off_chr_window <- tibble(
            window = seq_along(
              zoo::rollmean(
                seq(
                  min(df_cut_off_chr$x),
                  max(df_cut_off_chr$x),
                  100000
                ),
                2
              )
            ),
            mid_point = zoo::rollmean(
              seq(
                min(df_cut_off_chr$x),
                max(df_cut_off_chr$x),
                100000
              ),
              2
            )
          )
          df_cut_off_chr_window <- df_cut_off_chr %>%
            group_by(window) %>%
            summarize(
              y_cut_off = max(y),
              .groups = "drop"
            ) %>%
            left_join(
              x = df_cut_off_chr_window,
              y = .,
              by = "window"
            ) %>%
            mutate(
              y_cut_off = if_else(is.na(y_cut_off), 0, y_cut_off),
              chr_cat = first(pull(df_cut_off_chr, chr_cat))
            )
          df_cut_off_window <- df_cut_off_window %>%
            bind_rows(df_cut_off_chr_window)
        }
        fig <- fig +
          geom_ribbon(
            mapping = aes(
              x = mid_point, ymin = 0, ymax = y_cut_off,
              colour = chr_cat, fill = chr_cat
            ),
            data = df_cut_off_window,
            size = point_size,
            outline.type = "full",
            inherit.aes = FALSE
          )
      }
      fig <- fig +
        geom_point(
          mapping = aes(colour = chr_cat, fill = chr_cat),
          data = filter(df, y >= -log10(!!block_thresh) & highlight == 0),
          shape = 21, size = point_size, stroke = 0
        )
    } else {
      fig <- fig +
        geom_point(
          mapping = aes(colour = chr_cat, fill = chr_cat),
          shape = 21, size = point_size, stroke = point_size / 30
        )
    }
  }

  ## Colours
  if (!is.null(colours)) {
    if (length(colours) == 2) {
      colours <- rep(colours, ceiling(length(unique(pull(df, chr_int))) / 2))
      colours <- colours[seq_along(unique(pull(df, chr_int)))]
    }
    names(colours) <- levels(df$chr_cat)
    fig <- fig +
      scale_colour_manual(values = colours, drop = FALSE) +
      scale_fill_manual(values = colours, drop = FALSE)
  }

  ## Highlights
  if (any(pull(df, highlight) != 0)) {
    fig <- fig +
      new_scale_fill() +
      new_scale_colour()
    if (interactive == TRUE) {
      fig <- fig +
        geom_point_interactive(
          mapping = aes(
            colour = factor(highlight, levels = sort(unique(highlight))),
            fill = factor(highlight, levels = sort(unique(highlight))),
            tooltip = text
          ),
          data = filter(df, highlight != 0 & highlight_shape == 0),
          shape = 21, size = point_size, stroke = 0
        ) +
        geom_point_interactive(
          mapping = aes(
            fill = factor(
              highlight,
              levels = sort(unique(highlight))
            ),
            tooltip = text
          ),
          data = filter(df, highlight != 0 & highlight_shape == 1),
          colour = "black", shape = 21, size = 1.1 * point_size
        ) +
        geom_point_interactive(
          mapping = aes(
            colour = factor(
              highlight,
              levels = sort(unique(highlight))
            ),
            fill = factor(
              highlight,
              levels = sort(unique(highlight))
            ),
            tooltip = text
          ),
          data = filter(df, highlight != 0 & highlight_shape == 2),
          shape = 22, size = 1.1 * point_size, stroke = 0
        ) +
        geom_point_interactive(
          mapping = aes(
            fill = factor(
              highlight,
              levels = sort(unique(highlight))
            ),
            tooltip = text
          ),
          data = filter(df, highlight != 0 & highlight_shape == 3),
          colour = "black", shape = 22, size = 1.1 * point_size
        ) +
        geom_point_interactive(
          mapping = aes(
            colour = factor(highlight, levels = sort(unique(highlight))),
            fill = factor(
              highlight,
              levels = sort(unique(highlight))
            ),
            tooltip = text
          ),
          data = filter(df, highlight != 0 & highlight_shape == 4),
          shape = 23, size = 1.1 * point_size, stroke = 0
        ) +
        geom_point_interactive(
          mapping = aes(
            fill = factor(
              highlight,
              levels = sort(unique(highlight))
            ),
            tooltip = text
          ),
          data = filter(df, highlight != 0 & highlight_shape == 5),
          colour = "black", shape = 23, size = 1.1 * point_size
        )
    } else {
      fig <- fig +
        geom_point(
          mapping = aes(
            colour = factor(highlight, levels = sort(unique(highlight))),
            fill = factor(highlight, levels = sort(unique(highlight)))
          ),
          data = filter(df, highlight != 0 & highlight_shape == 0),
          shape = 21, size = point_size, stroke = 0
        ) +
        geom_point(
          mapping = aes(
            fill = factor(highlight, levels = sort(unique(highlight)))
          ),
          data = filter(df, highlight != 0 & highlight_shape == 1),
          colour = "black", shape = 21, size = 1.05 * point_size
        ) +
        geom_point(
          mapping = aes(
            colour = factor(highlight, levels = sort(unique(highlight))),
            fill = factor(highlight, levels = sort(unique(highlight)))
          ),
          data = filter(df, highlight != 0 & highlight_shape == 2),
          shape = 22, size = 1.05 * point_size, stroke = 0
        ) +
        geom_point(
          mapping = aes(
            fill = factor(highlight, levels = sort(unique(highlight)))
          ),
          data = filter(df, highlight != 0 & highlight_shape == 3),
          colour = "black", shape = 22, size = 1.05 * point_size
        ) +
        geom_point(
          mapping = aes(
            colour = factor(highlight, levels = sort(unique(highlight))),
            fill = factor(highlight, levels = sort(unique(highlight)))
          ),
          data = filter(df, highlight != 0 & highlight_shape == 4),
          shape = 23, size = 1.05 * point_size, stroke = 0
        ) +
        geom_point(
          mapping = aes(
            fill = factor(highlight, levels = sort(unique(highlight)))
          ),
          data = filter(df, highlight != 0 & highlight_shape == 5),
          colour = "black", shape = 23, size = 1.05 * point_size
        )
    }
  }

  ### Highlight colours
  if (any(pull(df, highlight) != 0) && !is.null(highlight_colours)) {
    names(highlight_colours) <- sort(unique(df$highlight[df$highlight != 0]))
    if (any(pull(df, highlight_shape) %in% c(0, 2, 4))) {
      fig <- fig +
        scale_colour_manual(values = highlight_colours, drop = FALSE)
    }
    fig <- fig +
      scale_fill_manual(values = highlight_colours, drop = FALSE)
  }

  ## Threshold
  if (!is.null(thresh)) {
    for (i in seq_along(thresh)) {
      if (length(thresh_colours) != 1) {
        fig <- fig +
          geom_hline(
            yintercept = -log10(thresh[i]),
            color = thresh_colours[i], linetype = "dashed", size = thresh_size
          )
      } else {
        fig <- fig +
          geom_hline(
            yintercept = -log10(thresh[i]),
            color = thresh_colours, linetype = "dashed", size = thresh_size
          )
      }
    }
  }

  ## Axes breaks & limits
  fig <- fig +
    scale_x_continuous(
      breaks = df_plot$mid_point, labels = paste0(" ", df_plot$chr),
      limits = c(-1000, (max(df$x) + 1000)), expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, (max(df$y) + 0.1 * max(df$y))), expand = c(0.01, 0.01)
    )

  ## Axes labels
  fig <- fig +
    xlab("Chromosome") +
    ylab(expression("-log"["10"] * paste("(", italic("p"), ")")))

  ## Title
  if (!is.null(title)) {
    title_hjust <- 0
    if (title_center == TRUE)
      title_hjust <- 0.5
    fig <- fig +
      ggtitle(title) +
      theme(plot.title = element_text(size = title_size, hjust = title_hjust))
  }

  ## Theme
  fig <- fig +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = unit(c(0.5, 0.5, 0, 0), "cm"),
      legend.position = "none"
    )
  if (!is.null(axis_text_size)) {
    fig <- fig +
      theme(axis.text = element_text(size = axis_text_size))
  }
  if (!is.null(axis_title_size)) {
    fig <- fig +
      theme(axis.title = element_text(size = axis_title_size))
  }
  if (x_labels == FALSE) {
    fig <- fig +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  ## Labels
  if (label_top == TRUE) {
    if (label_box == FALSE) {
      fig <- fig +
        geom_text_repel(
          mapping = aes(label = label),
          data = filter(df, label != ""),
          size = label_size, segment.size = 0.5, segment.color = "grey",
          ylim = c(label_ylim, NA), point.padding = 0, nudge_y = label_nudge_y
        )
    } else {
      fig <- fig +
        geom_label_repel(
          mapping = aes(label = label),
          data = filter(df, label != ""),
          size = label_size, segment.size = 0.5, segment.color = "grey",
          ylim = c(label_ylim, NA), point.padding = 0, nudge_y = label_nudge_y
        )
    }
  }

  # Return
  return(fig)

}
