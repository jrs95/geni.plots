#' @title QQ plot
#'
#' @description `fig_qq` creates a quantile-quantile (QQ) plot.
#'
#' @details This plotting function plots a quantile-quantile plot of
#'   -log10(p-values). Observations can be divided into groups and can include
#'   corresponding confidence intervals. This plot is based on the
#'   [QQ plot](https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_QQ_Plots_in_R)
#'   by Matthew Flickinger.
#'
#' @name fig_qq
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
#' @param pvalues the association p-value for each marker
#'   (default: `NULL`)
#'
#' @param group a `character` `vector` determining the observation group for
#'   each marker (default: `NULL`)
#'
#' @param data a `data.frame` containing the association statistics for each
#'   marker with the following columns:
#'   \itemize{
#'     \item{\code{pvalue}} {
#'       the association p-value for each marker
#'     }
#'     \item{\code{group}} {
#'       the optional grouping variable for each marker
#'     }
#'     \item{\code{label}} {
#'       the optional point labelling variable (e.g. genomic marker),
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
#' @param groups a `character` `vector` of groups describing the grouping
#'   variable in `data` (default: `NULL`)
#'
#' @param colours a `character` `vector` of colours corresponding to
#'   defined groups (default: `NULL`)
#'
#' @param interactive a `logical` value indicating whether the plot
#'   should be interactive (default: `FALSE`)
#'
#' @param thresh a `numeric` value providing the p-value threshold to be
#'   plotted (default: `NULL`)
#'
#' @param sample a `logical` value indicating whether a random subset of
#'   p-values above the plotting threshold should be plotted, the number of
#'   which is controlled by `sample_prop` (default = `FALSE`)
#'
#' @param sample_thresh a `numeric` value indicating the p-value threshold
#'   defining the sample from which additional observations are selected
#'   (default: `0.1`)
#'
#' @param sample_prop a `numeric` value indicating the proportion of
#'   sampled observations to be plotted (default: `0.1`)
#'
#' @param ci a `logical` value indicating whether confidence intervals
#'   should be displayed (default: `TRUE`)
#'
#' @param ci_alpha a `numeric` value providing the threshold defining the
#'   plotted confidence interval (default: `0.05`)
#'
#' @param ci_print a `logical` value indicating whether the proportion of points
#'   contained within the confidence interval band should be printed
#'   (default: `FALSE`)
#'
#' @param inf_factor a `logical` value indicating whether the inflation factor
#'   should be added to the plot (default: `FALSE`)
#'
#' @param point_size a `numeric` value indicating the size of each point
#'   (default: `3`)
#'
#' @param label_top a `logical` value whether the top associated points should
#'   be labelled (default: `TRUE`)
#'
#' @param label_thresh a `numeric` value providing a p-value threshold for
#'   labelling points (default: `1e-5`)
#'
#' @param label_n an `integer` value providing a limit on the number of top
#'   associations to label (default: `NULL`)
#'
#' @param label_size a `numeric` value indicating the size of each label
#'   (default: `3`)
#'
#' @param label_xlim a `numeric` value indicating maximum x-axis value at
#'   which labels can be displayed (default: `NULL`)
#'
#' @param label_box a `logical` value indicating whether labels should be
#'   surrounded by a box (default: `FALSE`)
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
#'   text for the plot (default: `NULL`)
#'
#' @param axis_title_size a `numeric` value indicating the size of the axis
#'   title text for the plot (default: `NULL`)
#'
#' @param legend a `logical`  value indicating whether a legend corresponding
#'   to the displayed groups should be included (default: `FALSE`)
#'
#' @param legend_title a `character` string providing a title for the legend
#'   (default: `"Group"`)
#'
#' @param legend_text_size a `numeric` value indicating the size of the legend
#'   text (default: `NULL`)
#'
#' @param legend_title_size a `numeric` value indicating the size of the legend
#'   title (default: `NULL`)
#'
#' @param legend_point_size a `numeric` value indicating the size of each point
#'   within the legend (default: `NULL`)
#'
#' @param legend_spacing_size a `numeric` value indicating spacing of points
#'   present in the legend (default: `NULL`)
#'
#' @param ymax a `numeric` value defining the maximum value of the y-axis
#'   (default: `NULL`)
#'
#' @param plot_width a `numeric` value indicating the width of the plot
#'   (default: `6`)
#'
#' @param plot_height a `numeric` value indicating the height of the plot
#'   (default: `6`)
#'
#' @param girafe a `logical` value indicating whether an interactive plot
#'   should be turned into an interactive graphic using
#'   [girafe()][ggiraph::girafe()] (default = `TRUE`)
#'
#' @return `fig_qq` returns a quantile-quantile plot.
#'
#' @examples
#' fig_qq(
#'   pvalues = geni.plots::geni_test_phewas$pvalue
#' )
#' fig_qq(
#'   data = geni.plots::geni_test_phewas %>%
#'     select(pvalue, group, label, text),
#'   interactive = TRUE
#' )
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
fig_qq <- function(pvalues = NULL, group = NULL, data = NULL, groups = NULL,
  colours = NULL, interactive = FALSE, thresh = NULL, sample = FALSE,
  sample_thresh = 0.1, sample_prop = 0.1, ci = TRUE, ci_alpha = 0.05,
  ci_print = FALSE, inf_factor = FALSE, point_size = 3, label_top = FALSE,
  label_thresh = 1e-5, label_n = 10, label_size = 3.25, label_xlim = NULL,
  label_box = FALSE, title = NULL, title_size = NULL, title_center = FALSE,
  axis_text_size = NULL, axis_title_size = NULL, legend = TRUE,
  legend_title = "Group", legend_text_size = NULL, legend_title_size = NULL,
  legend_point_size = NULL, legend_spacing_size = NULL, ymax = NULL,
  plot_width = 6, plot_height = 6, girafe = TRUE) {

  # Errors
  if (is.null(pvalues) && is.null(data))
    stop("pvalues and data both cannot be missing")
  if (!is.null(pvalues)) {
    if (!is.vector(pvalues))
      stop("pvalues has to be a vector")
    if (any(is.na(as.numeric(pvalues))))
      stop("pvalues cannot contain missing values")
  }
  if (!is.null(group)) {
    if (any(is.na(group)))
      stop("groups cannot contain missing values")
  }
  if (!is.null(data)) {
    if (any(is.na(data)))
      stop("data cannot contain missing values")
    if (!("pvalue" %in% names(data)))
      stop("pvalue variable must be included in data")
    if (any(is.na(as.numeric(data$pvalue))))
      stop(
        "the pvalue variable in data cannot contain ",
        "missing values"
      )
  }
  if (!is.null(groups)) {
    if (any(is.na(groups)))
      stop("groups cannot contain missing values")
    if (any(duplicated(groups)))
      stop("groups cannot contain duplicates")
    if (!is.null(data)) {
      if (!("group" %in% names(data)))
        stop(
          "the group variable has to be included in data ",
          "if groups is specified"
        )
      if (all(!(as.character(data$group) %in% as.character(groups))))
        stop(
          "there are no values in the group variable in data ",
          "that are in the groups variable"
        )
    } else {
      if (is.null(group))
        stop("the group variable has to be specified if groups is specified")
      if (all(!(as.character(group) %in% as.character(groups))))
        stop(
          "there are no values in the group variable in data ",
          "that are in the groups variable"
        )
    }
  }
  if (!is.null(colours)) {
    if (any(is.na(colours)))
      stop("colours cannot contain missing values")
    if (any(duplicated(colours)))
      stop("colours cannot contain duplicates")
    if (!is.null(groups) && length(groups) != length(colours))
      stop("the number of colours must equal the number of groups")
    if (
      is.null(groups) && !is.null(group) &&
        length(unique(group)) != length(colours)
    )
      stop("the number of colours must equal the number of groups")
    if (
      is.null(groups) && !is.null(data$group) &&
        length(unique(data$group)) != length(colours)
    )
      stop("the number of colours must equal the number of groups")
  }
  if (interactive == TRUE && is.null(data))
    stop(
      "data must be provided if the plot is specified to be ",
      "interactive"
    )
  if (interactive == TRUE && !("text" %in% names(data)))
    stop(
      "the text variable must be provided in data if the plot ",
      "is specified to be interactive"
    )
  if (interactive == TRUE && nrow(data) > 100000)
    stop(
      "an interactive QQ plot can only be plotted for a maximum of 100,000 ",
      "points"
    )
  if (label_top == TRUE && is.null(data))
    stop(
      "data must be provided if the plot is specified to have ",
      "labels"
    )
  if (label_top == TRUE && !("label" %in% names(data)))
    stop(
      "the label variable must be provided in data if the plot ",
      "is specified to be interactive"
    )

  # Data
  df <- fig_qq_data(
    pvalues, group, data, groups, interactive, sample, sample_thresh,
    sample_prop, ci, ci_alpha, ci_print, inf_factor, label_top,
    label_thresh, label_n
  )
  groups_exist <- "group" %in% names(df$df)

  # Plot
  fig <- fig_qq_plot(
    df, groups_exist, colours, interactive, thresh, ci, inf_factor, point_size,
    label_top, label_size, label_xlim, label_box, title, title_size,
    title_center, axis_text_size, axis_title_size, legend, legend_title,
    legend_text_size, legend_title_size, ymax
  )

  ## Legend
  if (!is.null(legend_point_size))
    fig <- fig +
      guides(fill = guide_legend(override.aes = list(size = legend_point_size)))
  if (!is.null(legend_spacing_size))
    fig <- fig +
      theme(legend.key.size = unit(legend_spacing_size, "lines"))

  ## Interactive
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

#' @title Figure QQ data
#'
#' @description `fig_qq_data` creates the `data.frame` used to plot the
#'   quantile-quantile plot.
#'
#' @name fig_qq_data
#'
#' @import dplyr
#'
#' @import tidyr
#'
#' @param pvalues association p-values; default: `NULL`
#'
#' @param group grouping variable; default: `NULL`
#'
#' @param data data.frame containing the columns:
#'   `pvalue`,
#'   `group` (grouping variable),
#'   `label` (point labelling variable, if `label = ""` for a point then
#'     no label is presented for that point),
#'   `text` (hover text variable, if `text = ""` for a point then no hover
#'     is presented for that point); default: `NULL`
#'
#' @param groups groups of the grouping variables; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param sample sample the p-values above the threshold for plotting with
#'   `sample_prop` proportion; default: `FALSE`
#'
#' @param sample_thresh sample threshold; default: `0.1`
#'
#' @param sample_prop sample proportion; default: `0.1`
#'
#' @param ci display confidence interval band; default: `TRUE`
#'
#' @param ci_alpha alpha for confidence interval; default: `0.05`
#'
#' @param ci_print print proportion of points contained within the confidence
#'   interval band; default: `FALSE`
#'
#' @param inf_factor print inflation factor; default: `FALSE`
#'
#' @param label_top label the top associated points; default: `FALSE`
#'
#' @param label_thresh label all points with p-value less than this threshold;
#' default: `1e-5`
#'
#' @param label_n label the top n points; default: `10`
#'
#' @return `fig_qq_data` returns the `data.frame` used to plot the
#'   quantile-quantile plot.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_qq_data <- function(pvalues = NULL, group = NULL, data = NULL,
  groups = NULL, interactive = FALSE, sample = FALSE, sample_thresh = 0.1,
  sample_prop = 0.1, ci = TRUE, ci_alpha = 0.05, ci_print = FALSE,
  inf_factor = FALSE, label_top = FALSE, label_thresh = 1e-5, label_n = 10) {

  # Process data
  if (is.null(data)) {
    data <- tibble(
      x = round(
        -log10(
          (rank(as.numeric(pvalues), ties.method = "first") - 0.5) /
            length(pvalues)
        ),
        4
      ),
      y = round(-log10(as.numeric(pvalues)), 4)
    )
    if (!is.null(group)) {
      if (!is.null(groups)) {
        groups <- unique(as.character(groups))
        data <- data %>%
          mutate(
            group = factor(!!as.character(group), levels = !!groups)
          )
      } else {
        data <- data %>%
          mutate(
            group = factor(
              !!as.character(group),
              levels = !!sort(unique(as.character(group)))
            )
          )
      }
    }
  } else {
    if (("group" %in% names(data)) && !is.null(groups))
      data <- data %>%
        filter(as.character(group) %in% !!as.character(groups))
    data <- data %>%
      as_tibble() %>%
      mutate(
        x = round(
          -log10(
            (rank(as.numeric(pvalue), ties.method = "first") - 0.5) /
              length(pvalue)
          ),
          4
        ),
        y = round(-log10(as.numeric(pvalue)), 4)
      ) %>%
      select(-pvalue) %>%
      relocate(x, y)
    if ("group" %in% names(data)) {
      if (!is.null(groups)) {
        groups <- unique(as.character(groups))
        data <- data %>%
          mutate(
            group = factor(as.character(group), levels = !!groups)
          )
      } else {
        data <- data %>%
          mutate(
            group = factor(
              as.character(group),
              levels = sort(unique(as.character(group)))
            )
          )
      }
    }
    if (label_top == TRUE) {
      data <- data %>%
        mutate(label = as.character(label))
    }
    if (interactive == TRUE) {
      data <- data %>%
        mutate(text = as.character(text))
    }
  }

  # Inflation factor
  if (inf_factor == TRUE) {
    inf_factor <- round(
      median(qchisq(10^(-pull(df, y)), 1, lower.tail = FALSE)) /
        qchisq(0.5, 1),
      4
    )
    print(paste0("Inflation factor: ", inf_factor))
  }

  # Confidence intervals
  if (ci == TRUE) {
    data <- data %>%
      arrange(desc(x))
    mpts <- tibble(
      x1 = round(-log10(ppoints(nrow(data), a = 0.5)), 4),
      z1 = round(
        -log10(qbeta((1 - ci_alpha / 2), 1:nrow(data), nrow(data):1)),
        4
      ),
      z2 = round(-log10(qbeta(ci_alpha / 2, 1:nrow(data), nrow(data):1)), 4)
    )
    if (nrow(data) != nrow(mpts))
      stop(
        "confidence interval data.frame is not the same size as the p-value ",
        "data.frame"
      )
    if (any(data$x != mpts$x1))
      stop("confidence interval x-axis different from p-value x-axis")
    data <- bind_cols(data, mpts)
    if (ci_print == TRUE) {
      print(
        paste0("Proportion of points within the CI: ",
          round(sum((data$y >= data$z1) & (data$y <= data$z2)) / nrow(data), 5))
      )
    }
    if (nrow(data) > 1000000) {
      ids <- seq_len(nrow(data))
      ids <- ids[
        -c(
          head(ids, n = floor(0.05 * length(ids))),
          tail(ids, n = floor(0.05 * length(ids)))
        )
      ]
      ids <- sample(ids, size = floor(0.9 * length(ids)))
      data$x1[ids] <- NA
      data$z1[ids] <- NA
      data$z2[ids] <- NA
    }
  }

  # Labels
  if (label_top == TRUE) {
    if (
      !is.null(label_thresh) && length(label_thresh) == 1 &&
        is.numeric(label_thresh)
    )
      data <- data %>%
        mutate(label = if_else(y < -log10(!!label_thresh), "", label))
    if (!is.null(label_n)) {
      if (sum(data$label != "") > label_n) {
        data <- data %>%
          arrange(desc(y)) %>%
          mutate(label = if_else(row_number(.) > !!label_n, "", label)) %>%
          arrange(x)
      }
    }
  }

  # Sample
  if (sample == TRUE) {
    df1 <- data %>%
      filter(y >= -log10(!!sample_thresh))
    df2 <- data %>%
      filter(y < -log10(!!sample_thresh)) %>%
      slice(sample(nrow(df2), ceiling(sample_prop * nrow(df2))))
    data <- bind_rows(df1, df2)
  }

  # Shuffle
  data <- data %>%
    slice(sample(nrow(.), replace = FALSE))

  # Output
  output <- list(df = data)
  if (inf_factor == TRUE) {
    output$inf_factor <- inf_factor
  }

  # Return
  return(output)

}

#' @title Figue QQ plot
#'
#' @description `fig_qq_plot` plots a quantile-quantile plot.
#'
#' @title fig_qq_plot
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
#' @param df `data.frame` containing the plotting variables:
#'   `x` (x-axis),
#'   `y` (y-axis),
#'   `group` (grouping variable),
#'   `label` (point labels,
#'     if `label = ""` for a point then no label is presented for that point),
#'   `text`` (hover text,
#'     if `text = ""` for a point then no hover is presented for that point),
#'   `x1` (x-axis for confidence interval),
#'   `z1` (y-axis lower confidence interval),
#'   `z2` (y-axis upper confidence interval)
#'
#' @param groups_exist display groups; default: `FALSE`
#'
#' @param colours colours of the groups; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param thresh p-value threshold; default: `NULL`
#'
#' @param ci display confidence interval band; default: `TRUE`
#'
#' @param inf_factor print inflation factor and add to plot; default: `FALSE`
#'
#' @param point_size point size; default: `3`
#'
#' @param label_top label the top associated points; default: `FALSE`
#'
#' @param label_size label text size; default: `3.25`
#'
#' @param label_xlim x-axis label limit maximum; default: `NULL`
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
#' @param legend display legend for groups; default: `TRUE`
#'
#' @param legend_title legend title; default: `"Group"`
#'
#' @param legend_text_size legend text size; default: `NULL`
#'
#' @param legend_title_size legend title size; default: `NULL`
#'
#' @param ymax maximum of the y-axis; default: `NULL`
#'
#' @return `fig_qq_plot` returns a quantile-quantile plot.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_qq_plot <- function(df, groups_exist = FALSE, colours = NULL,
  interactive = FALSE, thresh = NULL, ci = TRUE, inf_factor = FALSE,
  point_size = 3, label_top = FALSE, label_size = 3.25, label_xlim = NULL,
  label_box = FALSE, title = NULL, title_size = NULL, title_center = FALSE,
  axis_text_size = NULL, axis_title_size = NULL, legend = TRUE,
  legend_title = "Group", legend_text_size = NULL, legend_title_size = NULL,
  ymax = NULL) {

  # Dataset
  if (inf_factor == TRUE)
    inf_factor <- df$inf_factor
  df <- df$df

  # Plot

  ## Base plot
  if (label_top == TRUE) {
    fig <- ggplot(
      data = df,
      mapping = aes(x = x, y = y, label = label)
    )
  } else {
    fig <- ggplot(
      data = df,
      mapping = aes(x = x, y = y)
    )
  }

  ## Add confidence intervals
  if (ci == TRUE) {
    fig <- fig +
      geom_ribbon(
        mapping = aes(x = x1, ymin = z1, ymax = z2),
        colour = "lightblue1", fill = "lightblue1"
      ) +
      theme(legend.position = "none")
  }

  ## Identity line
  fig <- fig +
    geom_abline(slope = 1, intercept = 0, colour = "grey50")

  ## Threshold
  if (!is.null(thresh)) {
    fig <- fig +
      geom_hline(
        yintercept = -log10(thresh),
        color = "grey50", linetype = "dashed"
      )
  }

  ## Points
  if (interactive == TRUE) {
    if (groups_exist == TRUE) {
      fig <- fig +
        geom_point_interactive(
          mapping = aes(fill = group, tooltip = text),
          size = point_size, stroke = point_size / 50, shape = 21
        )
    } else {
      fig <- fig +
        geom_point_interactive(
          mapping = aes(tooltip = text), fill = "red",
          size = point_size, stroke = point_size / 50, shape = 21
        )
    }
  } else {
    if (groups_exist == TRUE) {
      fig <- fig +
        geom_point(
          mapping = aes(fill = group),
          size = point_size, stroke = point_size / 50, shape = 21
        )
    } else {
      fig <- fig +
        geom_point(
          fill = "red",
          size = point_size, stroke = point_size / 50, shape = 21
        )
    }
  }

  ## Colours
  if (!is.null(colours) && groups_exist == TRUE) {
    names(colours) <- levels(df$group)
    fig <- fig +
      scale_fill_manual(values = colours, drop = FALSE)
  }

  ## Axes labels
  fig <- fig +
    xlab(expression("Expected -log"["10"] * paste("(", italic("p"), ")"))) +
    ylab(expression("Observed -log"["10"] * paste("(", italic("p"), ")")))

  ## Y-axis limit
  if (!is.null(ymax)) {
    fig <- fig + ylim(NA, ymax)
  }

  ## Title
  if (!is.null(title)) {
    title_hjust <- 0
    if (title_center == TRUE) {
      title_hjust <- 0.5
    }
    fig <- fig +
      ggtitle(title) +
      theme(
        plot.title = element_text(
          size = title_size, face = "bold", hjust = title_hjust
        )
      )
  }

  ## Theme
  fig <- fig +
    theme_bw() +
    theme(
      axis.text = element_text(size = axis_text_size),
      axis.title = element_text(size = axis_title_size),
      plot.margin = unit(c(0.5, 0.5, 0, 0), "cm")
    )

  ### Legend
  if (legend == TRUE) {
    fig <- fig +
      labs(fill = legend_title) +
      theme(
        legend.text = element_text(size = legend_text_size),
        legend.title = element_text(size = legend_title_size)
      )
  } else {
    fig <- fig +
      theme(legend.position = "none")
  }

  ## Labels
  if (label_top == TRUE && any(df$label != "")) {
    if (label_box == FALSE) {
      if (!is.null(label_xlim)) {
        fig <- fig +
          geom_text_repel(
            data = filter(df, label != ""),
            size = label_size, segment.size = 0.5, segment.color = "grey",
            xlim = c(NA, label_xlim),
            nudge_x = (-1 * 0.1 * max(df$x)), nudge_y = (0.1 * max(df$y))
          )
      } else {
        fig <- fig +
          geom_text_repel(
            data = filter(df, label != ""), size = label_size,
            segment.size = 0.5, segment.color = "grey",
            nudge_x = (-1 * 0.1 * max(df$x)), nudge_y = (0.1 * max(df$y))
          )
      }
    } else {
      if (!is.null(label_xlim)) {
        fig <- fig +
          geom_label_repel(
            data = filter(df, label != ""),
            size = label_size, segment.size = 0.5, segment.color = "grey",
            xlim = c(NA, label_xlim),
            nudge_x = (-1 * 0.1 * max(df$x)), nudge_y = (0.1 * max(df$y))
          )
      } else {
        fig <- fig +
          geom_label_repel(
            data = filter(df, label != ""), size = label_size,
            segment.size = 0.5, segment.color = "grey",
            nudge_x = (-1 * 0.1 * max(df$x)), nudge_y = (0.1 * max(df$y))
          )
      }
    }
  }

  ## Inflation factor
  if (inf_factor == TRUE) {
    if (is.null(ymax) && ci == TRUE) {
      ymax <- max(df$z2, na.rm = TRUE)
    }
    if (is.null(ymax) && ci == FALSE) {
      ymax <- max(df$y, na.rm = TRUE)
    }
    df_inf <- tibble(
      x = (min(df$x, na.rm = TRUE) + 0.05 * min(df$x, na.rm = TRUE)),
      y = (ymax - 0.05 * ymax)
    )
    fig <- fig +
      geom_label(
        mapping = aes(x = x, y = y),
        data = df_inf,
        label = deparse(bquote(lambda[inf] ~ "=" ~ .(inf_factor))),
        parse = TRUE, hjust = 0, size = label_size
      )
  }

  # Return
  return(fig)

}
