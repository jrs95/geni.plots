#' @title PheWAS plot
#'
#' @description `fig_phewas` creates a plot visualising results from
#'   phenome-wide association studies (PheWAS).
#'
#' @details This plotting function visualises results from phenome-wide
#'   association studies (PheWAS) in the form of a Manhattan style plot.
#'   Associations are grouped into phenotype categories. By default the
#'   results are truncated using a p-value cut-off of `1e-30`.
#'
#' @name fig_phewas
#'
#' @import dplyr
#'
#' @import ggplot2
#'
#' @import ggiraph
#'
#' @import ggrepel
#'
#' @param data a `data.frame` containing the association statistics for each
#'   phenotype with the following columns:
#'   \itemize{
#'     \item{\code{pvalue}} {
#'       the association p-value for each phenotype
#'     }
#'     \item{\code{sign}} {
#'       the direction of the association with the phenotype, where
#'       `0` = missing,
#'       `1` positive association,
#'       `-1` negative association
#'     }
#'     \item{\code{group}} {
#'       the phenotype group for each phenotype
#'     }
#'     \item{\code{label}} {
#'       the optional point labelling variable (e.g. phenotype name),
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
#' @param colours a `character` `vector` of colours corresponding to defined
#'   groups (default: `NULL`)
#'
#' @param interactive a `logical` value indicating whether the plot
#'   should be interactive (default: `FALSE`)
#'
#' @param thresh a `numeric` value providing the p-value threshold to be
#'   plotted (default: `NULL`)
#'
#' @param thresh_size a `numeric` value indicating the width of the lines
#'   indicating the p-value thresholds (default: `0.5`)
#'
#' @param trunc a `numeric` value representing the maximum p-value for which
#'   results are displayed (default: `1e-30`)
#'
#' @param point_size a `numeric` value indicating the size of each point
#'   (default: `2`)
#'
#' @param group_dist a `numeric` value indicating the gap between different
#'   groups (default: `0.05`)
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
#' @param label_n an `integer` value providing a limit on the number of top
#'   associations to label (default: `NULL`)
#'
#' @param label_size a `numeric` value indicating the size of each label
#'   (default: `3`)
#'
#' @param label_ymax a `numeric` value indicating the p-value threshold for the
#'   maximum y-axis value at which labels can be displayed (default: `1e-5`)
#'
#' @param label_box a `logical` value indicating whether labels should be
#'   surrounded by a box (default: `FALSE`)
#'
#' @param label_nudge_x a `numeric` value indicating the degree to which label
#'   placement on the x-axis should be adjusted (default: `0`)
#'
#' @param label_nudge_y a `numeric` value indicating the degree to which label
#'   placement on the y-axis should be adjusted (default: `0`)`
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
#' @param axis_text_angle a `numeric` value indicating the angle of the text
#'   on the x-axis of the plot (default: `-60`)`
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
#' @param limit_padding a `numeric` value indicating the relative distance of
#'   plotted points from x-axis extremes (default: `20`)
#'
#' @param plot_width a `numeric` value indicating the width of the PheWAS
#'   plot (default: `9`)
#'
#' @param plot_height a `numeric` value indicating the height of the PheWAS
#'   plot (default: `6`)
#'
#' @param girafe a `logical` value indicating whether an interactive plot
#'   should be turned into an interactive graphic using
#'   [girafe()][ggiraph::girafe()] (default = `TRUE`)
#'
#' @return `fig_phewas` returns a PheWAS plot for phenome-wide association
#'   studies.
#'
#' @examples
#' fig_phewas(
#'   data = geni.plots::geni_test_phewas,
#'   axis_text_angle = -85,
#'   axis_text_size = 8
#' )
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @export
#' @md
fig_phewas <- function(data, groups = NULL, colours = NULL, interactive = FALSE,
  thresh = 1e-5, thresh_size = 0.5, trunc = 1e-30, point_size = 2,
  group_dist = 0.05, x_labels = TRUE, label_top = TRUE, label_thresh = 1e-5,
  label_n = NULL, label_size = 3, label_ymax = 1e-5, label_box = FALSE,
  label_nudge_x = 0, label_nudge_y = 0, title = NULL, title_size = NULL,
  title_center = FALSE, axis_text_size = NULL, axis_text_angle = -60,
  axis_title_size = NULL, legend = FALSE, legend_title = "Group",
  legend_text_size = NULL, legend_title_size = NULL, legend_point_size = NULL,
  legend_spacing_size = NULL, limit_padding = 20, plot_width = 9,
  plot_height = 6, girafe = TRUE) {

  # Errors
  if (any(is.na(data)))
    stop("data cannot contain missing values")
  if (!all(c("pvalue", "sign", "group") %in% names(data)))
    stop(
      "the variables pvalue, sign & group are required to be in the ",
      "data data.frame"
    )
  if (!is.null(groups)) {
    if (any(is.na(groups)))
      stop("groups cannot contain missing values")
    if (any(duplicated(groups)))
      stop("groups cannot contain duplicates")
    if (all(!(as.character(data$group) %in% as.character(groups))))
      stop(
        "there are no values in the group variable in the data data.frame ",
        "that are in the groups variable"
      )
  }
  if (!is.null(colours)) {
    if (any(is.na(colours)))
      stop("colours cannot contain missing values")
    if (any(duplicated(colours)))
      stop("colours cannot contain duplicates")
    if (!is.null(groups) && length(groups) != length(colours))
      stop("the number of colours must equal the number of groups")
    if (
      is.null(groups) && !is.null(data$group) &&
        length(unique(data$group)) != length(colours)
    )
      stop("the number of colours must equal the number of groups")
  }
  if (interactive == TRUE && !("text" %in% names(data)))
    stop(
      "the text variable must be provided in the data data.frame ",
      "if the plot is specified to be interactive"
    )
  if (label_top == TRUE && !("label" %in% names(data)))
    stop(
      "the label variable must be provided in the data data.frame ",
      "if the plot is specified to have labels"
    )

  # Data
  data <- fig_phewas_data(
    data, groups, interactive, thresh, trunc, group_dist,
    label_top, label_thresh, label_n
  )

  # Plot
  fig <- fig_phewas_plot(
    data$df, data$df_plot, colours, interactive, thresh, thresh_size,
    point_size, x_labels, label_top, label_size, label_ymax, label_box,
    label_nudge_x, label_nudge_y, title, title_size, title_center,
    axis_text_size, axis_text_angle, axis_title_size, legend, legend_title,
    legend_text_size, legend_title_size, limit_padding
  )

  ## Legend
  if (!is.null(legend_point_size))
    fig <- fig +
      guides(
        fill = guide_legend(override.aes = list(size = legend_point_size))
      )
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

#' @title Figure PheWAS data
#'
#' @description `fig_phewas_data` creates the data.frames used to plot the
#'   PheWAS plot.
#'
#' @title fig_phewas_data
#'
#' @import dplyr
#'
#' @param data `data.frame` containing the columns:
#'   `pvalue`,
#'   `sign` (sign of association, 0 = missing, -1 = negative, 1 = positive),
#'   `group` (grouping variable),
#'   `label` (point labelling variable, if `label = ""` for a point then no
#'     label is presented for that point),
#'   `text` (hover text variable, if `text = ""` for a point then no hover is
#'     presented for that point)
#'
#' @param groups groups of the grouping variables; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param thresh p-value threshold; default: `1e-5`
#'
#' @param trunc truncate the association results to this values;
#'   default: `1e-30`
#'
#' @param group_dist group distance measure to multiply the number of rows by;
#'   default: `0.05`
#'
#' @param label_top label the top associated points; default: `TRUE`
#'
#' @param label_thresh label all points with p-value less than this threshold;
#'   default: `1e-5`
#'
#' @param label_n label the top n points; default: `NULL`
#'
#' @return `fig_phewas_data` returns a `list` with the `data.frames` used to
#'   plot the PheWAS plot.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_phewas_data <- function(data, groups = NULL, interactive = FALSE,
  thresh = 1e-5, trunc = 1e-30, group_dist = 0.05, label_top = TRUE,
  label_thresh = 1e-5, label_n = NULL) {

  # Dataset
  data <- data %>%
    as_tibble()

  # Removals
  data <- data %>%
    filter(!is.na(as.numeric(pvalue)))
  if (!is.null(groups))
    data <- data %>%
      filter(as.character(group) %in% !!as.character(groups))

  # Variables

  ## mlog10p
  data <- data %>%
    mutate(y = round(-log10(as.numeric(pvalue)), 4)) %>%
    relocate(y) %>%
    select(-pvalue)
  if (!is.null(trunc))
    data <- data %>%
      mutate(y = if_else(y > -log10(!!trunc), round(-log10(!!trunc), 4), y, y))

  ## Groups
  data <- data %>%
    mutate(group = as.character(group))
  if (!is.null(groups)) {
    groups <- as.character(groups)
    data <- data %>%
      mutate(
        group_cat = factor(group, levels = !!groups),
        group_int = as.integer(
          as.character(
            factor(
              group,
              levels = groups[groups %in% group],
              labels = 0:(length(groups[groups %in% group]) - 1)
            )
          )
        )
      )
  } else {
    data <- data %>%
      mutate(
        group_cat = factor(group, levels = sort(unique(group))),
        group_int = as.integer(
          as.character(
            factor(
              group,
              levels = sort(unique(group)),
              labels = 0:(length(sort(unique(group))) - 1)
            )
          )
        )
      )
  }

  ## Shape
  data <- data %>%
    mutate(
      sign = as.integer(sign),
      sign = if_else(is.na(sign), as.integer(0), sign),
      shape = as.integer(0),
      shape = if_else(y >= -log10(!!thresh) & sign == -1, as.integer(1), shape),
      shape = if_else(y >= -log10(!!thresh) & sign == 1, as.integer(2), shape),
      shape = factor(shape, levels = c(0, 1, 2))
    ) %>%
    select(-sign)

  ## Labels & hover text
  if (label_top == TRUE)
    data <- data %>%
      mutate(label = as.character(label))
  if (interactive == TRUE)
    data <- data %>%
      mutate(text = as.character(text))

  ## Position
  data <- data %>%
    arrange(group_int, stats::runif(nrow(.))) %>%
    mutate(
      x = seq_len(nrow(.)) + (group_dist * nrow(.)) * group_int
    ) %>%
    relocate(x, .before = y)

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
          mutate(label = if_else(row_number() > !!label_n, "", label)) %>%
          arrange(x)
      }
    }
  }

  # Datasets
  data_plot <- data %>%
    group_by(group_int) %>%
    summarize(
      mid_point = (min(x) + (max(x) - min(x)) / 2),
      group = first(group),
      .groups = "drop"
    )

  # Output
  out <- list(df = data, df_plot = data_plot)

  # Return
  return(out)

}

#' @title Figure PheWAS plot
#'
#' @description `fig_phewas_plot` plots a PheWAS plot for phenome-wide
#'   association studies.
#'
#' @name fig_phewas_plot
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
#'   `shape` (sign shape; `0` = missing, `-1` = negative, `1` = positive),
#'   `group` (grouping variable), `group_cat` (grouping category),
#'   `group_int` (grouping plotting integer category),
#'   `label` (point labels, if `label = ""` for a point then no label is
#'     presented for that point),
#'   `text` (hover text, if `text = ""` for a point then no hover is
#'     presented for that point)
#'
#' @param df_plot `data.frame` containing the x-axis labelling:
#'   `group_int` (grouping plotting integer category),
#'   `mid_point` (x-axis mid-point),
#'   `group` (grouping variable)
#'
#' @param colours colours of the groups; default: `NULL`
#'
#' @param interactive make the plot interactive; default: `FALSE`
#'
#' @param thresh p-value threshold; default: `1e-5`
#'
#' @param thresh_size threshold size; default: `0.5`
#'
#' @param point_size point size; default: `2`
#'
#' @param x_labels label the x-axis; default: `TRUE`
#'
#' @param label_top label the top associated points; default: `TRUE`
#'
#' @param label_size label text size; default: `3`
#'
#' @param label_ymax label p-value position threshold; default: `1e-5`
#'
#' @param label_box put the label in a box; default: `FALSE`
#'
#' @param label_nudge_x nudge x-axis of labels; default: `0`
#'
#' @param label_nudge_y nudge y-axis of labels; default: `0`
#'
#' @param title figure title; default: `NULL`
#'
#' @param title_size title text size; default: `NULL`
#'
#' @param title_center center title; default: `FALSE`
#'
#' @param axis_text_size axis text size; default: `NULL`
#'
#' @param axis_text_angle angle of x-axis text; default: `-60`
#'
#' @param axis_title_size axis title size; default: `NULL`
#'
#' @param legend display legend for groups; default: `FALSE`
#'
#' @param legend_title legend title; default: `"Group"`
#'
#' @param legend_text_size legend text size; default: `NULL`
#'
#' @param legend_title_size legend title size; default: `NULL`
#'
#' @param limit_padding limit padding; default: `20`
#'
#' @return `fig_phewas_plot` returns a PheWAS plot for phenome-wide
#'   association studies.
#'
#' @author James Staley <jrstaley95@gmail.com>
#'
#' @noRd
#' @md
fig_phewas_plot <- function(df, df_plot, colours, interactive = FALSE,
  thresh = 1e-5, thresh_size = 0.5, point_size = 2, x_labels = TRUE,
  label_top = TRUE, label_size = 3, label_ymax = 1e-5, label_box = FALSE,
  label_nudge_x = 0, label_nudge_y = 0, title = NULL, title_size = NULL,
  title_center = FALSE, axis_text_size = NULL, axis_text_angle = -60,
  axis_title_size = NULL, legend = FALSE, legend_title = "Group",
  legend_text_size = NULL, legend_title_size = NULL, limit_padding = 20) {

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

  ## Threshold
  if (max(df$y) > -log10(thresh))
    fig <- fig +
      geom_hline(
        yintercept = -log10(thresh),
        color = "grey50",
        linetype = "dashed",
        size = thresh_size
      )

  ## Points
  if (interactive == TRUE) {
    fig <- fig +
      geom_point_interactive(
        mapping = aes(fill = group_cat, shape = shape, tooltip = text),
        size = point_size, stroke = point_size / 30
      )
  } else {
    fig <- fig +
      geom_point(
        mapping = aes(fill = group_cat, shape = shape),
        size = point_size, stroke = point_size / 30
      )
  }

  ## Colours
  if (!is.null(colours)) {
    names(colours) <- levels(df$group_cat)
    fig <- fig +
      scale_fill_manual(values = colours, drop = FALSE)
  }

  ## Shapes for directions
  ### 21 = circle, 25 = downwards triangle, 24 = upwards triangle
  shapes <- c(21, 25, 24)
  names(shapes) <- c(0, 1, 2)
  fig <- fig +
    scale_shape_manual(values = shapes, guide = "none")

  ## Axes breaks & limits
  fig <- fig +
    scale_x_continuous(
      breaks = df_plot$mid_point,
      labels = paste0(" ", df_plot$group),
      limits = c(-limit_padding, (max(df$x) + limit_padding)),
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, (max(df$y) + 0.1 * max(df$y))),
      expand = c(0.01, 0.01)
    )

  ## Axes labels
  fig <- fig +
    xlab(NULL) +
    ylab(expression("-log"["10"] * paste("(", italic("p"), ")")))

  ## Title
  if (!is.null(title)) {
    title_hjust <- 0
    if (title_center == TRUE) {
      title_hjust <- 0.5
    }
    fig <- fig +
      ggtitle(title) +
      theme(
        plot.title = element_text(size = title_size, hjust = title_hjust)
      )
  }

  ## Theme
  fig <- fig +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(
        hjust = 0, vjust = 0.25, angle = axis_text_angle
      ),
      plot.margin = unit(c(1, 1, 0, 0), "cm")
    )

  ### Axes
  if (!is.null(axis_text_size))
    fig <- fig +
      theme(axis.text = element_text(size = axis_text_size))
  if (!is.null(axis_title_size))
    fig <- fig +
      theme(axis.title.y = element_text(size = axis_title_size))
  if (x_labels == FALSE)
    fig <- fig +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )

  ### Legend
  if (legend == TRUE) {
    fig <- fig +
      labs(fill = legend_title) +
      theme(
        legend.text = element_text(size = legend_text_size),
        legend.title = element_text(size = legend_title_size)
      ) +
      guides(fill = guide_legend(override.aes = list(shape = 21)))
  } else {
    fig <- fig +
      theme(legend.position = "none")
  }

  ## Labels
  if (label_top == TRUE) {
    if (label_box == FALSE) {
      fig <- fig +
        geom_text_repel(
          data = filter(df, label != ""),
          size = label_size, segment.size = 0.5, segment.color = "grey",
          ylim = c(-log10(label_ymax), NA), point.padding = 0,
          nudge_x = label_nudge_x, nudge_y = label_nudge_y
        )
    } else {
      fig <- fig +
        geom_label_repel(
          data = filter(df, label != ""),
          size = label_size, segment.size = 0.5, segment.color = "grey",
          ylim = c(-log10(label_ymax), NA), point.padding = 0,
          nudge_x = label_nudge_x, nudge_y = label_nudge_y
        )
    }
  }

  # Return
  return(fig)

}
