# ===============================================================================
# SHARED DAG SPECIFICATIONS AND PLOTTING HELPERS
# ===============================================================================

DAG_SPEC <- 'dag {
bb="0,0,1,1"
Age [pos="0.100,0.100"]
Sex [pos="0.280,0.100"]
U [latent,pos="0.750,0.080"]
TMB_BRAF [pos="0.200,0.380"]
CRP [pos="0.720,0.320"]
T [exposure,pos="0.050,0.680"]
TxTMB [pos="0.350,0.560"]
TxCRP [pos="0.760,0.560"]
TLR [pos="0.560,0.560"]
PFS [outcome,pos="0.720,0.760"]
OS [outcome,pos="0.900,0.920"]

Age -> TMB_BRAF
Age -> OS
Sex -> TMB_BRAF
U -> TMB_BRAF
U -> CRP
U -> PFS
U -> OS
TMB_BRAF -> PFS
TMB_BRAF -> OS
TMB_BRAF -> TLR
TMB_BRAF -> TxTMB
CRP -> PFS
CRP -> OS
CRP -> TLR
CRP -> TxCRP
T -> PFS
T -> OS
T -> TLR
T -> TxTMB
T -> TxCRP
TxTMB -> PFS
TxTMB -> OS
TxCRP -> PFS
TxCRP -> OS
TLR -> PFS
PFS -> OS
}'

DAG_SENS_SPEC <- 'dag {
bb="0,0,1,1"
Age [pos="0.100,0.100"]
Sex [pos="0.280,0.100"]
U [latent,pos="0.750,0.080"]
TMB_BRAF [pos="0.200,0.380"]
CRP [pos="0.720,0.320"]
T [exposure,pos="0.050,0.680"]
TxTMB [pos="0.350,0.560"]
TxCRP [pos="0.760,0.560"]
TLR [pos="0.560,0.560"]
PFS [outcome,pos="0.720,0.760"]
OS [outcome,pos="0.900,0.920"]

Age -> OS
U -> TMB_BRAF
U -> CRP
U -> PFS
U -> OS
TMB_BRAF -> PFS
TMB_BRAF -> OS
TMB_BRAF -> TLR
TMB_BRAF -> TxTMB
CRP -> PFS
CRP -> OS
CRP -> TLR
CRP -> TxCRP
T -> PFS
T -> OS
T -> TLR
T -> TxTMB
T -> TxCRP
TxTMB -> PFS
TxTMB -> OS
TxCRP -> PFS
TxCRP -> OS
TLR -> PFS
PFS -> OS
}'

# Simplified DAG for a clinical audience. Age and Sex are separate nodes so
# that the graph matches the full DAG (issue #184): both point to TMB/BRAF,
# only Age points to survival (Age -> OS in the full DAG; there is no
# Sex -> outcome edge). Interaction and latent nodes are omitted.
DAG_SIMPLE_SPEC <- 'dag {
bb="0,0,1,1"
Sex      [pos="0.050,0.100"]
Age      [pos="0.400,0.100"]
TMB_BRAF [pos="0.180,0.380"]
CRP      [pos="0.180,0.780"]
T        [exposure,pos="0.050,0.580"]
Survival [outcome,pos="0.750,0.580"]

Age -> TMB_BRAF
Age -> Survival
Sex -> TMB_BRAF
TMB_BRAF -> Survival
CRP -> Survival
T -> Survival
}'

# Construct the canonical objects when this helper is sourced by a report.
dag <- dagitty::dagitty(DAG_SPEC)
dag_sens <- dagitty::dagitty(DAG_SENS_SPEC)
dag_simple <- dagitty::dagitty(DAG_SIMPLE_SPEC)

dag_palette <- c(
  exposure = "#2166ac",
  outcome = "#b2182b",
  latent = "#1b7837"
)

#' Prepare a full DAG for plotting
#'
#' Flips the y-axis coordinates used by DAGitty and marks the hypothesised
#' treatment-by-biomarker interaction edges.
#'
#' @param dag_obj A dagitty object using the canonical interaction-node names.
prepare_full_dag <- function(dag_obj) {
  tidy_dag_obj <- ggdag::tidy_dagitty(dag_obj)
  tidy_dag_obj$data <- tidy_dag_obj$data |>
    dplyr::mutate(
      y = 1 - y,
      yend = dplyr::if_else(is.na(yend), NA_real_, 1 - yend),
      hyp = !is.na(to) &
        name %in% c("TxTMB", "TxCRP") &
        to %in% c("PFS", "OS")
    )
  tidy_dag_obj
}

#' Plot a prepared full or sensitivity DAG
#'
#' @param tidy_dag_obj A ggdag tidy object with a logical `hyp` edge column.
#' @param caption_pal Named node-status colour palette.
plot_dag <- function(tidy_dag_obj, caption_pal = dag_palette) {
  node_data <- ggdag::node_status(tidy_dag_obj)$data |>
    dplyr::distinct(name, x, y, status)
  x_range <- range(node_data$x, na.rm = TRUE)
  y_range <- range(node_data$y, na.rm = TRUE)
  x_pad <- diff(x_range) * 0.10
  y_pad <- diff(y_range) * 0.10

  ggplot2::ggplot(
    tidy_dag_obj$data,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
  ) +
    ggdag::geom_dag_edges_link(
      data = function(d) dplyr::filter(d, !is.na(to), !hyp)
    ) +
    ggdag::geom_dag_edges_link(
      data = function(d) dplyr::filter(d, !is.na(to), hyp),
      edge_linetype = "dashed"
    ) +
    ggdag::geom_dag_node(
      data = node_data,
      mapping = ggplot2::aes(x = x, y = y, colour = status, fill = status),
      size = 24,
      inherit.aes = FALSE
    ) +
    ggdag::geom_dag_text(
      data = node_data,
      mapping = ggplot2::aes(x = x, y = y, label = name),
      colour = "white",
      size = 3.0,
      inherit.aes = FALSE
    ) +
    ggplot2::coord_equal(
      xlim = c(x_range[1] - x_pad, x_range[2] + x_pad),
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad),
      clip = "off"
    ) +
    ggdag::theme_dag() +
    ggplot2::scale_color_manual(
      name = "Node type",
      values = caption_pal,
      na.value = "grey55",
      breaks = c("exposure", "latent", "outcome"),
      labels = c(
        exposure = "Exposure (T)",
        latent = "Latent (U)",
        outcome = "Outcome (PFS, OS)"
      )
    ) +
    ggplot2::scale_fill_manual(
      name = "Node type",
      values = caption_pal,
      na.value = "grey55",
      breaks = c("exposure", "latent", "outcome"),
      labels = c(
        exposure = "Exposure (T)",
        latent = "Latent (U)",
        outcome = "Outcome (PFS, OS)"
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

#' Prepare and plot a canonical full or sensitivity DAG
#'
#' @param dag_obj A dagitty object using the canonical interaction-node names.
#' @param caption_pal Named node-status colour palette.
plot_full_dag <- function(dag_obj, caption_pal = dag_palette) {
  plot_dag(prepare_full_dag(dag_obj), caption_pal)
}

#' Plot the simplified clinical-effectiveness DAG
#'
#' @param dag_obj The simplified dagitty object.
plot_simple_dag <- function(dag_obj = dag_simple) {
  tidy_dag_obj <- ggdag::tidy_dagitty(dag_obj)
  tidy_dag_obj$data <- tidy_dag_obj$data |>
    dplyr::mutate(
      y = 1 - y,
      yend = dplyr::if_else(is.na(yend), NA_real_, 1 - yend),
      hyp = !is.na(to) & name == "T" & to == "Survival"
    )

  node_data <- ggdag::node_status(tidy_dag_obj)$data |>
    dplyr::distinct(name, x, y, status) |>
    dplyr::mutate(
      status = dplyr::if_else(is.na(status), "covariate", status),
      label = dplyr::case_when(
        name == "Survival" ~ "PFS/<br>OS",
        name == "TMB_BRAF" ~ "TMB/<br><i>BRAF</i>",
        TRUE ~ name
      )
    )
  x_range <- range(node_data$x, na.rm = TRUE)
  y_range <- range(node_data$y, na.rm = TRUE)
  palette <- c(
    exposure = "#2166ac",
    outcome = "#b2182b",
    covariate = "grey55"
  )

  ggplot2::ggplot(
    tidy_dag_obj$data,
    ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
  ) +
    ggdag::geom_dag_edges_link(
      data = function(d) dplyr::filter(d, !is.na(to), !hyp)
    ) +
    ggdag::geom_dag_edges_link(
      data = function(d) dplyr::filter(d, !is.na(to), hyp),
      edge_linetype = "dashed"
    ) +
    ggdag::geom_dag_node(
      data = node_data,
      mapping = ggplot2::aes(x = x, y = y, colour = status, fill = status),
      size = 24,
      inherit.aes = FALSE
    ) +
    ggtext::geom_richtext(
      data = node_data,
      mapping = ggplot2::aes(x = x, y = y, label = label),
      colour = "white",
      size = 3.5,
      fill = NA,
      label.color = NA,
      label.padding = grid::unit(c(0, 0, 0, 0), "pt"),
      inherit.aes = FALSE
    ) +
    ggplot2::coord_equal(
      xlim = c(x_range[1] - diff(x_range) * 0.10, x_range[2] + diff(x_range) * 0.10),
      ylim = c(y_range[1] - diff(y_range) * 0.10, y_range[2] + diff(y_range) * 0.10),
      clip = "off"
    ) +
    ggdag::theme_dag() +
    ggplot2::scale_color_manual(
      name = NULL,
      values = palette,
      breaks = c("exposure", "outcome", "covariate"),
      labels = c(
        exposure = "Treatment",
        outcome = "Outcome",
        covariate = "Covariates and biomarkers"
      )
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = palette,
      breaks = c("exposure", "outcome", "covariate"),
      labels = c(
        exposure = "Treatment",
        outcome = "Outcome",
        covariate = "Covariates and biomarkers"
      )
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

#' Format a DAGitty conditional-independence statement
format_ci_statement <- function(ci_obj) {
  z <- ci_obj$Z
  if (length(z) == 0) {
    sprintf("%s _||_ %s", ci_obj$X, ci_obj$Y)
  } else {
    sprintf(
      "%s _||_ %s | {%s}",
      ci_obj$X, ci_obj$Y, paste(sort(z), collapse = ", ")
    )
  }
}
