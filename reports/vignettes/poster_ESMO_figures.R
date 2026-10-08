# =============================================================================
# ESMO 2026 poster 4355P: poster-sized figures and number cross-check
# -----------------------------------------------------------------------------
# Regenerates the three figures embedded in the ESMO 2026 poster (poster_ESMO.pptx)
# (Kaplan-Meier grid, interaction forest plot, simplified DAG) at poster font
# sizes, and recomputes the hazard ratios typed into the poster's Table 1 from
# the trial data with the same fits as reports/clinical_effectiveness.qmd and
# R/dag_association_tests.R (Firth Cox, profile-likelihood CIs). Outputs go to
# data/output/poster_ESMO/ (ignored by git; override with the POSTER_OUT
# environment variable). The ridge point estimates are typed from the rendered
# clinical effectiveness report because the ridge fit lives in that document.
#
#   Rscript reports/vignettes/poster_ESMO_figures.R      # from the repo root
#   python docs/posters/esmo_2026/poster_ESMO_build.py   # then rebuild the .pptx (private repository)
# =============================================================================
out_dir <- Sys.getenv("POSTER_OUT", unset = here::here("data", "output", "poster_ESMO"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source(here::here("R/report_setup.R"))
setup_report(
  sources = c("02", "03"),
  funs = c("tlr_landmark", "dag_helpers", "cox_extract"),
  packages = c("here", "ggplot2", "dplyr", "survival", "survminer", "cowplot",
               "coxphf", "logistf", "dagitty", "ggdag", "ggtext", "ggraph"),
  set_knitr = FALSE, set_theme = FALSE, quiet_sources = TRUE
)

# ---- Cohorts (same rules as reports/clinical_effectiveness.qmd) ---------------
data_complete <- data[complete.cases(data[, c("Age", "sex", "Rx", "crp",
                                               "tmb_braf", "OSwk", "Death",
                                               "PFSwk", "Progression")]), ]
data_complete$sex <- as.factor(data_complete$sex)
data_complete$Rx  <- as.factor(data_complete$Rx)
data_complete$crp_num      <- as.numeric(as.character(data_complete$crp))
data_complete$tlr_num      <- as.numeric(as.character(data_complete$tlr))
data_complete$tmb_braf_num <- as.numeric(as.character(data_complete$tmb_braf))
data_complete$Rx_num <- as.numeric(data_complete$Rx == "Experimental arm")
# Shared TLR landmark base cohort (issue #184): the TLR-classified patients of
# the complete-case cohort, as in the clinical report, Figure 1 and Table S2.
data_tlr <- tlr_landmark_base_cohort(data_complete)
stopifnot(identical(sort(data_tlr$ID), sort(tlr_landmark_base_cohort(data)$ID)))
lm9 <- build_tlr_landmark_cohorts(data_tlr, landmark = "week9")

cat("n complete =", nrow(data_complete), "; n TLR =", nrow(data_tlr),
    "; landmark OS n =", nrow(lm9$os), "; landmark PFS n =", nrow(lm9$pfs), "\n")

# ---- Numbers -----------------------------------------------------------------
firth <- function(formula, d, term) {
  fit <- coxphf::coxphf(formula, data = d, firth = TRUE, pl = TRUE,
                        maxit = 100, maxstep = 0.1)
  data.frame(term = term, HR = exp(unname(fit$coefficients[term])),
             lo = unname(fit$ci.lower[term]), hi = unname(fit$ci.upper[term]),
             p = unname(fit$prob[term]))
}
rows <- list()
add <- function(label, endpoint, df) {
  rows[[length(rows) + 1]] <<- cbind(label = label, endpoint = endpoint, df)
}
# marginal (one term per model, as in the DAG edge tests)
for (ep in c("OS", "PFS")) {
  f <- if (ep == "OS") function(rhs) as.formula(paste("Surv(OSwk, Death) ~", rhs))
       else function(rhs) as.formula(paste("Surv(PFSwk, Progression) ~", rhs))
  add("Treatment (marginal)", ep, firth(f("Rx_num"), data_complete, "Rx_num"))
  add("CRP (marginal)", ep, firth(f("crp_num"), data_complete, "crp_num"))
  add("TMB/BRAF (marginal)", ep, firth(f("tmb_braf_num"), data_complete, "tmb_braf_num"))
}
# TLR, week-9 landmark, unadjusted, 65-patient TLR cohort
add("TLR (landmark, unadjusted)", "OS",
    firth(Surv(OSwk_lm, Death) ~ tlr_num, lm9$os, "tlr_num"))
add("TLR (landmark, unadjusted)", "PFS",
    firth(Surv(PFSwk_lm, Progression) ~ tlr_num, lm9$pfs, "tlr_num"))

# unified model (report formulas)
formula_os  <- Surv(OSwk, Death) ~ Age + sex + Rx + crp_num + tmb_braf_num +
  crp_num:Rx + tmb_braf_num:Rx
formula_pfs <- Surv(PFSwk, Progression) ~ Age + sex + Rx + crp_num + tmb_braf_num +
  crp_num:Rx + tmb_braf_num:Rx
model_os  <- fit_firth_cox(formula_os,  data_complete, "Unified OS")
model_pfs <- fit_firth_cox(formula_pfs, data_complete, "Unified PFS")
res_os  <- extract_coxphf_results(model_os)
res_pfs <- extract_coxphf_results(model_pfs)
write.csv(res_os,  file.path(out_dir, "unified_os.csv"),  row.names = FALSE)
write.csv(res_pfs, file.path(out_dir, "unified_pfs.csv"), row.names = FALSE)

# TLR landmark responder model (report formula, 65-patient cohort)
formula_lm_os  <- Surv(OSwk_lm, Death) ~ Age + sex + Rx + crp_num + tmb_braf_num +
  tlr_num + tlr_num:Rx
formula_lm_pfs <- Surv(PFSwk_lm, Progression) ~ Age + sex + Rx + crp_num + tmb_braf_num +
  tlr_num + tlr_num:Rx
lm_os  <- extract_coxphf_results(fit_firth_cox(formula_lm_os,  lm9$os,  "LM OS"))
lm_pfs <- extract_coxphf_results(fit_firth_cox(formula_lm_pfs, lm9$pfs, "LM PFS"))
write.csv(lm_os,  file.path(out_dir, "landmark_os.csv"),  row.names = FALSE)
write.csv(lm_pfs, file.path(out_dir, "landmark_pfs.csv"), row.names = FALSE)

marginal <- do.call(rbind, rows)
write.csv(marginal, file.path(out_dir, "marginal.csv"), row.names = FALSE)

# TLR prevalence by arm in the 65-patient cohort
tlr_tab <- table(data_tlr$Rx, data_tlr$tlr_num)
fisher_p <- fisher.test(tlr_tab)$p.value
lf <- logistf::logistf(tlr_num ~ Rx_num, data = data_tlr, firth = TRUE, pl = TRUE)
tlr_summary <- data.frame(
  arm = rownames(tlr_tab), n = rowSums(tlr_tab), tlr_pos = tlr_tab[, "1"],
  pct = 100 * tlr_tab[, "1"] / rowSums(tlr_tab)
)
tlr_summary$fisher_p <- fisher_p
tlr_summary$firth_or <- exp(unname(lf$coefficients["Rx_num"]))
tlr_summary$firth_lo <- exp(unname(lf$ci.lower["Rx_num"]))
tlr_summary$firth_hi <- exp(unname(lf$ci.upper["Rx_num"]))
tlr_summary$firth_p  <- unname(lf$prob["Rx_num"])
write.csv(tlr_summary, file.path(out_dir, "tlr_by_arm.csv"), row.names = FALSE)

# CRP positivity by arm (68-patient cohort) and demographics
crp_tab <- table(data_complete$Rx, data_complete$crp_num)
demo <- data.frame(
  n = nrow(data_complete),
  age_mean = mean(data_complete$Age), age_sd = sd(data_complete$Age),
  female_n = sum(data_complete$sex == "0"),
  female_pct = 100 * mean(data_complete$sex == "0"),
  deaths = sum(data_complete$Death), pfs_events = sum(data_complete$Progression),
  n_control = sum(data_complete$Rx == "Control arm"),
  n_exp = sum(data_complete$Rx == "Experimental arm"),
  crp_pos_control = crp_tab["Control arm", "1"],
  crp_pos_exp = crp_tab["Experimental arm", "1"],
  crp_fisher_p = fisher.test(crp_tab)$p.value,
  tmb_pos = sum(data_complete$tmb_braf_num), crp_pos = sum(data_complete$crp_num)
)
write.csv(demo, file.path(out_dir, "demo.csv"), row.names = FALSE)

# ---- Figure 2: Kaplan-Meier, 3 columns x 2 rows -----------------------------
make_group <- function(rx_num, biomarker, pos_label, neg_label) {
  arm <- ifelse(rx_num == 1, "FLOX/nivolumab", "FLOX alone")
  bm  <- ifelse(biomarker == 1, pos_label, neg_label)
  factor(paste0(arm, ", ", bm),
         levels = c(paste0("FLOX/nivolumab, ", pos_label),
                    paste0("FLOX/nivolumab, ", neg_label),
                    paste0("FLOX alone, ", pos_label),
                    paste0("FLOX alone, ", neg_label)))
}
data_complete$grp_crp <- make_group(data_complete$Rx_num, data_complete$crp_num, "CRP+", "CRP-")
data_complete$grp_tmb <- make_group(data_complete$Rx_num, data_complete$tmb_braf_num, "TMB/BRAF+", "TMB/BRAF-")
lm9$os$grp_tlr  <- make_group(lm9$os$Rx_num,  lm9$os$tlr_num,  "TLR+", "TLR-")
lm9$pfs$grp_tlr <- make_group(lm9$pfs$Rx_num, lm9$pfs$tlr_num, "TLR+", "TLR-")

km_palette   <- c("#B2182B", "#F4A582", "#2166AC", "#92C5DE")
km_linetypes <- c("solid", "dashed", "solid", "dashed")
km_theme <- theme_minimal(base_size = 18) + theme(
  axis.title       = element_text(size = 22),
  axis.text        = element_text(size = 20),
  legend.title     = element_blank(),
  legend.text      = element_text(size = 19),
  legend.position  = "bottom",
  legend.key.width = unit(2.4, "lines"),
  plot.margin      = margin(8, 14, 4, 6)
)
km_plot <- function(d, time_col, event_col, group_col, show_ylab, tag, xlab_text) {
  full_lvls <- levels(d[[group_col]])
  present   <- full_lvls[full_lvls %in% unique(as.character(d[[group_col]]))]
  idx       <- match(present, full_lvls)
  formula <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", group_col))
  km_fit  <- eval(bquote(survfit(.(formula), data = d)))
  sp <- ggsurvplot(
    km_fit, data = d, pval = FALSE, conf.int = FALSE,
    palette = km_palette[idx], linetype = km_linetypes[idx],
    legend = "bottom", legend.labs = present,
    xlab = xlab_text, ylab = if (show_ylab) "Survival probability" else "",
    ggtheme = km_theme, font.x = 22, font.y = 22, font.tickslab = 20,
    size = 1.6, censor.size = 7
  )
  p <- sp$plot + theme(legend.position = "none")  # one shared legend, see below
  if (!show_ylab) p <- p + theme(axis.title.y = element_blank())
  p + labs(tag = tag) +
    theme(plot.tag = element_text(size = 24, face = "bold"), plot.tag.position = "topleft")
}
x_rand <- "Weeks since randomization"
p_os_crp  <- km_plot(data_complete, "OSwk",  "Death",       "grp_crp", TRUE,  "a)", x_rand)
p_os_tmb  <- km_plot(data_complete, "OSwk",  "Death",       "grp_tmb", FALSE, "b)", x_rand)
p_os_tlr  <- km_plot(lm9$os,  "OSwk_lm",  "Death",       "grp_tlr", FALSE, "c)",
                     sprintf("Weeks since week-9 landmark (n = %d)", nrow(lm9$os)))
p_pfs_crp <- km_plot(data_complete, "PFSwk", "Progression", "grp_crp", TRUE,  "d)", x_rand)
p_pfs_tmb <- km_plot(data_complete, "PFSwk", "Progression", "grp_tmb", FALSE, "e)", x_rand)
p_pfs_tlr <- km_plot(lm9$pfs, "PFSwk_lm", "Progression", "grp_tlr", FALSE, "f)",
                     sprintf("Weeks since week-9 landmark (n = %d)", nrow(lm9$pfs)))

col_lab <- function(t) ggdraw() + draw_label(t, fontface = "bold", size = 26)
row_lab <- function(t) ggdraw() + draw_label(t, fontface = "bold", size = 26, angle = 90)
rw <- c(0.05, 1, 1, 1)
header  <- plot_grid(NULL, col_lab(sprintf("CRP <5 mg/L at week 4 (n = %d)", nrow(data_complete))),
                     col_lab(sprintf("TMB ≥9 mut/Mb or BRAF V600E (n = %d)", nrow(data_complete))),
                     col_lab(sprintf("TLR ≥10%%, week-9 landmark (n = %d)", nrow(data_tlr))),
                     ncol = 4, rel_widths = rw)
row_os  <- plot_grid(row_lab("Overall survival"), p_os_crp, p_os_tmb, p_os_tlr, ncol = 4, rel_widths = rw)
row_pfs <- plot_grid(row_lab("Progression-free survival"), p_pfs_crp, p_pfs_tmb, p_pfs_tlr, ncol = 4, rel_widths = rw)
fig1_title <- ggdraw() + draw_label(
  "Figure 2. Kaplan-Meier survival by biomarker status and treatment arm",
  fontface = "bold", size = 30, x = 0.005, hjust = 0)
# One shared legend for all six panels: the colour/linetype mapping (arm x
# biomarker status) is identical in every panel. cowplot::get_legend() returns
# an empty grob under ggplot2 3.5, so the guide box is extracted by name.
extract_legend <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  idx <- which(grepl("guide-box", g$layout$name))
  for (i in idx) if (!inherits(g$grobs[[i]], "zeroGrob")) return(g$grobs[[i]])
  stop("No legend found in the legend panel")
}
leg_labels <- c("FLOX/nivolumab, biomarker-positive", "FLOX/nivolumab, biomarker-negative",
                "FLOX alone, biomarker-positive", "FLOX alone, biomarker-negative")
leg_df <- data.frame(x = rep(1:2, 4), y = rep(1:2, 4),
                     g = factor(rep(leg_labels, each = 2), levels = leg_labels))
p_leg <- ggplot(leg_df, aes(x, y, colour = g, linetype = g)) +
  geom_line(linewidth = 1.6) +
  scale_colour_manual(values = km_palette, name = NULL) +
  scale_linetype_manual(values = km_linetypes, name = NULL) +
  guides(colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1)) +
  theme_minimal() +
  theme(legend.position = "bottom", legend.text = element_text(size = 22),
        legend.key.width = unit(3.2, "lines"), legend.spacing.x = unit(1.2, "lines"))
shared_legend <- extract_legend(p_leg)
fig1 <- plot_grid(fig1_title, header, row_os, row_pfs, shared_legend, ncol = 1,
                  rel_heights = c(0.10, 0.07, 1, 1, 0.11))
ggsave(file.path(out_dir, "poster_esmo_km.png"), fig1, width = 24, height = 13.5,
       dpi = 150, bg = "white")

# ---- Figure 3: interaction forest plot (Firth + ridge) -----------------------
gi <- function(res, pat) get_interaction_hr(res, pat)
hr_crp_os  <- gi(res_os,  "crp_num.*Rx|Rx.*crp_num")
hr_tmb_os  <- gi(res_os,  "tmb_braf_num.*Rx|Rx.*tmb_braf_num")
hr_crp_pfs <- gi(res_pfs, "crp_num.*Rx|Rx.*crp_num")
hr_tmb_pfs <- gi(res_pfs, "tmb_braf_num.*Rx|Rx.*tmb_braf_num")
# Ridge point estimates from reports/clinical_effectiveness.md (v3.11, issue #183:
# glmnet Cox, biomarker main effects and interactions penalized, CV lambda)
ridge <- c(crp_os = 0.80, tmb_os = 0.91, crp_pfs = 0.66, tmb_pfs = 0.74)
lab_crp <- "CRP × treatment"; lab_tmb <- "TMB/BRAF × treatment"
forest <- data.frame(
  Biomarker = factor(rep(c(lab_crp, lab_tmb), 2), levels = c(lab_tmb, lab_crp)),
  Outcome = factor(rep(c("Overall survival", "Progression-free survival"), each = 2),
                   levels = c("Overall survival", "Progression-free survival")),
  HR = c(hr_crp_os$HR, hr_tmb_os$HR, hr_crp_pfs$HR, hr_tmb_pfs$HR),
  lo = c(hr_crp_os$CI_Lower, hr_tmb_os$CI_Lower, hr_crp_pfs$CI_Lower, hr_tmb_pfs$CI_Lower),
  hi = c(hr_crp_os$CI_Upper, hr_tmb_os$CI_Upper, hr_crp_pfs$CI_Upper, hr_tmb_pfs$CI_Upper),
  ridge = c(ridge["crp_os"], ridge["tmb_os"], ridge["crp_pfs"], ridge["tmb_pfs"])
)
stopifnot(min(forest$lo) > 0.045, max(forest$hi) < 9)  # every interval inside the axis
forest$label <- sprintf("%.2f (%.2f–%.2f)", forest$HR, forest$lo, forest$hi)
write.csv(forest, file.path(out_dir, "forest.csv"), row.names = FALSE)
est_levels <- c("Firth Cox (95% profile-likelihood CI)", "Ridge Cox (point estimate)")
pts <- rbind(
  data.frame(forest[, c("Biomarker", "Outcome")], x = forest$HR, est = est_levels[1]),
  data.frame(forest[, c("Biomarker", "Outcome")], x = forest$ridge, est = est_levels[2])
)
pts$est <- factor(pts$est, levels = est_levels)
fig2 <- ggplot(forest, aes(y = Biomarker)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45", linewidth = 1) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, colour = "#0072B2", linewidth = 1.4) +
  geom_point(data = pts, aes(x = x, shape = est, fill = est), size = 7, colour = "#0072B2", stroke = 1.6) +
  geom_label(aes(x = HR, label = label), vjust = -0.9, size = 6.3, fill = "white",
             label.size = 0, label.padding = unit(0.18, "lines")) +
  scale_shape_manual(values = c(21, 23), name = NULL) +
  scale_fill_manual(values = c("#0072B2", "white"), name = NULL) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4, 8),
                labels = c("0.05", "0.1", "0.25", "0.5", "1", "2", "4", "8"),
                limits = c(0.045, 9)) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  facet_wrap(~ Outcome) +
  labs(x = "Hazard ratio for the treatment × biomarker interaction (log scale)", y = NULL,
       title = "Figure 3. Treatment × biomarker interactions",
       subtitle = "Unified model; HR < 1: biomarker-positive patients benefit more from nivolumab") +
  theme_minimal(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 22, face = "bold"),
        plot.title = element_text(size = 26, face = "bold"),
        plot.subtitle = element_text(size = 19),
        axis.text = element_text(size = 20), axis.title.x = element_text(size = 21),
        legend.position = "bottom", legend.text = element_text(size = 19),
        legend.key.size = unit(1.8, "lines"), plot.margin = margin(8, 16, 4, 8))
ggsave(file.path(out_dir, "poster_esmo_forest.png"), fig2, width = 13.1, height = 7.0,
       dpi = 150, bg = "white")

# ---- Figure 1: simplified DAG at poster sizes --------------------------------
# Same nodes and edges as DAG_SIMPLE_SPEC (dag_helpers.R); coordinates spread
# out so that poster-sized nodes do not overlap.
poster_dag <- dagitty::dagitty('dag {
bb="0,0,1,1"
Sex      [pos="0.02,0.06"]
Age      [pos="0.46,0.06"]
TMB_BRAF [pos="0.22,0.36"]
CRP      [pos="0.22,0.84"]
T        [exposure,pos="0.02,0.60"]
Survival [outcome,pos="0.90,0.60"]
Age -> TMB_BRAF
Age -> Survival
Sex -> TMB_BRAF
TMB_BRAF -> Survival
CRP -> Survival
T -> Survival
}')
stopifnot(identical(sort(dagitty::edges(poster_dag)$v), sort(dagitty::edges(dag_simple)$v)),
          identical(sort(dagitty::edges(poster_dag)$w), sort(dagitty::edges(dag_simple)$w)))
tidy_dag_obj <- ggdag::tidy_dagitty(poster_dag)
tidy_dag_obj$data <- tidy_dag_obj$data |>
  dplyr::mutate(y = 1 - y,
                yend = dplyr::if_else(is.na(yend), NA_real_, 1 - yend),
                hyp = !is.na(to) & name == "T" & to == "Survival")
node_data <- ggdag::node_status(tidy_dag_obj)$data |>
  dplyr::distinct(name, x, y, status) |>
  dplyr::mutate(status = dplyr::if_else(is.na(status), "covariate", status),
                label = dplyr::case_when(name == "Survival" ~ "PFS/<br>OS",
                                         name == "TMB_BRAF" ~ "TMB/<br><i>BRAF</i>",
                                         TRUE ~ name))
x_range <- range(node_data$x); y_range <- range(node_data$y)
palette <- c(exposure = "#2166ac", outcome = "#b2182b", covariate = "grey55")
cap <- ggraph::circle(12.5, "mm")
leg_labels <- c("Treatment", "Outcome", "Covariates / biomarkers")
fig3 <- ggplot(tidy_dag_obj$data, aes(x = x, y = y, xend = xend, yend = yend)) +
  ggdag::geom_dag_edges_link(data = function(d) dplyr::filter(d, !is.na(to), !hyp),
                             edge_width = 1.5, start_cap = cap, end_cap = cap,
                             arrow = grid::arrow(length = grid::unit(6, "mm"), type = "closed")) +
  ggdag::geom_dag_edges_link(data = function(d) dplyr::filter(d, !is.na(to), hyp),
                             edge_linetype = "dashed", edge_width = 1.5, start_cap = cap, end_cap = cap,
                             arrow = grid::arrow(length = grid::unit(6, "mm"), type = "closed")) +
  ggdag::geom_dag_node(data = node_data, aes(x = x, y = y, colour = status, fill = status),
                       size = 32, inherit.aes = FALSE) +
  ggtext::geom_richtext(data = node_data, aes(x = x, y = y, label = label), colour = "white",
                        size = 6.5, fill = NA, label.color = NA,
                        label.padding = grid::unit(c(0, 0, 0, 0), "pt"), inherit.aes = FALSE) +
  coord_equal(xlim = c(x_range[1] - diff(x_range) * 0.12, x_range[2] + diff(x_range) * 0.12),
              ylim = c(y_range[1] - diff(y_range) * 0.16, y_range[2] + diff(y_range) * 0.16),
              clip = "off") +
  ggdag::theme_dag() +
  scale_color_manual(name = NULL, values = palette, breaks = c("exposure", "outcome", "covariate"),
                     labels = leg_labels) +
  scale_fill_manual(name = NULL, values = palette, breaks = c("exposure", "outcome", "covariate"),
                    labels = leg_labels) +
  guides(colour = guide_legend(override.aes = list(size = 9)),
         fill = guide_legend(override.aes = list(size = 9))) +
  # The "Figure 1." title sits in the poster caption box under the image.
  theme(legend.position = "bottom", legend.text = element_text(size = 18),
        legend.key.size = unit(1.3, "lines"), plot.margin = margin(2, 2, 2, 2))
ggsave(file.path(out_dir, "poster_esmo_dag.png"), fig3, width = 7.0, height = 4.8,
       dpi = 200, bg = "white")

cat("done\n")
