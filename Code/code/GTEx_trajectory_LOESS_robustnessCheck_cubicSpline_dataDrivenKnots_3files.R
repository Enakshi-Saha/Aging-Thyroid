# Sensitivity analyses for aging pathway trajectory plots
# Addresses reviewer comments on robustness of non-linear transitions
# in females around ages ~40 and ~60

library(data.table)
library(ggplot2)
library(dplyr)
library(fgsea)
library(splines)
library(boot)
library(stringr)

set.seed(1)


# HELPER FUNCTIONS

# Bootstrap confidence bands: resample individuals, aggregate, then fit LOESS
bootstrap_smooth <- function(df, n_boot = 500, span_val = 0.75) {
  age_seq <- seq(min(df$age_group), max(df$age_group), length.out = 100)
  
  boot_preds <- lapply(1:n_boot, function(b) {
    df_boot <- df %>%
      group_by(sex) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      ungroup() %>%
      group_by(age_group, sex) %>%
      summarise(indegree = mean(indegree), .groups = "drop")
    
    df_boot %>%
      group_by(sex) %>%
      group_modify(~ {
        fit <- loess(indegree ~ age_group, data = .x, span = span_val)
        data.frame(age_seq = age_seq,
                   pred    = predict(fit, newdata = data.frame(age_group = age_seq)))
      }) %>%
      ungroup()
  })
  
  boot_df <- bind_rows(boot_preds)
  
  ci_df <- boot_df %>%
    group_by(sex, age_seq) %>%
    summarise(
      lower = quantile(pred, 0.025, na.rm = TRUE),
      upper = quantile(pred, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  return(ci_df)
}

# Subsampling stability: resample individuals, aggregate, then fit LOESS
subsample_smooth <- function(df, frac = 0.8, n_iter = 100, span_val = 0.75) {
  age_seq <- seq(min(df$age_group), max(df$age_group), length.out = 100)
  
  sub_preds <- lapply(1:n_iter, function(b) {
    df_sub <- df %>%
      group_by(sex) %>%
      slice_sample(prop = frac, replace = FALSE) %>%
      ungroup() %>%
      group_by(age_group, sex) %>%
      summarise(indegree = mean(indegree), .groups = "drop")
    
    df_sub %>%
      group_by(sex) %>%
      group_modify(~ {
        fit <- loess(indegree ~ age_group, data = .x, span = span_val)
        data.frame(age_seq  = age_seq,
                   pred     = predict(fit, newdata = data.frame(age_group = age_seq)),
                   iter     = b)
      }) %>%
      ungroup()
  })
  bind_rows(sub_preds)
}


# PLOT FUNCTIONS

plot_loess_spans_with_CI <- function(fulldata, pathway_class, span_values = c(0.6, 0.75, 0.9)) {
  
  fulldata_aggregated <- fulldata %>%
    group_by(age_group, sex) %>%
    summarise(indegree = mean(indegree), mean_pathScore = mean(indegree), .groups = "drop") %>%
    mutate(sex = factor(sex, levels = c("FEMALE", "MALE")))
  
  ci_df <- bootstrap_smooth(fulldata, n_boot = 500, span_val = 0.75)
  ci_df$sex <- factor(ci_df$sex, levels = c("FEMALE", "MALE"))
  
  age_seq <- seq(min(fulldata_aggregated$age_group), max(fulldata_aggregated$age_group), length.out = 100)
  span_lines <- lapply(span_values, function(sp) {
    fulldata_aggregated %>%
      group_by(sex) %>%
      group_modify(~ {
        fit <- loess(indegree ~ age_group, data = .x, span = sp)
        data.frame(age_seq = age_seq,
                   pred    = predict(fit, newdata = data.frame(age_group = age_seq)),
                   span    = as.character(sp))
      }) %>%
      ungroup()
  })
  span_df <- bind_rows(span_lines)
  span_df$sex <- factor(span_df$sex, levels = c("FEMALE", "MALE"))
  
  p <- ggplot() +
    geom_ribbon(data = ci_df %>% filter(sex == "FEMALE"),
                aes(x = age_seq, ymin = lower, ymax = upper),
                fill = "#F8766D", alpha = 0.15) +
    geom_ribbon(data = ci_df %>% filter(sex == "MALE"),
                aes(x = age_seq, ymin = lower, ymax = upper),
                fill = "#00BFC4", alpha = 0.15) +
    geom_line(data = span_df,
              aes(x = age_seq, y = pred, color = sex, linetype = span),
              size = 0.8) +
    geom_point(data = fulldata_aggregated,
               aes(x = age_group, y = mean_pathScore, color = sex),
               size = 1.2, alpha = 0.5) +
    geom_vline(xintercept = c(40, 60), linetype = "dotted",
               color = "grey40", size = 0.5) +
    annotate("text", x = 40, y = Inf, label = "~40y", vjust = 1.5,
             hjust = -0.1, size = 3, color = "grey30") +
    annotate("text", x = 60, y = Inf, label = "~60y", vjust = 1.5,
             hjust = -0.1, size = 3, color = "grey30") +
    scale_linetype_manual(name = "LOESS span",
                          values = c("0.6" = "dotdash", "0.75" = "solid", "0.9" = "dashed")) +
    labs(title = str_to_title(pathway_class),
         x = "Age", y = "score", color = "Sex") +
    theme_bw(base_size = 25) +
    theme(plot.title    = element_text(size = 20, face = "bold"),
          legend.position = "right")
  
  return(p)
}

plot_subsample_stability <- function(fulldata, pathway_class, frac = 0.8, n_iter = 100) {
  
  sub_df <- subsample_smooth(fulldata, frac = frac, n_iter = n_iter)
  sub_df$sex <- factor(sub_df$sex, levels = c("FEMALE", "MALE"))
  
  median_df <- sub_df %>%
    group_by(sex, age_seq) %>%
    summarise(median_pred = median(pred, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot() +
    geom_line(data = sub_df,
              aes(x = age_seq, y = pred, group = interaction(iter, sex), color = sex),
              alpha = 0.07, size = 0.4) +
    geom_line(data = median_df,
              aes(x = age_seq, y = median_pred, color = sex),
              size = 1.2) +
    geom_vline(xintercept = c(40, 60), linetype = "dotted",
               color = "grey40", size = 0.5) +
    annotate("text", x = 40, y = Inf, label = "~40y", vjust = 1.5,
             hjust = -0.1, size = 3, color = "grey30") +
    annotate("text", x = 60, y = Inf, label = "~60y", vjust = 1.5,
             hjust = -0.1, size = 3, color = "grey30") +
    labs(title = str_to_title(pathway_class),
         x = "Age", y = "score", color = "Sex") +
    theme_bw(base_size = 25) +
    theme(plot.title    = element_text(size = 20, face = "bold"),
          legend.position = "right")
  
  return(p)
}

plot_cubic_spline <- function(fulldata, pathway_class) {
  
  fulldata_mod <- fulldata %>%
    group_by(age_group, sex) %>%
    summarise(indegree = mean(indegree), .groups = "drop") %>%
    mutate(sex = factor(sex, levels = c("FEMALE", "MALE")))
  
  age_seq  <- seq(min(fulldata_mod$age_group), max(fulldata_mod$age_group), length.out = 200)
  newdata  <- data.frame(age_group = age_seq)
  
  # Knots at quantiles of observed age
  age_knots <- quantile(fulldata_mod$age_group, probs = c(0.25, 0.5, 0.75))
  
  spline_preds <- fulldata_mod %>%
    group_by(sex) %>%
    group_modify(~ {
      fit      <- lm(indegree ~ ns(age_group, knots = age_knots), data = .x)
      pred_obj <- predict(fit, newdata = newdata, interval = "confidence", level = 0.95)
      data.frame(age_seq = age_seq,
                 fit     = pred_obj[, "fit"],
                 lower   = pred_obj[, "lwr"],
                 upper   = pred_obj[, "upr"])
    }) %>%
    ungroup()
  spline_preds$sex <- factor(spline_preds$sex, levels = c("FEMALE", "MALE"))
  
  loess_lines <- fulldata_mod %>%
    group_by(sex) %>%
    group_modify(~ {
      fit <- loess(indegree ~ age_group, data = .x, span = 0.75)
      data.frame(age_seq = age_seq,
                 pred    = predict(fit, newdata = newdata))
    }) %>%
    ungroup()
  loess_lines$sex <- factor(loess_lines$sex, levels = c("FEMALE", "MALE"))
  
  p <- ggplot() +
    geom_ribbon(data = spline_preds,
                aes(x = age_seq, ymin = lower, ymax = upper, fill = sex),
                alpha = 0.15) +
    geom_line(data = spline_preds,
              aes(x = age_seq, y = fit, color = sex,
                  linetype = "Natural cubic spline"),
              size = 1.1) +
    geom_line(data = loess_lines,
              aes(x = age_seq, y = pred, color = sex,
                  linetype = "LOESS (span 0.75)"),
              size = 0.8, alpha = 0.7) +
    geom_vline(xintercept = age_knots, linetype = "dotted",
               color = "grey50", size = 0.4) +
    annotate("text", x = age_knots, y = Inf,
             label = paste0("Q", c(25, 50, 75), "\n(", round(age_knots, 1), "y)"),
             vjust = 1.5, hjust = -0.05, size = 2.8, color = "grey40") +
    scale_color_manual(values = c("FEMALE" = "#F8766D", "MALE" = "#00BFC4")) +
    scale_fill_manual(values  = c("FEMALE" = "#F8766D", "MALE" = "#00BFC4")) +
    scale_linetype_manual(name   = "Method",
                          values = c("Natural cubic spline" = "solid",
                                     "LOESS (span 0.75)"    = "dashed")) +
    labs(title = str_to_title(pathway_class),
         x = "Age", y = "score",
         color = "Sex", fill = "Sex") +
    theme_bw(base_size = 25) +
    theme(plot.title    = element_text(size = 20, face = "bold"),
          legend.position = "right")
  
  return(p)
}


print_cover_page <- function(title_text, description_lines) {
  grid::grid.newpage()
  y <- 0.88
  grid::grid.text(title_text,
                  x = 0.05, y = y, just = "left",
                  gp = grid::gpar(fontsize = 16, fontface = "bold"))
  y <- y - 0.07
  for (line in description_lines) {
    grid::grid.text(line,
                    x = 0.05, y = y, just = "left",
                    gp = grid::gpar(fontsize = 11))
    y <- y - 0.055
  }
}


# PATHWAY SUPERCLASS
pathway_superclass <- c(
  "metabolism",
  "genetic information processing",
  "cell signaling, cell growth and death",
  "cell communication",
  "immune system and immune diseases",
  "endocrine, circulatory, nervous system & disease",
  "metabolic diseases",
  "infectious diseases"
)

pathway_superclass_index <- c(
  "c(2,4)", "c(14,15,17)", "c(20,23)", "c(21, 24)",
  "c(25, 34)", "c(26, 27, 29, 33, 35)", "37", "38"
)
path_superclass <- cbind(pathway_superclass, pathway_superclass_index)


# PDF 1: LOESS alternative spans + bootstrap CI
pdf("/home/esaha/thyroid_revision_results/revision_plots/sensitivity_1_LOESS_spans_bootstrapCI.pdf",
    width = 9, height = 6)

print_cover_page(
  title_text = "Sensitivity Analysis 1: LOESS with Alternative Smoothing Spans + Bootstrap CI",
  description_lines = c(
    "Method: LOESS fitted on age-aggregated mean pathway scores per sex.",
    "Three smoothing spans shown: 0.6 (dot-dash), 0.75 (solid, default), 0.9 (dashed).",
    "Shaded ribbon: 95% bootstrap CI around the default span (0.75),",
    "  derived from 500 bootstrap resamples of individuals, aggregated before smoothing.",
    "Vertical dotted lines mark the inferred hormonal transition ages (~40y, ~60y).",
    "",
    "Rationale: Demonstrates that the observed non-linear transitions in females are not",
    "  an artifact of the chosen smoothing parameter. Consistent inflection patterns",
    "  across all three spans support the robustness of the result.",
    "",
    "Note: Trajectories represent inferred network-level trends, not precise thresholds."
  )
)

for (i in 1:length(pathway_superclass)) {
  pathway_class <- pathway_superclass[i]
  j        <- eval(parse(text = path_superclass[i, 2]))
  pathList <- intersect(unlist(pathway_hierarchy[j]), sig_pathways)
  if (length(pathList) == 0) next
  fulldata$indegree <- apply(
    indegree[which(geneNames %in% unlist(pathways[which(names(pathways) %in% pathList)])), ],
    MARGIN = 2, mean)
  print(plot_loess_spans_with_CI(fulldata, pathway_class))
}
dev.off()

# PDF 2: Subsampling stability
pdf("/home/esaha/thyroid_revision_results/revision_plots/sensitivity_2_subsampling_stability.pdf",
    width = 9, height = 6)

print_cover_page(
  title_text = "Sensitivity Analysis 2: Subsampling Stability",
  description_lines = c(
    "Method: LOESS fitted on age-aggregated mean pathway scores per sex.",
    "100 iterations of 80% random subsampling (without replacement) of individuals.",
    "Each faint line = one subsampled run; bold line = median across all 100 runs.",
    "Vertical dotted lines mark the inferred hormonal transition ages (~40y, ~60y).",
    "",
    "Rationale: Demonstrates that the inferred transitions are not driven by any",
    "  particular subset of samples or by uneven sample density across age groups.",
    "  Consistent trajectory shapes across all subsamples support robustness.",
    "",
    "Note: Trajectories represent inferred network-level trends, not precise thresholds."
  )
)

for (i in 1:length(pathway_superclass)) {
  pathway_class <- pathway_superclass[i]
  j        <- eval(parse(text = path_superclass[i, 2]))
  pathList <- intersect(unlist(pathway_hierarchy[j]), sig_pathways)
  if (length(pathList) == 0) next
  fulldata$indegree <- apply(
    indegree[which(geneNames %in% unlist(pathways[which(names(pathways) %in% pathList)])), ],
    MARGIN = 2, mean)
  print(plot_subsample_stability(fulldata, pathway_class, frac = 0.8, n_iter = 100))
}
dev.off()

# PDF 3: Natural cubic spline vs LOESS
pdf("/home/esaha/thyroid_revision_results/revision_plots/sensitivity_3_cubic_spline_vs_LOESS.pdf",
    width = 9, height = 6)

print_cover_page(
  title_text = "Sensitivity Analysis 3: Natural Cubic Spline (data-driven knots) vs. LOESS",
  description_lines = c(
    "Method: Natural cubic spline regression (solid, +/-95% CI ribbon) vs. LOESS span 0.75 (dashed).",
    "Knots placed at the 25th, 50th, and 75th percentiles of the observed age distribution",
    "  — purely data-driven, with no reference to ages 40 or 60.",
    "Vertical dotted lines show where the data-driven knots actually fell.",
    "",
    "Rationale: Addresses potential circularity. Because knot locations are not pre-specified",
    "  at biologically motivated ages, any inflections near ~40y and ~60y recovered by the",
    "  spline represent independent confirmation from a parametric model with fundamentally",
    "  different assumptions than LOESS — not a built-in artifact of the smoothing choice.",
    "  Convergence between the two methods strengthens the biological conclusion.",
    "",
    "Note: Trajectories represent inferred network-level trends, not precise thresholds."
  )
)

for (i in 1:length(pathway_superclass)) {
  pathway_class <- pathway_superclass[i]
  j        <- eval(parse(text = path_superclass[i, 2]))
  pathList <- intersect(unlist(pathway_hierarchy[j]), sig_pathways)
  if (length(pathList) == 0) next
  fulldata$indegree <- apply(
    indegree[which(geneNames %in% unlist(pathways[which(names(pathways) %in% pathList)])), ],
    MARGIN = 2, mean)
  print(plot_cubic_spline(fulldata, pathway_class))
}
dev.off()

cat("Done. Three PDFs saved:\n")
