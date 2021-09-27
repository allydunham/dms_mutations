#!/us/bin/env Rscript
# Generate figure on DMS benchmark for thesis
source("src/config.R")
library(multipanelfigure)
library(ggtext)

theme_set(theme_pubclean() + theme(legend.position = 'right',
                                   plot.title = element_text(hjust = 0.5),
                                   plot.subtitle = element_text(hjust = 0.5),
                                   strip.background = element_blank(),
                                   legend.key = element_blank()))

tool_colours <- c(Envision="#984ea3", EVCouplings="#4daf4a", FoldX="#377eb8", PolyPhen2="#ff7f00", SIFT4G="#e41a1c")

study_meta <- readRDS('data/rdata/study_meta_data.RDS')

combine_study <- function(x) {
  tbl <- select(x$single_variants, any_of(c("variants", "score", "norm_score", "evcoup_epistatic",
                                            "pph2_prob", "envision_prediction", "sift_score")),
                ends_with("_ddG")) %>%
    filter(!str_ends(variants, "\\*"))
  
  if (any(str_detect(names(tbl), "_ddG"))) {
    tbl <- pivot_longer(tbl, starts_with("foldx_"), names_to = "pdb_id", values_to = "ddg", names_prefix = "foldx_") %>%
      group_by(variants) %>%
      mutate(foldx_ddg = mean(ddg)) %>%
      select(-pdb_id, -ddg) %>%
      distinct() %>%
      ungroup()
  }
  return(tbl)
}

# New SIFT4G results have more sig figs
new_preds <- read_tsv("data/long_combined_mutational_scans.tsv") %>%
  select(gene, position, wt, mut, new_sift_score = sift, new_foldx = total_energy)

dms <- readRDS('data/rdata/processed_variant_data.RDS') %>%
  lapply(combine_study) %>%
  bind_rows(.id = "study") %>%
  tidyr::extract(variants, into = c("wt", "position", "mut"), regex = "([A-Z])([0-9]*)([A-Z])", convert = TRUE) %>%
  left_join(select(study_meta, study, gene=gene_name, uniprot_id, species, authour, thresh, norm_thresh), ., by = "study") %>%
  mutate(pretty_study = map_chr(study, format_study)) %>%
  left_join(new_preds, by = c("gene", "position", "wt", "mut")) %>%
  mutate(sift_score = ifelse(is.na(new_sift_score), sift_score, new_sift_score),
         foldx_ddg = ifelse(is.na(new_foldx), foldx_ddg, new_foldx)) %>%
  select(-new_sift_score, -new_foldx)

#### Panel - Dataset Description ####
variant_counts <- group_by(dms, study = pretty_study) %>%
  summarise(variants = n(),
            SIFT4G = sum(!is.na(sift_score)),
            FoldX = sum(!is.na(foldx_ddg)),
            EVCouplings = sum(!is.na(evcoup_epistatic)),
            PolyPhen2 = sum(!is.na(envision_prediction)),
            Envision = sum(!is.na(pph2_prob))) %>%
  pivot_longer(c(-study, -variants), names_to = "tool", values_to = "count")

p_data <- ggplot(distinct(variant_counts, study, variants), aes(x = study, y = variants, label = variants)) +
  geom_col(fill = "#a6cee3") +
  geom_text(hjust = -0.25, size = 3) +
  geom_point(data = variant_counts, mapping = aes(y = count, fill = tool, shape = tool),
             position = position_dodge(0.9), size = 1.5, colour = "transparent") +
  labs(x = "", y = "Variants Measured") +
  coord_flip() +
  scale_fill_manual(name = "", values = tool_colours) +
  scale_shape_manual(name = "", values = c(SIFT4G = 21, FoldX = 22, EVCouplings = 23, PolyPhen2 = 24, Envision = 25)) + 
  scale_y_continuous(expand = expansion(c(0.01, 0.18))) +
  guides(fill = guide_legend(reverse = TRUE, override.aes = list(size = 3)),
         shape = guide_legend(reverse = TRUE)) +
  theme(panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.y = element_blank(),
        legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, "mm"),
        legend.key.size = unit(2, "mm"),
        legend.text = element_text(size = 8.5),
        axis.ticks.y = element_blank())

#### Panel - Score Correlation per study ####
correlation <- select(dms, study = pretty_study, score, SIFT4G = sift_score, FoldX = foldx_ddg, EVCouplings = evcoup_epistatic,
                      PolyPhen2 = pph2_prob, Envision = envision_prediction) %>%
  mutate(SIFT4G = -log10(SIFT4G + 0.005),
         PolyPhen2 = -log10(PolyPhen2),
         Envision = log2(Envision),
         EVCouplings = -EVCouplings) %>%
  pivot_longer(c(-study, -score), names_to = "tool", values_to = "value") %>%
  mutate(score = ifelse(tool %in% c("Envision", "PolyPhen2"), score, abs(score))) %>%
  filter(is.finite(value)) %>%
  group_by(study, tool) %>%
  group_modify(~broom::tidy(cor.test(.x$score, .x$value))) %>%
  ungroup() %>%
  mutate(p_cat = pretty_p_values(p.value, breaks = c(1e-12, 1e-6, 1e-3, 0.01, 0.05), markdown_exp = TRUE, prefix_p = TRUE))

study_order <- filter(correlation, tool == "SIFT4G") %>%
  arrange(estimate) %>%
  pull(study) %>%
  c("Ashenberg et al. 2016 (NP)", .)

p_cor <- ggplot(correlation, aes(x = factor(study, levels = study_order), y = estimate, ymin = conf.low, ymax = conf.high, fill = p_cat)) +
  facet_wrap(~tool, nrow = 1) +
  geom_col(width = 0.7) +
  geom_errorbar(width = 0.4, size = 0.25) +
  geom_hline(yintercept = 0, size = 0.3) +
  coord_flip() +
  scale_fill_viridis_d(direction = -1) +
  scale_y_continuous(breaks = seq(-0.2, 0.6, 0.2), labels = c("-0.2", "0", "0.2", "0.4", "0.6"), expand = expansion(0)) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 1)) +
  labs(x = "", y = "Pearson Correlation Coefficient") +
  theme(panel.grid.major.x = element_line(linetype = "dotted", colour = "grey", size = 0.3),
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.text = element_markdown(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, "mm"),
        legend.key.size = unit(2, "mm"))

#### Panel - ROC Curves ####
step_auc <- function(x, y){
  sum(diff(x) * y[-length(y)])
}

calc_roc <- function(tbl, true_col, var_col, greater=TRUE, max_steps = 500){
  true_col <- enquo(true_col)
  var_col <- enquo(var_col)
  
  tbl <- select(tbl, !!true_col, !!var_col) %>%
    drop_na()
  if (nrow(tbl) == 0){
    return(tibble(thresh=NA, tp=NA, tn=NA, fp=NA, fn=NA, tpr=NA, tnr=NA, fpr=NA, precision=NA, auc=NA))
  }
  
  true <- pull(tbl, !!true_col)
  var <- pull(tbl, !!var_col)
  
  unique_values <- unique(var)
  if (length(unique_values) > max_steps) {
    steps <- c(-Inf, seq(from = min(unique_values), to = max(unique_values), length.out = max_steps), Inf)
  } else {
    steps <- c(-Inf, sort(unique_values), Inf)
  }
  
  true_mat <- matrix(true, nrow = length(true), ncol = length(steps))
  var_mat <- matrix(var, nrow = length(var), ncol = length(steps))
  thresh_mat <- matrix(steps, nrow = length(var), ncol = length(steps), byrow = TRUE)
  
  if (greater){
    preds <- var_mat >= thresh_mat
  } else {
    preds <- var_mat <= thresh_mat
  }
  
  tp <- colSums(preds & true_mat)
  tn <- colSums(!preds & !true_mat)
  fp <- colSums(preds & !true_mat)
  fn <- colSums(!preds & true_mat)
  tbl <- tibble(thresh = steps, tp = tp, tn = tn, fp = fp, fn = fn,
                tpr = tp / (tp + fn),
                tnr = tn / (tn + fp),
                fpr = fp / (tn + fp),
                precision = tp / (tp + fp))
  
  if (greater){
    tbl <- arrange(tbl, desc(thresh))
  } else {
    tbl <- arrange(tbl, thresh)
  }
  
  tbl$auc <- step_auc(tbl$fpr, tbl$tpr)
  return(tbl)
}

greater <- c(SIFT4G = FALSE, FoldX = TRUE, EVCouplings = FALSE, PolyPhen2 = TRUE, Envision = FALSE)
if (file.exists("data/tool_roc.tsv")) {
  roc <- read_tsv("data/tool_roc.tsv") # Load cached result by default
} else {
  roc <- select(dms, study = pretty_study, score, thresh, SIFT4G = sift_score, FoldX = foldx_ddg, EVCouplings = evcoup_epistatic,
         PolyPhen2 = pph2_prob, Envision = envision_prediction) %>%
    mutate(del = score < thresh) %>%
    select(study, del, SIFT4G, FoldX, EVCouplings, PolyPhen2, Envision) %>%
    pivot_longer(c(-study, -del), names_to = "tool", values_to = "value") %>%
    drop_na() %>%
    {
      bind_rows(group_by(., study, tool) %>%
                  group_modify(~calc_roc(.x, del, value, greater = greater[.y$tool], max_steps = 6000)) %>%
                  ungroup(),
                group_by(., tool) %>%
                  group_modify(~calc_roc(.x, del, value, greater = greater[.y$tool], max_steps = 6000)) %>%
                  ungroup() %>%
                  mutate(study = "All"))
    }
  write_tsv(roc, "data/tool_roc.tsv")
}

p_roc <- ggplot(mapping = aes(x = fpr, y = tpr)) +
  facet_wrap(~tool, nrow = 1, scales = "free") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
  geom_step(data = filter(roc, study != "All"), mapping = aes(group = study), linetype = "dotted", colour = "darkgrey", size = 0.4) +
  geom_step(data = filter(roc, study == "All"), mapping = aes(colour = tool), show.legend = FALSE) +
  scale_colour_manual(values = tool_colours) + 
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = "FPR", y = "TPR") +
  theme(panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 8))

p_pr <- ggplot(mapping = aes(x = tpr, y = precision)) +
  facet_wrap(~tool, nrow = 1, scales = "free") +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black") +
  geom_step(data = filter(roc, study != "All"), mapping = aes(group = study), linetype = "dotted", colour = "darkgrey", size = 0.4) +
  geom_step(data = filter(roc, study == "All"), mapping = aes(colour = tool), show.legend = FALSE) +
  scale_colour_manual(values = tool_colours) + 
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  labs(x = "Recall", y = "Precision") +
  theme(panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "grey"),
        axis.text.x = element_text(size = 8))

### Panel Overall Performance ###
auc <- select(roc, study, tool, value = auc) %>%
  distinct() %>%
  mutate(stat = "auc", study = ifelse(study == "All", "Overall", study))

stats <- select(dms, study = pretty_study, thresh, score, evcoup_epistatic:foldx_ddg) %>%
  mutate(true = score < thresh,
         SIFT4G = sift_score < 0.05,
         FoldX = foldx_ddg > 1,
         PolyPhen2 = pph2_prob > 0.5,
         Envision = envision_prediction < 0.75,
         EVCouplings = evcoup_epistatic < -5) %>%
  select(study, true, SIFT4G:EVCouplings) %>%
  bind_rows(mutate(., study = "Overall"), .) %>%
  pivot_longer(SIFT4G:EVCouplings, names_to = "tool", values_to = "pred") %>%
  drop_na() %>%
  group_by(study, tool) %>%
  summarise(tp = sum(pred & true),
            tn = sum(!pred & !true),
            fp = sum(pred & !true),
            fn = sum(!pred & true),
            .groups = "drop") %>%
  mutate(accuracy = (tp + tn) / (tp + tn + fp + fn),
         precision = tp / (tp + fp),
         recall = tp / (tp + fn),
         f1 = 2 * tp / (2 * tp + fp + fn),
         kappa = 2 * (tp * tn - fn * fp) / (((as.numeric(tp) + fp) * (fp + tn)) + (tp + fn) * (fn + tn))) %>%
  select(study, tool, accuracy, precision, recall, f1, kappa) %>%
  pivot_longer(accuracy:kappa, names_to = "stat") %>%
  bind_rows(., auc) %>%
  mutate(stat = factor(stat, levels = c("accuracy", "precision", "recall", "f1", "kappa", "auc")))

stat_names <- c(accuracy = '"Accuracy"', precision = '"Precision"', recall = '"Recall"',
                f1 = '"F1 Score"', kappa = '"Cohen\'s"~kappa', auc = '"ROC AUC"')
stat_lims <- tibble(stat = c("accuracy", "accuracy", "precision", "precision", "recall", "recall", "f1", "f1", "kappa", "kappa", "auc", "auc"),
                    value = c(0, 1, 0, 1, 0, 1, 0, 1, -0.04, 1, 0, 1),
                    tool = NA) %>%
  mutate(stat = factor(stat, levels = c("accuracy", "precision", "recall", "f1", "kappa", "auc")))

p_stats <- ggplot(mapping = aes(x = tool, y = value)) +
  facet_wrap(vars(stat), nrow = 2, scales = "fixed", strip.position = "left",
             labeller = labeller(stat = as_labeller(stat_names, default = label_parsed))) +
  geom_boxplot(data = filter(stats, study != "Overall"), mapping = aes(fill = tool), width = 0.75, show.legend = FALSE, outlier.shape = 20) +
  geom_point(data = filter(stats, study == "Overall"), mapping = aes(shape = "All Studies"), size = 2.5) +
  geom_blank(data = stat_lims) +
  scale_fill_manual(values = tool_colours) + 
  scale_shape_manual(name = "", values = c(`All Studies` = 10)) +
  labs(x = "", y = "") + 
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(2, "mm"),
        strip.placement = "outside",
        strip.text = element_markdown(margin = margin(0, 0, 0, 0, unit = "mm")))

### ROC Figure Assembly
size <- theme(text = element_text(size = 12))
pr1 <- p_stats + labs(tag = 'A') + size
pr2 <- p_roc + labs(tag = 'B') + size
pr3 <- p_pr + labs(tag = 'C') + size

figure2 <- multi_panel_figure(width = 180, height = c(100, 40, 40), columns = 1,
                              panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(pr1, row = 1, column = 1) %>%
  fill_panel(pr2, row = 2, column = 1) %>%
  fill_panel(pr3, row = 3, column = 1)

ggsave('figures/thesis_figure_roc.pdf', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_roc.tiff', figure2, width = figure_width(figure2), height = figure_height(figure2), units = 'mm')

# Individual figures
ggsave('figures/thesis_figure_data.pdf', p_data + size, width = 180, height = 180, units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure_cor.pdf', p_cor + size, width = 220, height = 220, units = 'mm', device = cairo_pdf)

