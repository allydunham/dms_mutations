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
new_sift <- read_tsv("data/long_combined_mutational_scans.tsv") %>%
  select(gene, position, wt, mut, new_sift_score = sift)

dms <- readRDS('data/rdata/processed_variant_data.RDS') %>%
  lapply(combine_study) %>%
  bind_rows(.id = "study") %>%
  tidyr::extract(variants, into = c("wt", "position", "mut"), regex = "([A-Z])([0-9]*)([A-Z])", convert = TRUE) %>%
  left_join(select(study_meta, study, gene=gene_name, uniprot_id, species, authour, thresh, norm_thresh), ., by = "study") %>%
  mutate(pretty_study = map_chr(study, format_study)) %>%
  left_join(new_sift, by = c("gene", "position", "wt", "mut")) %>%
  mutate(sift_score = ifelse(is.na(new_sift_score), sift_score, new_sift_score)) %>%
  select(-new_sift_score)

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
  geom_text(hjust = -0.25, size = 2) +
  geom_point(data = variant_counts, mapping = aes(y = count, colour = tool),
             position = position_dodge(0.9), shape = 20, size = 0.75) +
  labs(x = "", y = "Variants") +
  coord_flip() +
  scale_colour_brewer(name = "", type = "qual", palette = "Set1") +
  scale_y_continuous(expand = expansion(c(0.01, 0.18))) +
  guides(colour = guide_legend(reverse = TRUE, override.aes = list(size = 2))) +
  theme(panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.y = element_blank(),
        legend.position = "right",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, "mm"),
        legend.key.size = unit(2, "mm"),
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

p_cor <- ggplot(correlation, aes(x = study, y = estimate, ymin = conf.low, ymax = conf.high, fill = p_cat)) +
  facet_wrap(~tool, nrow = 1) +
  geom_col(width = 0.8) +
  geom_errorbar(width = 0.6) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_fill_viridis_d() +
  guides(fill = guide_legend(reverse = TRUE, nrow = 1)) +
  labs(x = "", y = "Pearson Correlation Coefficient") +
  theme(panel.grid.major.x = element_line(linetype = "dotted", colour = "grey"),
        panel.grid.major.y = element_blank(),
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

auc <- select(roc, study, tool, auc) %>%
  distinct()

p_roc <- ggplot(mapping = aes(x = fpr, y = tpr)) +
  facet_wrap(~tool, nrow = 2, scales = "free") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "black") +
  geom_step(data = filter(roc, study != "All"), mapping = aes(group = study), linetype = "dotted", colour = "darkgrey") +
  geom_step(data = filter(roc, study == "All"), mapping = aes(colour = tool), show.legend = FALSE) +
  scale_colour_brewer(palette = "Dark2") + 
  labs(x = "FPR", y = "TPR") +
  theme(panel.grid.major.y = element_blank(),
        axis.line = element_line(colour = "grey"))

p_roc_summary <- ggplot(mapping = aes(x = tool, y = auc)) +
  geom_boxplot(data = filter(auc, study != "All"), mapping = aes(fill = tool), show.legend = FALSE, outlier.shape = 20) +
  geom_point(data = filter(auc, study == "All"), mapping = aes(shape = "All Studies"), colour = "black", size = 3) +
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(name = "", values = c(`All Studies` = 10)) +
  labs(x = "", y = "ROC AUC") + 
  lims(y = c(0, 1)) +
  theme(legend.position = "top",
        legend.margin = margin(0,0,0,0),
        legend.box.spacing = unit(1, "mm"),
        legend.key.size = unit(2, "mm"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#### Figure Assembly ####
size <- theme(text = element_text(size = 9))
p1 <- p_data + labs(tag = 'A') + size
p2 <- p_cor + labs(tag = 'B') + size
p3 <- p_roc + labs(tag = 'C') + size
p4 <- p_roc_summary + labs(tag = 'D') + size

figure <- multi_panel_figure(width = 180, height = 240, columns = 3, rows = 3,
                             panel_label_type = 'none', row_spacing = 0, column_spacing = 0) %>%
  fill_panel(p1, row = 1, column = 1:3) %>%
  fill_panel(p2, row = 2, column = 1:3) %>%
  fill_panel(p3, row = 3, column = 1:2) %>%
  fill_panel(p4, row = 3, column = 3)

ggsave('figures/thesis_figure.pdf', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm', device = cairo_pdf)
ggsave('figures/thesis_figure.tiff', figure, width = figure_width(figure), height = figure_height(figure), units = 'mm')
