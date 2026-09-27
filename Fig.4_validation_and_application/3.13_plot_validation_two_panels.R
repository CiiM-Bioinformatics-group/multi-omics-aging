# 3.13 Model comparison in 500FG (training / internal test) and 300BCG (external validation, from 3.12)
#     Paired t-tests (paired by iteration): TabPFN vs each other model, BH-adjusted per panel.
#      Run from the repository root: Rscript validation/<this script>

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggpubr)
})

# ---- Paths ----
model_dir <- "results/tabpfn_300bcg"
fig_dir   <- "results/figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

tabpfn <- read.csv(file.path(model_dir, "model_results_300bcg_external_with_methy_share_feattop235.csv"))

# ---- Settings ----
model_levels   <- c("TabPFN", "ElasticNet", "RandomForest", "XGBoost", "CatBoost")
panel_levels   <- c("500FG", "300BCG")
panel_f_levels <- c("500FG training", "300BCG external validation")
panel_label    <- function(p) {
  factor(case_when(p == "500FG" ~ "500FG training",
                   p == "300BCG" ~ "300BCG external validation"), levels = panel_f_levels)
}

# ---- Long table ----
plot_df_2panel <- bind_rows(
  tabpfn %>% transmute(Panel = "500FG",  Model, Iteration, R2 = Train_R2,    SetType = "Training"),
  tabpfn %>% transmute(Panel = "500FG",  Model, Iteration, R2 = Test_R2,     SetType = "Internal test"),
  tabpfn %>% transmute(Panel = "300BCG", Model, Iteration, R2 = External_R2, SetType = "External test")
) %>%
  mutate(
    Model   = factor(as.character(Model), levels = model_levels),
    Panel   = factor(Panel, levels = panel_levels),
    SetType = factor(SetType, levels = c("Training", "Internal test", "External test")),
    Panel_f = panel_label(as.character(Panel)),
    FillKey = factor(case_when(
      Panel == "500FG" & SetType == "Training"      ~ "Training",
      Panel == "500FG" & SetType == "Internal test" ~ "Internal test",
      TRUE                                          ~ "External test"
    ), levels = c("Training", "Internal test", "External test"))
  ) %>%
  filter(!is.na(R2), Model %in% model_levels, !is.na(Panel))

# ---- Paired t-tests (500FG internal test and 300BCG external test) ----
wide_by_panel <- plot_df_2panel %>%
  filter((Panel == "500FG" & SetType == "Internal test") | Panel == "300BCG") %>%
  distinct(Panel, Iteration, Model, R2) %>%
  pivot_wider(names_from = Model, values_from = R2)

pval_df_2panel <- wide_by_panel %>%
  group_split(Panel) %>%
  map_dfr(function(w) {
    pn <- as.character(unique(w$Panel))
    others <- setdiff(names(w), c("Panel", "Iteration", "TabPFN"))
    map_dfr(others, function(m) {
      d <- w %>% select(Iteration, TabPFN, all_of(m)) %>% drop_na()
      p <- if (nrow(d) < 2) NA_real_ else t.test(d$TabPFN, d[[m]], paired = TRUE)$p.value
      tibble(Panel = pn, group1 = "TabPFN", group2 = m, p = p)
    })
  }) %>%
  group_by(Panel) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    Panel  = factor(Panel, levels = panel_levels),
    group1 = factor(group1, levels = model_levels),
    group2 = factor(group2, levels = model_levels),
    p.adj.signif = symnum(p.adj, corr = FALSE, na = FALSE,
                          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                          symbols = c("***", "**", "*", ".", "ns")),
    Panel_f = panel_label(as.character(Panel))
  )

# Bracket positions above each panel's data
panel_rng <- plot_df_2panel %>%
  group_by(Panel) %>%
  summarise(ymax = max(R2, na.rm = TRUE),
            span = pmax(max(R2, na.rm = TRUE) - min(R2, na.rm = TRUE), 1e-6),
            .groups = "drop")

pval_df_2panel <- pval_df_2panel %>%
  left_join(panel_rng, by = "Panel") %>%
  group_by(Panel) %>%
  arrange(group2, .by_group = TRUE) %>%
  mutate(y.position = ymax + span * (0.05 + 0.04 * row_number())) %>%
  ungroup()
print(pval_df_2panel)

# ---- Plot (Nature style, 1 x 2) ----
axis_text_pt  <- 5
axis_title_pt <- 5
strip_text_pt <- 5

fill_manual <- c("Training" = "#7EB6D9", "Internal test" = "#82C09A", "External test" = "#E8A8A8")
y_top <- max(c(plot_df_2panel$R2, pval_df_2panel$y.position), na.rm = TRUE) * 1.18

p2 <- ggplot(plot_df_2panel, aes(x = Model, y = R2, fill = FillKey)) +
  geom_boxplot(outlier.shape = NA, colour = "black", linewidth = 0.35, alpha = 0.85,
               position = position_dodge(width = 0.78)) +
  geom_jitter(aes(color = FillKey),
              position = position_jitterdodge(jitter.width = 0.12, jitter.height = 0, dodge.width = 0.78),
              alpha = 0.7, size = 0.65, show.legend = FALSE) +
  facet_wrap(~Panel_f, nrow = 1, ncol = 2, scales = "fixed") +
  scale_fill_manual(name = NULL, values = fill_manual) +
  scale_colour_manual(values = fill_manual, guide = "none") +
  scale_y_continuous(limits = c(0, y_top), breaks = seq(0, 1, 0.2),
                     labels = function(x) ifelse(x <= 1, x, ""),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(x = "Model", y = expression(Test ~ R^2)) +
  theme_minimal(base_family = "sans") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = axis_text_pt, colour = "black"),
    axis.text.y = element_text(size = axis_text_pt, colour = "black"),
    axis.title = element_text(size = axis_title_pt, colour = "black"),
    strip.text = element_text(size = strip_text_pt, colour = "black", lineheight = 0.95),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.box.spacing = unit(0.05, "cm"),
    legend.margin = margin(t = 0, b = 0),
    legend.key.size = unit(0.30, "cm"),
    legend.text = element_text(size = axis_text_pt, colour = "black"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35)
  ) +
  stat_pvalue_manual(data = pval_df_2panel, label = "p.adj.signif",
                     xmin = "group1", xmax = "group2", y.position = "y.position",
                     facet.by = "Panel_f", inherit.aes = FALSE,
                     tip.length = 0.007, bracket.size = 0.24, step.increase = 0.02,
                     hide.ns = FALSE, size = 2.2)

# ---- Save (95 mm x 60 mm) ----
ggsave(file.path(fig_dir, "validation_two_panels_nature.pdf"), p2,
       width = 95 / 25.4, height = 60 / 25.4, device = "pdf", useDingbats = FALSE)
write.csv(pval_df_2panel, file.path(fig_dir, "validation_two_panels_paired_ttests.csv"), row.names = FALSE)
