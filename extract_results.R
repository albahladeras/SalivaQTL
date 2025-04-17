library(dplyr)
library(readr)
library(stringr)

results_dir <- "results/outputs"
share_files <- list.files(results_dir, pattern = "\\.share$", full.names = TRUE)
traits <- str_replace(basename(share_files), "\\.share$", "")

ci_bounds <- function(est, se) {
  lower <- est - 1.96 * se
  upper <- est + 1.96 * se
  c(lower = round(lower, 5), upper = round(upper, 5))
}

results_table <- data.frame()

for (trait in traits) {
  share_path <- file.path(results_dir, paste0(trait, ".share"))
  hers_path <- file.path(results_dir, paste0(trait, ".hers.liab"))
  overlap_path <- file.path(results_dir, paste0(trait, ".overlap"))

  if (!file.exists(share_path) | !file.exists(hers_path) | !file.exists(overlap_path)) next

  share_df <- read_table(share_path, col_names = c("Component", "Share", "SE"), skip = 1)
  hers_df <- read_table(hers_path, col_names = c("Component", "Heritability", "SE", "Influence", "Influence_SE"), skip = 1)
  predictors <- as.numeric(str_split(readLines(overlap_path)[2], "\\s+")[[1]][2])

  h2_total <- hers_df$Heritability[hers_df$Component == "Her_All"]
  h2_mqtl  <- hers_df$Heritability[hers_df$Component == "Her_A1"]
  h2_eqtl  <- hers_df$Heritability[hers_df$Component == "Her_A2"]

  se_total <- hers_df$SE[hers_df$Component == "Her_All"]
  se_mqtl  <- hers_df$SE[hers_df$Component == "Her_A1"]
  se_eqtl  <- hers_df$SE[hers_df$Component == "Her_A2"]

  share_mqtl <- share_df$Share[share_df$Component == "Share_A1"]
  share_eqtl <- share_df$Share[share_df$Component == "Share_A2"]
  share_combined <- share_mqtl + share_eqtl

  se_share_mqtl <- share_df$SE[share_df$Component == "Share_A1"]
  se_share_eqtl <- share_df$SE[share_df$Component == "Share_A2"]

  results_table <- bind_rows(results_table, tibble(
    Trait = trait,
    predictors = predictors,
    h2_liability = round(h2_total, 5),
    h2_liab_CI_lower = ci_bounds(h2_total, se_total)["lower"],
    h2_liab_CI_upper = ci_bounds(h2_total, se_total)["upper"],
    h2_mQTL = round(h2_mqtl, 5),
    h2_mQTL_CI_lower = ci_bounds(h2_mqtl, se_mqtl)["lower"],
    h2_mQTL_CI_upper = ci_bounds(h2_mqtl, se_mqtl)["upper"],
    h2_eQTL = round(h2_eqtl, 5),
    h2_eQTL_CI_lower = ci_bounds(h2_eqtl, se_eqtl)["lower"],
    h2_eQTL_CI_upper = ci_bounds(h2_eqtl, se_eqtl)["upper"],
    pct_h2_mQTL = round(100 * share_mqtl, 1),
    pct_mQTL_CI_lower = round(100 * ci_bounds(share_mqtl, se_share_mqtl)["lower"], 1),
    pct_mQTL_CI_upper = round(100 * ci_bounds(share_mqtl, se_share_mqtl)["upper"], 1),
    pct_h2_eQTL = round(100 * share_eqtl, 1),
    pct_eQTL_CI_lower = round(100 * ci_bounds(share_eqtl, se_share_eqtl)["lower"], 1),
    pct_eQTL_CI_upper = round(100 * ci_bounds(share_eqtl, se_share_eqtl)["upper"], 1),
    pct_combined = round(100 * share_combined, 1)
  ))
}

write_csv(results_table, "ldak_qtl_heritability_summary.csv")
