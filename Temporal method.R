# ========== Load and install necessary packages ==========
packages <- c("hdf5r", "dplyr", "tidyr", "data.table", "ggplot2", "RColorBrewer", "patchwork", "readr")
invisible(lapply(packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}))

# ========== Load metadata ==========
metadata <- read.csv('/Users/kunyaoxu/Desktop/BFgam_metadata.csv', header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(time_point = paste0(year, "-", sprintf("%02d", month)))

# ========== Function: Process HDF5 and estimate Ne ==========
process_chromosome <- function(hdf5_path, output_path) {
  dat <- H5File$new(hdf5_path, 'r')
  genotypes <- dat[["genotype"]]$read()
  sample_id <- dat[["sample_id"]]$read()
  genotypes_merged <- genotypes[1,,] + genotypes[2,,]
  rownames(genotypes_merged) <- sample_id
  positions <- dat[["POS"]]$read()
  time_points <- unique(metadata$time_point)
  time_pairs <- combn(time_points, 2, simplify = FALSE)
  results <- list()
  
  filter_by_distance <- function(positions, min_dist = 1000) {
    keep <- logical(length(positions))
    keep[1] <- TRUE
    last_pos <- positions[1]
    for (i in 2:length(positions)) {
      if ((positions[i] - last_pos) >= min_dist) {
        keep[i] <- TRUE
        last_pos <- positions[i]
      }
    }
    return(keep)
  }
  
  for (pair in time_pairs) {
    time0 <- pair[1]; time1 <- pair[2]
    samples0 <- metadata %>% filter(time_point == time0) %>% pull(sample_id)
    samples1 <- metadata %>% filter(time_point == time1) %>% pull(sample_id)
    S0 <- length(samples0); St <- length(samples1)
    if (S0 < 5 || St < 5) next
    
    geno0 <- genotypes_merged[samples0, , drop = FALSE]
    geno1 <- genotypes_merged[samples1, , drop = FALSE]
    
    allele_freq0 <- colMeans(geno0, na.rm = TRUE) / 2
    allele_freq1 <- colMeans(geno1, na.rm = TRUE) / 2
    
    maf0 <- pmin(allele_freq0, 1 - allele_freq0)
    maf1 <- pmin(allele_freq1, 1 - allele_freq1)
    valid_snps <- which(maf0 >= 0.05 & maf1 >= 0.05)
    if (length(valid_snps) < 100) next
    
    allele_freq0 <- allele_freq0[valid_snps]
    allele_freq1 <- allele_freq1[valid_snps]
    pos_valid <- positions[valid_snps]
    
    ord <- order(pos_valid)
    pos_sorted <- pos_valid[ord]
    keep_flags <- filter_by_distance(pos_sorted)
    keep_idx <- order(valid_snps)[ord[keep_flags]]
    allele_freq0 <- allele_freq0[keep_idx]
    allele_freq1 <- allele_freq1[keep_idx]
    
    year0 <- as.numeric(substr(time0, 1, 4))
    month0 <- as.numeric(substr(time0, 6, 7))
    year1 <- as.numeric(substr(time1, 1, 4))
    month1 <- as.numeric(substr(time1, 6, 7))
    t <- (year1 - year0) * 12 + (month1 - month0)
    if (t <= 0) next
    
    numerator <- (allele_freq0 - allele_freq1)^2
    denominator <- ((allele_freq0 + allele_freq1)/2) - (allele_freq0 * allele_freq1)
    Fc_values <- numerator / denominator
    Fc_hat <- mean(Fc_values, na.rm = TRUE)
    K <- length(Fc_values)
    chi_upper <- qchisq(0.975, df = K)
    chi_lower <- qchisq(0.025, df = K)
    Fc_lower <- (K * Fc_hat) / chi_upper
    Fc_upper <- (K * Fc_hat) / chi_lower
    
    denom_ne <- Fc_hat - 1/(2*S0) - 1/(2*St)
    denom_lower <- Fc_upper - 1/(2*S0) - 1/(2*St)
    denom_upper <- Fc_lower - 1/(2*S0) - 1/(2*St)
    Ne <- ifelse(denom_ne > 0, 1 / (2 * denom_ne), NA)
    Ne_lower <- ifelse(denom_upper > 0, 1 / (2 * denom_upper), NA)
    Ne_upper <- ifelse(denom_lower > 0, 1 / (2 * denom_lower), NA)
    
    results[[length(results) + 1]] <- data.frame(
      time0 = time0, timet = time1,
      S0 = S0, St = St, t = t, K = K,
      Fc = round(Fc_hat, 4),
      Fc_lower = round(Fc_lower, 4),
      Fc_upper = round(Fc_upper, 4),
      Ne = Ne, Ne_lower = Ne_lower, Ne_upper = Ne_upper
    )
  }
  
  final_results <- bind_rows(results) %>% arrange(time0, timet)
  
  # ========= 99% CI Calculation =========
  K_eff <- round(final_results$K / 6)
  Fc_lower99 <- K_eff * final_results$Fc / qchisq(0.995, df = K_eff)
  Fc_upper99 <- K_eff * final_results$Fc / qchisq(0.005, df = K_eff)
  Ne_point <- final_results$t / (2 * (final_results$Fc - 0.5 / final_results$S0 - 0.5 / final_results$St))
  Ne_upper99 <- final_results$t / (2 * (Fc_lower99 - 0.5 / final_results$S0 - 0.5 / final_results$St))
  Ne_lower99 <- final_results$t / (2 * (Fc_upper99 - 0.5 / final_results$S0 - 0.5 / final_results$St))
  
  dat_out <- data.frame(
    time0 = final_results$time0,
    timet = final_results$timet,
    S0 = final_results$S0,
    St = final_results$St,
    t = final_results$t,
    K = final_results$K,
    Fc = final_results$Fc,
    Fc_lower = round(Fc_lower99, 4),
    Fc_upper = round(Fc_upper99, 4),
    Ne = Ne_point,
    Ne_lower = Ne_lower99,
    Ne_upper = Ne_upper99
  )
  
  fwrite(dat_out, output_path, na = "NA")
  cat("\u2705 Results with 99% CI saved to", output_path, "\n")
  
  if (!is.null(dat) && dat$is_valid) dat$close_all()
}

# ========== Run for both 2R and 3R ==========
process_chromosome(
  hdf5_path = "/Users/kunyaoxu/Desktop/BFgam_2R_biallelic_CDS4.hdf5",
  output_path = "/Users/kunyaoxu/Desktop/Ne_estimates_2R_with99CI.csv"
)

process_chromosome(
  hdf5_path = "/Users/kunyaoxu/Desktop/BFgam_3R_biallelic_CDS4.hdf5",
  output_path = "/Users/kunyaoxu/Desktop/Ne_estimates_3R_with99CI.csv"
)

# ========== Cross-validation 2R vs 3R ==========
dat2R <- read.csv("/Users/kunyaoxu/Desktop/Ne_estimates_2R_with99CI.csv")
dat3R <- read.csv("/Users/kunyaoxu/Desktop/Ne_estimates_3R_with99CI.csv")
plot(dat2R$Fc - 0.5 / dat2R$S0 - 0.5 / dat2R$St, dat3R$Fc - 0.5 / dat3R$S0 - 0.5 / dat3R$St,
     xlab = "Fc minus sampling (2R)", ylab = "Fc minus sampling (3R)")
abline(a = 0, b = 1, col = "red", lty = 2)
print(cor(dat2R$Fc - 0.5 / dat2R$S0 - 0.5 / dat2R$St, dat3R$Fc - 0.5 / dat3R$S0 - 0.5 / dat3R$St))

get_season <- function(month) {
  if (month %in% 6:10) return("Wet")
  if (month %in% c(11, 12, 1, 2, 3, 4)) return("Dry")
  return(NA)
}
dat1$month_timet <- as.numeric(substr(dat1$timet, 6, 7))
dat1$season <- sapply(dat1$month_timet, get_season)
season_colors <- c(Wet = "blue", Dry = "red")
df_plot <- dat1
df_plot$season <- factor(df_plot$season, levels = names(season_colors))

# ========== Overall regression plot ==========
lm_model <- lm(Fc ~ t, data = df_plot)
p_val <- summary(lm_model)$coefficients["t", "Pr(>|t|)"]
p_text <- ifelse(is.na(p_val), "p = NA", ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", formatC(p_val, format = "e", digits = 2))))

K_eff <- round(df_plot$K / 6)
df_plot$Fc_lower99 <- K_eff * df_plot$Fc / qchisq(0.995, df = K_eff)
df_plot$Fc_upper99 <- K_eff * df_plot$Fc / qchisq(0.005, df = K_eff)

overall_plot <- ggplot(df_plot, aes(x = t, y = Fc, color = season)) +
  geom_point(alpha = 0.8, size = 2.5) +
  geom_errorbar(aes(ymin = Fc_lower99, ymax = Fc_upper99), width = 0.5, alpha = 0.3) +
  scale_color_manual(values = season_colors) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 1) +  # <- changed
  annotate("text", x = max(df_plot$t, na.rm = TRUE) * 0.7, y = max(df_plot$Fc, na.rm = TRUE) * 0.95,
           label = paste("Regression:", p_text), size = 5, hjust = 0, fontface = "italic") +
  labs(title = "Fc vs. t (2R)", x = "t (months)", y = "Fc", color = "Season") +
  theme_minimal(base_size = 14)

print(overall_plot)
# ========== Per-time0 subplots ==========
x_max <- max(df_plot$t, na.rm = TRUE)
y_max <- max(df_plot$Fc_upper99, na.rm = TRUE)
group_levels <- levels(as.factor(df_plot$time0))
plot_list <- lapply(group_levels, function(g) {
  sub_df <- subset(df_plot, time0 == g)
  sub_df <- sub_df[complete.cases(sub_df$t, sub_df$Fc), ]
  if (nrow(sub_df) < 2 || length(unique(sub_df$t)) < 2) return(NULL)
  lm_model <- lm(Fc ~ t, data = sub_df)
  pval <- summary(lm_model)$coefficients["t", "Pr(>|t|)"]
  p_text <- if (is.na(pval)) "p = NA" else if (pval < 0.001) "p < 0.001" else paste0("p = ", formatC(pval, format = "e", digits = 2))
  ggplot(sub_df, aes(x = t, y = Fc, color = season)) +
    geom_point(size = 2.5, alpha = 0.8) +
    geom_errorbar(aes(ymin = Fc_lower99, ymax = Fc_upper99), width = 0.5, alpha = 0.4) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 1) +
    scale_color_manual(values = season_colors) +
    annotate("text", x = x_max * 0.65, y = y_max * 0.92, label = paste("Regression:", p_text), size = 4.2, hjust = 0, fontface = "italic") +
    labs(title = paste("Fc vs. t | time0 =", g), x = "t (months)", y = "Fc") +
    xlim(0, x_max) +
    ylim(0, y_max) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
})
plot_list <- Filter(Negate(is.null), plot_list)
if (length(plot_list) > 0) {
  combined_plot <- wrap_plots(plot_list, ncol = 2)
  print(combined_plot)
} else {
  message("No valid time0 groups for plotting.")
}