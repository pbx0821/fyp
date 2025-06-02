# Required libraries
if (!requireNamespace("hdf5r", quietly = TRUE)) install.packages("hdf5r")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(hdf5r)
library(dplyr)
library(data.table)
library(ggplot2)

# Clear memory
rm(list = ls())
gc()

# Output path
output_path <- "/Users/kunyaoxu/Desktop/LD_summary_2R_3R_before.csv"
if (file.exists(output_path)) file.remove(output_path)

# Read data
dat_2R <- H5File$new('/Users/kunyaoxu/Desktop/BFgam_2R_biallelic_CDS4.hdf5', 'r')
dat_3R <- H5File$new('/Users/kunyaoxu/Desktop/BFgam_3R_biallelic_CDS4.hdf5', 'r')
metadata <- read.csv('/Users/kunyaoxu/Desktop/BFgam_metadata.csv', header = TRUE, stringsAsFactors = FALSE)
metadata$year_month <- paste0(metadata$year, "_", sprintf("%02d", metadata$month))
gt_sample_id <- dat_2R[['sample_id']][]

# Merge months
merge_groups <- list("2015_01/02" = c("2015_01", "2015_02"))
merged_lookup <- unlist(lapply(names(merge_groups), function(name) {
  setNames(rep(name, length(merge_groups[[name]])), merge_groups[[name]])
}))
metadata$group_label <- ifelse(
  metadata$year_month %in% names(merged_lookup),
  merged_lookup[metadata$year_month],
  metadata$year_month
)
group_list <- unique(metadata$group_label)

summary_table <- data.frame(
  group_label = character(),
  sample_size = integer(),
  raw_r2 = numeric(),
  adjusted_r2 = numeric(),
  num_poly_2R = numeric(),
  stringsAsFactors = FALSE
)

for (grp in group_list) {
  cat("Processing", grp, "\n")
  subset_sample_id <- metadata[metadata$group_label == grp, "sample_id"]
  matched_idx <- gt_sample_id %in% subset_sample_id
  sample_size <- sum(matched_idx)
  
  if (sample_size < 10) {
    summary_table <- rbind(summary_table, data.frame(
      group_label = grp, sample_size = sample_size,
      raw_r2 = NA, adjusted_r2 = NA, num_poly_2R = NA
    ))
    next
  }
  
  gt_2R_raw <- dat_2R[['genotype']][, matched_idx, , drop = FALSE]
  gt_3R_raw <- dat_3R[['genotype']][, matched_idx, , drop = FALSE]
  
  gt_2R <- gt_2R_raw[1,,] + gt_2R_raw[2,,]
  gt_3R <- gt_3R_raw[1,,] + gt_3R_raw[2,,]
  
  non_missing_2R <- apply(gt_2R, 2, function(x) all(x >= 0))
  non_missing_3R <- apply(gt_3R, 2, function(x) all(x >= 0))
  
  maf_fun <- function(x) {
    af <- mean(x)
    p <- af / 2
    if (p > 0.5) p <- 1 - p
    return(p)
  }
  
  maf_2R <- apply(gt_2R, 2, maf_fun)
  maf_3R <- apply(gt_3R, 2, maf_fun)
  
  keep_2R <- which(non_missing_2R & maf_2R >= 0.2)
  keep_3R <- which(non_missing_3R & maf_3R >= 0.2)
  
  gt_2R <- gt_2R[, keep_2R, drop = FALSE]
  gt_3R <- gt_3R[, keep_3R, drop = FALSE]
  
  num_poly_2R <- ncol(gt_2R)
  if (num_poly_2R == 0 || ncol(gt_3R) == 0) {
    summary_table <- rbind(summary_table, data.frame(
      group_label = grp, sample_size = sample_size,
      raw_r2 = NA, adjusted_r2 = NA, num_poly_2R = num_poly_2R
    ))
    next
  }
  
  gc()
  batch_size <- 100
  n_batches <- ceiling(ncol(gt_2R) / batch_size)
  s <- nrow(gt_2R)
  result_LD <- numeric(1e6)
  valid_idx <- 0
  
  for (batch in 1:n_batches) {
    idx_start <- (batch - 1) * batch_size + 1
    idx_end <- min(batch * batch_size, ncol(gt_2R))
    gt_2R_batch <- gt_2R[, idx_start:idx_end, drop = FALSE]
    
    dpA <- colSums(gt_2R_batch) / (2 * s)
    dpB <- colSums(gt_3R) / (2 * s)
    dhA <- colSums(gt_2R_batch == 2) / s
    dhB <- colSums(gt_3R == 2) / s
    
    delta_AB_mat <- crossprod(gt_2R_batch, gt_3R) / (2 * s)
    val_matrix <- delta_AB_mat - 2 * outer(dpA, dpB)
    val_matrix <- val_matrix * s / (s - 1)
    denom_matrix <- outer(dpA - 2 * dpA^2 + dhA, dpB - 2 * dpB^2 + dhB)
    valid_mask <- denom_matrix > 0
    
    r2_matrix <- matrix(NA_real_, nrow = nrow(val_matrix), ncol = ncol(val_matrix))
    r2_matrix[valid_mask] <- (val_matrix[valid_mask]^2) / denom_matrix[valid_mask]
    
    result_LD_batch <- as.vector(r2_matrix)
    result_LD_batch <- result_LD_batch[!is.na(result_LD_batch) & !is.nan(result_LD_batch)]
    
    if (length(result_LD_batch) > 0) {
      result_LD[(valid_idx + 1):(valid_idx + length(result_LD_batch))] <- result_LD_batch
      valid_idx <- valid_idx + length(result_LD_batch)
    }
    gc()
  }
  
  result_LD <- result_LD[1:valid_idx]
  if (length(result_LD) == 0) {
    summary_table <- rbind(summary_table, data.frame(
      group_label = grp, sample_size = sample_size,
      raw_r2 = NA, adjusted_r2 = NA, num_poly_2R = num_poly_2R
    ))
    next
  }
  
  raw_r2 <- round(mean(result_LD), 4)
  adjusted_r2 <- round(raw_r2 - (1 / s) - (3.19 / s^2), 4)
  
  summary_table <- rbind(summary_table, data.frame(
    group_label = grp, sample_size = sample_size,
    raw_r2 = raw_r2, adjusted_r2 = adjusted_r2, num_poly_2R = num_poly_2R
  ))
}

# Add season info
get_year <- function(x) as.numeric(sub("_.*", "", x))
get_first_month <- function(x) {
  m <- sub(".*_", "", x)
  m <- sub("[/-].*", "", m)
  as.numeric(m)
}

summary_table <- summary_table %>%
  mutate(
    year_part = get_year(group_label),
    month_part = get_first_month(group_label),
    sort_key = year_part * 100 + month_part,
    season = case_when(
      month_part %in% c(6, 7, 8, 9) ~ "Wet",
      month_part %in% c(11, 12, 1, 2, 3, 4) ~ "Dry",
      TRUE ~ "Transition"
    )
  ) %>%
  arrange(sort_key) %>%
  select(group_label, season, sample_size, num_poly_2R, raw_r2, adjusted_r2)

# Save CSV
fwrite(summary_table, output_path)
cat("\u2705 Done! Summary saved to:", output_path, "\n")

# Plotting
season_colors <- c("Wet" = "blue", "Dry" = "red", "Transition" = "black")

p1 <- ggplot(summary_table %>% filter(!is.na(adjusted_r2)),
             aes(x = group_label, y = adjusted_r2, color = season)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_size), vjust = -1.2, size = 3.2, color = "black") +
  scale_color_manual(values = season_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Year-Month",
    y = expression(Adjusted~r^2),
    color = "Season",
    title = "Adjusted r² Before Filtering"
  )

p2 <- ggplot(summary_table %>% filter(!is.na(raw_r2)),
             aes(x = group_label, y = raw_r2, color = season)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_size), vjust = -1.2, size = 3.2, color = "black") +
  scale_color_manual(values = season_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Year-Month",
    y = expression(Raw~r^2),
    color = "Season",
    title = "Raw r² Before Filtering"
  )

print(p1)
print(p2)

# LD after filter
# ============= SETUP =============
rm(list = ls())
gc()

library(hdf5r)
library(matrixStats)
library(inline)
library(data.table)
library(dplyr)
library(ggplot2)

output_path <- "/Users/kunyaoxu/Desktop/LD_summary_2R_3R_after.csv"
if (file.exists(output_path)) file.remove(output_path)

dat_2R <- H5File$new('/Users/kunyaoxu/Desktop/BFgam_2R_biallelic_CDS4.hdf5', 'r')
dat_3R <- H5File$new('/Users/kunyaoxu/Desktop/BFgam_3R_biallelic_CDS4.hdf5', 'r')
metadata <- read.csv('/Users/kunyaoxu/Desktop/BFgam_metadata.csv', header = TRUE, stringsAsFactors = FALSE)
metadata$year_month <- paste0(metadata$year, "_", sprintf("%02d", metadata$month))
gt_sample_id <- dat_2R[['sample_id']][]

merge_groups <- list("2015_01/02" = c("2015_01", "2015_02"))
merged_lookup <- unlist(lapply(names(merge_groups), function(name) {
  setNames(rep(name, length(merge_groups[[name]])), merge_groups[[name]])
}))
metadata$group_label <- ifelse(
  metadata$year_month %in% names(merged_lookup),
  merged_lookup[metadata$year_month],
  metadata$year_month
)
group_list <- unique(metadata$group_label)

body_jackknife <- "
    int len = LENGTH(x);
    SEXP output = PROTECT(Rf_allocVector(REALSXP, len));
    double *d_x = REAL(x);
    double *d_output = REAL(output);
    double ss = 0; 
    double s = (double) len - 1;
    for (int i = 0; i < len; i++) {
        ss = ss + d_x[i];
    }
    for (int i = 0; i < len; i++) {
        d_output[i] = (ss - d_x[i]) / s;
    }
    UNPROTECT(1);
    return output;
"
jackknife_C <- cfunction(sig = c(x = 'array'), body = body_jackknife, language = 'C', convention = '.Call')

filter_snps_by_distance <- function(positions, min_dist = 10000) {
  keep <- logical(length(positions))
  keep[1] <- TRUE
  last_pos <- positions[1]
  for (i in 2:length(positions)) {
    if ((positions[i] - last_pos) >= min_dist) {
      keep[i] = TRUE
      last_pos = positions[i]
    }
  }
  return(keep)
}

summary_table <- data.frame(
  group_label = character(),
  sample_size = integer(),
  raw_r2 = numeric(),
  adjusted_r2 = numeric(),
  jackknife_mean = numeric(),
  jackknife_variance = numeric(),
  Ne = numeric(),
  num_poly_2R = numeric(),
  stringsAsFactors = FALSE
)

total_snp_2R <- 0
total_snp_3R <- 0
keep_2R <- filter_snps_by_distance(dat_2R[["POS"]][])
keep_3R <- filter_snps_by_distance(dat_3R[["POS"]][])
cat("2R SNPs before filtering:", length(dat_2R[["POS"]][]), 
    " | after filtering (≥10kb apart):", sum(keep_2R), "\n")
cat("3R SNPs before filtering:", length(dat_3R[["POS"]][]), 
    " | after filtering (≥10kb apart):", sum(keep_3R), "\n")

for (grp in group_list) {
  cat("Processing", grp, "\n")
  subset_sample_id <- metadata[metadata$group_label == grp, "sample_id"]
  matched_idx <- gt_sample_id %in% subset_sample_id
  sample_size <- sum(matched_idx)
  
  if (sample_size < 10) {
    summary_table <- rbind(summary_table, data.frame(
      group_label = grp, sample_size = sample_size,
      raw_r2 = NA, adjusted_r2 = NA, jackknife_mean = NA,
      jackknife_variance = NA, Ne = NA, num_poly_2R = NA
    ))
    next
  }
  
  gt_2R_raw <- dat_2R[['genotype']][, matched_idx, , drop = FALSE]
  gt_3R_raw <- dat_3R[['genotype']][, matched_idx, , drop = FALSE]
  pos_2R <- dat_2R[["POS"]][]
  pos_3R <- dat_3R[["POS"]][]
  
  geno_2R <- gt_2R_raw[1,,] + gt_2R_raw[2,,]
  geno_3R <- gt_3R_raw[1,,] + gt_3R_raw[2,,]
  
  non_missing_2R <- apply(geno_2R, 2, function(x) all(x >= 0))
  non_missing_3R <- apply(geno_3R, 2, function(x) all(x >= 0))
  
  maf_fun <- function(x) {
    af <- mean(x)
    p <- af / 2
    if (p > 0.5) p <- 1 - p
    return(p)
  }
  
  maf_2R <- apply(geno_2R, 2, maf_fun)
  maf_3R <- apply(geno_3R, 2, maf_fun)
  
  maf_ok_2R <- maf_2R >= 0.2
  maf_ok_3R <- maf_3R >= 0.2
  
  valid_2R <- which(non_missing_2R & maf_ok_2R)
  valid_3R <- which(non_missing_3R & maf_ok_3R)
  
  geno_2R_filt <- geno_2R[, valid_2R]
  geno_3R_filt <- geno_3R[, valid_3R]
  pos_2R_filt <- pos_2R[valid_2R]
  pos_3R_filt <- pos_3R[valid_3R]
  
  keep_2R <- filter_snps_by_distance(pos_2R_filt)
  keep_3R <- filter_snps_by_distance(pos_3R_filt)
  
  gt_2R <- geno_2R_filt[, keep_2R, drop = FALSE]
  gt_3R <- geno_3R_filt[, keep_3R, drop = FALSE]
  
  num_poly_2R <- ncol(gt_2R)
  
  if (num_poly_2R == 0 || ncol(gt_3R) == 0) {
    summary_table <- rbind(summary_table, data.frame(
      group_label = grp, sample_size = sample_size,
      raw_r2 = NA, adjusted_r2 = NA, jackknife_mean = NA,
      jackknife_variance = NA, Ne = NA, num_poly_2R = num_poly_2R
    ))
    next
  }
  
  total_snp_2R <- total_snp_2R + ncol(gt_2R)
  total_snp_3R <- total_snp_3R + ncol(gt_3R)
  
  gc()
  batch_size <- 100
  n_batches <- ceiling(ncol(gt_2R) / batch_size)
  s <- nrow(gt_2R)
  result_LD <- numeric(1e6)
  valid_idx <- 0
  
  for (batch in 1:n_batches) {
    idx_start <- (batch - 1) * batch_size + 1
    idx_end <- min(batch * batch_size, ncol(gt_2R))
    gt_2R_batch <- gt_2R[, idx_start:idx_end, drop = FALSE]
    dpA <- colSums(gt_2R_batch) / (2 * s)
    dpB <- colSums(gt_3R) / (2 * s)
    dhA <- colSums(gt_2R_batch == 2) / s
    dhB <- colSums(gt_3R == 2) / s
    delta_AB_mat <- crossprod(gt_2R_batch, gt_3R) / (2 * s)
    val_matrix <- delta_AB_mat - 2 * outer(dpA, dpB)
    val_matrix <- val_matrix * s / (s - 1)
    denom_matrix <- outer(dpA - 2 * dpA^2 + dhA, dpB - 2 * dpB^2 + dhB)
    valid_mask <- denom_matrix > 0
    r2_matrix <- matrix(NA_real_, nrow = nrow(val_matrix), ncol = ncol(val_matrix))
    r2_matrix[valid_mask] <- (val_matrix[valid_mask]^2) / denom_matrix[valid_mask]
    result_LD_batch <- as.vector(r2_matrix)
    result_LD_batch <- result_LD_batch[!is.na(result_LD_batch) & !is.nan(result_LD_batch)]
    if (length(result_LD_batch) > 0) {
      result_LD[(valid_idx + 1):(valid_idx + length(result_LD_batch))] <- result_LD_batch
      valid_idx <- valid_idx + length(result_LD_batch)
    }
    gc()
  }
  
  result_LD <- result_LD[1:valid_idx]
  if (length(result_LD) == 0) {
    summary_table <- rbind(summary_table, data.frame(
      group_label = grp, sample_size = sample_size,
      raw_r2 = NA, adjusted_r2 = NA, jackknife_mean = NA,
      jackknife_variance = NA, Ne = NA, num_poly_2R = num_poly_2R
    ))
    next
  }
  
  raw_r2 <- mean(result_LD)
  adjusted_r2 <- round(raw_r2 - (1 / s) - (3.19 / s^2), 4)
  
  discriminant <- (1 / 9) - 2.76 * adjusted_r2
  Ne <- if (discriminant >= 0 && adjusted_r2 > 0) {
    (1 / 3 + sqrt(discriminant)) / (2 * adjusted_r2)
  } else { NA }
  
  jackknife_means <- jackknife_C(rep(adjusted_r2, length(result_LD)))
  jackknife_mean <- mean(jackknife_means)
  n <- length(result_LD)
  jackknife_var <- (n - 1) / n * sum((jackknife_means - jackknife_mean)^2)
  
  summary_table <- rbind(summary_table, data.frame(
    group_label = grp,
    sample_size = sample_size,
    raw_r2 = raw_r2,
    adjusted_r2 = adjusted_r2,
    jackknife_mean = jackknife_mean,
    jackknife_variance = jackknife_var,
    Ne = Ne,
    num_poly_2R = num_poly_2R
  ))
}

# =========== TIME INFO + SEASON =============
get_year <- function(x) as.numeric(sub("_.*", "", x))
get_first_month <- function(x) {
  m <- sub(".*_", "", x)
  m <- sub("[/-].*", "", m)
  as.numeric(m)
}

summary_table <- summary_table %>%
  mutate(
    year_part = get_year(group_label),
    month_part = get_first_month(group_label),
    sort_key = year_part * 100 + month_part,
    season = case_when(
      month_part %in% c(6, 7, 8, 9) ~ "Wet",
      month_part %in% c(11, 12, 1, 2, 3, 4) ~ "Dry",
      TRUE ~ "Transition"
    )
  ) %>%
  arrange(sort_key) %>%
  select(group_label, season, sample_size, num_poly_2R,
         raw_r2, adjusted_r2, Ne, sort_key) %>%
  mutate(
    raw_r2 = round(raw_r2, 4),
    adjusted_r2 = round(adjusted_r2, 4),
    Ne = round(Ne, 1)
  )

fwrite(summary_table, output_path, na = "NA")
cat("✅ Done! Summary saved to:", output_path, "\n")

# =========== PLOTTING =============
season_colors <- c("Wet" = "blue", "Dry" = "red", "Transition" = "black")

p1 <- ggplot(summary_table %>% filter(!is.na(adjusted_r2)),
             aes(x = group_label, y = adjusted_r2, color = season)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_size), vjust = -1.2, size = 3.2, color = "black") +
  scale_color_manual(values = season_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year-Month", y = expression(Adjusted~r^2), color = "Season", title = "Adjusted r²")

p2 <- ggplot(summary_table %>% filter(!is.na(raw_r2)),
             aes(x = group_label, y = raw_r2, color = season)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_size), vjust = -1.2, size = 3.2, color = "black") +
  scale_color_manual(values = season_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year-Month", y = expression(Raw~r^2), color = "Season", title = "Raw r²")

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

y_limits <- range(c(summary_table$raw_r2, summary_table$adjusted_r2), na.rm = TRUE)
p1 <- p1 + ylim(y_limits)
p2 <- p2 + ylim(y_limits)
p3 <- (p2 + p1) + plot_annotation(tag_levels = 'a', tag_prefix = "(", tag_suffix = ") ")
print(p3)

# New Ne vs Time plot
p4 <- ggplot(summary_table %>% filter(!is.na(Ne)),
             aes(x = reorder(group_label, sort_key), y = Ne, color = season)) +
  geom_point(size = 3) +
  geom_text(aes(label = sample_size), vjust = -1.2, size = 3.2, color = "black") +
  scale_color_manual(values = season_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Year-Month", y = expression(N[e]), color = "Season",
       title = expression(paste("Effective Population Size ", N[e], " by Season")))

print(p4)