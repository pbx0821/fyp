require(hdf5r)

### === 3R ===
dat_3R <- H5File$new('/Users/kunyaoxu/Desktop/BFgam_3R_biallelic_CDS4.hdf5', 'r')
gt_3R <- dat_3R[['genotype']]
gt_3R <- gt_3R[1,,] + gt_3R[2,,]
POS_3R <- dat_3R[['POS']][]
gt_sample_id_3R <- dat_3R[['sample_id']][]
metadata <- read.csv('/Users/kunyaoxu/Desktop/BFgam_metadata.csv', header=T, stringsAsFactors=F)

# Filter missing and MAF
missing_3R <- apply(gt_3R, 2, function(x) sum(x < 0))
gt_3R <- gt_3R[, missing_3R == 0]
POS_3R <- POS_3R[missing_3R == 0]
af_3R <- apply(gt_3R, 2, function(x) mean(x)/2)
maf_3R <- ifelse(af_3R < 0.5, af_3R, 1 - af_3R)
gt_3R <- gt_3R[, maf_3R >= 0.05]
POS_3R <- POS_3R[maf_3R >= 0.05]

# Subsample one SNP per 1000 bp
temp_POS <- ceiling(POS_3R / 1000)
temp_gt_3R <- NULL
for (i in 1:max(temp_POS)) {
  temp <- as.matrix(gt_3R[, temp_POS == i], nr = nrow(gt_3R))
  if (ncol(temp) > 0) {
    temp_gt_3R <- cbind(temp_gt_3R, temp[, 1])
  }
}

# Allele count function
f <- function(month, year, temp_gt, gt_sample_id) {
  subset_sample_id <- metadata[metadata$month == month & metadata$year == year, ]$sample_id
  temp <- temp_gt[gt_sample_id %in% subset_sample_id, ]
  s <- nrow(temp)
  count <- apply(temp, 2, function(x) sum(x))
  return(cbind(count, 2 * s - count))
}

# Sample sizes and generation times
s <- c(99, 74, 94, 113, 91, 37)
gen <- c(0, 24, 27, 36, 39, 45)

# Get counts and write to file for 3R
count_list_3R <- list(
  f(7, 2012, temp_gt_3R, gt_sample_id_3R),
  f(7, 2014, temp_gt_3R, gt_sample_id_3R),
  f(10, 2014, temp_gt_3R, gt_sample_id_3R),
  f(7, 2015, temp_gt_3R, gt_sample_id_3R),
  f(10, 2016, temp_gt_3R, gt_sample_id_3R),
  f(4, 2017, temp_gt_3R, gt_sample_id_3R)
)

write.table(count_list_3R[[1]], col.names=F, row.names=F, file='multigenerational_3R.txt')
for (i in 2:length(count_list_3R)) {
  write('', file='multigenerational_3R.txt', append=T)
  write.table(count_list_3R[[i]], col.names=F, row.names=F, file='multigenerational_3R.txt', append=T)
}

### === 2R ===
dat_2R <- H5File$new('/Users/kunyaoxu/Desktop/BFgam_2R_biallelic_CDS4.hdf5', 'r')
gt_2R <- dat_2R[['genotype']]
gt_2R <- gt_2R[1,,] + gt_2R[2,,]
POS_2R <- dat_2R[['POS']][]
gt_sample_id_2R <- dat_2R[['sample_id']][]

missing_2R <- apply(gt_2R, 2, function(x) sum(x < 0))
gt_2R <- gt_2R[, missing_2R == 0]
POS_2R <- POS_2R[missing_2R == 0]
af_2R <- apply(gt_2R, 2, function(x) mean(x)/2)
maf_2R <- ifelse(af_2R < 0.5, af_2R, 1 - af_2R)
gt_2R <- gt_2R[, maf_2R >= 0.05]
POS_2R <- POS_2R[maf_2R >= 0.05]

temp_POS <- ceiling(POS_2R / 1000)
temp_gt_2R <- NULL
for (i in 1:max(temp_POS)) {
  temp <- as.matrix(gt_2R[, temp_POS == i], nr = nrow(gt_2R))
  if (ncol(temp) > 0) {
    temp_gt_2R <- cbind(temp_gt_2R, temp[, 1])
  }
}

# Get counts and write to file for 2R
count_list_2R <- list(
  f(7, 2012, temp_gt_2R, gt_sample_id_2R),
  f(7, 2014, temp_gt_2R, gt_sample_id_2R),
  f(10, 2014, temp_gt_2R, gt_sample_id_2R),
  f(7, 2015, temp_gt_2R, gt_sample_id_2R),
  f(10, 2016, temp_gt_2R, gt_sample_id_2R),
  f(4, 2017, temp_gt_2R, gt_sample_id_2R)
)

write.table(count_list_2R[[1]], col.names=F, row.names=F, file='multigenerational_2R.txt')
for (i in 2:length(count_list_2R)) {
  write('', file='multigenerational_2R.txt', append=T)
  write.table(count_list_2R[[i]], col.names=F, row.names=F, file='multigenerational_2R.txt', append=T)
}

### === Run Ne estimation for both ===
result_3R <- NB.estimator('multigenerational_3R.txt', alleles=rep(2, ncol(temp_gt_3R)), bound=c(1e4, 1e8), sample.interval=gen)
result_2R <- NB.estimator('multigenerational_2R.txt', alleles=rep(2, ncol(temp_gt_2R)), bound=c(1e4, 1e8), sample.interval=gen)

print(result_3R$N)
print(result_2R$N)