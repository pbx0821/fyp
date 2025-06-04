require(hdf5r)
# 2R DATA
dat_2R <- H5File$new('/Users/kunyaoxu/Desktop/BFgam_2R_biallelic_CDS4.hdf5', 'r')
dat_2R

# METADATA
metadata <- read.csv('/Users/kunyaoxu/Desktop/BFgam_metadata.csv', header=T, stringsAsFactors=F)
dim(metadata)
# SAMPLE SIZE PER TIME POINT TABLE
monthyear <- interaction(metadata$month, metadata$year, drop=T)
by(metadata$sample_id, monthyear, length)

# USE ONLY POPULATIONS WITH >=37 SAMPLES
# WHICH ARE 072012, 072014, 102014, 072015, 102016, 042017
# SAMPLE SIZE AND NUMBER OF GENERATIONS VECTORS
s <- c(99, 74, 94, 113, 91, 37)
gen <- c(0, 24, 27, 36, 39, 45)

# EXTRACT GENOTYPE, POS, sample_id FROM hdf5
gt_2R <- dat_2R[['genotype']]
gt_2R <- gt_2R[1,,] + gt_2R[2,,]
POS <- dat_2R[['POS']][] 
gt_sample_id <- dat_2R[['sample_id']][]

# REMOVE MISSING ALLELES
missing_2R <- apply(gt_2R, 2, function(x) { sum(x < 0) })
gt_2R <- gt_2R[, missing_2R == 0]
POS <- POS[missing_2R == 0]
dim(gt_2R)
length(POS)

# CALCULATE MAF
af_2R <- apply(gt_2R, 2, function(x) { mean(x)/2 })
maf_2R <- ifelse(af_2R < 0.5, af_2R, 1 - af_2R)

# MAF CUTOFF AND FILTER
maf_cutoff <- 0.05
gt_2R <- gt_2R[, maf_2R >= maf_cutoff]
POS <- POS[maf_2R >= maf_cutoff]
dim(gt_2R)
length(POS)

# CHOOSE ONE SNP PER ~1000BP
temp_POS <- ceiling(POS / 1000)
temp_gt <- NULL
for (i in 1:max(temp_POS)) {
  temp <- as.matrix(gt_2R[, temp_POS == i], nr = nrow(gt_2R))
  if (ncol(temp) > 0) {
    temp_gt <- cbind(temp_gt, temp[,1])
  }
}
rm(temp)
dim(temp_gt)

# FIND ALLELE COUNTS PER POPULATION
f <- function(month, year) {
  subset_sample_id <- metadata[metadata$month == month & metadata$year == year, ]$sample_id
  temp <- temp_gt[gt_sample_id %in% subset_sample_id, ]
  s <- nrow(temp)
  count <- apply(temp, 2, function(x) { sum(x) })
  print(dim(temp))
  return(cbind(count, 2 * s - count))
}
count1 <- f(month=7, year=2012)
count2 <- f(month=7, year=2014)
count3 <- f(month=10, year=2014)
count4 <- f(month=7, year=2015)
count5 <- f(month=10, year=2016)
count6 <- f(month=4, year=2017)

# WRITE TO FILE
write.table(count1, col.name=F, row.name=F, file='multigenerational.txt')
write('', file='multigenerational.txt', append=T)
write.table(count2, col.name=F, row.name=F, file='multigenerational.txt', append=T)
write('', file='multigenerational.txt', append=T)
write.table(count3, col.name=F, row.name=F, file='multigenerational.txt', append=T)
write('', file='multigenerational.txt', append=T)
write.table(count4, col.name=F, row.name=F, file='multigenerational.txt', append=T)
write('', file='multigenerational.txt', append=T)
write.table(count5, col.name=F, row.name=F, file='multigenerational.txt', append=T)
write('', file='multigenerational.txt', append=T)
write.table(count6, col.name=F, row.name=F, file='multigenerational.txt', append=T)

# RUN MULTIGENERATIONAL NB
result <- NB.estimator('multigenerational.txt', alleles=rep(2, ncol(temp_gt)), bound=c(1e4, 1e8), sample.interval=gen)
# POINT ESTIMATE NE AND THIS IS OUR AVERAGE NE BETWEEN 072012 TILL 042017
result$N
