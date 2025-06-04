require(hdf5r)
# 3R DATA
dat_3R<-H5File$new('/Users/kunyaoxu/Desktop/BFgam_3R_biallelic_CDS4.hdf5', 'r')
dat_3R

# METADATA
metadata<-read.csv('/Users/kunyaoxu/Desktop/BFgam_metadata.csv', header=T, stringsAsFactors=F)
dim(metadata)
# SAMPLE SIZE PER TIME POINT TABLE
monthyear<-interaction(metadata$month, metadata$year, drop=T)
by(metadata$sample_id, monthyear, length)

# USE ONLY POPULATIONS WITH >=37 SAMPLES
# WHICH ARE 072012, 072014, 102014, 072015, 102016, 042017
# SAMPLE SIZE AND NUMBER OF GENERATIONS VECTORS
s<-c(99, 74, 94, 113, 91, 37)
gen<-c(0, 24, 27, 36, 39, 45)

# EXTRACT GENOTYPE, POS, sample_id FROM hdf5
gt_3R<-dat_3R[['genotype']]
gt_3R<-gt_3R[1,,]+gt_3R[2,,]
POS<-dat_3R[['POS']][]
gt_sample_id<-dat_3R[['sample_id']][]

# REMOVE MISSING ALLELES
missing_3R<-apply(gt_3R, 2, function(x) {sum(x<0)})
gt_3R<-gt_3R[,missing_3R==0]
POS<-POS[missing_3R==0]
dim(gt_3R)
length(POS)

# CALCULATE MAF
af_3R<-apply(gt_3R, 2, function(x) {mean(x)/2})
maf_3R<-ifelse(af_3R<0.5, af_3R, 1-af_3R)
# MAF CUTOFF AND FILTER
maf_cutoff<-0.05
gt_3R<-gt_3R[,maf_3R>=maf_cutoff]
POS<-POS[maf_3R>=maf_cutoff]
dim(gt_3R)
length(POS)

# CHOOSE ONE SNP PER ~1000BP
temp_POS<-ceiling(POS/1000)
temp_gt<-NULL
for (i in 1:max(temp_POS))
{
  temp<-as.matrix(gt_3R[, temp_POS==i], nr=nrow(gt_3R))
  if (ncol(temp)>0)
  {
    temp_gt<-cbind(temp_gt, temp[,1])
  }
}
rm(temp)
dim(temp_gt)

# FIND ALLELE COUNTS PER POPULATION?
f<-function(month, year)
{
  subset_sample_id<-metadata[metadata$month==month & metadata$year==year,]$sample_id
  temp<-temp_gt[gt_sample_id %in% subset_sample_id,]
  s<-nrow(temp)
  count<-apply(temp, 2, function(x) {sum(x)})
  print(dim(temp))
  return(cbind(count, 2*s-count))
}
count1<-f(month=7, year=2012)
count2<-f(month=7, year=2014)
count3<-f(month=10, year=2014)
count4<-f(month=7, year=2015)
count5<-f(month=10, year=2016)
count6<-f(month=4, year=2017)
# WRITE TO FILE?
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

# ORIGINAL SOURCE CODE FOR NB (HUI AND BURT 2015)
# DON'T WORRY ABOUT IT

NB.estimator <-
  function(infile, alleles, sample.interval, bound=c(50,1e7), profile.likelihood=FALSE)
  {
    # INPUT CHECKING
    if (length(sample.interval)<2)
    {stop('sample.interval must be a vector of length 2 or more. See help for details.')}
    if (any(sample.interval<0))
    {stop('sample.interval must be a non-negative vector. See help for details.')}
    if (any(sample.interval%%1!=0))
    {stop('sample.interval must contain non-negative integers. See help for details.')}
    if (any(alleles%%1!=0))
    {stop('alleles must be positive integers. See help for details.')}
    
    #####
    # INFILE TOOL
    dirmulti.infile<-function(infile, alleles, sample.interval)
    {
      f<-function(dat.line)
      {
        temp<-strsplit(dat.line, ' ')[[1]]
        temp<-as.numeric(temp[temp!=''])
        return(temp)
      }
      
      dat<-readLines(infile)
      dat<-unname(sapply(dat, f))
      dat.length<-sapply(dat, length)
      dat<-dat[dat.length!=0]
      temporal.samples<-length(sample.interval)
      
      # INFILE DIMENSION CHECKING
      if (length(dat)!=temporal.samples*length(alleles)) 
      {stop('infile error. Number of loci per generation does not match. See help for details.')}
      
      id<-rep(1:length(alleles), temporal.samples)
      
      # CHECK IF EVERY LOCI HAS THE SAME NUMBER OF ALLELES AS INPUTED
      test.statistic<-sum((sapply(dat, length)-alleles)^2)
      if (test.statistic>0)
      {stop('infile error. Number of alleles within locus does not match. See help for details.')}
      
      # OUTPUT IS A SEPARATE MATRIX FOR EACH LOCUS, REMOVE ZERO COUNTS
      output<-list()
      for (i in 1:length(alleles))
      {
        output[[i]]<-matrix(unlist(dat[id==i]), ncol=alleles[i], byrow=T)
        temp<-apply(output[[i]], 2, sum)
        output[[i]]<-output[[i]][,temp!=0]
        if (any(temp==0)) 
        {
          mess<-paste(c('In locus ', i, ', some alleles are removed because of zero count. '), collapse='')
          print(mess)
        }
      }
      return(output)
    }
    #####
    # DIRICHLET-MULTINOMIAL PMF
    ddirmulti<-function(x, alpha, log=T)
    {
      if (any(x%%1!=0))
      {stop('allele counts must be non-negative integers')}
      if (any(x<0))
      {stop('allele counts must be non-negative integers')}
      temp<-lgamma(sum(x)+1)-sum(lgamma(x+1))+lgamma(sum(alpha))-lgamma(sum(alpha)+sum(x))+sum(lgamma(alpha+x)-lgamma(alpha))
      if (log==T) 
      {return(temp)}
      else 
      {return(exp(temp))}
    }
    #####
    # PARAMETER UPDATING TOOL
    dirichlet.updating<-function(dat.list, N.dip, sample.interval=sample.interval)
    {
      dirichlet.parms<-dat.list*0
      dirichlet.parms[1,]<-1
      for (i in 2:temporal.samples)
      {
        time.diff<-sample.interval[i]-sample.interval[i-1]
        kt<-(1-1/(2*N.dip))^time.diff
        drift.parms<-kt/(1-kt)
        temp<-dirichlet.parms[i-1,]+dat.list[i-1,]
        dirichlet.parms[i,]<-temp*drift.parms/(1+sum(temp)+drift.parms)
      }
      return(dirichlet.parms)
    }
    #####
    # LIKELIHOOD FOR EACH LOCUS
    dirichlet.log.likelihood<-function(dat.list, N.dip, sample.interval)
    {
      dirichlet.parms<-dirichlet.updating(dat.list=dat.list, N.dip=N.dip, sample.interval=sample.interval)
      likelihood.value<-rep(0, nrow(dat.list))
      for (i in 2:temporal.samples) 
      {likelihood.value[i]<-ddirmulti(x=dat.list[i,], alpha=dirichlet.parms[i,])}
      return(sum(likelihood.value))
    }
    #####
    # NEED TO WRAP IT ACROSS LOCI
    lapply.wrapper<-function(N.dip, dat, z=0)
    {
      log.likelihood.overall<-lapply(dat, dirichlet.log.likelihood, N.dip=N.dip, sample.interval)
      return(sum(unlist(log.likelihood.overall))-z)
    }
    
    # RUN HERE, READ IN DATA FILE
    temporal.samples<-length(sample.interval)
    dat<-dirmulti.infile(infile=infile, alleles=alleles, sample.interval=sample.interval)
    
    # MAXIMISE THE LIKELIHOOD
    result<-optimize(lapply.wrapper, interval=bound, dat=dat, maximum=T, tol= .Machine$double.eps^0.1)
    N.point<-result$maximum
    log.likelihood<-result$objective
    
    # 95% CONFIDENCE INTERVAL
    N.lb<-min(bound)
    try(N.lb<-uniroot(lapply.wrapper, c(min(bound),N.point), dat=dat, z=log.likelihood-2, tol= .Machine$double.eps^0.15)$root, silent=T)
    N.ub<-max(bound)
    try(N.ub<-uniroot(lapply.wrapper, c(N.point,max(bound)), dat=dat, z=log.likelihood-2, tol= .Machine$double.eps^0.15)$root, silent=T)
    
    # IF YOU WANT profile.likelihood
    if (profile.likelihood==TRUE)
    {
      N.value<-seq(N.lb, N.ub, length.out=100)
      profile.CI<-cbind(N.value, sapply(N.value, lapply.wrapper, dat=dat))
      colnames(profile.CI)<-c('log.like', 'N')
      return(list('N'=N.point, 'CI'=c(N.lb, N.ub), 'log.like'=log.likelihood, 'profile.CI'=profile.CI))
    }
    
    # OUTPUT LIST
    return(list('N'=N.point, 'CI'=c(N.lb, N.ub), 'log.like'=log.likelihood))
  }

# RUN MULTIGENERATIONAL NB
result<-NB.estimator('multigenerational.txt', alleles=rep(2, ncol(temp_gt)), bound=c(1e4, 1e8), sample.interval=gen)
# POINT ESTIMATE NE AND THIS IS OUR AVERAGE NE BETWEEN 072012 TILL 042017
result$N