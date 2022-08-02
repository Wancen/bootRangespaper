### Comparsion with gsc on small data
library(vroom)
library(GenomicRanges)
library(plyranges)
library(dplyr)
library(nullranges)
y <- vroom("/proj/milovelab/mu/nullranges/encodegsc/test_data/conserved_sequences.bed",col_names = c("chrom","Start","End"))
x <- vroom("/proj/milovelab/mu/nullranges/encodegsc/test_data/ENCODE_annotations.bed",col_names = c("chrom","Start","End"))

y_gr <- makeGRangesFromDataFrame(y, keep.extra.columns=TRUE,
                                 start.field="Start",end.field="End")
x_gr <- makeGRangesFromDataFrame(x, keep.extra.columns=TRUE,
                                 start.field="Start",end.field="End")
seqlengths(x_gr) <- 1877426+1
seqlengths(y_gr) <- 1877426+1
blockLength <- 0.041*seqlengths(y_gr) %>% ceiling()

basepair_overlap <- function(y_gr,x_gr,R) { 
  boots <- bootRanges(y_gr, blockLength=blockLength,R=R,withinChrom=F)
  combined <- y_gr %>% mutate(block=0,iter=0) %>%
    bind_ranges(boots) %>% select(iter,block) 
  
  null <- combined %>% join_overlap_intersect(x_gr) %>%
    group_by(iter,block) %>%
    summarize(overlaps=sum(width))
  
  y2 <- combined %>% 
    group_by(iter,block) %>%
    summarize(widthdenom=sum(width))
  ## Generate basepair overlap fraction
  dat_merge<-merge(x=null,y=y2,by=c("iter","block"),all.x=TRUE)
  frac <- dat_merge$overlaps/dat_merge$widthdenom
  obs <-frac[1] # notate obs data
  frac <- frac[-1] 
  z = (mean(frac)-obs)/(sd(frac)*sqrt(0.0411))
  p = pnorm(z) 
} 

system.time({basepair_overlap(y_gr,x_gr,R=100)})

microbenchmark(
  basepair_overlap(y_gr,x_gr,R=100), times=10)


## expectation under the null follow gsc formula Jn = sum(feature covering length)/n
# Nmean <- x_gr %>% width() %>% sum()/(1877426+1)
# ## mean of bootstrap
# obsmean <- mean(frac)
# dist <- frac-Nmean
# dist_sd = sqrt(0.041)*sd(dist)
# -0.23933637022 +Nmean
# var(frac)
# 
# RegionLength = (1877426+1)
# 
# I = y_gr %>% width() %>% sum()
# J = x_gr %>% width() %>% sum()
# I_num = 0.041*I*J/(RegionLength^2) 
# Nmean <- I_num/(I*0.041/RegionLength) # expected
# real_stat_expectation = (IJ*0.041/RegionLength)/(I*0.041/RegionLength) # observed
# IJ = GenomicRanges::intersect(y_gr,x_gr) %>% width() %>% sum()
# obs_mean = real_stat_expectation - Nmean
# 
# obs <- GenomicRanges::intersect(y_gr,x_gr) %>% width() %>% sum()/ y_gr %>% width() %>% sum()
# z <- obs_mean/dist_sd

## Comparison with bootRanges on genome-wise analysis
# ywidth2 <- yprime %>% join_overlap_intersect(x_gr) %>%
#   group_by(iter,block) %>%
#   summarize(overlaps=sum(width))
# yprimew <- yprime %>% 
#   group_by(iter,block) %>%
#   summarize(widthdenom=sum(width))
# dat_merge<-merge(x=ywidth2,y=yprimew,by=c("iter"),all.x=TRUE)
# frac2 <- dat_merge$overlaps/dat_merge$widthdenom
# hist(frac2)
# obsmean2 <- mean(frac2)
# sd <- sd(frac2)
# z2 <-(obsmean2-obs)/sd

## read in GSC output file
gsc <- read_table("/proj/milovelab/mu/nullranges/encodegsc/statistics.txt",col_names = "stat",skip=1)
hist(gsc$stat)
mean(gsc$stat)
sd(gsc$stat)*0.041*2

## compute region overlap according to blocks
len <- yprime %>% group_by(iter,block) %>% summarise(n=dplyr::n())
dat <- yprime %>% mutate(count = countOverlaps(.,x_gr), rate=ifelse(count>0,1,0)) %>% 
  mcols() %>% as.data.frame() %>%
  group_by(iter,block) %>% 
  summarise(region = sum(rate),total = n()) %>% 
  mutate(region_overlap = region/total)
hist(dat$region_overlap)
mean(dat$region_overlap)
dist_sd<-sqrt(0.041)*sd(dat$region_overlap)

datobs <- y_gr %>% mutate(count = countOverlaps(.,x_gr), rate=ifelse(count>0,1,0)) %>% 
  mcols() %>% as.data.frame() %>%
  summarise(region = sum(rate),total = n()) %>% 
  mutate(region_overlap = region/total)
z = (mean(dat$region_overlap)-datobs$region_overlap)/dist_sd

## compute region overlap statistics according to genome-wide data
dat2 <- yprime %>% mutate(count = countOverlaps(.,x_gr), rate=ifelse(count>0,1,0)) %>% 
  mcols() %>% as.data.frame() %>%
  group_by(iter,block) %>% 
  summarise(region = sum(rate),total = n()) %>% 
  mutate(region_overlap = region/total)
hist(dat2$region_overlap)
mean(dat2$region_overlap)
dist_sd2<-sd(dat2$region_overlap)
z2 = (mean(dat2$region_overlap)-datobs$region_overlap)/dist_sd2


## Compare gsc on genome-wise
library(nullrangesData)
dhs <- DHSA549Hg38()
dhs
library(rtracklayer)
export.bed(dhs,con='/proj/milovelab/mu/nullranges/encodegsc/dhs.bed')
set.seed(5) # reproducibility
library(microbenchmark)
blockLength <- 5e5
microbenchmark(
    b_across=bootRanges(dhs, blockLength=blockLength,
                        type="bootstrap", withinChrom=FALSE), times=10)
dhs <- dhs %>% plyranges::filter(signalValue > 100) %>%
  mutate(id = seq_along(.)) %>%
  plyranges::select(id)
length(dhs)
load('/proj/milovelab/mu/nullranges/encodegsc/hg38data.rda')
x <- GRanges("chr2", IRanges(1 + 50:99 * 1e6, width=1e6), x_id=1:50)
export.bed(x,con='/proj/milovelab/mu/nullranges/encodegsc/x.bed')
length<- data.frame(chr = names(seqlengths(dhs)), start=rep(1,length(seqlengths(dhs))),end=seqlengths(dhs))
length <- makeGRangesFromDataFrame(length)
export.bed(length,con='/proj/milovelab/mu/nullranges/encodegsc/length.bed')

library(AnnotationHub)
ah <- AnnotationHub()
kidney_pks <- ah[["AH43443"]]
bladder_pks <- ah[["AH44180"]]
library(GenomeInfoDb)
x <- kidney_pks
y <- bladder_pks
x <- keepStandardChromosomes(x)
seqlevels(x, pruning.mode="coarse") <- setdiff(seqlevels(x), "chrM")
seqlevels(y, pruning.mode="coarse") <- seqlevels(x)
seqlevels(exclude, pruning.mode="coarse") <- seqlevels(x)
q_thr <- 3
s_thr <- 9
x <- x %>% filter(qValue > q_thr & signalValue > s_thr) %>%
  sort()
y <- y %>% filter(qValue > q_thr & signalValue > s_thr) %>%
  sort()
obs <- x %>% mutate(n_overlaps = count_overlaps(., y))
obs %>% summarize(total = sum(n_overlaps))
pks_to_boot <- y %>%
  mutate(id = seq_along(y)) %>%
  select(id, signal = signalValue)
R=1
exclude <- ah[["AH95912"]]
exclude <- exclude %>%
  filter(width(exclude) >= 500)

basepair_overlap <- function(pks_to_boot,x,R) { 
  boots <- bootRanges(pks_to_boot, blockLength=5e5, R=R)
  combined <- pks_to_boot %>% mutate(block=0) %>%
    bind_ranges(boots) %>% select(block) %>% 
    mutate(block = Rle(factor(block)))
  null <- combined %>% join_overlap_intersect(x) %>%
    group_by(block) %>%
    summarize(overlaps=sum(width)) %>% 
    as.data.frame() %>% 
    complete(block, fill = list(overlaps=0))

  y2 <- combined %>% 
    group_by(block) %>%
    summarize(widthdenom=sum(width))
  ## Generate basepair overlap fraction
  dat_merge<-merge(x=null,y=y2,by=c("block"),all.x=TRUE)
  frac <- dat_merge$overlaps/dat_merge$widthdenom
  obs <-frac[1] # notate obs data
  frac <- frac[-1] 
  z = (mean(frac)-obs)/(sd(frac)*sqrt(0.0021))
  p = pnorm(z) 
} 

system.time({basepair_overlap(pks_to_boot,x,R)})

microbenchmark(
  basepair_overlap(pks_to_boot,x,R), times=10)

export.bed(pks_to_boot,con='/proj/milovelab/mu/nullranges/encodegsc/bladder.bed')
export.bed(x,con='/proj/milovelab/mu/nullranges/encodegsc/kidney.bed')
## report length
length<- data.frame(chr = names(seqlengths(x)), start=rep(1,length(seqlengths(x))),end=seqlengths(x))
length <- makeGRangesFromDataFrame(length)
export.bed(length,con='/proj/milovelab/mu/nullranges/encodegsc/length.bed')


