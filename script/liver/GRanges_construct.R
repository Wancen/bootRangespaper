library(readxl)
library(plyranges)
library(dplyr)
library(vroom)
library(GenomeInfoDb)
gwas <- vroom("/proj/milovelab/mu/nullranges/data/total_cholesterol-2021-09-22.csv")
# read consensus peaks
peaks<-vroom("/proj/milovelab/mu/nullranges/data/GSE164870_consensusNarrowPeaksWithPeakIDs-1bpOverlap-min3samples_liver-firstPass-20samples.bed")
# read each sample peaks
l <- lapply(list.files(path="/proj/milovelab/mu/nullranges/data/",pattern = ".narrowPeak.gz$", full=TRUE), function(f) {
  vroom(f, col_names = c("peak_chromosome", "peak_start", "peak_stop", "peak_name", "score", "strand", "fold enrichment", "-log10(p_value)", "-log10(q_value)",  "peak summit location relative to peak_start"))
})
lapply(1:20, function(t) quantile(l[[t]]$`fold enrichment`))

b<-import("/proj/milovelab/mu/nullranges/data/GSM5021283_469_peaks_blacklistFilter.narrowPeak.gz",format = "narrowpeak")
## concatenate all samples together
dat <- do.call(rbind, l)

### Construct peaks_gr #######################################################
## overlap merged peaks to consenses peaks to derive mean fc
library(GenomicRanges)
library(GenomeInfoDb)
si <- Seqinfo(genome="hg19")
si <- keepSeqlevels(si, value=paste0("chr",1:22))
peaks_gr <- makeGRangesFromDataFrame(peaks,seqinfo = si, keep.extra.columns = T,seqnames.field = "#chr")
table(seqnames(peaks_gr))
dat_gr <- makeGRangesFromDataFrame(dat,seqinfo = si, keep.extra.columns = T,seqnames.field = "peak_chromosome",
                                   start.field = "peak_start",end.field = "peak_stop")
dat_gr$lfc <- log2(dat_gr$`fold enrichment`)
olap <- peaks_gr %>% join_overlap_left(dat_gr)
olap$peak_ID <-factor(olap$peak_ID, levels = paste0("peak",1:length(peaks_gr)))
peaks_fc <- olap %>%
  group_by(peak_ID) %>%
  summarize(
    lfc=mean(`lfc`),
    `-log10qvalue`=mean(`-log10(q_value)`))
peaks_gr$lfc <- peaks_fc$lfc
# peaks_gr$lfc=log2(peaks_gr$fc)
peaks_gr$log10qvalue <- peaks_fc$`X.log10qvalue`

## liftover peaks from hg19 to hg38
library(rtracklayer)
ch = import.chain("/proj/milovelab/mu/nullranges/GWAS/hg18ToHg38.over.chain")
peaks_gr = liftOver(peaks_gr, ch) %>% unlist()
genome(peaks_gr) = "hg38"
peaks_gr <- keepSeqlevels(peaks_gr, value=paste0("chr",1:22),pruning.mode = "coarse")

### Construct gwas_gr #######################################################
## extract gwas seqnames and location
library(stringr)
idx<-which(gwas$Location=="Mapping not available")
gwas$Location[idx] <- gwas$`Variant and risk allele`[idx]  ## only mapping "^rs" allele. When risk allele name is "chr" start, it's not

temp<-str_split_fixed(gwas$Location, ":", 2) %>% `colnames<-`(c("seqnames","start")) %>% as.data.frame()
idx2<-setdiff(seq_len(nrow(temp)),idx)
temp$seqnames[idx2] <- paste0("chr", temp$seqnames[idx2])
temp$start <- str_split(temp$start, "-", simplify = T)[,1]
temp$start <- as.numeric(temp$start)
temp$end <- temp$start+1
gwas_ranges <-cbind(gwas,temp)
gwas_ranges <-gwas_ranges[!is.na(gwas_ranges$start)&gwas_ranges$seqnames%in%paste0("chr",1:22),]
gwas_ranges <- gwas_ranges[!duplicated(gwas_ranges$`Variant and risk allele`),]


si <- Seqinfo(genome="hg38")
si <- keepSeqlevels(si, value=paste0("chr",1:22))
gwas_gr <- makeGRangesFromDataFrame(gwas_ranges,seqinfo = si, keep.extra.columns = T)
gwas_gr <- sortSeqlevels(gwas_gr)
gwas_gr<- sort(gwas_gr)

save(peaks_gr, file="/proj/milovelab/mu/nullranges/data/liver_peakshg38.rda")
save(gwas_gr, file="/proj/milovelab/mu/nullranges/data/liver_gwashg38.rda")
