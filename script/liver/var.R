cmd_args=commandArgs(TRUE)

# blockLength <- as.numeric(cmd_args[1]) # blockLength
seg <- cmd_args[1]
Lb <- cmd_args[2]

library(GenomeInfoDb)
library(matrixStats)
library(plyranges)
library(nullranges)
library(dplyr)
library(tidyr)

load("/proj/milovelab/mu/nullranges/data/liver_peakshg38.rda")
load("/proj/milovelab/mu/nullranges/data/liver_gwashg38.rda")

obs <-gwas_gr %>% join_overlap_left(peaks_gr,maxgap=10e3) %>%
  summarize(countOverlaps=sum(!is.na(`peak_ID`))/length(gwas_gr), rateOverlaps = 1-sum(is.na(`peak_ID`))/length(gwas_gr)) %>% 
  as.data.frame() %>% as.numeric()

R=1000

  load(sprintf("/proj/milovelab/mu/nullranges/liver/bootranges_%s/%s_Lb%s.rda",seg,seg,format(Lb,scientific = F)))
  
  stack <- bind_ranges(as.list(gwas_prime), .id="iter") %>%
    mutate(iter=factor(iter, levels=seq_len(R))) 
  len <- do.call(c,lapply(gwas_prime, length))
  
  Overlaps <-stack %>% join_overlap_left(peaks_gr,maxgap=10e3) %>%
    group_by(iter) %>%
    summarize(countOverlaps=sum(!is.na(`peak_ID`))/len, rateOverlaps = 1-sum(is.na(`peak_ID`))/len) %>%
    as.data.frame() %>% select(countOverlaps,rateOverlaps) %>% 
    summarise_all(funs(mean,var,IQR)) %>% 
    mutate(countOverlaps_Z = abs(obs[1]-countOverlaps_mean)/sqrt(countOverlaps_var),
           rateOverlaps_z = abs(obs[2]-rateOverlaps_mean)/sqrt(rateOverlaps_var))
Overlaps$method = seg
Overlaps$Lb = Lb
# 
# write out as a table
write.table(Overlaps[1,3:10], file=sprintf("/proj/milovelab/mu/nullranges/liver/csv/%s_Lb%s_var.csv",seg,format(Lb,scientific = F)), row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
