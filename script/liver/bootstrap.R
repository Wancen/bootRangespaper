library(readxl)
library(dplyr)
library(vroom)
library(GenomeInfoDb)
library(nullranges)
load("/proj/milovelab/mu/nullranges/data/liver_gwashg38.rda")
# load("/proj/milovelab/mu/nullranges/data/liver_peakshg38.rda")
library(nullrangesData)
### Based on peaks density
data("deny")
deny <- deny %>% plyranges::filter(seqnames(deny)!="chrY"&seqnames(deny)!="chrX")
seqlevels(deny) <- paste0("chr",1:22)
L_chr = seqlengths(gwas_gr)
seqlengths(deny)<-L_chr
accept <- GenomicRanges::gaps(deny) %>% plyranges::filter(strand=="*")%>% as.data.frame()
accept$cum <- cumsum(accept$width %>% as.numeric())
sequence <- do.call(c,sapply(1:nrow(accept),function(i){
  seq(accept[i,2],accept[i,3])
}))
# chr_name= do.call(c,sapply(1:nrow(accept),function(i){
#   rep(accept[i,1],accept[i,4])
# }))
library(plyranges)
library(pbapply)
R=1000
N= length(gwas_gr)
total =sum(L_chr)
n=sample(1:22,N,replace = T,prob = L_chr) %>% table

# deny_split<-split(deny,seqnames(deny))
gwas_prime <- replicate(R,{
  # generate random postion with prob according to observsed gwas count
  random <- unlist(lapply(seq_len(length(L_chr)), function(s){
    range<-which(accept$seqnames==paste0("chr",s))
    sample(sequence[accept$cum[range[1]]:accept$cum[range[length(range)]]],n[s],replace = T)
  }))
  
  # disrupt order
  index <- sample(seq_len(N),N)
  x_prime <- gwas_gr
  seqnames(x_prime) <- seqnames(gwas_gr)[index]
  ranges(x_prime) <- IRanges(start=random[index],width = 1)
  x_prime <- sort(x_prime)
})

gwas_prime <- GRangesList(gwas_prime)
mcols(gwas_prime)$iter <- seq_len(R)

save(gwas_prime, file="/proj/milovelab/mu/nullranges/liver/bootranges/bootranges.rda")
