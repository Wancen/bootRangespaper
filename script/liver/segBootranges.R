cmd_args=commandArgs(TRUE)

blockLength <- as.numeric(cmd_args[1]) # blockLength
out <- cmd_args[2]

library(readxl)
library(plyranges)
library(vroom)
library(GenomeInfoDb)
library(nullranges)
load("/proj/milovelab/mu/nullranges/data/liver_gwashg38.rda")
load("/proj/milovelab/mu/nullranges/data/liver_peakshg38.rda")
library(nullrangesData)

## Segmentation on gene density
ah <- AnnotationHub()
query(ah, c("excluderanges","hg38"))
exclude <- ah[["AH95917"]]
exclude <- exclude %>% plyranges::filter(seqnames(exclude)!="chrY"&seqnames(exclude)!="chrX")
exclude2 <- exclude %>% plyranges::filter(width(exclude) >= 500)
L_s <- 2e6

suppressPackageStartupMessages(library(EnsDb.Hsapiens.v86))
edb <- EnsDb.Hsapiens.v86
filt <- AnnotationFilterList(GeneIdFilter("ENSG", "startsWith"))
g <- genes(edb, filter = filt)

library(GenomeInfoDb)
g <- keepStandardChromosomes(g, pruning.mode = "coarse")
seqlevels(g, pruning.mode="coarse") <- setdiff(seqlevels(g), c("MT","X","Y"))
# normally we would assign a new style, but for recent host issues
## seqlevelsStyle(g) <- "UCSC"
seqlevels(g) <- paste0("chr", seqlevels(g))
genome(g) <- "hg38"

g <- sortSeqlevels(g)
g <- sort(g)

seg <- segmentDensity(g, n = 2, L_s = L_s,
                       exclude = exclude2, type = "cbs")

# bootstrap on GWAS
# blockLength <- 0.7^seq(13,22,0.6)*1e8
R=1000
gwas_prime = bootRanges(gwas_gr, seg, blockLength = blockLength, R = R, exclude = exclude, withinChrom=FALSE)
save(gwas_prime, file=out)

# prop = bootRanges(gwas_gr, seg, blockLength, 1, exclude, proportionLength = TRUE)
# no_prop = bootRanges(gwas_gr, seg, blockLength, 1, exclude, proportionLength = FALSE)

