library(data.table)
library(ggplot2)
library(Seurat)
library(plyranges)

## dowonload it from url("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/seurat.rds")
seurat <- readRDS("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/seurat.rds")
passIndex <- seurat@meta.data$pass_accQC==TRUE & seurat@meta.data$pass_rnaQC==TRUE
## use cells pass both assay
seurat <- seurat %>%
  .[,seurat@meta.data$pass_accQC==TRUE & seurat@meta.data$pass_rnaQC==TRUE]


table(seurat@meta.data$celltype)
table(seurat@meta.data$broad_celltype,seurat@meta.data$celltype)
## get chr postion information for each assay
load("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/feature_metadata.rds")

# feature_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/filtered_feature_bc_matrix/features.tsv.gz") %>%
#   setnames(c("ens_id","gene","view","chr","start","end"))
feature_metadata.rna <- feature_metadata[view=="Gene Expression"]
head(feature_metadata.rna,n=3)

feature_metadata.atac <- feature_metadata[view=="Peaks"] %>%
  .[,ens_id:=NULL] %>% setnames("gene","peak")
head(feature_metadata.atac,n=3)

## classify peaks into distal and promoter, select promoter
# foo <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/10x_rna_atac_vignette/atac_peak_annotation.tsv") %>%
#   .[,c("peak","peak_type")] %>%
#   .[peak_type%in%c("distal", "promoter")]
feature_metadata.atac <- feature_metadata.atac %>%
  merge(foo,by="peak",all.x=TRUE)

## construct RNA Granges with counts as CompressedNumericList
feature_metadata.rna2 <-feature_metadata.rna[feature_metadata.rna$gene %in% (seurat@assays$RNA@counts %>% rownames())]
rna_name <- feature_metadata.rna2$gene[feature_metadata.rna2$chr%in% paste0("chr",1:22)]
rna_Granges<-feature_metadata.rna2[which(chr %in% paste0("chr",1:22))] %>% .[,c("chr","start","end","gene")] %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)
# seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", assay = "RNA")
# seurat <- ScaleData(seurat, do.center = TRUE, do.scale = FALSE)


## construct promoter Granges with counts as compressedNumericList
promoter_name<-feature_metadata.atac[which(chr %in% paste0("chr",1:22) & peak_type=="promoter")]$peak
promoter_Granges <-feature_metadata.atac[which(chr %in% paste0("chr",1:22)& peak_type=="promoter")] %>% .[,c("chr","start","end","peak")] %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)

groups <- seurat@meta.data$celltype
levels(groups)<-setdiff(unique(groups),NA)
# Aggregate across cluster-sample groups
getPseudobulk <- function(mat.sparse, celltype) {
  mat.summary <- do.call(cbind, lapply(levels(celltype), function(ct) {
    cells <- which(celltype==ct)
    pseudobulk <- rowSums(mat.sparse[, cells])
    return(pseudobulk)
  }))
  colnames(mat.summary) <- levels(celltype)
  return(mat.summary)
}

## call function
rna.summary <- getPseudobulk(GetAssayData(seurat)[rna_name,], groups)
promoter.summary <- getPseudobulk(seurat@assays$ATAC[promoter_name,], groups)

rna.sd <-apply(rna.summary, 1, sd)
rna.summary <- rna.summary[-which(rna.sd==0),]
promoter.sd <-apply(promoter.summary, 1, sd)
promoter.summary <- promoter.summary[-which(promoter.sd==0),]

library(edgeR)
rna.scaled <-edgeR::cpm(rna.summary,log = T)
promoter.scaled <-edgeR::cpm(promoter.summary,log = T)

#######################
## plyranges example ##
#######################
## split sparse count matrix into NumericList
rna <- rna_Granges[-which(rna.sd==0)] %>%
  mutate(counts1 = NumericList(asplit(rna.scaled, 1)))%>% sort()
promoter <- promoter_Granges[-which(promoter.sd==0)] %>%
  mutate(counts2 = NumericList(asplit(promoter.scaled, 1))) %>% sort()

## use join_overlap_inner to separate read counts all=0 from non-overlaping
rna$counts1 <- NumericList(lapply(rna$counts1, function(z) (z-mean(z))/sd(z)))
rna$counts1 <- NumericList(lapply(rna$counts1, scale))
promoter$counts2 <- NumericList(lapply(promoter$counts2, function(z) (z-mean(z))/sd(z)))

load(file = "C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/sc_omics.rda")

cor_gr<- rna %>% join_overlap_inner(promoter, maxgap=1000) %>%
  mutate(rho = 1/13 * sum(counts1 * counts2)) %>%
  select(gene, peak, rho)
cor_df <- mcols(cor_gr)
mean(cor_gr$rho,na.rm=T)

library(nullranges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

seqlengths(promoter) <- seqlengths(genome)[1:22] # pull chrom lens from USCS
R=100
bootranges <- bootRanges(promoter,blockLength = 5e5,R=R,type = "bootstrap", withinChrom = F)

## draw mean correlation distribution plot
cor_whole<-rna %>% join_overlap_inner(bootranges, maxgap=1000) %>%
  mutate(rho = 1/13 * sum(counts1 * counts2)) %>%
  select(rho,iter) %>%
  group_by(iter) %>%
  summarise(meanCor = mean(rho)) %>%
  as.data.frame()

s2.3<-cor_whole %>% ggplot(aes(meanCor))+
  geom_histogram(bins=15)+
  theme_classic(base_size = 16)+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))
# h<-hist(x, breaks=10, col="red", xlab="Mean correlation per iteration",
#         main="Histogram with Normal Curve")
# xfit<-seq(min(x),max(x),length=40)
# yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
# yfit <- yfit*diff(h$mids[1:2])*length(x)
# lines(xfit, yfit, col="blue", lwd=2)

## Find confidence interval bound per gene
boot_cor_df <-rna %>% join_overlap_inner(bootranges, maxgap=1000) %>%
  mutate(cor = 1/13 * sum(counts1 * counts2)) %>%
  select(gene, peak, cor) %>%
  group_by(gene) %>%
  summarise(lowTile = quantile(cor,0.05),highTile = quantile(cor,0.95)) %>%
  as.data.frame()

res<- cor_df %>% as.data.frame() %>% dplyr::left_join(boot_cor_df,by="gene") %>%
  filter(rho<lowTile.5. | rho>highTile.95.)

## draw specific gene's bootstrap correlation distribution plot
cor_df<- rna[which(rna$gene=="DPM1")] %>% join_overlap_inner(bootranges, maxgap=1000) %>%
  mutate(cor=cor(counts1, counts2)) %>% filter(!is.na(cor)) %>%
  select(gene, peak, cor) %>%
  as.data.frame()
x <- cor_df$cor
h<-hist(x, breaks=10, col="red", xlab="Mean correlation per iteration",
        main="Histogram with Normal Curve")
xfit<-seq(min(x),max(x),length=40)
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x))
yfit <- yfit*diff(h$mids[1:2])*length(x)
lines(xfit, yfit, col="blue", lwd=2)

####################
## tidySE example ##
####################
# from Stefano Mangiola on #tidiness_in_bioc branch:adapt_to_rangedSE
# devtools::install_github("stemangiola/tidySummarizedExperiment@adapt_to_rangedSE")
library(tidySummarizedExperiment)
library(purrr)
# make an SE where each row is an overlap
se_rna <- SummarizedExperiment(
  assays=list(rna=rna.scaled),
  rowRanges=rna_Granges)
se_promoter <- SummarizedExperiment(
  assays=list(promoter=promoter.scaled),
  rowRanges=promoter_Granges)

# make an SE where each row is an overlap
makeOverlapSE <- function(se_rna, se_promoter) {
  idx <- rowRanges(se_rna) %>% join_overlap_inner(rowRanges(se_promoter),maxgap = 1000)
  assay_x <- assay(se_rna, "rna")[ idx$gene, ]
  assay_y <- assay(se_promoter, "promoter")[ idx$peak, ]
  # this is needed to build the SE
  rownames(assay_x) <- rownames(assay_y) <- seq_along( idx$gene )
  new_ranges <- rowRanges(se_rna)[ idx$gene ]
  names(new_ranges) <- seq_along( idx$gene )
  SummarizedExperiment(
    assays=list(x=assay_x, y=assay_y),
    rowRanges=new_ranges
  )
}
se <- makeOverlapSE(se_rna, se_promoter)
### two ways ############
## code shorter but time longer
se <- se %>%
  nest(data = -.feature) %>%
  mutate(rho = map_dbl(data, ~cor(pull(.x, x), pull(.x, y)))) %>%
  unnest(data)

## code longer but time shorter
se %>%
  as_tibble() %>%
  nest(data = -.feature) %>%
  mutate(rho = map(data,
                   function(data) data %>% summarize(rho = cor(x, y))
  )) %>%
  unnest(rho) %>%
  select(-data)
mcols(se)

## Draw specific gene observed correlation plot
df <- data.frame(rna = rna$counts1[which(rna$gene=="TET3")],peak= promoter$counts2[which(promoter$peak=="chr2:74000098-74003475")])
library(ggplot2)
library(ggsci)
s2.1 <- df %>% mutate(celltype = rownames(df)) %>%  ggplot(df,mapping=aes(y = rna.value,x = peak.value,col=celltype))+
  geom_point()+theme_classic(base_size = 16)+
  ggtitle("TET3")+
  # geom_smooth(se = FALSE, method = lm,formula = y ~ x)+
  scale_color_igv()+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))
cor(df$rna,df$peak)

df <- data.frame(rna = rna$counts1[which(rna$gene=="CD83")],peak= promoter$counts2[which(promoter$peak=="chr6:14116971-14139988")])
s2.2 <- df %>% mutate(celltype = rownames(df)) %>%  ggplot(df,mapping=aes(y = rna.value,x = peak.value,col=celltype))+
  geom_point()+theme_classic(base_size = 16)+
  ggtitle("CD83")+
  # geom_smooth(se = FALSE, method = lm,formula = y ~ x)+
  scale_color_igv()+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))
jpeg(file="C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/plots/sfig2.jpeg",width = 18, height = 6,units = "in",res=450)
s2.3 + s2.1+s2.2
dev.off()

cor(df$rna.value,df$peak.value)
