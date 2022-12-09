library(nullranges)
library(patchwork)
suppressPackageStartupMessages(library(plyranges))
library(AnnotationHub)
library(readr)

library(macrophage)
data(gse)
# model assys
dds_macrophage <- DESeqDataSet(gse, design = ~line + condition)
keep_macrophage <- rowSums(counts(dds_macrophage) >= 10) >= 6
table(keep_macrophage)
dds_macrophage <- dds_macrophage[keep_macrophage,]
dds_macrophage <- DESeq(dds_macrophage)

res <- results(dds_macrophage,contrast=c("condition","IFNg","naive"),
               lfcThreshold=1, alpha=0.01)
# use cauchy prior shrink
genes_macrophage <-lfcShrink(dds_macrophage,
                             coef="condition_IFNg_vs_naive", type="apeglm", format="GRanges")%>%
  names_to_column("gene_id")
#
# ## De genes
degenes_macrophage <- genes_macrophage %>% filter(padj < 0.01) %>% plyranges::select(gene_id, de_log2FC = log2FoldChange,de_padj=padj)
# degenes_macrophage
# ## define promoter and enhancer regions
degenes_macrophage2<-promoters(degenes_macrophage, 10e3, 10e3)
# library(fluentGenomics)
# peaks
# peaks <- sortSeqlevels(peaks)
# peaks <- sort(peaks)
# table(seqnames(peaks))
# da_peaks <- peaks %>%
#   filter(da_padj < 0.01)
# seqlengths(da_peaks) <- seqlengths(genes_macrophage)

## load da_peaks and de_genes directly without preprocessing
load(file = "~/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/macrophage.rda")

## Segmentation on gene density
ah <- AnnotationHub()
query(ah, c("excluderanges","hg38"))
exclude <- ah[["AH95917"]]

exclude <- exclude %>% plyranges::filter(seqnames(exclude)!="chrY")
## Bootstraps on da peaks
R=100
set.seed(5) # set seed for reproducibility
peaks_prime = bootRanges(da_peaks, blockLength = 5e5, R=R,exclude = exclude,
                         type="bootstrap", withinChrom=FALSE)

## Summary statistics - Use with plyranges ##############################
## First whether there is an enrichment
len <- peaks_prime %>% group_by(iter) %>% summarise(n=n())
boot_stat <-degenes_macrophage2 %>% join_overlap_inner(peaks_prime,maxgap=10e3) %>%
  group_by(iter) %>%
  summarize(counts = n()) %>% as.data.frame() %>%
  mutate(avg = counts/4546)

res <- t.test(boot_stat$avg, mu = mean(data2$countOverlaps))
res
# helper function
find_fit<-function(data,logFC){
  index <- mapply(function(x) which.min(abs(data$de_log2FC - x)),logFC)
  return(data$fit[index])
}

library(tidymv)
library(mgcv)
y <- degenes_macrophage2 %>%
  join_overlap_inner(peaks_prime) %>%
  group_by(gene_id,iter) %>%
  summarize(countOverlaps = plyranges::n(), de_log2FC = max(de_log2FC)) %>%
  as.data.frame() %>%
  complete(gene_id, iter, fill=list(countOverlaps = 0)) %>%
  select(iter, countOverlaps, de_log2FC) %>%
  nest(-iter) %>%
  mutate(fit= map(data, ~gam(countOverlaps ~ s(de_log2FC), data = ., method="REML", family=poisson)),
         pred  = map(fit, ~predict_gam(model = ., length_out = 2000)),
         fitted = map(pred,~find_fit(data=.,logFC = seq(-8,10,1))))

######################
fitted <- exp(do.call(rbind, y$fitted))
mean <- colMeans(fitted)
sd <- apply(fitted, 2, sd)

## draw black line for observed data
data2 <-degenes_macrophage2 %>% mutate( countOverlaps= count_overlaps(., da_peaks)) %>%
  mcols() %>% as.data.frame()
fit <- mgcv::gam(countOverlaps ~ s(de_log2FC) , data = data2, method="REML", family=poisson)
fit_data <- predict_gam(fit,length_out = 2000)
fitted_obs <- exp(find_fit(fit_data,logFC = seq(-8,10,1)) )
z <- abs(fitted_obs-mean)/sd

## conditional density plot
library(ggridges)
dataPlot<- fitted %>% as.data.frame() %>% `colnames<-`(seq(-8,10,1)) %>%
  pivot_longer(everything(), names_to = "de_logFC", values_to = "fitted_rate")
dataPlot_num <- transform(dataPlot, de_logFC = as.numeric(de_logFC))
dataPlot_num$z <- rep(z,R)

data_obs = fit_data %>%
  mutate(fit_format = exp(fit),lwr=exp(fit-1.96*se.fit),upr=exp(fit+1.96*se.fit))
p4<- data_obs %>% ggplot(aes(y=de_log2FC, x=fit_format)) +
  geom_jitter()+
  # geom_line(size=2)+
  # geom_ribbon(aes(xmin = lwr, xmax = upr), color = "grey70", alpha = 0.25, fill = "grey50")+
  labs(y = "logFC(degenes)", x = "fitted overlap count")+
  geom_density_ridges(data = dataPlot_num, aes(x = fitted_rate, y = de_logFC,group = de_logFC,fill=z,
                                               point_color = z),rel_min_height = 0.01,
                      alpha = .6,scale=1.2, point_alpha = 0.2, jittered_points = T,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 2)+
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_viridis_c(name = "z score", option = "C") +
  coord_flip()

## draw z plot
s2 <- data.frame(z=z,logFC=seq(-8,10,1)) %>%
  ggplot(aes(x=logFC,y=z))+
  geom_line()+
  labs(y = "z score", x = "logFC(degenes)")+
  geom_point()+theme_classic(base_size = 16)+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))
save(y,p4,z,s2, file = "~/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/macrophage_summary.rda")
load(file = "~/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/macrophage_summary.rda")
