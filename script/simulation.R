##########################################
## visualizing genes and other features ##
##########################################

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
g <- genes(txdb)

# add symbols to genes
suppressPackageStartupMessages(library(plyranges))

g <- g %>%
  mutate(symbol = mapIds(org.Hs.eg.db, gene_id,
                         "SYMBOL", "ENTREZID"))

# for visualizing, restrict to a range of chr4
chrom <- "chr4"
rng <- c(98.8e6, 99.8e6) # where we will zoom into, 1 Mb
rng_big <- c(90e6, 110e6) # where features live, 20 Mb

# filtering the genes to this range
r <- data.frame(seqnames=chrom, start=rng[1]+1, end=rng[2]) %>%
  as_granges()

# just look at the genes in this small range
g <- g %>%
  filter_by_overlaps(r) %>%
  sort() %>%
  arrange(strand)

source("boot_and_match_script.R")

suppressPackageStartupMessages(library(plotgardener))

plotSomeGenes(chrom, rng, showGuides=FALSE)

# make n features in clumps of ~lambda
seqlens <- seqlengths(g)[chrom]
gr <- makeClusterRanges(chrom, rng_big, n=300, lambda=5, seqlens)

# define some plotting parameters for plotgardener,
# e.g. a palette for feature 'score':
pal <- colorRampPalette(c("dodgerblue2", "firebrick2"))

# shared genomic location, width & height, x position, fill, etc.
p <- pgParams(
  chrom=chrom, chromstart=rng[1], chromend=rng[2],
  width=5.5, height=1, x=.25,
  fill=colorby("score", palette=pal),
  order="random", baseline=TRUE,
)

# shared parameters for text labels
textp <- pgParams(x=.1, rot=90, just="left")

# plot the original GRanges
plotRanges(gr, params=p, y=2)
plotText("original", params=textp, y=3)

# uniform shuffling
shuf <- shuffle(gr, rng_big)

# plot shuffled ranges
plotRanges(shuf, params=p, y=1)
plotText("shuffled", params=textp, y=2)

# segmented block bootstrapping
# blocks 100kb, not proportion to segment length
library(nullranges)
seg <- makeSegmentation(chrom, rng_big, seqlens)
boot <- bootRanges(gr, blockLength=1e5, R=1,
                   seg=seg, proportionLength=FALSE)

# plot bootstrapped ranges
plotRanges(boot, params=p, y=0)
plotText("boot", params=textp, y=1)

# for genome-wide analysis, consider excluding gaps, repeats, etc.
# see https://dozmorovlab.github.io/excluderanges for details

#library(AnnotationHub)
#ah <- AnnotationHub()
#query(ah, "excluderanges")

###########################
## bootstrapping example ##
###########################

# first just counts as statistic

g %>%
  mutate(n_overlaps = count_overlaps(., gr))

g %>% join_overlap_left(gr) %>%
  group_by(symbol) %>%
  summarize(n_overlaps = sum(!is.na(id))) %>%
  as_tibble()


# working with metadata

gr_res <- g %>% join_overlap_left(gr) %>%
  group_by(symbol) %>%
  summarize(ave_score = sum(score))
res_true <- sum(gr_res$ave_score, na.rm = TRUE)
# simple violin plot

library(tibble)
library(ggplot2)
# inner instead of left: leaves out no-overlap genes
g %>% join_overlap_inner(gr) %>%
  mutate(type = "original") %>%
  group_by(symbol, type) %>%
  summarize(ave_score = mean(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, ave_score)) +
  geom_violin() +
  geom_jitter()

# adding more draws from the distribution for simulated features

niter <- 50
sim_list <- replicate(niter, {
  makeClusterRanges(chrom, rng_big, n=300, lambda=5, seqlens)
})
sim_long <- bind_ranges(sim_list, .id="iter")

g %>% join_overlap_inner(sim_long) %>%
  mutate(type = "original") %>%
  group_by(symbol, iter, type) %>%
  summarize(ave_score = mean(score)) %>%
  as_tibble() %>%
  ggplot(aes(type, ave_score)) +
  geom_violin() +
  geom_jitter()

# shuffling and bootstrapping multiple times

shuf_list <- replicate(niter, shuffle(gr, rng_big))
shuf_long <- bind_ranges(shuf_list, .id="iter")

boot_long <- bootRanges(gr, blockLength=5e5, R=niter,
                   seg=seg, proportionLength=FALSE)

# bind together

lvls <- c("sim","shuffle","block_bootstrap")
all <- bind_ranges(sim=sim_long, shuffle=shuf_long,
                   block_bootstrap=boot_long, .id="type") %>%
  mutate(type = factor(type, levels=lvls))

# total overlaps peaks
res <- g %>% join_overlap_left(all) %>%
  group_by(symbol, iter, type) %>%
  summarize(n_overlaps = sum(!is.na(id))) %>%
  as_tibble()%>%
  tidyr::complete(iter, symbol,type, fill=list(n_overlaps = 0)) %>%
  group_by(iter, type) %>%
  summarize(sumOverlaps = sum(n_overlaps))

## mean overlap peaks per gene
res <- g %>% join_overlap_left(all) %>%
  group_by(symbol, iter, type) %>%
  summarize(sumOverlaps = sum(!is.na(id))) %>%
  as_tibble()%>%
  tidyr::complete(iter, symbol,type, fill=list(sumOverlaps = 0))

res_wide <- spread(res,type,sumOverlaps)
res %>%
  ggplot(aes(type, sumOverlaps)) +
  geom_violin() +
  geom_boxplot() +
  geom_jitter(width=.25, alpha=.15)

p1 <- ggplot(res)+
  geom_density(aes(sumOverlaps, fill = type),alpha=0.4)+
  # labs(x = "-log10(pvalue)")+
  scale_fill_igv()+
  scale_color_igv()+
  theme_minimal(base_size = 16)+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))

shuffle_0.95 <- quantile(res_wide$shuffle, 0.95)
blockboot_0.95 <- quantile(res_wide$block_bootstrap, 0.95)
percentile <- ecdf(res_wide$sim)
percentile(blockboot_0.95)
percentile(shuffle_0.95)

# show table of features per iteration
head(table(all$iter, all$type))

# look at ave_score
res <- g %>% join_overlap_inner(all) %>%
  group_by(symbol, iter, type) %>%
  summarize(sum_score = sum(score)) %>%
  as_tibble()

# final plot of distributions:
# multiple draws, shuffling one instance, bootstrapping one instances

res_long <- g %>% join_overlap_inner(all) %>%
  group_by(symbol, iter, type) %>%
  summarize(sumScore = sum(score)) %>%
  as_tibble()

res_long %>%
  ggplot(aes(type, sum_score)) +
  geom_violin() +
  geom_boxplot() +
  geom_jitter(width=.25, alpha=.15)

res_wide <- spread(res_long,type,sum_score)

p2 <- ggplot(res_long)+
  geom_density(aes(sumScore, fill = type),alpha=0.4)+
  # labs(x = "-log10(pvalue)")+
  scale_fill_igv()+
  scale_color_igv()+
  theme_minimal(base_size = 16)+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))

library(patchwork)
jpeg(file="plots/simulation.jpeg",width = 10, height = 4,units = "in",res=450)
p1+p2
dev.off()

