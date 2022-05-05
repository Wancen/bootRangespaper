load("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/cbs2e6_Lb600000.rda")
load(sprintf("/proj/milovelab/mu/nullranges/liver/bootranges_cbs2e6/cbs2e6_Lb%s.rda",format(blockLength,scientific = F)))
load("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/liver_peakshg38.rda")
load("C:/Users/wancen/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/liver_gwashg38.rda")

library(tidyverse)
library(tidymv)
library(plyranges)
library(pbapply)
library(mgcv)

R=1000

find_fit<-function(data,pvalue){
  index <- mapply(function(x) which.min(abs(data$log10p - x)),pvalue)
  return(data$fit[index])
}
stack <- bind_ranges(as.list(gwas_prime), .id="iter") %>%
  mutate(iter=factor(iter, levels=seq_len(R)))

y <-stack %>% mutate( countOverlaps= count_overlaps(., peaks_gr, maxgap=10e3)) %>%
  mcols() %>% as.data.frame() %>%
  mutate(rateOverlaps = ifelse(countOverlaps>0,1,0)) %>%
  filter(prank<1800) %>%
  select(iter, rateOverlaps, log10p) %>%
  nest(-iter) %>%
  mutate(fit= map(data, ~gam(rateOverlaps ~ s(log10p), data = ., method="REML", family=binomial)),
         pred  = map(fit, ~predict_gam(model = ., length_out = 2000)),
         fitted = map(pred,~find_fit(data=.,pvalue = seq(5,20,1))))

fitted <- plogis(do.call(rbind, y$fitted))
mean <- colMeans(fitted)
sd <- apply(fitted, 2, sd)

## draw black line for observed data
data2 <-gwas_gr %>% mutate( countOverlaps= count_overlaps(., peaks_gr, maxgap=10e3)) %>%
  mcols() %>% as.data.frame() %>%
  mutate(rateOverlaps = ifelse(countOverlaps>0,1,0)) %>% filter(prank<1800)
fit <- mgcv::gam(rateOverlaps ~ s(log10p) , data = data2, method="REML", family=binomial)
library(tidymv)
fit_data <- predict_gam(fit,length_out = 2000)
fitted_obs <- plogis(find_fit(fit_data,pvalue = seq(5,20,1)) )
z <- abs(fitted_obs-mean)/sd
# plot(fit, pages = 1, trans = plogis,shift = coef(fit)[1],seWithMean = TRUE,ylab="Prob of Overlap",xlab="Log10(pvalue)")
# summary(fit)

library(ggsci)
x = seq(5,20,1)
curves <- lapply(seq_len(length(z)), function(i) {
  mu <- mean[i]
  sd <- sd[i]
  seq <- seq(mu+2*sd, mu-2*sd, length.out = 100)
  data.frame(
    x = -0.1 * dnorm(seq, mean = mu, sd=sd) + x[i],
    y = seq,
    grp = i
  )
})
# Combine above densities in one data.frame
curves <- do.call(rbind, curves)
curves$z <- rep(z,each=100)
p1 <- fit_data %>%
  mutate(fit_format = plogis(fit),lwr=plogis(fit-1.96*se.fit),upr=plogis(fit+1.96*se.fit)) %>% filter(log10p<30) %>%
  ggplot(aes(x=log10p, y=fit_format)) +
  geom_jitter()+
  geom_line(aes(col=pal_igv("default")(2)[1]),size=2)+
  geom_ribbon(aes(ymin = lwr, ymax = upr), color = "grey70", alpha = 0.25, fill = "grey50")+
  labs(x = "-log10(pvalue)", y = "fitted overlap rate")+
  theme_classic(base_size = 16)+
  scale_color_igv()+
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))+
  geom_line(data = curves, aes(x=x, y= y, group = grp))+
  geom_vline(xintercept=x, lty=2)

ggplot()+
  geom_path(data = curves, aes(x=x, y= y, group = grp))

for (i in 1:1000) {
  p <- p + geom_line(data =y$pred[[i]], aes(x = log10p,y=plogis(fit),col=pal_igv("default")(2)[2]))
}
p+ geom_vline(xintercept = seq(from=5, to=10, by = 0.5),linetype='dashed', color='gray')
boot_meanrate<-do.call(c,res[seq(1,2000,2)])
hist(boot_meanrate)

## conditional density plot
library(ggridges)
dataPlot<- fitted %>% as.data.frame() %>% `colnames<-`(seq(5,20,1)) %>%
  pivot_longer(everything(), names_to = "log10p", values_to = "fitted_rate")
dataPlot_num <- transform(dataPlot, log10p = as.numeric(log10p))
dataPlot_num$z <- rep(z,1000)

data_obs = fit_data %>%
  mutate(fit_format = plogis(fit),lwr=plogis(fit-1.96*se.fit),upr=plogis(fit+1.96*se.fit))%>% filter(log10p<30)
p2<- data_obs %>% ggplot(aes(y=log10p, x=fit_format)) +
  geom_jitter()+
  geom_line(size=2)+
  geom_ribbon(aes(xmin = lwr, xmax = upr), color = "grey70", alpha = 0.25, fill = "grey50")+
  labs(y = "-log10(pvalue)", x = "fitted overlap rate")+
  geom_density_ridges(data = dataPlot_num, aes(x = fitted_rate, y = log10p,group = log10p,fill=z,
                                           point_color = z),rel_min_height = 0.01,
                      alpha = .6,scale=1.2, point_alpha = 0.2, jittered_points = T,
                      position = position_points_jitter(width = 0.05, height = 0),
                      point_shape = '|', point_size = 2)+
  theme_ridges(center_axis_labels = TRUE)+
  scale_fill_viridis_c(name = "z score", option = "C") +
  coord_flip()
library(patchwork)
p2+p3
