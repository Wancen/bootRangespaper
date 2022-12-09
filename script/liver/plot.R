library(tidyverse)
l<- lapply(list.files(path="~/OneDrive - University of North Carolina at Chapel Hill/Lab/project2/data/csv/", full=TRUE), function(f) {
  read.csv(f, header=FALSE)
})
dat <- do.call(rbind, l[c(1:10,12:81)])

# blockvector <- c(1000000,800000,600000,500000,400000,300000,200000,100000,7000,4000)
colnames(dat) <- c("countOverlaps_var","rateOverlaps_var", "countOverlaps_z","rateOverlaps_z","method","Lb")

library(ggplot2)
library(ggsci)
### rateOverlap ######################
p1 <- ggplot(dat,aes(x=Lb,y=rateOverlaps_var,col=method))+
  geom_line()+
  labs(x = "Block Length", y = "Var(overlap rate)")+
  geom_point()+theme_classic(base_size = 16)+
  scale_color_igv()+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))

p2<-ggplot(dat,aes(x=Lb,y=rateOverlaps_z,col=method))+
  geom_line()+
  labs(x = "Block Length", y = "z(overlap rate)")+
  geom_point()+theme_classic(base_size = 16)+
  ## horizontal line represent convential shuffling derived z score
  geom_hline(yintercept = l[[11]][1,4])+
  scale_color_igv()+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))

p2 <- ggplot(dat %>% filter(Lb>250000),aes(x=Lb,y=countOverlaps_z,col=method))+
  geom_line()+
  labs(x = "Block Length", y = "z(overlap count)")+
  geom_point()+theme_grey(base_size = 16)+
  geom_hline(yintercept = l[[11]][1,3])+
  scale_color_igv()+
  theme(legend.position="right",
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = 'transparent'))

library(patchwork)
p1+p2
