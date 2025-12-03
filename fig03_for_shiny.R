source('Gini_combined_master.R')

g1 <- 0.4
g2 <- 0.3
corr <- 0.15
defrate <- 0.1

w2a <- function(w){w/(1-w)}
ws<-1:199/200
(gres1 = gini_combine_calculator(g1, g2, corr, defrate))
gresdf <- data.frame(rbind(gres1))

gini_from_aop <- function(a, corr=0.25){gini_from_r((a * g["rho1"] + g["rho2"]) / sqrt(a ^ 2 + 2 * a * corr +  1), defrate)}
gini_from_aop <- Vectorize(gini_from_aop)

gdf <- data.frame(ws, 
                  corr1 = gini_from_aop(w2a(ws), corr = corr))
#View(gdf)
gresdf$y30 <- as.numeric(gdf[gdf$ws==g2,2])

library(ggplot2)
library(dplyr)
library(ggrepel)

p3 <- gdf %>% ggplot()+ theme_bw()+ 
  xlab("Normalized weight of scorecard 1") + 
  ylab("Gini of the combined scorecard")+ 
  scale_x_continuous(breaks=0:5/5) +
  geom_line(aes(x=ws, y=corr1)) +
  geom_point(aes(x=score_1_weight, y = new_gini), data=gresdf) +
  geom_text_repel(aes(x=score_1_weight, y = new_gini, label = round(new_gini,3)), data=gresdf) +
  geom_label(aes(x=g2, y=y30, label=paste0("corr = ", corr)), data=gresdf) +
  geom_text_repel(aes(x=x, y=y, label=la), data=data.frame(x=-0.05, y=g["gini2"], la=paste0("Gini of scorecard 2 = ",g["gini2"]))) +
  geom_text_repel(aes(x=x, y=y, label=la), data=data.frame(x=1+0.1, y=g["gini1"]-0.005, la=paste0("Gini of scorecard 1 = ",g["gini1"])))

print(p3)


 
                   