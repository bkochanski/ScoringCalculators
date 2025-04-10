source('Gini_combined_master.R')

corrs = 0.25
defaultrate = .1
g<-gini_combine_calculator(.4, .3, corrs, defaultrate)
g

corr_opt<-(g["a_opt"] * g["rho1"] + g["rho2"]) / sqrt(g["a_opt"] ^ 2 + 2 * g["a_opt"] * corr +  1)
gini_from_r(corr_opt, defaultrate)

gini_from_aop <- function(a, corr=0.25){gini_from_r((a * g["rho1"] + g["rho2"]) / sqrt(a ^ 2 + 2 * a * corr +  1), defaultrate)}
#gini_from_aop(g["a_opt"])
gini_from_aop <- Vectorize(gini_from_aop)

#a_opt_range<-g["a_opt"]+(-5:5)/10
a2w <- function(a){a/(1+a)}
w2a <- function(w){w/(1-w)}

a2w(g["a_opt"])

ws<-1:199/200


plot(ws, gini_from_aop(w2a(ws), corr = 0.2), ylim=c(0,.5), type="l")
points(ws, gini_from_aop(w2a(ws), corr=0.35), type="l", lty=2)
points(ws, gini_from_aop(w2a(ws), corr=0.6), type="l", lty=3)

c1 <- 0.15
c2 <- 0.4
c3 <- 0.65

(gres1 = gini_combine_calculator(.4, .3, c1, defaultrate))
(gres2 = gini_combine_calculator(.4, .3, c2, defaultrate))
(gres3 = gini_combine_calculator(.4, .3, c3, defaultrate))

gresdf <- data.frame(rbind(gres1, gres2, gres3))

gdf <- data.frame(ws, 
                  corr1 = gini_from_aop(w2a(ws), corr = c1), 
                  corr2 = gini_from_aop(w2a(ws), corr = c2), 
                  corr3 = gini_from_aop(w2a(ws), corr = c3))
#View(gdf)
gresdf$y30 <- as.numeric(gdf[gdf$ws==0.3,2:4])

library(ggplot2)
library(dplyr)
library(ggrepel)

gdf %>% ggplot()+ theme_bw()+ 
  xlab("Normalized weight of scorecard 1") + 
  ylab("Gini of the combined scorecard")+ 
  scale_x_continuous(breaks=1:5/5) +
  geom_line(aes(x=ws, y=corr1)) +
  geom_line(aes(x=ws, y=corr2), lty=2) +
  geom_line(aes(x=ws, y=corr3), lty=3) +
  geom_point(aes(x=score_1_weight, y = new_gini), data=gresdf) +
  geom_text_repel(aes(x=score_1_weight, y = new_gini, label = round(new_gini,3)), data=gresdf) +
  geom_label(aes(x=0.30, y=y30, label=paste0("corr = ", corr)), data=gresdf) +
  geom_text_repel(aes(x=x, y=y, label=la), data=data.frame(x=-0.05, y=g["gini2"], la=paste0("Gini of scorecard 2 = ",g["gini2"]))) +
  geom_text_repel(aes(x=x, y=y, label=la), data=data.frame(x=1+0.1, y=g["gini1"]-0.005, la=paste0("Gini of scorecard 1 = ",g["gini1"])))
 
                   