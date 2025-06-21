#source('Gini_combined_master.R')
#r=0:99/100
#g=1:99/100
#e=expand.grid(g,r)
#tmp<-Vectorize(function(g, r){gini_combine_calculator(.6, g, r, .1)[[1]]})
#res<-tmp(e$Var1, e$Var2)
#e2<-cbind(e,res)
#write.csv2(e2,'e2_v2.csv')
#df<-e2
df<-read.csv2('e3_v2.csv')

library(ggplot2)
plot1<-ggplot(df, aes(x = Var1, y = Var2, z = res)) +
  stat_contour(aes(colour = ..level..)) + xlab('Gini2 (Gini of the additional scorecard)') + ylab('Correlation')+ scale_colour_continuous(low = "gray40", high = "black")

df_points=data.frame(x=c(.64, .41), y=c(.83, .26), z=c(0, 0))
df_points$label=paste0('(Gini1=0.60, Gini2=', df_points$x, ', correlaton=', df_points$y, ',\ncombined Gini = 0.65)')
plot2<-plot1+geom_point(data=df_points, aes(x=x, y=y, z=z))   + 
  geom_text(data=df_points, aes( x=x, y=y, z=z, label=label), nudge_y = -0.05) + theme_bw()
library(directlabels)
p4 <- direct.label(plot2, "get.means")
print(p4)

# Save as EPS
ggsave("fig04.eps", plot = p4, device = "eps", width = 7, height = 5)
