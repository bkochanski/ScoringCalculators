source('Gini_combined_master.R')

r=0:99/100
g=1:99/100
e=expand.grid(g,r)
tmp<-Vectorize(function(g, r){gini_combine_calculator(.6, g, r, .1)[[1]]})
res<-tmp(e$Var1, e$Var2)

e2<-cbind(e,res)
#write.csv2(e2_v1,'e2_v1.csv')
#write.csv2(e2,'e2_v2.csv')
#write.csv2(e2,'e2_v4.csv')

df<-e2
df<-read.csv2('e3_v2.csv')

library(ggplot2)
plot1<-ggplot(df, aes(x = Var1, y = Var2, z = res)) +
  stat_contour(aes(colour = ..level..)) + xlab('Gini 2') + ylab('Correlation') 

df_points=data.frame(x=c(.64, .41), y=c(.83, .26), z=c(0, 0))
df_points$label=paste0('(Gini=', df_points$x, ', correlaton=', df_points$y, ')')
plot2<-plot1+geom_point(data=df_points, aes(x=x, y=y, z=z))   + 
  geom_text(data=df_points, aes( x=x, y=y, z=z, label=label), nudge_y = -0.02)
library(directlabels)
direct.label(plot2, "get.means")



r=0:9/10
g=1:9/10
f=expand.grid(g,r)
tmp2<-Vectorize(function(g, r){gini_combine_calculator(.4, g, r, .1)[[1]]})
res<-tmp2(f$Var1, f$Var2)

f2<-cbind(f,res)
#write.csv2(f2,'f2_v1.csv')
f2<-read.csv2('f2_v1.csv')

library(ggplot2)
plot2<-ggplot(f2, aes(x = Var1, y = Var2, z = res)) +
  stat_contour(aes(colour = ..level..)) + xlab('Gini 2') + ylab('Correlation') 

direct.label(plot2, 'get.means')
