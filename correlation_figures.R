
gini_combine_calculator<-function(g1, g2, corr, defaultrate){
  ginic<-function(bc, gc){
    #function for gini when we have cumulative goods and bads vectors
    sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*
          (bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1
  } 
  
  gini_from_r<-function(rho=.5, defrate=.1){
    F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
    F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
    F1_<-Vectorize(function(x){F1(x, defrate, rho)})
    F2_<-Vectorize(function(x){F2(x, defrate, rho)})
    2*integrate(function(x){F1_(x)*F2_(x)}, 
                lower=-Inf, upper=Inf)$value/defrate/(1-defrate)-1
  }
  
  #rho_s1
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g1}
  rho_s1<-uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
  #rho_s2
  phi_s2<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g2}
  rho_s2<-uniroot(phi_s2,lower=0,upper=1,tol = .Machine$double.eps)$root
  
  (a_opt<-(corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2))
  corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
  
  g_result0<-if(abs(corr_opt)>1) {NaN} else {gini_from_r(corr_opt, defaultrate)}
  g_result1<-if(a_opt<0 | a_opt>1000) {NaN} else {g_result0}
  return(c(new_gini=g_result1, 
           #new_gini_2=g_result0, 
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt))
}

gv<-Vectorize(function(r){gini_combine_calculator(.6, .6, r, .1)[[1]]})
corrs<-0:9/10
ginis<-gv(corrs)
plot(corrs, ginis, xlab="correlation", ylab="new scorecard gini")

library(ggplot2)
ggplot(data.frame(corrs, ginis), aes(x=corrs, y=ginis)) + 
  geom_point() +
  geom_text(label=paste0('r=', round(corrs,2), ' Gini=', round(ginis,3)), hjust=0.4, vjust=1) +
  xlab("Correlation") + ylab("Gini") + ggtitle("Combined Gini of two scorecards with Gini=0.6 depending on the correlation between them")

(target1<-unname(gini_combine_calculator(.6, .6, .7, .1)[1]))
ginis1<-c(.6)
corrs1<-c(.7)

i=6

  corrs1[length(corrs1)+1]<-i
  search1<-function(x){gini_combine_calculator(.6, x, corrs1[length(corrs1)], .1)[1]-target1}
  ginis1[length(ginis1)+1]<-uniroot(search1,lower=ginis1[length(ginis1)]-.1,upper=ginis1[length(ginis1)],tol = .Machine$double.eps)$root

df<-data.frame(corrs1, ginis1)
library(ggplot2)
ggplot(df, aes(x=corrs1, y=ginis1)) + geom_point() +
  geom_text(label=paste0('r=', round(corrs1,2), ' Gini2=', round(ginis1,3)), hjust=0.4, vjust=1)+
  xlab("Correlation") + ylab("Gini") + ggtitle("Gini/correlation trade-off: all these combinations have the same effect (Combined Gini=.648) when combining with a scorecard of Gini1=.6") +
  theme_bw()

gini_combine_calculator(.6, .521, .7, .1)
gini_combine_calculator(.6, .6, .7, .1)
gini_combine_calculator(.6, .474, .6, .1)



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

gini_combine_calculator(.6, .64, .83, .1)
gini_combine_calculator(.6, .41, .26, .1)

direct.label(ggplot(df[round(df$Var1,2) %in% round(0:10/10,2),], aes(x=Var2, y=res)) + 
  geom_line(aes(group = factor(Var1), col=Var1)), "get.means") + xlab("Correlation") +
  ylab("New Gini")




# for (i in 1:dim(e)[1]){
#   print(e[i,])
#   print(gini_combine_calculator(.6, e[i,1], e[i,2], .1))
# }
# 
# gini_combine_calculator(.6, 0.9, 0, .1)

# Zadanie: znajdź korelację taką, żeby po doklejeniu karty g2=.8, do karty scoringowej g1=.6 otrzymać gnew=.65 (bad rate=.1)



r=0:9/10
g=1:9/10
f=expand.grid(g,r)
tmp2<-Vectorize(function(g, r){gini_combine_calculator(.4, g, r, .1)[[1]]})
res<-tmp2(f$Var1, f$Var2)

f2<-cbind(f,res)
#write.csv2(f2,'f2_v1.csv')

plot2<-ggplot(f2, aes(x = Var1, y = Var2, z = res)) +
  stat_contour(aes(colour = ..level..)) + xlab('Gini 2') + ylab('Correlation') 
direct.label(plot2, 'get.means')
