ginic<-function(bc, gc){
  #function for gini when we have cumulative goods and bads vectors
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*
        (bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1
} 

gini_crd<-function(rho=0.5, defrate=0.1, gran=10000) {
  #function translating rho into gini for a given defrate
  drates_i<-pnorm(qnorm(defrate), rho*(qnorm(1:(gran-1)/gran)),sqrt(1-rho^2))
  drates_2i<-(c(1, drates_i)+c(drates_i, 0))/2
  return(ginic(cumsum((drates_2i)/sum(drates_2i)), cumsum((1-drates_2i)/sum(1-drates_2i))))
}

gini_combine_calculator<-function(g1, g2, corr, defaultrate){
  #rho_s1
  phi_s1<-function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-g1}
  rho_s1<-uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
  #rho_s2
  phi_s2<-function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-g2}
  rho_s2<-uniroot(phi_s2,lower=0,upper=1,tol = .Machine$double.eps)$root
  
  (a_opt<-(corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2))
  corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
  g_result0<-gini_crd(corr_opt, defaultrate, gran=100000)
  g_result1<-if(a_opt<0 | a_opt>1000) {NaN} else {g_result0}
  return(c(new_gini=g_result1, 
           #new_gini_2=g_result0, 
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt))
}

gini_combine_calculator(.6, .4, .1, .04)

# generate the simulated data (... many ways..., assume multivariate distribution)
# build the logistic regression on the two scorecards and compare with the calculator


cors<-c((-2:9)/10,.99)
ginic4<-function(x){gini_combine_calculator(.6,.6,x,.1)[1]}
ginic4<-Vectorize(ginic4)
ginis<-ginic4(cors)
plot(cors, ginis)
library(ggplot2)
ggplot(data.frame(cors, ginis), aes(x=cors, y=ginis)) + 
  geom_point() +
  geom_text(label=paste0('r=', round(cors,2), ' Gini=', round(ginis,3)), hjust=0.4, vjust=1) +
  xlab("Correlation") + ylab("Gini") + ggtitle("Combined Gini of two scorecards with Gini=0.6 depending on the correlation between them")


(target1<-unname(gini_combine_calculator(.6, .6, .7, .1)[1]))
ginis1<-c(.6)
corrs1<-c(.7)

for (i in c(.6,.5,.4,.3,.2,.1, 0, -.1, -.2)){
corrs1[length(corrs1)+1]<-i
search1<-function(x){gini_combine_calculator(.6, x, corrs1[length(corrs1)], .1)[1]-target1}
ginis1[length(ginis1)+1]<-uniroot(search1,lower=ginis1[length(ginis1)]-.1,upper=ginis1[length(ginis1)],tol = .Machine$double.eps)$root
}

df<-data.frame(corrs1, ginis1)
library(ggplot2)
ggplot(df, aes(x=corrs1, y=ginis1)) + geom_point() +
  geom_text(label=paste0('r=', round(corrs1,2), ' Gini2=', round(ginis1,3)), hjust=0.4, vjust=1)+
  xlab("Correlation") + ylab("Gini") + ggtitle("Gini/correlation trade-off: all these combinations have the same effect (Combined Gini=.648) when combining with a scorecard of Gini1=.6") +
  theme_bw()

#View(df)
#plot(df)

#saul new
# Sample	N	Region	AUC_Bureau	AUC_Psych	AUC_Combined	Bad_Rate	r_Bureau_Psych	rho_Bureau_Psych
# 1	4166	S. America	0,678	0,664	0,729	4,25%	0,099	0,092
# 2	721	S. America	0,694	0,727	0,784	4,85%	0,104	0,088
# 3	1306	Europe	0,642	0,627	0,688	16,77%	0,142	0,135
# 4	3236	S. America	0,715	0,642	0,741	17,27%	0,135	0,129

# gini_from_auc<-function(x){(x-.5)*2}

#which is which

# # 1	4166	S. America	0,678	0,664	0,729	4,25%	0,099	0,092
# 
# #case 1
# (gini1<-gini_from_auc(.678))
# (gini2<-gini_from_auc(.664))
# (gini3<-gini_from_auc(.729))
# (res<-gini_combine_calculator(gini1, gini2, 0.099, .0425))
# res[1]/2+.5
# (res<-gini_combine_calculator(gini1, gini2, 0.092, .0425))
# res[1]/2+.5
# 
# # 2	721	S. America	0,694	0,727	0,784	4,85%	0,104	0,088
# 
# #case 3
# (gini1<-gini_from_auc(.694))
# (gini2<-gini_from_auc(.727))
# (res<-gini_combine_calculator(gini1, gini2, 0.104, .0485))
# res[1]/2+.5
# (res<-gini_combine_calculator(gini1, gini2, 0.088, .0485))
# res[1]/2+.5
# 
# # 3	1306	Europe	0,642	0,627	0,688	16,77%	0,142	0,135
# 
# #case 3
# (gini1<-gini_from_auc(.642))
# (gini2<-gini_from_auc(.627))
# (res<-gini_combine_calculator(gini1, gini2, 0.142, .1677))
# res[1]/2+.5
# (res<-gini_combine_calculator(gini1, gini2, 0.135, .1677))
# res[1]/2+.5
# 
# # 4	3236	S. America	0,715	0,642	0,741	17,27%	0,135	0,129
# 
# (gini1<-gini_from_auc(.715))
# (gini2<-gini_from_auc(.642))
# (res<-gini_combine_calculator(gini1, gini2, 0.135, .1727))
# res[1]/2+.5
# (res<-gini_combine_calculator(gini1, gini2, 0.129, .1727))
# res[1]/2+.5



