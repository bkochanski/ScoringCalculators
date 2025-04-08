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


#data<-readxl::read_excel('C:/Users/blaze/Downloads/data.xlsx')
#data<-readxl::read_excel('C:/Users/Błażej/Downloads/data.xlsx')
data <- readxl::read_excel('D:/BK_dokumenty_D/googlepgedu/Combining/data.xlsx')
#data<-readxl::read_excel('/Users/blakocha/Downloads/data.xlsx')

gini_from_auc<-function(x){(x-.5)*2}
gini_from_auc<-Vectorize(gini_from_auc)

auc_from_gini<-function(x){x/2+.5}
auc_from_gini<-Vectorize(auc_from_gini)

#gini_from_auc(.7)
#auc_from_gini(.4)

data$Gini_Bureau<-gini_from_auc(data$AUC_Bureau)
data$Gini_Psych<-gini_from_auc(data$AUC_Psych)
data$Gini_Combined<-gini_from_auc(data$AUC_Combined)

#View(data)

gini_predict<-function(g1, g2, corr, defaultrate){unname(gini_combine_calculator(g1, g2, corr, defaultrate)[1])}
gini_predict<-Vectorize(gini_predict)

data$Gini_Predict_r<-gini_predict(data$Gini_Bureau, data$Gini_Psych, data$r_Bureau_Psych, data$Bad_Rate)
data$Gini_Predict_rho<-gini_predict(data$Gini_Bureau, data$Gini_Psych, data$rho_Bureau_Psych, data$Bad_Rate)
rho2r <- function(rho){2*sin(rho*pi/6)}
data$Gini_Predict_rho2<-gini_predict(data$Gini_Bureau, data$Gini_Psych, rho2r(data$rho_Bureau_Psych), data$Bad_Rate)

data$AUC_predict_r<-auc_from_gini(data$Gini_Predict_r)
data$AUC_predict_rho<-auc_from_gini(data$Gini_Predict_rho)
data$AUC_predict_rho2<-auc_from_gini(data$Gini_Predict_rho2)

AUC_lwr<-function(auc, brate, n){presize::prec_auc(auc, brate, n)$lwr}
AUC_lwr<-Vectorize(AUC_lwr)
AUC_upr<-function(auc, brate, n){presize::prec_auc(auc, brate, n)$upr}
AUC_upr<-Vectorize(AUC_upr)

data$AUC_predict_r_lwr<-AUC_lwr(data$AUC_predict_r, data$Bad_Rate, data$N)
data$AUC_predict_r_upr<-AUC_upr(data$AUC_predict_r, data$Bad_Rate, data$N)

data$Gini_predict_r_lwr<-gini_from_auc(data$AUC_predict_r_lwr)
data$Gini_predict_r_upr<-gini_from_auc(data$AUC_predict_r_upr)

#presize::prec_auc(auc_from_gini(.376), .1677, 1306)$lwr

#data$AUC_predict_rlwr<-auc_from_gini(data$Gini_Predict_r)


gini_from_auc(presize::prec_auc(auc_from_gini(.376), .1677, 1306)$lwr)
gini_from_auc(presize::prec_auc(auc_from_gini(.376), .1677, 1306)$upr)


library(ggplot2)
ggplot(data, aes(y=paste0(Sample, ': ', Region, '\nN=', N, '\nbad rate=', Bad_Rate*100, '%\ncorr=', r_Bureau_Psych))) +
  geom_point(aes(x=Gini_Bureau), col='darkgreen') +
  geom_text(aes(x=Gini_Bureau, label=paste('Bureau:', Gini_Bureau)), hjust=0.5, vjust=1.5, col='darkgreen') +
  geom_point(aes(x=Gini_Psych), col='darkcyan')+
  geom_text(aes(x=Gini_Psych, label=paste('Psych:', Gini_Psych)), hjust=0.5, vjust=-1, col='darkcyan') +
  geom_point(aes(x=Gini_Combined), col='black', shape=17, size=2)+
  geom_text(aes(x=Gini_Combined, label=paste('Log. reg.:', Gini_Combined)), hjust=0.5, vjust=-1) +
  geom_point(aes(x=Gini_Predict_r), col='blue', shape=17, size=2)+
  geom_errorbar(aes(xmin = Gini_predict_r_lwr, xmax = Gini_predict_r_upr), col='blue', width=.1)+
  geom_text(aes(x=Gini_Predict_r, label=paste('MVN calc.:', round(Gini_Predict_r,3))), hjust=0, vjust=1.5, col='blue') +
  #  geom_point(aes(x=Gini_Predict_rho), col='dark blue')+
  #  geom_text(aes(x=Gini_Predict_rho, label=paste('Predicted (rho):', round(Gini_Predict_rho,3))), hjust=0, vjust=1.5, col='dark blue') +
  scale_y_discrete(limits=rev)+
  xlab('')+ylab('')+xlim(c(0.2, .75)) +
  theme_bw()


ggplot(data, aes(y=paste0(Sample, ': ', Region, '\nN=', N, '\nbad rate=', Bad_Rate*100, '%\ncorr=', r_Bureau_Psych, '\nrho=', rho_Bureau_Psych))) +
  geom_point(aes(x=Gini_Bureau), col='darkgreen') +
  geom_text(aes(x=Gini_Bureau, label=paste('Bureau:', Gini_Bureau)), hjust=0.5, vjust=1.5, col='darkgreen') +
  geom_point(aes(x=Gini_Psych), col='darkcyan')+
  geom_text(aes(x=Gini_Psych, label=paste('Psych:', Gini_Psych)), hjust=0.5, vjust=-1, col='darkcyan') +
  geom_point(aes(x=Gini_Combined), col='black', shape=17, size=2)+
  geom_text(aes(x=Gini_Combined, label=paste('Log. reg.:', Gini_Combined)), hjust=0.5, vjust=-1) +
  geom_point(aes(x=Gini_Predict_r), col='blue', shape=17, size=2)+
  geom_errorbar(aes(xmin = Gini_predict_r_lwr, xmax = Gini_predict_r_upr), col='blue', width=.1)+
  geom_text(aes(x=Gini_Predict_r, label=paste('MVN calc.:', round(Gini_Predict_rho,3))), hjust=0, vjust=1.5, col='red') +
  geom_point(aes(x=Gini_Predict_rho), col='red')+
#  geom_text(aes(x=Gini_Predict_rho, label=paste('Predicted (rho):', round(Gini_Predict_rho,3))), hjust=0, vjust=1.5, col='red') +
  scale_y_discrete(limits=rev)+
  xlab('')+ylab('')+xlim(c(0.2, .75)) +
  theme_bw()



ggplot(data, aes(y=paste0(Sample, ': ', Region, '\nN=', N, '\nbad rate=', Bad_Rate*100, '%\ncorr=', r_Bureau_Psych, '\nrho=', rho_Bureau_Psych))) +
  geom_point(aes(x=Gini_Bureau), col='darkgreen') +
  geom_text(aes(x=Gini_Bureau, label=paste('Bureau:', Gini_Bureau)), hjust=0.5, vjust=1.5, col='darkgreen') +
  geom_point(aes(x=Gini_Psych), col='darkcyan')+
  geom_text(aes(x=Gini_Psych, label=paste('Psych:', Gini_Psych)), hjust=0.5, vjust=-1, col='darkcyan') +
  geom_point(aes(x=Gini_Combined), col='black', shape=17, size=2)+
  geom_text(aes(x=Gini_Combined, label=paste('Log. reg.:', Gini_Combined)), hjust=0.5, vjust=-1) +
  geom_point(aes(x=Gini_Predict_r), col='blue', shape=17, size=2)+
  geom_errorbar(aes(xmin = Gini_predict_r_lwr, xmax = Gini_predict_r_upr), col='blue', width=.1)+
  geom_text(aes(x=Gini_Predict_r, label=paste('MVN calc.:', round(Gini_Predict_rho2,3))), hjust=0, vjust=1.5, col='lightgreen') +
  geom_point(aes(x=Gini_Predict_rho2), col='lightgreen')+
  #  geom_text(aes(x=Gini_Predict_rho, label=paste('Predicted (rho):', round(Gini_Predict_rho,3))), hjust=0, vjust=1.5, col='red') +
  scale_y_discrete(limits=rev)+
  xlab('')+ylab('')+xlim(c(0.2, .75)) +
  theme_bw()

source('Gini_combined_master.R')

gini_combine_calculator <- Vectorize(gini_combine_calculator)
new_gini_1 <-gini_combine_calculator(data$Gini_Bureau, data$Gini_Psych, data$r_Bureau_Psych, data$Bad_Rate)[1,]
new_gini_2 <-gini_combine_calculator(data$Gini_Bureau, data$Gini_Psych, data$rho_Bureau_Psych, data$Bad_Rate, method="spearman")[1,]
data.frame(new_gini_1, compare1 = data$Gini_Predict_r, new_gini_2, compare2 = data$Gini_Predict_rho2)
round(new_gini_1-new_gini_2,4)
round(data$Gini_Predict_r-data$Gini_Predict_rho2,4)
sd(round(new_gini_1-new_gini_2,4))
sd(round(data$Gini_Predict_r-data$Gini_Predict_rho2,4))
