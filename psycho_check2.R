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

gini_combine_calculator<-function(g1, g2, corr, defaultrate){
  #rho_s1
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g1}
  rho_s1<-uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
  #rho_s2
  phi_s2<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g2}
  rho_s2<-uniroot(phi_s2,lower=0,upper=1,tol = .Machine$double.eps)$root
  
  (a_opt<-(corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2))
  corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
  g_result0<-gini_from_r(corr_opt, defaultrate)
  g_result1<-if(a_opt<0 | a_opt>1000) {NaN} else {g_result0}
  return(c(new_gini=g_result1, 
           #new_gini_2=g_result0, 
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt))
}


data<-readxl::read_excel('C:/Users/blaze/Downloads/data.xlsx')

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

data$AUC_predict_r<-auc_from_gini(data$Gini_Predict_r)
data$AUC_predict_rho<-auc_from_gini(data$Gini_Predict_rho)

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
  geom_point(aes(x=Gini_Bureau), col='dark green') +
  geom_text(aes(x=Gini_Bureau, label=paste('Bureau:', Gini_Bureau)), hjust=0.5, vjust=1.5) +
  geom_point(aes(x=Gini_Psych), col='green')+
  geom_text(aes(x=Gini_Psych, label=paste('Psych:', Gini_Psych)), hjust=0.5, vjust=-1) +
  geom_point(aes(x=Gini_Combined), col='black')+
  geom_text(aes(x=Gini_Combined, label=paste('Combined:', Gini_Combined)), hjust=0.5, vjust=-1) +
  geom_point(aes(x=Gini_Predict_r), col='blue')+
  geom_errorbar(aes(xmin = Gini_predict_r_lwr, xmax = Gini_predict_r_upr), col='blue', width=.1)+
  geom_text(aes(x=Gini_Predict_r, label=paste('Predicted (r):', round(Gini_Predict_r,3))), hjust=0, vjust=1.5, col='blue') +
#  geom_point(aes(x=Gini_Predict_rho), col='dark blue')+
#  geom_text(aes(x=Gini_Predict_rho, label=paste('Predicted (rho):', round(Gini_Predict_rho,3))), hjust=0, vjust=1.5, col='dark blue') +
  scale_y_discrete(limits=rev)+
  xlab('')+ylab('')+xlim(c(0.2, .75))




