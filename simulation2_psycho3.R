#data<-readxl::read_excel('C:/Users/Błażej/Downloads/data.xlsx')
#data<-data[-c(1),]
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


gini_from_auc<-function(x){(x-.5)*2}
gini_from_auc<-Vectorize(gini_from_auc)

auc_from_gini<-function(x){x/2+.5}
auc_from_gini<-Vectorize(auc_from_gini)


data$Gini_Bureau<-gini_from_auc(data$AUC_Bureau)
data$Gini_Psych<-gini_from_auc(data$AUC_Psych)
data$Gini_Combined<-gini_from_auc(data$AUC_Combined)

gini_predict<-function(g1, g2, corr, defaultrate){unname(gini_combine_calculator(g1, g2, corr, defaultrate)[1])}
gini_predict<-Vectorize(gini_predict)


r_from_gini<-function(gini, defaultrate=.1){
  phi_s1<-function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-gini}
  uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
}

for (j in 1:1){
  print(j)
  n=data$N[j]
  sigma <- matrix(c(1, data$r_Bureau_Psych[j], r_from_gini(gini_from_auc(data$AUC_Bureau[j]), data$Bad_Rate[j]),
                    data$r_Bureau_Psych[j], 1, r_from_gini(gini_from_auc(data$AUC_Psych[j]), data$Bad_Rate[j]),
                    r_from_gini(gini_from_auc(data$AUC_Bureau[j]), data$Bad_Rate[j]), r_from_gini(gini_from_auc(data$AUC_Psych[j]), data$Bad_Rate[j]), 1), 
                      nrow=3)
  
  }


#auc_from_gini(gini_crd(sigma[1,3], defrate = data$Bad_Rate[1]))
#r_from_gini(gini_from_auc(data$AUC_Bureau[1]), data$Bad_Rate[1])

for (i in 1:100){
  z <- mvrnorm(n,mu=rep(0, 3),Sigma=sigma)
  df<-z
  df[,3]<-(z[,3]<qnorm(dfrate))*1
  df<-data.frame(df)
  names(df)<-c('s1', 's2', 'default_flag')
  
  #my_gini(df$default_flag, df$s1)
  #my_gini(df$default_flag, df$s2)
  #cor(df$s1, df$s2)
  #mean(df$default_flag)
  res1<-gini_combine_calculator(my_gini(df$default_flag, df$s1), my_gini(df$default_flag, df$s2), cor(df$s2, df$s1), mean(df$default_flag))
  
}