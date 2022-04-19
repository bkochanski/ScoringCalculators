data<-readxl::read_excel('C:/Users/Błażej/Downloads/data.xlsx')
data<-data[-c(1),]
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
              lower=-Inf, upper=Inf, subdivisions=500, rel.tol = .Machine$double.eps^.24)$value/defrate/(1-defrate)-1
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
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt, gini1=g1, gini2=g2, corr=corr, defrate=defaultrate))
}

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 


gini_from_auc<-function(x){(x-.5)*2}
gini_from_auc<-Vectorize(gini_from_auc)

auc_from_gini<-function(x){x/2+.5}
auc_from_gini<-Vectorize(auc_from_gini)


data$Gini_Bureau<-gini_from_auc(data$AUC_Bureau)
data$Gini_Psych<-gini_from_auc(data$AUC_Psych)
data$Gini_Combined<-gini_from_auc(data$AUC_Combined)

gini_predict<-function(g1, g2, corr, defaultrate){unname(gini_combine_calculator(g1, g2, corr, defaultrate)[1])}
gini_predict<-Vectorize(gini_predict)

data$Gini_Predicted<-gini_predict(data$Gini_Bureau, data$Gini_Psych,
                                  data$r_Bureau_Psych, data$Bad_Rate)


r_from_gini<-function(gini, defaultrate=.1){
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-gini}
  uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
}

tot_results<-list()

for (j in 1:5){
  print(j)
  n=data$N[j]
  dfrate=data$Bad_Rate[j]
  sigma <- matrix(c(1, data$r_Bureau_Psych[j], r_from_gini(gini_from_auc(data$AUC_Bureau[j]), dfrate),
                    data$r_Bureau_Psych[j], 1, r_from_gini(gini_from_auc(data$AUC_Psych[j]), dfrate),
                    r_from_gini(gini_from_auc(data$AUC_Bureau[j]), dfrate), r_from_gini(gini_from_auc(data$AUC_Psych[j]), dfrate), 1), 
                      nrow=3)
  
#auc_from_gini(gini_from_r(sigma[1,3], defrate = data$Bad_Rate[1]))
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
  mylogit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = "binomial")
  logit1_a<-mylogit$coefficients[2]/mylogit$coefficients[3]
  gini_logit1<-my_gini(df$default_flag, logit1_a*df$s1+df$s2)
  tot_results[[length(tot_results)+1]]<-c(j=j, 
                                          i=i, 
                                          n=n, 
                                          res1,
                                          gini_logit1=gini_logit1)
}}
  
bind_tot_results<-do.call("rbind", tot_results)
df_tot_results<-data.frame(matrix(unlist(bind_tot_results), ncol=length(tot_results[[1]])))
names(df_tot_results)<-names(tot_results[[1]])
df_tot_results$gini1_increase<-df_tot_results$new_gini-df_tot_results$gini1
df_tot_results$gini2_increase<-df_tot_results$new_gini-df_tot_results$gini2
df_tot_results$gini1_increl<-df_tot_results$new_gini/df_tot_results$gini1-1
df_tot_results$gini2_increl<-df_tot_results$new_gini/df_tot_results$gini2-1

#View(df_tot_results)
sd((df_tot_results$gini1_increase)[df_tot_results$j==4])
aggregate(df_tot_results$gini1_increase, sd, by=list(df_tot_results$j))
aggregate(df_tot_results$gini2_increase, sd, by=list(df_tot_results$j))
aggregate(df_tot_results$gini1_increl, sd, by=list(df_tot_results$j))
aggregate(df_tot_results$gini2_increl, sd, by=list(df_tot_results$j))

aggregate(df_tot_results$gini1, mean, by=list(df_tot_results$j))
data$Gini_Bureau
aggregate(df_tot_results$gini2, mean, by=list(df_tot_results$j))
data$Gini_Psych
aggregate(df_tot_results$new_gini, mean, by=list(df_tot_results$j))
aggregate(df_tot_results$new_gini, sd, by=list(df_tot_results$j))
aggregate(df_tot_results$gini1_increase, sd, by=list(df_tot_results$j))
data$Gini_Predicted

# for (i in 201:210){
# print(gini_predict(df_tot_results$gini1[i], df_tot_results$gini2[i],
#              df_tot_results$corr[i], df_tot_results$defrate[i])-df_tot_results$new_gini[i])
# }
aggregate((df_tot_results$gini_logit1-df_tot_results$new_gini), mean, by=list(df_tot_results$j))
aggregate((df_tot_results$gini_logit1-df_tot_results$new_gini), median, by=list(df_tot_results$j))
aggregate((df_tot_results$gini_logit1-df_tot_results$new_gini), sd, by=list(df_tot_results$j))
hist(I(df_tot_results$gini_logit1-df_tot_results$new_gini)[df_tot_results$j==3], breaks=30)
boxplot(I(df_tot_results$gini_logit1-df_tot_results$new_gini)~df_tot_results$j)
data$N
data$r_Bureau_Psych

df_tot_results[df_tot_results$j==1,]
