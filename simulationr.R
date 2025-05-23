# random  correlation 0-.99
# random gini1 0.1-.9
# random gini2 .1 - .9
# random bad rate from .01 - .5


# store: 
# input_corr
# input_gini1
# input_gini2
# input_badrate
# theor_gini_combined
# actual_gini1
# actual_gini2
# actual_corr
# actual_bad_rate
# gini_model
# gini_logistic
# gini_probit

library(MASS)

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

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
  tryCatch(2*integrate(function(x){F1_(x)*F2_(x)}, 
                       lower=-Inf, upper=Inf, subdivisions=500, rel.tol = .Machine$double.eps^.24)$value/defrate/(1-defrate)-1, error=function(err) NA)
}

r_from_gini<-function(gini, defaultrate=.1){
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-gini}
  uniroot(phi_s1,lower=0,upper=.9999,tol = .Machine$double.eps)$root
}

gini_combine_calculator<-function(g1, g2, corr, defaultrate){
  #rho_s1
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g1}
  rho_s1<-uniroot(phi_s1,lower=0,upper=.9999,tol = .Machine$double.eps)$root
  #rho_s2
  phi_s2<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g2}
  rho_s2<-uniroot(phi_s2,lower=0,upper=.9999,tol = .Machine$double.eps)$root
  
  (a_opt<-(corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2))
  corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
  g_result0<-if(corr_opt>1) {NaN} else {gini_from_r(corr_opt, defaultrate)}
  g_result1<-if(a_opt<0 | a_opt>1000 | corr_opt>1) {NaN} else {g_result0}
  return(c(new_gini=g_result1, 
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt, gini1=g1, gini2=g2, corr=corr, defrate=defaultrate))
}

gini_combine_1<-function(g1, g2, corr, defaultrate){
  tryCatch(unname(gini_combine_calculator(g1, g2, corr, defaultrate)[1]), error=function(err) NA)}

gini_combine_1<-Vectorize(gini_combine_1)

myseed <- 125
set.seed(myseed)

number_of_sim<-1000
sample_size<-1000

# random  correlation 0-.99
# random gini1 0.1-.9
# random gini2 .1 - .9
# random bad rate from .01 - .5
input_corr<-runif(number_of_sim, 0, .99)
input_gini1<-runif(number_of_sim, 0.1, .9)
input_gini2<-runif(number_of_sim, 0.1, .9)
input_badrate<-runif(number_of_sim, 0.01, .5)
theor_gini_combined<-gini_combine_1(input_gini1, input_gini2, input_corr, input_badrate)

#x=2
#gini_combine_1(input_gini1[x], input_gini2[x], input_corr[x], input_badrate[x])
#gini_combine_1(0.421377, 0.887435, 0.954183, 0.2187387)


actual_gini1<-c()
actual_gini2<-c()
actual_corr<-c()
actual_rho<-c()
actual_rho2<-c()
actual_bad_rate<-c()
gini_model<-c()
gini_model_rho<-c()
gini_model_rho2<-c()
gini_logistic<-c()
gini_probit<-c()


for (i in 1:number_of_sim){
  print(i)
  gini1<-input_gini1[i]
  gini2<-input_gini2[i]
  corrs<-input_corr[i]
  dfrate<-input_badrate[i]
  if (is.na(theor_gini_combined[i])){
    actual_gini1[i]<-NA
    actual_gini2[i]<-NA
    actual_corr[i]<-NA
    actual_rho[i]<-NA
    actual_rho2[i]<-NA
    actual_bad_rate[i]<-NA
    gini_model[i]<-NA
    gini_model_rho[i]<-NA
    gini_model_rho2[i]<-NA    
    gini_logistic[i]<-NA
    gini_probit[i]<-NA
    next
      }
  sigma <- matrix(c(1, corrs, r_from_gini(gini1, dfrate),
                    corrs, 1, r_from_gini(gini2, dfrate),
                    r_from_gini(gini1, dfrate), r_from_gini(gini2, dfrate), 1), 
                  nrow=3)
  
  if(any(eigen(sigma)$values<=0)){
    actual_gini1[i]<-NA
    actual_gini2[i]<-NA
    actual_corr[i]<-NA
    actual_bad_rate[i]<-NA
    gini_model[i]<-NA
    gini_logistic[i]<-NA
    gini_probit[i]<-NA
    next    
  }
  
  z <- mvrnorm(sample_size,mu=rep(0, 3),Sigma=sigma)
  df<-z
  df[,3]<-(z[,3]<qnorm(dfrate))*1
  df<-data.frame(df)
  names(df)<-c('s1', 's2', 'default_flag')
  actual_gini1[i]<-my_gini(df$default_flag, df$s1)
  actual_gini2[i]<-my_gini(df$default_flag, df$s2)
  actual_corr[i]<-cor(df$s1, df$s2)
  actual_rho[i]<-cor(df$s1, df$s2, method="spearman")
  actual_rho2[i]<-2*sin(actual_rho[i]*pi/6)
  actual_bad_rate[i]<-mean(df$default_flag)
  if (any(c(actual_gini1[i]<=0, actual_gini2[i]<=0, actual_corr[i]<=0, actual_rho[i]<=0, actual_bad_rate[i]<=0))){
    gini_model[i]<-NA
    gini_logistic[i]<-NA
    gini_probit[i]<-NA
    next
  }
  gini_model[i]<-gini_combine_1(actual_gini1[i], actual_gini2[i], actual_corr[i], actual_bad_rate[i])
  gini_model_rho[i]<-gini_combine_1(actual_gini1[i], actual_gini2[i], actual_rho[i], actual_bad_rate[i])
  gini_model_rho2[i]<-gini_combine_1(actual_gini1[i], actual_gini2[i], actual_rho2[i], actual_bad_rate[i])
  mylogit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = "binomial")
  myprobit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = binomial(link="probit"))
  gini_logistic[i]<-my_gini(mylogit$y, -mylogit$fitted.values)
  gini_probit[i]<-my_gini(myprobit$y, -myprobit$fitted.values)
}

diff.probit<-gini_probit-gini_model
diff.logistic<-gini_logistic-gini_model

results<-data.frame(input_corr
  , input_gini1
  , input_gini2
  , input_badrate
  , theor_gini_combined
  , actual_gini1
  , actual_gini2
  , actual_corr
  , actual_rho
  , actual_rho2
  , actual_bad_rate
  , gini_model
  , gini_model_rho
  , gini_model_rho2
  , gini_logistic
  , gini_probit
  , diff.logistic
  , diff.probit)

summary(results$theor_gini_combined-results$gini_model)
summary(results$theor_gini_combined-results$gini_logistic)

#View(results)
write.csv(results, paste0('simresults', '_', myseed, '_', sample_size, '_', number_of_sim, '_', format(Sys.time(), "%Y%m%d%H%M%S"), '.csv')
)
sd(diff.logistic, na.rm=TRUE)
sd(diff.probit, na.rm=TRUE)

res40000<-read.csv("simresults4000020220428233229.csv")
res10000<-read.csv("simresults1000020220428223014.csv")
res1000<-read.csv("simresults100020220428212247.csv")
res100000<-read.csv("simresults10000020220429000403.csv")
# res40000<-read.csv("simresults4000020220428233229.csv")
# res10000<-read.csv("simresults1000020220428223014.csv")
# res1000<-read.csv("simresults100020220428212247.csv")
# res100000<-read.csv("simresults10000020220429000403.csv")

res1000 <- read.csv("simresults_125_1000_1000_20250331193250.csv")
res10000 <- read.csv("simresults1000020250331191346.csv")
res100000 <- read.csv("simresults10000020250331185114.csv")

hist(res1000$diff.logistic, freq=FALSE, ylim=c(0,600), breaks=(-35:35)/1000)
hist(res10000$diff.logistic, col=rgb(1,0,0,.3), freq=FALSE, add=TRUE, breaks=(-35:35)/1000)
hist(res40000$diff.logistic, col=rgb(0,0,1,.3), freq=FALSE,add=TRUE, breaks=(-70:70)/2000)
hist(res100000$diff.logistic, col=rgb(0,1,0,.3), freq=FALSE,add=TRUE, breaks=(-70:70)/2000)

sd(res1000$diff.logistic, na.rm=TRUE)
sd(res10000$diff.logistic, na.rm=TRUE)
sd(res40000$diff.logistic, na.rm=TRUE)
sd(res100000$diff.logistic, na.rm=TRUE)

boxplot(list(res1000$diff.logistic[!is.na(res1000$diff.logistic)],
     res10000$diff.logistic[!is.na(res10000$diff.logistic)],
     res100000$diff.logistic[!is.na(res100000$diff.logistic)]),
     horizontal=TRUE
)

mean(res1000$theor_gini_combined-res1000$gini_model, na.rm=TRUE)
mean(res10000$theor_gini_combined-res10000$gini_model, na.rm=TRUE)
mean(res100000$theor_gini_combined-res100000$gini_model, na.rm=TRUE)

mean(res1000$theor_gini_combined-res1000$gini_logistic, na.rm=TRUE)
mean(res10000$theor_gini_combined-res10000$gini_logistic, na.rm=TRUE)
mean(res100000$theor_gini_combined-res100000$gini_logistic, na.rm=TRUE)

mean(res1000$gini_logistic-res1000$gini_model, na.rm=TRUE)
mean(res10000$gini_logistic-res10000$gini_model, na.rm=TRUE)
mean(res100000$gini_logistic-res100000$gini_model, na.rm=TRUE)


sd(res1000$theor_gini_combined-res1000$gini_model, na.rm=TRUE)
sd(res1000$theor_gini_combined-res1000$gini_logistic, na.rm=TRUE)
sd(res1000$gini_model-res1000$gini_logistic, na.rm=TRUE)

sd(res10000$theor_gini_combined-res10000$gini_model, na.rm=TRUE)
sd(res10000$theor_gini_combined-res10000$gini_logistic, na.rm=TRUE)
sd(res10000$gini_model-res10000$gini_logistic, na.rm=TRUE)

sd(res100000$theor_gini_combined-res100000$gini_model, na.rm=TRUE)
sd(res100000$theor_gini_combined-res100000$gini_logistic, na.rm=TRUE)
sd(res100000$gini_model-res100000$gini_logistic, na.rm=TRUE)

auc_from_gini<-function(x){x/2+.5}
auc_from_gini<-Vectorize(auc_from_gini)

gini_conf_width<-Vectorize(function(auc, brate, n){if(any(is.na(brate), is.na(auc))) {NA} else 2*{presize::prec_auc(auc, brate, n)$conf.width}})
mean(gini_conf_width(auc_from_gini(res1000$gini_logistic), res1000$actual_bad_rate, 1000), na.rm=TRUE)/(2*1.96)

sum(res1000$gini_logistic-res1000$gini_model, na.rm=TRUE)/
  mean(res1000$gini_logistic-res1000$gini_model, na.rm=TRUE)
sum(res10000$gini_logistic-res10000$gini_model, na.rm=TRUE)/
  mean(res10000$gini_logistic-res10000$gini_model, na.rm=TRUE)
sum(res100000$gini_logistic-res100000$gini_model, na.rm=TRUE)/
  mean(res100000$gini_logistic-res100000$gini_model, na.rm=TRUE)
dim(res100000)
dim(res10000)
dim(res1000)

summary_table<-data.frame(n = c(1e3, 1e4, 1e5), 
                          mnv_log_m = c(mean(res1000$gini_model-res1000$gini_logistic, na.rm=TRUE),
                                       mean(res10000$gini_model-res10000$gini_logistic, na.rm=TRUE),
                                       mean(res100000$gini_model-res100000$gini_logistic, na.rm=TRUE)),
                          mnv_log_sd = c(sd(res1000$gini_model-res1000$gini_logistic, na.rm=TRUE),
                                        sd(res10000$gini_model-res10000$gini_logistic, na.rm=TRUE),
                                        sd(res100000$gini_model-res100000$gini_logistic, na.rm=TRUE)),
                          log_the_m = c(mean(res1000$gini_logistic-res1000$theor_gini_combined, na.rm=TRUE),
                                       mean(res10000$gini_logistic-res10000$theor_gini_combined, na.rm=TRUE),
                                       mean(res100000$gini_logistic-res100000$theor_gini_combined, na.rm=TRUE)),
                          log_the_sd = c(sd(res1000$gini_logistic-res1000$theor_gini_combined, na.rm=TRUE),
                                        sd(res10000$gini_logistic-res10000$theor_gini_combined, na.rm=TRUE),
                                        sd(res100000$gini_logistic-res100000$theor_gini_combined, na.rm=TRUE)),
                          mnv_the_m = c(mean(res1000$gini_model-res1000$theor_gini_combined, na.rm=TRUE),
                                        mean(res10000$gini_model-res10000$theor_gini_combined, na.rm=TRUE),
                                        mean(res100000$gini_model-res100000$theor_gini_combined, na.rm=TRUE)),
                          mnv_the_sd = c(sd(res1000$gini_model-res1000$theor_gini_combined, na.rm=TRUE),
                                         sd(res10000$gini_model-res10000$theor_gini_combined, na.rm=TRUE),
                                         sd(res100000$gini_model-res100000$theor_gini_combined, na.rm=TRUE)),
                          mnvs_the_m = c(mean(res1000$gini_model_rho-res1000$theor_gini_combined, na.rm=TRUE),
                                         mean(res10000$gini_model_rho-res10000$theor_gini_combined, na.rm=TRUE),
                                         mean(res100000$gini_model_rho-res100000$theor_gini_combined, na.rm=TRUE)),
                          mnvs_the_sd = c(sd(res1000$gini_model_rho-res1000$theor_gini_combined, na.rm=TRUE),
                                          sd(res10000$gini_model_rho-res10000$theor_gini_combined, na.rm=TRUE),
                                          sd(res100000$gini_model_rho-res100000$theor_gini_combined, na.rm=TRUE)),
                          mnvs2_the_m = c(mean(res1000$gini_model_rho2-res1000$theor_gini_combined, na.rm=TRUE),
                                         mean(res10000$gini_model_rho2-res10000$theor_gini_combined, na.rm=TRUE),
                                         mean(res100000$gini_model_rho2-res100000$theor_gini_combined, na.rm=TRUE)),
                          mnvs2_the_sd = c(sd(res1000$gini_model_rho2-res1000$theor_gini_combined, na.rm=TRUE),
                                          sd(res10000$gini_model_rho2-res10000$theor_gini_combined, na.rm=TRUE),
                                          sd(res100000$gini_model_rho2-res100000$theor_gini_combined, na.rm=TRUE)),
                          
                          gini_cf_se = c(
                            mean(gini_conf_width(auc_from_gini(res1000$gini_logistic), res1000$actual_bad_rate, 1000), na.rm=TRUE)/(2* qnorm(.975)),
                            mean(gini_conf_width(auc_from_gini(res10000$gini_logistic), res10000$actual_bad_rate, 10000), na.rm=TRUE)/(2* qnorm(.975)),
                            mean(gini_conf_width(auc_from_gini(res100000$gini_logistic), res100000$actual_bad_rate, 100000), na.rm=TRUE)/(2* qnorm(.975))
                          ))
options(scipen=999)
View(summary_table)

library(knitr)
kable(summary_table)
print(summary_table)
