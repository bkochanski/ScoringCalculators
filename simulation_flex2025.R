library(MASS)
library(splines)
library(copula)

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

######## simul function #################

simul <- function(inp_seed = 125, inp_nsim = 999, inp_sample_size= 1e4, 
                  inp_mycop = "Normal", inp_mymarginals = "Normal"){
myseed <- inp_seed
set.seed(myseed)

number_of_sim <- inp_nsim

sample_size <- inp_sample_size #c(1e3, 1e4, 1e5)

mycop <- inp_mycop #"T5"; "T4"; "CGN"; "CBN"; "NNC"

mymarginals <- inp_mymarginals #"Normal" "Uniform"; "Skewed6"

input_corr<-runif(number_of_sim, 0, .99)
input_gini1<-runif(number_of_sim, 0.1, .9)
input_gini2<-runif(number_of_sim, 0.1, .9)
input_badrate<-runif(number_of_sim, 0.01, .5)
theor_gini_combined<-gini_combine_1(input_gini1, input_gini2, input_corr, input_badrate)

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
gini_logistic2<-c()
gini_logs<-c()
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
    gini_logistic2[i]<-NA
    gini_logs[i]<-NA
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
  
  if (mycop == "Normal") {
    norm_cop <- normalCopula(param = sigma[lower.tri(sigma)], dim = 3, dispstr = "un")
    z <- rCopula(sample_size, norm_cop)
  }
  
  if (mycop == "T5") {
    t_cop <- tCopula(param = sigma[lower.tri(sigma)], dim = 3, dispstr = "un", df=5)
    z <- rCopula(sample_size, t_cop)
    }

  if (mycop == "T4") {
    t_cop <- tCopula(param = sigma[lower.tri(sigma)], dim = 3, dispstr = "un", df=4)
    z <- rCopula(sample_size, t_cop)
  }
  
  
  df<-z
  df[,3]<-(z[,3]<dfrate)*1
  df<-data.frame(df)
  names(df)<-c('s1', 's2', 'default_flag')

  if (mymarginals == "Normal") {
    df$s1 <- qnorm(df$s1)
    df$s2 <- qnorm(df$s2)
  }
  
  if (mymarginals == "Skewed6") {
    df$s1 <- df$s1^6
    df$s2 <- df$s2^(1/6)
  }
  
  
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
  mylogs <- glm(I(1-default_flag) ~ ns(qnorm(pobs(s1)), df=2) + ns(qnorm(pobs(s2)), df=2), data = df, family = binomial)
  gini_logistic[i]<-my_gini(mylogit$y, -mylogit$fitted.values)
  gini_probit[i]<-my_gini(myprobit$y, -myprobit$fitted.values)
  gini_logs[i]<-my_gini(mylogs$y, -mylogs$fitted.values)
  
}

diff.probit<-gini_probit-gini_model
diff.logistic<-gini_logistic-gini_model
diff.logistic2<-gini_logistic-gini_model_rho2
diff.logs2<-gini_logs-gini_model_rho2

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
                    , gini_logs
                    , gini_probit
                    , diff.logistic
                    , diff.logistic2
                    , diff.logs2
                    , diff.probit)

#View(results)
write.csv(results, paste0('simresults_flex', '_', mycop, '_', mymarginals, '_', sample_size, '_', number_of_sim, '_', myseed, '_', format(Sys.time(), "%Y%m%d%H%M%S"), '.csv')
          , row.names = FALSE)

return(results)}

results_1e3_NS <- simul(125, 444, 1e3, "Normal", "Skewed6")
results_1e4_NS <- simul(125, 444, 1e4, "Normal", "Skewed6")
results_1e5_NS <- simul(125, 444, 1e5, "Normal", "Skewed6")
results_1e3_Nu <- simul(125, 444, 1e3, "Normal", "Uniform")
results_1e4_Nu <- simul(125, 444, 1e4, "Normal", "Uniform")
results_1e5_Nu <- simul(125, 444, 1e5, "Normal", "Uniform")
results_1e3_Nn <- simul(125, 444, 1e3, "Normal", "Normal")
results_1e4_Nn <- simul(125, 444, 1e4, "Normal", "Normal")
results_1e5_Nn <- simul(125, 444, 1e5, "Normal", "Normal")
results_1e3_Tn <- simul(125, 444, 1e3, "T4", "Normal")
results_1e4_Tn <- simul(125, 444, 1e4, "T4", "Normal")
results_1e5_Tn <- simul(125, 444, 1e5, "T4", "Normal")
results_1e3_Tu <- simul(125, 444, 1e3, "T4", "Uniform")
results_1e4_Tu <- simul(125, 444, 1e4, "T4", "Uniform")
results_1e5_Tu <- simul(125, 444, 1e5, "T4", "Uniform")
results_1e3_Ts <- simul(125, 444, 1e3, "Normal", "Skewed6")
results_1e4_Ts <- simul(125, 444, 1e4, "Normal", "Skewed6")
results_1e5_Ts <- simul(125, 444, 1e5, "Normal", "Skewed6")



res1000 <- read.csv("simresults_flex_Normal_Normal_1000_999_125_20250402213447.csv")
res10000 <- read.csv("simresults_flex_Normal_Normal_10000_999_125_20250402220915.csv")
res100000 <- read.csv("simresults_flex_Normal_Normal_100000_999_125_20250403010128.csv")  
res10000 <-results



gini_conf_width<-Vectorize(function(auc, brate, n){if(any(is.na(brate), is.na(auc))) {NA} else 2*{presize::prec_auc(auc, brate, n)$conf.width}})
auc_from_gini<-Vectorize(function(x){x/2+.5})
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
                            mnvs2_logs_m = c(mean(res1000$gini_model_rho2-res1000$gini_logs, na.rm=TRUE),
                                             mean(res10000$gini_model_rho2-res10000$gini_logs, na.rm=TRUE),
                                             mean(res100000$gini_model_rho2-res100000$gini_logs, na.rm=TRUE)),
                            mnvs2_logs_sd = c(sd(res1000$gini_model_rho2-res1000$gini_logs, na.rm=TRUE),
                                             sd(res10000$gini_model_rho2-res10000$gini_logs, na.rm=TRUE),
                                             sd(res100000$gini_model_rho2-res100000$gini_logs, na.rm=TRUE)),
                          logs_the_m = c(mean(res1000$gini_logs-res1000$theor_gini_combined, na.rm=TRUE),
                                          mean(res10000$gini_logs-res10000$theor_gini_combined, na.rm=TRUE),
                                          mean(res100000$gini_logs-res100000$theor_gini_combined, na.rm=TRUE)),
                          logs_the_sd = c(sd(res1000$gini_logs-res1000$theor_gini_combined, na.rm=TRUE),
                                           sd(res10000$gini_logs-res10000$theor_gini_combined, na.rm=TRUE),
                                           sd(res100000$gini_logs-res100000$theor_gini_combined, na.rm=TRUE)),
                          gini_cf_se = c(
                              mean(gini_conf_width(auc_from_gini(res1000$gini_logistic), res1000$actual_bad_rate, 1000), na.rm=TRUE)/(2* qnorm(.975)),
                              mean(gini_conf_width(auc_from_gini(res10000$gini_logistic), res10000$actual_bad_rate, 10000), na.rm=TRUE)/(2* qnorm(.975)),
                              mean(gini_conf_width(auc_from_gini(res100000$gini_logistic), res100000$actual_bad_rate, 100000), na.rm=TRUE)/(2* qnorm(.975))
                            ))
options(scipen=999)
View(summary_table)

