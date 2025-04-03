
library(copula)
library(MASS)
library(splines)


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

t_cop <- tCopula(param =c( 0.6,.6,.6), dim = 3, dispstr = "un", df=4)
z <- rCopula(1e5, t_cop)
cor(z)
cor(z, method="spearman")

df<-z
dfrate <- .2
df[,3]<-(z[,3]<dfrate)*1
df<-data.frame(df)
names(df)<-c('s1', 's2', 'default_flag')
df$s1 <- qnorm(df$s1)
df$s2 <- qnorm(df$s2)

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

actual_gini1<-my_gini(df$default_flag, df$s1)
#2*bigstatsr::AUC(-df$s1, df$default_flag)-1
actual_gini2<-my_gini(df$default_flag, df$s2)
(actual_corr<-cor(df$s1, df$s2))
(actual_rho<-cor(df$s1, df$s2, method="spearman"))
(actual_rho2<-2*sin(actual_rho*pi/6))
actual_bad_rate<-mean(df$default_flag)


gini_model<-gini_combine_1(actual_gini1, actual_gini2, actual_corr, actual_bad_rate)
gini_model_rho<-gini_combine_1(actual_gini1, actual_gini2, actual_rho, actual_bad_rate)
gini_model_rho2<-gini_combine_1(actual_gini1, actual_gini2, actual_rho2, actual_bad_rate)
mylogit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = "binomial")
mylogs <- glm(I(1-default_flag) ~ ns(qnorm(pobs(s1)), df=2) + ns(qnorm(pobs(s2)), df=2), data = df, family = binomial)
gini_logistic<-my_gini(mylogit$y, -mylogit$fitted.values)
gini_logs<-my_gini(mylogs$y, -mylogs$fitted.values)


#validation:
t_cop <- tCopula(param =c( 0.6,.6,.6), dim = 3, dispstr = "un", df=4)
z_val <- rCopula(1e5, t_cop)
cor(z_val)
cor(z_val, method="spearman")
df_val<-z_val
dfrate <- .2
df_val[,3]<-(z_val[,3]<dfrate)*1
df_val<-data.frame(df_val)
names(df_val)<-c('s1', 's2', 'default_flag')
df_val$s1 <- qnorm(df_val$s1)
df_val$s2 <- qnorm(df_val$s2)
gini_logistic_val<-my_gini(df_val$default_flag, predict(mylogit, newdata = df_val))
gini_logs_val<-my_gini(df_val$default_flag, predict(mylogs, newdata = df_val))

plot(cumsum(1-df$default_flag[order(df$s1)])/sum(1-df$default_flag), 
     cumsum(df$default_flag[order(df$s1)])/sum(df$default_flag), type="l")

points(cumsum(1-df$default_flag[order(df$s2)])/sum(1-df$default_flag), 
     cumsum(df$default_flag[order(df$s2)])/sum(df$default_flag), type="l", col='blue')


