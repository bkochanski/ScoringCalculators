library(mgcv) 
library(splines)

sample_size=100000
my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

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
gini1<-.6
gini2<-.6
corrs<-.5
dfrate<-.2
sigma <- matrix(c(1, corrs, r_from_gini(gini1, dfrate),
                  corrs, 1, r_from_gini(gini2, dfrate),
                  r_from_gini(gini1, dfrate), r_from_gini(gini2, dfrate), 1), 
                nrow=3)

z <- MASS::mvrnorm(sample_size,mu=rep(0, 3),Sigma=sigma)
df<-z
df[,3]<-(z[,3]<qnorm(dfrate))*1
df<-data.frame(df)
names(df)<-c('s1', 's2', 'default_flag')
df$s1 <- pnorm(df$s1)^6
df$s2 <- pnorm(df$s2)^(1/6)
actual_gini1<-my_gini(df$default_flag, df$s1)
actual_gini2<-my_gini(df$default_flag, df$s2)
actual_corr<-cor(df$s1, df$s2)
actual_rho<-cor(df$s1, df$s2, method="spearman")
actual_rho2<-2*sin(actual_rho*pi/6)
actual_bad_rate<-mean(df$default_flag)
gini_model<-gini_combine_1(actual_gini1, actual_gini2, actual_corr, actual_bad_rate)
gini_model_rho<-gini_combine_1(actual_gini1, actual_gini2, actual_rho, actual_bad_rate)
gini_model_rho2<-gini_combine_1(actual_gini1, actual_gini2, actual_rho2, actual_bad_rate)
mylogit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = "binomial")
#mytree <- rpart(I(1-default_flag) ~ s1 + s2, data = df, method = "class")
#mygam <- gam(I(1-default_flag) ~ s(s1)+ s(s2), data = df, family = binomial)
#mybam <- bam(I(1-default_flag) ~ s(s1, k=10)+ s(s2, k=10), data = df, family = binomial)

mylogs <- glm(I(1-default_flag) ~ ns(s1, df=8) + ns(s2, df=8), data = df, family = binomial)

gini_logistic<-my_gini(mylogit$y, -mylogit$fitted.values)
#gini_tree<-my_gini(df$default_flag, predict(mytree, type = "prob")[,2])
#gini_gam <- my_gini(df$default_flag, predict(mygam, type = "response"))
#gini_bam <- my_gini(df$default_flag, predict(mybam, type = "response"))
gini_logs<-my_gini(mylogs$y, -mylogs$fitted.values)


theor_gini_combined<-gini_combine_1(gini1, gini2, corrs, dfrate)

