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


number_of_sim<-300
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
actual_bad_rate<-c()
gini_model<-c()
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
    actual_bad_rate[i]<-NA
    gini_model[i]<-NA
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
  actual_bad_rate[i]<-mean(df$default_flag)
  if (any(c(actual_gini1[i]<=0, actual_gini2[i]<=0, actual_corr[i]<=0, actual_bad_rate[i]<=0))){
    gini_model[i]<-NA
    gini_logistic[i]<-NA
    gini_probit[i]<-NA
    next
  }
  gini_model[i]<-gini_combine_1(actual_gini1[i], actual_gini2[i], actual_corr[i], actual_bad_rate[i])
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
  , actual_bad_rate
  , gini_model
  , gini_logistic
  , gini_probit
  , diff.logistic
  , diff.probit)



#View(results)
write.csv(results, paste0('simresults', sample_size, format(Sys.time(), "%Y%m%d%H%M%S"), '.csv')
)
sd(diff.logistic, na.rm=TRUE)
sd(diff.probit, na.rm=TRUE)

#res40000<-read.csv("ScoringCalculators/simresults4000020220428233229.csv")
#res10000<-read.csv("ScoringCalculators/simresults1000020220428223014.csv")
#res1000<-read.csv("ScoringCalculators/simresults100020220428212247.csv")
#res100000<-read.csv("ScoringCalculators/simresults10000020220429000403.csv")
res40000<-read.csv("simresults4000020220428233229.csv")
res10000<-read.csv("simresults1000020220428223014.csv")
res1000<-read.csv("simresults100020220428212247.csv")
res100000<-read.csv("simresults10000020220429000403.csv")


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

