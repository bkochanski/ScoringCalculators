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
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf, subdivisions=500, rel.tol = .Machine$double.eps^.24)$value/defrate/(1-defrate)-1
}

r_from_gini<-function(gini, defaultrate=.1){
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-gini}
  uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
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

gini_combine_1<-function(g1, g2, corr, defaultrate){
  tryCatch(unname(gini_combine_calculator(g1, g2, corr, defaultrate)[1]), error=function(err) NA)}

gini_combine_1<-Vectorize(gini_combine_1)

gini1_here<-0.4
gini2_here<-0.6
corr_here<-0.3
badrate_here<-0.1

sigma <- matrix(c(1, corr_here, r_from_gini(gini1_here, dfrate),
                  corr_here, 1, r_from_gini(gini1_here, dfrate),
                  r_from_gini(gini1_here, badrate_here), r_from_gini(gini1_here, badrate_here), 1), 
                nrow=3)
z <- mvrnorm(100000,mu=rep(0, 3),Sigma=sigma)
df<-z
df[,3]<-(z[,3]<qnorm(badrate_here))*1
df<-data.frame(df)
names(df)<-c('s1', 's2', 'default_flag')
df$ls1<-log(df$s1-min(df$s1)+1)
#sum(is.na(df$ls1))
hist(df$ls1, breaks=100)
df$ls2<-log(df$s2-min(df$s2)+1)
#sum(is.na(df$ls1))
hist(df$ls2, breaks=100)

gini_combine_1(gini1_here, gini2_here, corr_here, badrate_here)

gini_combine_1(my_gini(df$default_flag, df$s1), my_gini(df$default_flag, df$s2), 
               cor(df$s1, df$s2), mean(df$default_flag))

gini_combine_1(my_gini(df$default_flag, df$ls1), my_gini(df$default_flag, df$ls2), 
               cor(df$ls1, df$ls2), mean(df$default_flag))

sum(is.na(df$ls1))
sum(is.na(df$ls2))
mean(df$ls1)
mean(df$ls2)
