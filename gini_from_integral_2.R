gini_combine_calculator<-function(g1, g2, corr, defaultrate){
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
                lower=-Inf, upper=Inf, subdivisions = 200)$value/defrate/(1-defrate)-1
  }
  
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

gini_combine_calculator(0.5008577, 
                        0.5888478,
                        0.5560054,
                        0.1017)

g1= 0.5008577
g2 = 0.5888478
corr = 0.5560054
defaultrate = 0.1017



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
              lower=-Inf, upper=Inf, subdivisions=200)$value/defrate/(1-defrate)-1
}

#rho_s1
phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g1}
rho_s1<-uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
#rho_s2
phi_s2<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g2}
rho_s2<-uniroot(phi_s2,lower=0,upper=1,tol = .Machine$double.eps)$root
