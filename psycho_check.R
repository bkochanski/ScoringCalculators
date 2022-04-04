# from point-biserial to rho
# -.35 i -.36 & bad rate = 23.1
#https://www.jstor.org/stable/pdf/2333725.pdf?casa_token=TzSP-FbzfPQAAAAA:BTauIBYixcwCg33pL_oWkdVODKC2hGlJA4dQcos8u5UjZUnXliIIDl5NSqPxApj44DndIyc9c_qqmqQt7OiFsbgMv7NesdJNU1hDYoiAeznJUQ-_Atk


pbs_from_rho<-function(rho=0.5, defrate=0.1) {
  return(-dnorm(qnorm(defrate))*rho/sqrt(defrate-defrate^2))
}

rho_from_pbs<-function(pbs=0.5, defrate=.1){
  return(-sqrt(defrate-defrate^2)/dnorm(qnorm(defrate))*pbs)
}

pbs_from_rho(.5,.2)
rho_from_pbs(pbs_from_rho(.5,.2), .2)

# Sample 1
# drate<-.231
# rho_1<-rho_from_pbs(-.35,drate)
# rho_2<-rho_from_pbs(-.36,drate)
# rho_s<-.37

drate<-.0415
rho_1<-rho_from_pbs(-.13,drate)
rho_2<-rho_from_pbs(-.12,drate)
rho_s<-.09

# defrate=0.1
# rho=0.32
# sigma <- matrix(c(1, rho,
#                   rho, 1), 
#                 nrow=2)
# z <- mvrnorm(1000000,mu=rep(0, 2),Sigma=sigma)
# cor(z[,1],z[,2]<qnorm(defrate))
# cor(z[,1],z[,2])
# 
# dnorm(qnorm(defrate))*rho/sqrt(defrate-defrate^2)

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

gini1<-gini_crd(rho=rho_1, drate)
gini2<-gini_crd(rho=rho_2, drate)


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

(results<-gini_combine_calculator(gini1, gini2, rho_s, drate))

gini1
gini2
results[1]


#new_score= results[3]*score_1 + results[4]*score_2
#new_score= results[3]*(score_1o-65.39)/12.12 + results[4]*(score_2o-570)/73.69
results[3]/12.12
results[4]/73.69


gini1/2+.5
gini2/2+.5
results[1]/2+.5

#sample3
gini_crd(rho_from_pbs(-.35,.291), .291)


#Sample 1
drate<-.231
rho_1<-rho_from_pbs(-.35,drate)
rho_2<-rho_from_pbs(-.36,drate)
rho_s<-.37
(gini1<-gini_crd(rho=rho_1, drate))
(gini2<-gini_crd(rho=rho_2, drate))
gini1/2+.5
gini2/2+.5
(results<-gini_combine_calculator(gini1, gini2, rho_s, drate))
results[1]/2+.5

(w1<-results[3]/12.8)
(w2<-results[4]/110.2)
(w1/w2)

(-.061/-.005)

#Sample 2
drate<-.0415
rho_1<-rho_from_pbs(-.13,drate)
rho_2<-rho_from_pbs(-.12,drate)
rho_s<-.09
(gini1<-gini_crd(rho=rho_1, drate))
(gini2<-gini_crd(rho=rho_2, drate))
gini1/2+.5
gini2/2+.5
(results<-gini_combine_calculator(gini1, gini2, rho_s, drate))
results[1]/2+.5

(w1<-results[3]/12.12)
(w2<-results[4]/73.69)
(w1/w2)

(-.05/-.01)

#sample3
gini_crd(rho_from_pbs(-.35,.291), .291)
gini_crd(rho_from_pbs(-.35,.291), .291)/2+.5
