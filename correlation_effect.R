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

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

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

rescor<-rep(0,20)
length(rescor)
for (i in 0:(length(rescor)-1)){
rescor[i+1]<-gini_combine_calculator(.6, .6, i/length(rescor), .1)[1]
}
plot(rescor/.6-1, ylim=c(0,.4))
rescor/.6-1

(0:(length(rescor)-1))/length(rescor)

gini_combine_calculator(.6, .6, .96, .1)[1]
gini_combine_calculator(.6, .6, .97, .1)[1]
gini_combine_calculator(.6, .6, .98, .1)[1]
gini_combine_calculator(.6, .6, .99, .1)[1]

gini_combine_calculator(.6, .6, .8, .1)[1]/.6

