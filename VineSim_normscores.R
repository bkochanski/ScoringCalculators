# define R-vine tree structure matrix
Matrix <- c(3, 2, 1,
            0, 2, 1,
            0, 0, 1)
Matrix <- matrix(Matrix, 3, 3)
# define R-vine pair-copula family matrix
# https://tnagler.github.io/VineCopula/reference/BiCopSelect.html
# ?RVineMatrix
family <- c(0, 1, 3,
            0, 0, 10,
            0, 0, 0)
family <- matrix(family, 3, 3)
# define R-vine pair-copula parameter matrix

par <- c(0, .5, 0.7,
         0, 0, 3.32,
         0, 0, 0)
par <- matrix(par, 3, 3)
# define second R-vine pair-copula parameter matrix
par2 <- c(0, 0, 0,
          0, 0, 0.88,
          0, 0, 0)
par2 <- matrix(par2, 3, 3)



# define RVineMatrix object
library(VineCopula)
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3"))
#RVM

# plot(RVM,1)
# plot(RVM,2)
# contour(RVM,1)
# contour(RVM,2)

# simulate a sample of size 300 from the R-vine copula model
#set.seed(123)
sizesim<-10000
simdata <- RVineSim(sizesim, RVM)
#psych::pairs.panels(simdata)

cor(simdata)
#cor(simdata, method="kendall")
copula::rho(copula::claytonCopula(.7))
copula::tau(VC2copula::BB8Copula(param=c(3.32, 0.88)))
copula::tau(copula::claytonCopula(.7))



# BiCopSelect(pobs(simdata[1:1000,1]), pobs(simdata[1:1000,2]), familyset = NA) #BB8 (3.32, 0.88)
# BiCopSelect(pobs(simdata[1:2000,1]), pobs(simdata[1:2000,3]), familyset = NA) #Clayton (0.7)
# BiCopSelect(pobs(simdata[1:1000,2]), pobs(simdata[1:1000,3]), familyset = NA) #Gaussian ?

simdf<-data.frame(simdata)
names(simdf)<-c("latentunif", "s1", "s2")
simdf$latent<-(qnorm(simdf$latentunif))
mean(simdf$latent)
sd(simdf$latent)
mean(simdf$latentunif<.1)
simdf$bad<-1*(simdf$latentunif<.1)

simdf$s1n<-qnorm(simdf$s1)
simdf$s2n<-qnorm(simdf$s2)

my_gini_roc<-function(resp, pred){
  #c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
#  plot(gc, bc)
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

my_gini<-function(resp, pred){
  abs(2*bigstatsr::AUC(pred, resp)-1)}


#my_gini_roc(simdf$bad, simdf$s1)
my_gini_roc(simdf$bad, simdf$s1n)
#my_gini_roc(simdf$bad, simdf$s2)
my_gini_roc(simdf$bad, simdf$s2n)
#cor(simdf$s1, simdf$s2)
cor(simdf$s1n, simdf$s2n)

# Two ROC curves...
d1<-simdf$bad[order(simdf$s1n)]
bc1<-c(0, cumsum(d1)/sum(d1))
gc1<-c(0, cumsum(1-d1)/sum(1-d1))
plot(gc1, bc1, type='l')
d2<-simdf$bad[order(simdf$s2n)]
bc2<-c(0, cumsum(d2)/sum(d2))
gc2<-c(0, cumsum(1-d2)/sum(1-d2))
lines(gc2, bc2, col='blue')



my_gini_roc(simdf$bad, simdf$s1n+simdf$s2n)
w1=0.3424
my_gini_roc(simdf$bad, w1*simdf$s1+(1-w1)*simdf$s2)

gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf)$value/defrate/(1-defrate)-1
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

gini_combine_calculator(0.5008577, 
                        0.5888478,
                        0.5560054,
                        0.1017)

gini_combine_calculator(my_gini_roc(simdf$bad, simdf$s1n), 
                        my_gini_roc(simdf$bad, simdf$s2n),
                        cor(simdf$s1n, simdf$s2n),
                        mean(simdf$bad))



gini_profile<-Vectorize(function(w){my_gini_roc(simdf$bad, w*simdf$s1+(1-w)*simdf$s2)})
gp<-gini_profile(20:40/100)
plot(20:40/100, gp)

