FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}

f<-function(x){FuncBiNormal(x, .6, .7)}

curve(f)

size<-100
prct<-1:(size-1)/size

B=.1

optfV<-Vectorize(function(p)
{
optf<-function(x){x*(1-B)+f(x)*B-p}
uniroot(optf, lower=0, upper=1)$root}
)

prctx<-optfV(prct)

library(numDeriv)
options(scipen=10)
oddsx<-grad(f, prctx)*B/(1-B)
badratesx<-oddsx/(1+oddsx)

# prctx[1]
# badratesx[1]

# #should be B
# mean(badratesx) 

# sum(prctx*(1-B)+f(prctx)*B-prct)

cutoff<-qnorm(B)
#(cutoff-qnorm(badratesx))
#pnorm(cutoff, mean=(cutoff-qnorm(badratesx)))

sims<-10000
scores<-sample(1:(size-1),sims, replace=TRUE)
#hist(scores)
latent<-rnorm(sims, mean=(cutoff-qnorm(badratesx[scores])))

plot(scores, latent)
psych::pairs.panels(data.frame(scores, latent))

my_gini_roc<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  plot(gc, bc)
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

my_gini_roc(latent<cutoff, scores)
curve(f, add=TRUE, col='red', lwd=4)

library(VineCopula)
selectedCopula <- BiCopSelect(pobs(scores), pobs(latent), familyset = NA)
selectedCopula

# Bivariate copula: BB8 (par = 3.46, par2 = 0.86, tau = 0.45) 
# Bivariate copula: BB8 (par = 3.74, par2 = 0.82, tau = 0.45) 
# Bivariate copula: BB8 (par = 3.66, par2 = 0.83, tau = 0.45) 

# mean(latent)
# sd(latent)
# 
# hist(scores)
# mean(scores)
# sd(scores)


library(copula)
library(VC2copula)
# set.seed(100)
myCop <- claytonCopula(par=0.73, dim=2)
myMvd <- mvdc(copula=myCop, margins=c("unif", "norm"),
              paramMargins=list(list(min=0, max=1),
                                list(mean=mean(latent), sd=sd(latent))))

Z2 <- rMvdc(10000, myMvd)
colnames(Z2) <- c("x1", "x2")

psych::pairs.panels(Z2)

S<-data.frame(score=Z2[,1], bad=1*(Z2[,2]<qnorm(B)))

my_gini_roc(S$bad, S$score)
curve(f, add=TRUE, col='red', lwd=4)

# cor(scores, latent)
# cor(Z2)
# mean(Z2[,1])
# mean(Z2[,2])
# sd(Z2[,1])
# sd(Z2[,2])



