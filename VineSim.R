# define R-vine tree structure matrix
Matrix <- c(3, 2, 1,
            0, 2, 1,
            0, 0, 1)
Matrix <- matrix(Matrix, 3, 3)
# define R-vine pair-copula family matrix
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

# plot(RVM,1)
# plot(RVM,2)
# contour(RVM,1)
# contour(RVM,2)

# simulate a sample of size 300 from the R-vine copula model
#set.seed(123)
simdata <- RVineSim(10000, RVM)
#psych::pairs.panels(simdata)

cor(simdata)

# BiCopSelect(pobs(simdata[,1]), pobs(simdata[,2]), familyset = NA)
# BiCopSelect(pobs(simdata[,1]), pobs(simdata[,3]), familyset = NA)
# BiCopSelect(pobs(simdata[,2]), pobs(simdata[,3]), familyset = NA)

simdf<-data.frame(simdata)
names(simdf)<-c("latentunif", "s1", "s2")
simdf$latent<-(qnorm(simdf$latentunif))
mean(simdf$latent)
sd(simdf$latent)
mean(simdf$latentunif<.1)
simdf$bad<-1*(simdf$latentunif<.1)

my_gini_roc<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  plot(gc, bc)
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

my_gini<-function(resp, pred){
  abs(2*bigstatsr::AUC(pred, resp)-1)}


my_gini_roc(simdf$bad, simdf$s1)
#my_gini(simdf$bad, simdf$s1)
my_gini_roc(simdf$bad, simdf$s2)
cor(simdf$s1, simdf$s2)
my_gini_roc(simdf$bad, simdf$s1+simdf$s2)
w1=0.335797663146443
my_gini_roc(simdf$bad, w1*simdf$s1+(1-w1)*simdf$s2)

gini_combine_calculator(0.5008577, 
                        0.5888478,
                        0.5560054,
                        0.1017)

gini_combine_calculator(my_gini_roc(simdf$bad, simdf$s1), 
                        my_gini_roc(simdf$bad, simdf$s2),
                        cor(simdf$s1, simdf$s2),
                        mean(simdf$bad))



gini_profile<-Vectorize(function(w){my_gini_roc(simdf$bad, w*simdf$s1+(1-w)*simdf$s2)})
gp<-gini_profile(20:40/100)


