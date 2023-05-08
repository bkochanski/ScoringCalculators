library(copula)
set.seed(100)
myCop <- normalCopula(param=c(0.4,0.2,-0.8), dim = 3, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("gamma", "beta", "t"),
              paramMargins=list(list(shape=2, scale=1),
                                list(shape1=2, shape2=2), 
                                list(df=5)) )
Z2 <- rMvdc(2000, myMvd)
colnames(Z2) <- c("x1", "x2", "x3")
psych::pairs.panels(Z2)




my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 
my_gini2<-function(resp, pred){1-2*bigstatsr::AUC(pred, resp)}

#my_gini2(S$bad, S$score)

my_gini_roc<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  plot(gc, bc)
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 



library(copula)
# set.seed(100)
myCop <- normalCopula(param=c(0.6), dim = 2, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("norm", "norm"),
              paramMargins=list(list(mean=0, sd=1),
                                list(mean=0, sd=1)))

myMvd2 <- mvdc(copula=myCop, margins=c("gamma", "norm"),
              paramMargins=list(list(shape=2, scale=1),
                                list(mean=0, sd=1)))


Z2 <- rMvdc(10000, myMvd)
G2 <- rMvdc(10000, myMvd2)
colnames(Z2) <- c("x1", "x2")
colnames(G2) <- c("x1", "x2")

psych::pairs.panels(Z2)
psych::pairs.panels(G2)

S<-data.frame(score=Z2[,1], bad=1*(Z2[,2]<qnorm(.1)))

SG2<-data.frame(score=G2[,1], bad=1*(G2[,2]<qnorm(.1)))

my_gini_roc(SG2$bad, SG2$score)




