copula::rho(copula::claytonCopula(.7))
#copula::tau(VC2copula::BB8Copula(param=c(3.32, 0.88)))
copula::rho(copula::gumbelCopula(1.4))

clayton2rho <- Vectorize(function(x){copula::rho(copula::claytonCopula(x))})
plot(0:1000/20, clayton2rho(0:1000/20))
rho2clayton <- splinefun(clayton2rho(0:1000/20), 0:1000/20)
rho2clayton(.5)

gumbel2rho <- Vectorize(function(x){copula::rho(copula::gumbelCopula(x))})
plot(20:1000/20, gumbel2rho(20:1000/20))
rho2gumbel <- splinefun(gumbel2rho(20:1000/20), 20:1000/20)
rho2gumbel(.5)

rho2r <- function(rho){2*sin(rho*pi/6)}
r2rho <- function(r){6/pi*asin(r/2)}


copsimn<-100000

library(copula)
clayton <- claytonCopula(param = rho2clayton(r2rho(.6)), dim = 2)
sim <- rCopula(copsimn, clayton) 
sim2<-data.frame(sim)
#plot(sim2, pch='.')
cor(sim2)
names(sim2)<-c("score", "bad")
sim2$bad<-1*(sim2$bad<.2)

plot(cumsum(1-sim2$bad[order(sim2$score)])/sum(1-sim2$bad), 
     cumsum(sim2$bad[order(sim2$score)])/sum(sim2$bad), type="l")

library(copula)
gumbel <- gumbelCopula(param = rho2gumbel(r2rho(.6)), dim = 2)
simg <- rCopula(copsimn, gumbel) 
simg2<-data.frame(simg)
#plot(simg2, pch='.')
names(simg2)<-c("score", "bad")
cor(simg2)
simg2$bad<-1*(simg2$bad<.2)
# mean(simg2$bad)
# hist(simg[,1])

points(cumsum(1-simg2$bad[order(simg2$score)])/sum(1-simg2$bad), 
     cumsum(simg2$bad[order(simg2$score)])/sum(simg2$bad), type="l", col='blue')


library(copula)
nor <- normalCopula(param = 0.65, dim = 2)
simn <- rCopula(copsimn, nor) 
simn2<-data.frame(simn)
names(simn2)<-c("score", "bad")
simn2$bad<-1*(simn2$bad<.2)

points(cumsum(1-simn2$bad[order(simn2$score)])/sum(1-simn2$bad), 
       cumsum(simn2$bad[order(simn2$score)])/sum(simn2$bad), type="l", col='green')

cor(simn)


