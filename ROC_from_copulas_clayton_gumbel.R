copula::rho(copula::claytonCopula(.7))
#copula::tau(VC2copula::BB8Copula(param=c(3.32, 0.88)))
copula::rho(copula::gumbelCopula(1.4))

clayton2rho <- Vectorize(function(x){copula::rho(copula::claytonCopula(x))})
#plot(0:1000/20, clayton2rho(0:1000/20))
rho2clayton <- splinefun(clayton2rho(0:1000/20), 0:1000/20)
rho2clayton(.5)

gumbel2rho <- Vectorize(function(x){copula::rho(copula::gumbelCopula(x))})
#plot(20:1000/20, gumbel2rho(20:1000/20))
rho2gumbel <- splinefun(gumbel2rho(20:1000/20), 20:1000/20)
rho2gumbel(.5)

rho2r <- function(rho){2*sin(rho*pi/6)}
r2rho <- function(r){6/pi*asin(r/2)}



#### drawing

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

#
# for gumbel
# {sim <- copula::rCopula(1e8,copula::claytonCopula(param = x), dim = 2)
#   2*bigstatsr::AUC(-sim[,1], 1*(sim[,2]<.1))-1}


simulcgn <- function(){

gumbelseq <- rho2gumbel(1:90/100)
claytonseq <- rho2clayton(1:90/100)

gpar<-sample(gumbelseq, 1)
cpar<-sample(claytonseq, 1)
npar<- sample(1:80/100, 1)

Matrix <- c(3, 2, 1,
            0, 2, 1,
            0, 0, 1)
Matrix <- matrix(Matrix, 3, 3)
family <- c(0, 1, 3,
            0, 0, 4,
            0, 0, 0)
family <- matrix(family, 3, 3)
par <- c(0, npar, cpar,
         0, 0, gpar,
         0, 0, 0)
par <- matrix(par, 3, 3)
par2 <- c(0, 0, 0,
          0, 0, 0,
          0, 0, 0)
par2 <- matrix(par2, 3, 3)
library(VineCopula)
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3"))
sizesim<-10000
simdata <- RVineSim(sizesim, RVM)
cm<-cor(simdata)

return(data.frame(gpar, cpar, npar, gcor=cm[1,2], ccor=cm[1,3], ncor=cm[2,3]))
}

set.seed(123)  # For reproducibility
n_simulations <- 1000  # Number of simulations

results <- do.call(rbind, replicate(n_simulations, simulcgn(), simplify = FALSE))

# View the first few rows
head(results)

plot(rho2r(gumbel2rho(results$gpar)), results$gcor)
plot(rho2r(clayton2rho(results$cpar)), results$ccor)
plot(results$npar, results$ncor)

ncorf <- lm(ncor ~ gpar * cpar * npar, results)
summary(ncorf)

plot(ncorf$fitted.values, results$ncor)


plot(rho2r(gumbel2rho(results$gpar)), results$gcor)
plot(rho2r(clayton2rho(results$cpar)), results$ccor)
plot(results$npar, results$ncor)

ncorf <- lm(npar ~ gpar * cpar * ncor, results)
summary(ncorf)

plot(ncorf$fitted.values, results$npar)

predict(ncorf, newdata = data.frame(gpar = 2, cpar = 2, ncor = .5))

saveRDS(ncorf, "ncorf.RDS")

library(splines)
ncorf2 <- lm(npar ~ ns(gpar, df=4) * ns(gpar, df=4)  * ns(ncor, df=4), results)
summary(ncorf2)
plot(ncorf2$fitted.values, results$ncor)


