library(copula)

# clayton gumbel check
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

sample_size=1e3

cpar <- rho2clayton(r2rho(0.6469349))
gpar <- rho2gumbel(r2rho(0.5301213))

ncorf <- readRDS("ncorf.RDS")
npar <- predict(ncorf, newdata = data.frame(gpar, cpar, ncor = 0.5132728 ))
npar <- min(max(npar, 0.01),0.99)
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
simdata <- RVineSim(sample_size, RVM)
z <- simdata[,c(3, 1:2)]

cor(z)

df<-z
df[,3]<-(z[,3]<0.1646967)*1
df<-data.frame(df)
names(df)<-c('s1', 's2', 'default_flag')

  df$s1 <- qnorm(df$s1)
  df$s2 <- qnorm(df$s2)


plot(cumsum(1-df$default_flag[order(df$s1)])/sum(1-df$default_flag), 
         cumsum(df$default_flag[order(df$s1)])/sum(df$default_flag), type="l", col='blue')
points(cumsum(1-df$default_flag[order(df$s2)])/sum(1-df$default_flag), 
     cumsum(df$default_flag[order(df$s2)])/sum(df$default_flag), type="l", col='green', add=TRUE)

2*bigstatsr::AUC(-df$s1, df$default_flag)-1
2*bigstatsr::AUC(-df$s2, df$default_flag)-1
