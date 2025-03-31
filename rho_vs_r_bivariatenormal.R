library(pracma)
fun <- function(x,y) sqrt(2.5^2 * (1 - x^2/5^2 - y^2/8.5^2))
xmin <- 0; xmax <- 5
ymin <- 0; ymax <- function(x) sqrt(8.5^2 * (1-x^2/5^2))
integral2(fun, xmin, xmax, ymin, ymax)

library(mvtnorm)
mu = c(0, 0)
Sigma = matrix (c(1, 0.5, 0.5, 1), 2, 2) # var-covariance matrix
# Generate a sample from N(mu, Sigma)

fxy <- function(x, y) {
  x1<-x
  y1<-y
  return(dmvnorm(c(x1, y1), mean = mu, sigma = Sigma))}

fxy<-Vectorize(fxy)
fxy(0,0)

integral2(fxy, -10, 10, -10, 10)

fxy2 <- Vectorize(function(x,y){x*y*fxy(x,y)})
integral2(fxy2, -20, 20, -20, 20)

fxy3 <- Vectorize(function(x,y){pnorm(x)*pnorm(y)*fxy(x,y)})
12*(integral2(fxy3, -10, 10, -10, 10)$Q-0.5*0.5)

#https://www.scirp.org/pdf/ojs_2016111714440618.pdf
r2rho <- function(r){6/pi*asin(r/2)}
curve(r2rho)
abline(0,1, col="blue", lty=3)

#https://stats.stackexchange.com/questions/262061/relation-between-spearman-and-pearson-correlation-in-gaussian-copula
rho2r <- function(rho){2*sin(rho*pi/6)}
rho2r(0.482583740)

#https://pmc.ncbi.nlm.nih.gov/articles/PMC3558980/