# define R-vine tree structure matrix
Matrix <- c(3, 2, 1,
            0, 2, 1,
            0, 0, 1)
Matrix <- matrix(Matrix, 3, 3)
# define R-vine pair-copula family matrix
family <- c(0, 1, 1,
            0, 0, 1,
            0, 0, 0)
family <- matrix(family, 3, 3)
# define R-vine pair-copula parameter matrix

par <- c(0, 0.4, 0.5,
         0, 0, 0.5,
         0, 0, 0)
par <- matrix(par, 3, 3)
# define second R-vine pair-copula parameter matrix
par2 <- c(0, 0, 0,
         0, 0, 0,
         0, 0, 0)
par2 <- matrix(par2, 3, 3)
# define RVineMatrix object
library(VineCopula)
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3"))

plot(RVM,1)
plot(RVM,2)
contour(RVM,1)
contour(RVM,2)

# simulate a sample of size 300 from the R-vine copula model
#set.seed(123)
simdata <- RVineSim(300, RVM)
psych::pairs.panels(simdata)
