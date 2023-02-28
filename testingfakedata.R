library(SimMultiCorrData)

# Testing ideas from https://rviews.rstudio.com/2020/09/09/fake-data-with-r/

# Binary, Ordinal, Continuous, Poisson, and Negative Binomial Variables
options(scipen = 999) # decides when to switch between fixed or exponential notation
seed <- 1234
n <- 10000 # number of random draws to simulate
Dist <- c("Logistic", "Weibull") # the 2 continuous rvs to be simulated 
#Params list: first element location and scale for logistic rv
#             second element shape and scale for Weibull rv
Params <- list(c(0, 1), c(3, 5))
# Calculate theoretical parameters for Logistic rv including 4th, 5th and 6th moments
Stcum1 <- calc_theory(Dist[1], Params[[1]]) 
# Calculate theoretical parameters for a Weibull rv 
Stcum2 <- calc_theory(Dist[2], Params[[2]])
Stcum <- rbind(Stcum1, Stcum2)
rownames(Stcum) <- Dist
colnames(Stcum) <- c("mean", "sd", "skew", "skurtosis", "fifth", "sixth")
Stcum # matrix of parameters for continuous random variables

##           mean    sd   skew skurtosis  fifth sixth
## Logistic 0.000 1.814 0.0000    1.2000  0.000 6.857
## Weibull  4.465 1.623 0.1681   -0.2705 -0.105 0.595

# Six is a list of vectors of correction values to add to the sixth cumulants
# if no valid pdf constants are found
Six <- list(seq(1.7, 1.8, 0.01), seq(0.10, 0.25, 0.01))
marginal <- list(0.3) # cum prob for Weibull rv, Logistic assumed = 1
lam <- 0.5 # Constant for Poisson rv
size <- 2 # Size parameter for Negative Binomial RV
prob <- 0.75 # prob of success
Rey <- matrix(0.4, 5, 5)  # target correlation matrix: the order is important
diag(Rey) <- 1
# Make sure Rey is within upper and lower correlation limits
# and that a valid power method exists
valid <- valid_corr(k_cat = 1, k_cont = 2, k_pois = 1, k_nb = 1,
                    method = "Polynomial", means = Stcum[, 1],
                    vars = Stcum[, 2]^2, skews = Stcum[, 3],
                    skurts = Stcum[, 4], fifths = Stcum[, 5],
                    sixths = Stcum[, 6], Six = Six, marginal = marginal,
                    lam = lam, size = size, prob = prob, rho = Rey,
                    seed = seed)

## 
##  Constants: Distribution  1  
## 
##  Constants: Distribution  2  
## All correlations are in feasible range!

# Simulate the data

Sim1 <- rcorrvar(n = n, k_cat = 1, k_cont = 2, k_pois = 1, k_nb = 1,
                 method = "Polynomial", means = Stcum[, 1],
                 vars = Stcum[, 2]^2, skews = Stcum[, 3],
                 skurts = Stcum[, 4], fifths = Stcum[, 5],
                 sixths = Stcum[, 6], Six = Six, marginal = marginal,
                 lam = lam, size = size, prob = prob, rho = Rey,
                 seed = seed)

## 
##  Constants: Distribution  1  
## 
##  Constants: Distribution  2  
## 
## Constants calculation time: 0.141 minutes 
## Intercorrelation calculation time: 0.008 minutes 
## Error loop calculation time: 0 minutes 
## Total Simulation time: 0.149 minutes

