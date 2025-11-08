# Required packages
# install.packages("rootSolve") # if not already installed

library(rootSolve)

# Define parameters
c <- 0       # correlation-like parameter
d_a <- Gini2Da(.7, c)     # separation parameter
x0 <- 0.5      # point at which to compute dy/dx

# Derived quantity
alpha <- sqrt(1 + c^2)
Heaviside <- function(x){1*(x>=0)}

# Define x(s)
x_s <- function(s) {
  term1 <- pnorm(-(1 - c) * s - 0.5 * d_a * alpha)
  term2 <- pnorm(-(1 - c) * s + 0.5 * d_a * alpha / c) - Heaviside(c)
  return(term1 + term2)
}

# Define y(s)
y_s <- function(s) {
  term1 <- pnorm(-(1 + c) * s + 0.5 * d_a * alpha)
  term2 <- pnorm(-(1 + c) * s + 0.5 * d_a * alpha / c) - Heaviside(c)
  return(term1 + term2)
}

# Define derivative dy/dx
dy_dx <- function(s) {
  z1 <- -(1 - c) * s - 0.5 * d_a * alpha
  z2 <- -(1 - c) * s + 0.5 * d_a * alpha / c
  z3 <- -(1 + c) * s + 0.5 * d_a * alpha
  z4 <- -(1 + c) * s + 0.5 * d_a * alpha / c
  
  dxds <- -(1 - c) * (dnorm(z1) + dnorm(z2))
  dyds <- -(1 + c) * (dnorm(z3) + dnorm(z4))
  
  return(dyds / dxds)
}

# Find s such that x(s) = x0
root_fun <- function(s) x_s(s) - x0

s_solution <- uniroot(root_fun, lower = -10, upper = 10)$root

# Compute dy/dx at x0
s_solution

v_c_from_tpf(function(x_){tpf(x_,c,d_a)},x0)
uniroot(function(s) x_s(s) - x0, lower = -10, upper = 10)$root
uniroot(function(s) tpf(s,c,d_a) - x0, lower = -10, upper = 10)$root
tpf(.7328,c,d_a)
x_s(.7328)

x_s(0.73)
tpf(0.73,c,d_a)


slope <- dy_dx(s_solution)
slope


x_s_check <- function(s) {
  #term1 <- pnorm(-(1 - c) * s - 0.5 * d_a * alpha)
  #term2 <- pnorm(-(1 - c) * s + 0.5 * d_a * alpha / c) - Heaviside(c)
  return(pnorm(-(1 - c) * s - 0.5 * d_a * alpha) + pnorm(-(1 - c) * s + 0.5 * d_a * alpha / c) - Heaviside(c))
}

tpf_check <- function(s) {
  return(pnorm(-(1+c)*s+d_a/2*sqrt(1+c^2)) + ( pnorm(-(1+c)*s+d_a/2/c*sqrt(1+c^2)) - Heaviside(c)))
}

x_s_check(.73)
tpf_check(.73)
