# Required package
library(rootSolve)

# Define parameters
c <- 0       # correlation-like parameter
d_a <- Gini2Da(.7, c)     # separation parameter
x0 <- 0.36      # point at which to compute dy/dx

# Derived quantity
alpha <- sqrt(1 + c^2)
Heaviside <- function(x){1*(x>=0)}

# Define x(s)
x_s <- function(s) {
  term1 <- pnorm(-(1 + c) * s + 0.5 * d_a * alpha)
  term2 <- pnorm(-(1 + c) * s + 0.5 * d_a * alpha / c) - Heaviside(c)
  return(term1 + term2)
}

# Define y(s)
y_s <- function(s) {
  term1 <- pnorm(-(1 - c) * s - 0.5 * d_a * alpha)
  term2 <- pnorm(-(1 - c) * s + 0.5 * d_a * alpha / c) - Heaviside(c)
  return(term1 + term2)
}

# Define dy/dx via chain rule
dy_dx <- function(s) {
  z1 <- -(1 + c) * s + 0.5 * d_a * alpha
  z2 <- -(1 + c) * s + 0.5 * d_a * alpha / c
  z3 <- -(1 - c) * s - 0.5 * d_a * alpha
  z4 <- -(1 - c) * s + 0.5 * d_a * alpha / c
  
  dxds <- -(1 + c) * (dnorm(z1) + dnorm(z2))
  dyds <- -(1 - c) * (dnorm(z3) + dnorm(z4))
  
  return(dyds / dxds)
}

# Find s such that x(s) = x0
root_fun <- function(s) x_s(s) - x0
s_solution <- uniroot(root_fun, lower = -10, upper = 10)$root

# Compute dy/dx at x0
s_solution  # s corresponding to x0
slope <- dy_dx(s_solution)
slope       # dy/dx at x0


v_c_from_tpf(function(x_){tpf(x_,c,d_a)},x0)
uniroot(function(s) x_s(s) - x0, lower = -10, upper = 10)$root
uniroot(function(s) tpf(s,c,d_a) - x0, lower = -10, upper = 10)$root
tpf(.7328,c,d_a)
x_s(.7328)

x_s(0.73)
tpf(0.73,c,d_a)

ProperBiNormalDeriv1 <- Vectorize(function(x, c1, d_a1){
  s = v_c_from_tpf(function(x_){tpf(x_,c1,d_a1)},x)
  # res = 1/((1+c1)*(dnorm(-(1+c1)*s+d_a1/2*sqrt(1+c1^2)) + dnorm(-(1+c1)*s+d_a1/2/c1*sqrt(1+c1^2)) ))*
  #   ((1-c1)*(dnorm(-(1-c1)*s-d_a1/2*sqrt(1+c1^2)) + dnorm(-(1-c1)*s+d_a1/2/c1*sqrt(1+c1^2))))
  res <- ((1-c1)*(dnorm(-(1-c1)*s-d_a1*sqrt(1+c1^2)/2)+dnorm(-(1-c1)*s+d_a1*sqrt(1+c1^2)/2/c1)))/((1+c1)*(dnorm(-(1+c1)*s+d_a1*sqrt(1+c1^2)/2)+dnorm(-(1+c1)*s+d_a1*sqrt(1+c1^2)/2/c1)))
  return(res)
})

x0 <- .5
numDeriv::grad(function(x){FuncBiNormal(x, .7, 1)}, .36)
#numDeriv::grad(function(x){FuncBiProper(x, 0, Gini2Da(.7, 0))}, x0)
FuncBiProperDeriv(.36, 0, Gini2Da(.7,0))
ProperBiNormalDeriv1(.36, 0, Gini2Da(.7,0))

FuncBiNormal(.36, .7, 1)
FuncBiProper(.36, 0, Gini2Da(.7,0))

FuncBiNormal(.37, .7, 1)
FuncBiProper(.37, 0, Gini2Da(.7,0))

mala <- 0.0001
(FuncBiProper(.36+mala, 0, Gini2Da(.7,0)) - FuncBiProper(.36, 0, Gini2Da(.7,0)))/mala
