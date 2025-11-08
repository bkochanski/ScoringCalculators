# Proper

c = +.1
A = 0.8
ffind <- function(x){pnorm(x/sqrt(2)) + 2 * VGAM::pbinorm(-x/sqrt(2), 0, cov12 = -(1-c^2)/(1+c^2))-A}
#d_a = 1.2
d_a <- uniroot(ffind, lower = 0, upper = 10)$root
2*(AUCexp = pnorm(d_a/sqrt(2)) + 2 * VGAM::pbinorm(-d_a/sqrt(2), 0, cov12 = -(1-c^2)/(1+c^2)))-1

lower = -Inf
upper = Inf
if (c<0) {lower = d_a/4/c*(sqrt(1+c^2))}
if (c>0) {upper = d_a/4/c*(sqrt(1+c^2))}

H <- function(x){1*(x>=0)}
v_c = qlogis(seq(plogis(lower), plogis(upper), length.out=100))
options(scipen=20)
FPF <- pnorm(-(1-c)*v_c-d_a/2*sqrt(1+c^2)) + ( pnorm(-(1-c)*v_c+d_a/2/c*sqrt(1+c^2)) - H(c))
TPF <- pnorm(-(1+c)*v_c+d_a/2*sqrt(1+c^2)) + ( pnorm(-(1+c)*v_c+d_a/2/c*sqrt(1+c^2)) - H(c))
plot(FPF, TPF)

#pnorm(-(1-c)*v_c-d_a/2*sqrt(1+c^2)) + ( pnorm(-(1-c)*v_c+d_a/2/c*sqrt(1+c^2)) - H(c)) -  FPF
# FPF = pnorm(-(1-c)*v_c-d_a/2*sqrt(1+c^2)) + pnorm(-(1-c)*v_c+d_a/2/c*sqrt(1+c^2)) 

fpf <- function(v_c) {
  return(pnorm(-(1-c)*v_c-d_a/2*sqrt(1+c^2)) + ( pnorm(-(1-c)*v_c+d_a/2/c*sqrt(1+c^2)) - H(c)))
}

tpf <- function(v_c) {
  return(pnorm(-(1+c)*v_c+d_a/2*sqrt(1+c^2)) + ( pnorm(-(1+c)*v_c+d_a/2/c*sqrt(1+c^2)) - H(c)))
}

# fpf_ <- 0.5
# f <- function(x){fpf(x)-fpf_}

v_c_from_fpf <- function(fpf_){
  uniroot(function(x){fpf(x)-fpf_}, lower = -100, upper = 100)$root
}

roc <- function(x){tpf(v_c_from_fpf(x))}
Vroc <- Vectorize(roc)
curve(Vroc)

v_c_from_tpf <- function(tpf_){
  uniroot(function(x){tpf(x)-tpf_}, lower = -100, upper = 100)$root
}
revroc <- function(x){fpf(v_c_from_tpf(x))}
Vrevroc <- Vectorize(revroc)
curve(Vrevroc, add=TRUE, col="blue")

#points <- c(0.1, .2, .3, .4)
#numDeriv::grad(Vroc, x = points)

ProperBiNormalFunc <- function(x){
  Vroc(x)
}
ProperBiNormalDerivNum <- function(x){
  numDeriv::grad(Vroc, x=x)
}
ProperBiNormalDeriv <- Vectorize(function(x){
  s = v_c_from_fpf(x)
  res = ((1+c)*(dnorm(-(1+c)*s+d_a/2*sqrt(1+c^2)) + dnorm(-(1+c)*s+d_a/2/c*sqrt(1+c^2)) ))/
    ((1-c)*(dnorm(-(1-c)*s-d_a/2*sqrt(1+c^2)) + dnorm(-(1-c)*s+d_a/2/c*sqrt(1+c^2))))
  return(res)
})

# ProperBiNormalDerivNum(.1)
# ProperBiNormalDeriv(.1)

ProperBiNormalInverse <- function(x){
  Vrevroc(x) 
}


# Wykresy:
## ROC
library(ggplot2)
ggplot(data.frame(x = c(0, 1)), aes(x = x)) +
  stat_function(fun = ProperBiNormalFunc, color = "darkblue") +
  labs(x = "TPR",
       y = "FPR") +
  theme_minimal()

## Dwie gęstości
goodnsim <- 1e4
BadRate <- 0.1
BadOdds <- BadRate / (1 + BadRate)
GoodLogits <- -log(ProperBiNormalDeriv(randtoolbox::sobol(goodnsim))*BadOdds)
BadLogits <- -log(ProperBiNormalDeriv(ProperBiNormalInverse(randtoolbox::sobol(goodnsim*BadOdds)))*BadOdds)

break_probs <- c(.99, .9, .8, .65, .5, .2, .1, .05, .02, .01, .005, .002, .001, .0001)
break_points <- -qlogis(break_probs)
break_labels <- paste0(break_probs*100,"%")

target <- c(rep(0, length(GoodLogits)), rep(1, length(BadLogits)))
score <- c(GoodLogits, BadLogits)
bigstatsr::AUC(-score, target)

dff <- data.frame(score, target=as.factor(target))

p1 <- ggplot(dff, aes(x = score)) +
  #  geom_density(alpha = 0.5) +
  labs(title = "",
       x = "Predicted Bad rate (logit scale)",
       y = "") +
  scale_x_continuous(
    limits = c(min(break_points), max(break_points)),
    breaks = break_points,
    labels = break_labels
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("skyblue", "salmon"))

print(p1 + geom_density(aes(fill=target), alpha = 0.5))

## Dwie gęstości prop
print(p1 + geom_density(alpha = 0.5, position = "identity", aes(fill= target, y = ..count..)))

## Dwie gęstości stack
print(p1 + geom_density(alpha = 0.5, position = "stack", aes(fill= target, y = ..count..)))

## Jedna gęstość
print(p1 + geom_density(alpha = 0.5))

## Fit the 
