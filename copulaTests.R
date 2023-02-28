library(copula)
set.seed(100)
myCop <- normalCopula(param=c(0.4,0.2,-0.8), dim = 3, dispstr = "un")
myMvd <- mvdc(copula=myCop, margins=c("gamma", "beta", "t"),
              paramMargins=list(list(shape=2, scale=1),
                                list(shape1=2, shape2=2), 
                                list(df=5)) )
Z2 <- rMvdc(2000, myMvd)
colnames(Z2) <- c("x1", "x2", "x3")
psych::pairs.panels(Z2, method="spearman")
psych::pairs.panels(Z2, method="pearson")

#?normalCopula



(Xtras <- copula:::doExtras())
n <- if(Xtras) 200 else 64

## A Gumbel copula
set.seed(7) # for reproducibility
gumbel.cop <- gumbelCopula(2, dim=3)
x <- rCopula(n, gumbel.cop) # "true" observations (simulated)
u <- pobs(x)                # pseudo-observations
## Inverting Kendall's tau
fit.tau <- fitCopula(gumbelCopula(), u, method="itau")
fit.tau
confint(fit.tau) # work fine !
confint(fit.tau, level = 0.98)
summary(fit.tau) # a bit more, notably "Std. Error"s
coef(fit.tau)# named vector
coef(fit.tau, SE = TRUE)# matrix

## Inverting Spearman's rho
fit.rho <- fitCopula(gumbelCopula(), u, method="irho")
summary(fit.rho)
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(gumbelCopula(), u, method="mpl")
fit.mpl
## Maximum likelihood -- use 'x', not 'u' ! --
fit.ml <- fitCopula(gumbelCopula(), x, method="ml")
summary(fit.ml) # now prints a bit more than simple 'fit.ml'
## ... and what's the log likelihood (in two different ways):
(ll. <- logLik(fit.ml))
stopifnot(all.equal(as.numeric(ll.),
                    loglikCopula(coef(fit.ml), u=x, copula=gumbel.cop)))

## A Gauss/normal copula

## With multiple/*un*constrained parameters
set.seed(6) # for reproducibility
normal.cop <- normalCopula(c(0.6, 0.36, 0.6), dim=3, dispstr="un")
x <- rCopula(n, normal.cop) # "true" observations (simulated)
u <- pobs(x)                # pseudo-observations
## Inverting Kendall's tau
fit.tau <- fitCopula(normalCopula(dim=3, dispstr="un"), u, method="itau")
fit.tau
## Inverting Spearman's rho
fit.rho <- fitCopula(normalCopula(dim=3, dispstr="un"), u, method="irho")
fit.rho
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(normalCopula(dim=3, dispstr="un"), u, method="mpl")
summary(fit.mpl)
coef(fit.mpl) # named vector
coef(fit.mpl, SE = TRUE) # the matrix, with SE
## Maximum likelihood (use 'x', not 'u' !)
fit.ml <- fitCopula(normalCopula(dim=3, dispstr="un"), x, method="ml", traceOpt=TRUE)
summary(fit.ml)
confint(fit.ml)
confint(fit.ml, level = 0.999) # clearly non-0

## Fix some of the parameters
param <- c(.6, .3, NA_real_)
attr(param, "fixed") <- c(TRUE, FALSE, FALSE)
ncp <- normalCopula(param = param, dim = 3, dispstr = "un")
fixedParam(ncp) <- c(TRUE, TRUE, FALSE)
## 'traceOpt = 5': showing every 5-th log likelihood evaluation:
summary(Fxf.mpl <- fitCopula(ncp, u, method = "mpl", traceOpt = 5))
Fxf.mpl@copula # reminding of the fixed param. values

## With dispstr = "toep" :
normal.cop.toep <- normalCopula(c(0, 0), dim=3, dispstr="toep")
## Inverting Kendall's tau
fit.tau <- fitCopula(normalCopula(dim=3, dispstr="toep"), u, method="itau")
fit.tau
## Inverting Spearman's rho
fit.rho <- fitCopula(normalCopula(dim=3, dispstr="toep"), u, method="irho")
summary(fit.rho)
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(normalCopula(dim=3, dispstr="toep"), u, method="mpl")
fit.mpl
## Maximum likelihood (use 'x', not 'u' !)
fit.ml <- fitCopula(normalCopula(dim=3, dispstr="toep"), x, method="ml")
summary(fit.ml)

## With dispstr = "ar1"
normal.cop.ar1 <- normalCopula(c(0), dim=3, dispstr="ar1")
## Inverting Kendall's tau
summary(fit.tau <- fitCopula(normalCopula(dim=3, dispstr="ar1"), u, method="itau"))
## Inverting Spearman's rho
summary(fit.rho <- fitCopula(normalCopula(dim=3, dispstr="ar1"), u, method="irho"))
## Maximum pseudo-likelihood
summary(fit.mpl <- fitCopula(normalCopula(dim=3, dispstr="ar1"), u, method="mpl"))
## Maximum likelihood (use 'x', not 'u' !)
fit.ml <- fitCopula(normalCopula(dim=3, dispstr="ar1"), x, method="ml")
summary(fit.ml)

## A t copula with variable df (df.fixed=FALSE)
(tCop <- tCopula(c(0.2,0.4,0.6), dim=3, dispstr="un", df=5))
set.seed(101)
x <- rCopula(n, tCop) # "true" observations (simulated)
## Maximum likelihood (start = (rho[1:3], df))
summary(tc.ml <- fitCopula(tCopula(dim=3, dispstr="un"), x, method="ml",
                           start = c(0,0,0, 10)))
## Maximum pseudo-likelihood (the asymptotic variance cannot be estimated)
u <- pobs(x)          # pseudo-observations
tc.mpl <- fitCopula(tCopula(dim=3, dispstr="un"),
                    u, method="mpl", estimate.variance=FALSE,
                    start= c(0,0,0, 10))
summary(tc.mpl)

library(VineCopula)
gumbel.cop <- gumbelCopula(2, dim=2)
x <- rCopula(2000, gumbel.cop) # "true" observations (simulated)
u <- pobs(x)                # pseudo-observations
u1 <- u[,1]
v1 <- u[,2]
selectedCopula <- BiCopSelect(u1,v1,familyset=NA)
summary(selectedCopula)

cor(x, method="spearman")
cor(x, method="pearson")
cor(x, method="kendall")


