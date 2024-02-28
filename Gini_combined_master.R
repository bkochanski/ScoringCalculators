gini_from_r <- function(rho = .5, defrate = .1) {
  F1 <- function(s, d, rho) {
      integrate(function(x) {
        pnorm((qnorm(d) - rho * x) / sqrt(1 - rho ^ 2)) * dnorm(x)
      }, lower = -Inf, upper = s)$value
    }
  F2 <- function(s, d, rho) {
      pnorm((-qnorm(d) + rho * s) / sqrt(1 - rho ^ 2)) * dnorm(s)
    }
  F1_ <- Vectorize(function(x) {
    F1(x, defrate, rho)
  })
  F2_ <- Vectorize(function(x) {
    F2(x, defrate, rho)
  })
  2 * integrate(
    function(x) {F1_(x) * F2_(x)},
    lower = -Inf, upper = Inf,
    subdivisions = 500,
    rel.tol = .Machine$double.eps ^ .24
  )$value / defrate / (1 - defrate) - 1
}

r_from_gini <- function(gini, defaultrate = .1) {
  phi_s1 <- function(x) {gini_from_r(rho = x, defrate = defaultrate) - gini}
  uniroot(phi_s1,
          lower = 0,
          upper = 1,
          tol = .Machine$double.eps)$root
}

gini_combine_calculator <- function(g1, g2, corr, defaultrate) {
  phi_s1 <- function(x) {gini_from_r(rho = x, defrate = defaultrate) - g1}
  rho_s1 <- uniroot(phi_s1,
            lower = 0,
            upper = 1,
            tol = .Machine$double.eps)$root
  phi_s2 <- function(x) { gini_from_r(rho = x, defrate = defaultrate) - g2}
  rho_s2 <-
    uniroot(phi_s2,
            lower = 0,
            upper = 1,
            tol = .Machine$double.eps)$root
  
  (a_opt <- (corr * rho_s2 - rho_s1) / (corr * rho_s1 - rho_s2))
  corr_opt <- (a_opt * rho_s1 + rho_s2) / 
    sqrt(a_opt ^ 2 + 2 * a_opt * corr +  1)
  g_result0 <- if (abs(corr_opt) > 1) {
    NaN
  } else {
    gini_from_r(corr_opt, defaultrate)
  }
  
  g_result1 <- if (a_opt < 0 | a_opt > 1000) {
    NaN
  } else {
    g_result0
  }
  return(
    c(
      new_gini = g_result1,
      a_opt = a_opt,
      score_1_weight = a_opt / (1 + a_opt),
      score_2_weight = 1 / (1 + a_opt),
      rho1 = rho_s1,
      rho2 = rho_s2,
      new_corr = corr_opt,
      gini1 = g1,
      gini2 = g2,
      corr = corr,
      defrate = defaultrate
    )
  )
}

#Example:
gini_combine_calculator(.6, .6, .5, .1)
