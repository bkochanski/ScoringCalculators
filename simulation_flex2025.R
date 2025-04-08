library(MASS)
library(splines)
library(copula)

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

ginic<-function(bc, gc){
  #function for gini when we have cumulative goods and bads vectors
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*
        (bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1
} 

gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  tryCatch(2*integrate(function(x){F1_(x)*F2_(x)}, 
                       lower=-Inf, upper=Inf, subdivisions=500, rel.tol = .Machine$double.eps^.24)$value/defrate/(1-defrate)-1, error=function(err) NA)
}

r_from_gini<-function(gini, defaultrate=.1){
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-gini}
  uniroot(phi_s1,lower=0,upper=.9999,tol = .Machine$double.eps)$root
}

gini_combine_calculator<-function(g1, g2, corr, defaultrate){
  #rho_s1
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g1}
  rho_s1<-uniroot(phi_s1,lower=0,upper=.9999,tol = .Machine$double.eps)$root
  #rho_s2
  phi_s2<-function(x){gini_from_r(rho=x, defrate=defaultrate)-g2}
  rho_s2<-uniroot(phi_s2,lower=0,upper=.9999,tol = .Machine$double.eps)$root
  
  (a_opt<-(corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2))
  corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
  g_result0<-if(corr_opt>1) {NaN} else {gini_from_r(corr_opt, defaultrate)}
  g_result1<-if(a_opt<0 | a_opt>1000 | corr_opt>1) {NaN} else {g_result0}
  return(c(new_gini=g_result1, 
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt, gini1=g1, gini2=g2, corr=corr, defrate=defaultrate))
}

gini_combine_1<-function(g1, g2, corr, defaultrate){
  tryCatch(unname(gini_combine_calculator(g1, g2, corr, defaultrate)[1]), error=function(err) NA)}

gini_combine_1<-Vectorize(gini_combine_1)

ncorf <- readRDS("ncorf.RDS")

clayton2rho <- Vectorize(function(x){copula::rho(copula::claytonCopula(x))})
#plot(0:1000/20, clayton2rho(0:1000/20))
rho2clayton <- splinefun(clayton2rho(0:1000/20), 0:1000/20)

gumbel2rho <- Vectorize(function(x){copula::rho(copula::gumbelCopula(x))})
#plot(20:1000/20, gumbel2rho(20:1000/20))
rho2gumbel <- splinefun(gumbel2rho(20:1000/20), 20:1000/20)

rho2r <- function(rho){2*sin(rho*pi/6)}
r2rho <- function(r){6/pi*asin(r/2)}


######## simul function #################

simul <- function(inp_seed = 125, inp_nsim = 999, inp_sample_size= 1e4, 
                  inp_mycop = "Normal", inp_mymarginals = "Normal"){
myseed <- inp_seed
set.seed(myseed)

number_of_sim <- inp_nsim

sample_size <- inp_sample_size #c(1e3, 1e4, 1e5)

mycop <- inp_mycop #"T5"; "T4"; "CGN"; "CBN"; "NNC"

mymarginals <- inp_mymarginals #"Normal" "Uniform"; "Skewed6"

input_corr<-runif(number_of_sim, 0, .99)
input_gini1<-runif(number_of_sim, 0.1, .9)
input_gini2<-runif(number_of_sim, 0.1, .9)
input_badrate<-runif(number_of_sim, 0.01, .5)
theor_gini_combined<-gini_combine_1(input_gini1, input_gini2, input_corr, input_badrate)

actual_gini1<-c()
actual_gini2<-c()
actual_corr<-c()
actual_rho<-c()
actual_rho2<-c()
actual_bad_rate<-c()
gini_model<-c()
gini_model_rho<-c()
gini_model_rho2<-c()
gini_logistic<-c()
gini_logistic2<-c()
gini_logs<-c()
gini_probit<-c()

for (i in 1:number_of_sim){
  print(i)
  gini1<-input_gini1[i]
  gini2<-input_gini2[i]
  corrs<-input_corr[i]
  dfrate<-input_badrate[i]
  if (is.na(theor_gini_combined[i])){
    actual_gini1[i]<-NA
    actual_gini2[i]<-NA
    actual_corr[i]<-NA
    actual_rho[i]<-NA
    actual_rho2[i]<-NA
    actual_bad_rate[i]<-NA
    gini_model[i]<-NA
    gini_model_rho[i]<-NA
    gini_model_rho2[i]<-NA    
    gini_logistic[i]<-NA
    gini_logistic2[i]<-NA
    gini_logs[i]<-NA
    gini_probit[i]<-NA
    next
  }
  sigma <- matrix(c(1, corrs, r_from_gini(gini1, dfrate),
                    corrs, 1, r_from_gini(gini2, dfrate),
                    r_from_gini(gini1, dfrate), r_from_gini(gini2, dfrate), 1), 
                  nrow=3)
  
  if(any(eigen(sigma)$values<=0)){
    actual_gini1[i]<-NA
    actual_gini2[i]<-NA
    actual_corr[i]<-NA
    actual_bad_rate[i]<-NA
    gini_model[i]<-NA
    gini_logistic[i]<-NA
    gini_probit[i]<-NA
    next    
  }
  
  if (mycop == "Normal") {
    norm_cop <- normalCopula(param = sigma[lower.tri(sigma)], dim = 3, dispstr = "un")
    z <- rCopula(sample_size, norm_cop)
  }
  
  if (mycop == "T5") {
    t_cop <- tCopula(param = sigma[lower.tri(sigma)], dim = 3, dispstr = "un", df=5)
    z <- rCopula(sample_size, t_cop)
    }

  if (mycop == "T4") {
    t_cop <- tCopula(param = sigma[lower.tri(sigma)], dim = 3, dispstr = "un", df=4)
    z <- rCopula(sample_size, t_cop)
  }
  
  if (mycop == "CGN") {
    cpar <- rho2clayton(r2rho(sigma[1,3]))
    gpar <- rho2gumbel(r2rho(sigma[2,3]))
    npar <- predict(ncorf, newdata = data.frame(gpar, cpar, ncor = corrs))
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
  }
  
  if (mycop == "NNC") {
    #print(sigma[1,3])
    npar1 <- sigma[1,3]
    npar2 <- sigma[2,3]
    cpar <- rho2clayton(r2rho(sigma[1,3]))
    Matrix <- c(3, 2, 1,
                0, 2, 1,
                0, 0, 1)
    Matrix <- matrix(Matrix, 3, 3)
    family <- c(0, 3, 1,
                0, 0, 1,
                0, 0, 0)
    family <- matrix(family, 3, 3)
    par <- c(0, cpar, npar1,
             0, 0, npar2,
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
  }
  
    
  df<-z
  df[,3]<-(z[,3]<dfrate)*1
  df<-data.frame(df)
  names(df)<-c('s1', 's2', 'default_flag')

  if (mymarginals == "Normal") {
    df$s1 <- qnorm(df$s1)
    df$s2 <- qnorm(df$s2)
  }
  
  if (mymarginals == "Skewed6") {
    df$s1 <- df$s1^6
    df$s2 <- df$s2^(1/6)
  }
  
  
  actual_gini1[i]<-my_gini(df$default_flag, df$s1)
  actual_gini2[i]<-my_gini(df$default_flag, df$s2)
  actual_corr[i]<-cor(df$s1, df$s2)
  actual_rho[i]<-cor(df$s1, df$s2, method="spearman")
  actual_rho2[i]<-2*sin(actual_rho[i]*pi/6)
  actual_bad_rate[i]<-mean(df$default_flag)
  if (any(c(actual_gini1[i]<=0, actual_gini2[i]<=0, actual_corr[i]<=0, actual_rho[i]<=0, actual_bad_rate[i]<=0))){
    gini_model[i]<-NA
    gini_logistic[i]<-NA
    gini_probit[i]<-NA
    next
  }
  gini_model[i]<-gini_combine_1(actual_gini1[i], actual_gini2[i], actual_corr[i], actual_bad_rate[i])
  gini_model_rho[i]<-gini_combine_1(actual_gini1[i], actual_gini2[i], actual_rho[i], actual_bad_rate[i])
  gini_model_rho2[i]<-gini_combine_1(actual_gini1[i], actual_gini2[i], actual_rho2[i], actual_bad_rate[i])
  mylogit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = "binomial")
  myprobit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = binomial(link="probit"))
  mylogs <- glm(I(1-default_flag) ~ ns(qnorm(pobs(s1)), df=2) + ns(qnorm(pobs(s2)), df=2), data = df, family = binomial)
  gini_logistic[i]<-my_gini(mylogit$y, -mylogit$fitted.values)
  gini_probit[i]<-my_gini(myprobit$y, -myprobit$fitted.values)
  gini_logs[i]<-my_gini(mylogs$y, -mylogs$fitted.values)
  
}

diff.probit<-gini_probit-gini_model
diff.logistic<-gini_logistic-gini_model
diff.logistic2<-gini_logistic-gini_model_rho2
diff.logs2<-gini_logs-gini_model_rho2

results<-data.frame(input_corr
                    , input_gini1
                    , input_gini2
                    , input_badrate
                    , theor_gini_combined
                    , actual_gini1
                    , actual_gini2
                    , actual_corr
                    , actual_rho
                    , actual_rho2
                    , actual_bad_rate
                    , gini_model
                    , gini_model_rho
                    , gini_model_rho2
                    , gini_logistic
                    , gini_logs
                    , gini_probit
                    , diff.logistic
                    , diff.logistic2
                    , diff.logs2
                    , diff.probit)

#View(results)
write.csv(results, paste0('simresults_flex', '_', mycop, '_', mymarginals, '_', sample_size, '_', number_of_sim, '_', myseed, '_', format(Sys.time(), "%Y%m%d%H%M%S"), '.csv')
          , row.names = FALSE)

return(results)}

#results_1e3_CGNn <- simul(125, 444, 1e3, "CGN", "Normal")
# results_1e4_CGNn <- simul(125, 444, 1e4, "CGN", "Normal")
# results_1e5_CGNn <- simul(125, 444, 1e5, "CGN", "Normal")
# results_1e3_CGNu <- simul(125, 444, 1e3, "CGN", "Uniform")
# results_1e4_CGNu <- simul(125, 444, 1e4, "CGN", "Uniform")
# results_1e5_CGNu <- simul(125, 444, 1e5, "CGN", "Uniform")
# results_1e3_NNCu <- simul(125, 444, 1e3, "NNC", "Uniform")
# results_1e4_NNCu <- simul(125, 444, 1e4, "NNC", "Uniform")
# results_1e5_NNCu <- simul(125, 444, 1e5, "NNC", "Uniform")
# results_1e3_NNCn <- simul(125, 444, 1e3, "NNC", "Normal")
# results_1e4_NNCn <- simul(125, 444, 1e4, "NNC", "Normal")
# results_1e5_NNCn <- simul(125, 444, 1e5, "NNC", "Normal")

# results_1e3_NS <- simul(125, 444, 1e3, "Normal", "Skewed6")
# results_1e4_NS <- simul(125, 444, 1e4, "Normal", "Skewed6")
# results_1e5_NS <- simul(125, 444, 1e5, "Normal", "Skewed6")
# results_1e3_Nu <- simul(125, 444, 1e3, "Normal", "Uniform")
# results_1e4_Nu <- simul(125, 444, 1e4, "Normal", "Uniform")
# results_1e5_Nu <- simul(125, 444, 1e5, "Normal", "Uniform")
# results_1e3_Nn <- simul(125, 444, 1e3, "Normal", "Normal")
# results_1e4_Nn <- simul(125, 444, 1e4, "Normal", "Normal")
# results_1e5_Nn <- simul(125, 444, 1e5, "Normal", "Normal")
# results_1e3_Tn <- simul(125, 444, 1e3, "T4", "Normal")
# results_1e4_Tn <- simul(125, 444, 1e4, "T4", "Normal")
# results_1e5_Tn <- simul(125, 444, 1e5, "T4", "Normal")
# results_1e3_Tu <- simul(125, 444, 1e3, "T4", "Uniform")
# results_1e4_Tu <- simul(125, 444, 1e4, "T4", "Uniform")
# results_1e5_Tu <- simul(125, 444, 1e5, "T4", "Uniform")
# results_1e3_Ts <- simul(125, 444, 1e3, "T4", "Skewed6")
# results_1e4_Ts <- simul(125, 444, 1e4, "T4", "Skewed6")
# results_1e5_Ts <- simul(125, 444, 1e5, "T4", "Skewed6")

library(dplyr)

results_1e3_Tu <- read.csv("simresults_flex_T4_Uniform_1000_444_125_20250403034027.csv")
results_1e4_Tu <- read.csv("simresults_flex_T4_Uniform_10000_444_125_20250403034723.csv")
results_1e5_Tu <- read.csv("simresults_flex_T4_Uniform_100000_444_125_20250403035905.csv")
results_1e3_Ts <- read.csv("simresults_flex_T4_Skewed6_1000_444_125_20250403072637.csv")
results_1e4_Ts <- read.csv("simresults_flex_T4_Skewed6_10000_444_125_20250403073318.csv")
results_1e5_Ts <- read.csv("simresults_flex_T4_Skewed6_100000_444_125_20250403074441.csv")
results_1e3_Tn <- read.csv("simresults_flex_T4_Normal_1000_444_125_20250403031517.csv")
results_1e4_Tn <- read.csv("simresults_flex_T4_Normal_10000_444_125_20250403032212.csv")
results_1e5_Tn <- read.csv("simresults_flex_T4_Normal_100000_444_125_20250403033358.csv")
results_1e3_Nu <- read.csv("simresults_flex_Normal_Uniform_1000_444_125_20250403022535.csv")
results_1e4_Nu <- read.csv("simresults_flex_Normal_Uniform_10000_444_125_20250403023229.csv")
results_1e5_Nu <- read.csv("simresults_flex_Normal_Uniform_100000_999_125_20250402232018.csv")
results_1e3_NS <- read.csv("simresults_flex_Normal_Skewed6_1000_444_125_20250403040536.csv")
results_1e4_NS <- read.csv("simresults_flex_Normal_Skewed6_10000_444_125_20250403041231.csv")
results_1e5_NS <- read.csv("simresults_flex_Normal_Skewed6_100000_444_125_20250403042448.csv")
results_1e3_Nn <- read.csv("simresults_flex_Normal_Normal_1000_444_125_20250403025019.csv")
results_1e4_Nn <- read.csv("simresults_flex_Normal_Normal_10000_444_125_20250403025712.csv")
results_1e5_Nn <- read.csv("simresults_flex_Normal_Normal_100000_444_125_20250403030850.csv")
results_1e3_CGNn <- read.csv("simresults_flex_CGN_Normal_1000_444_125_20250403162307.csv")
results_1e4_CGNn <- read.csv("simresults_flex_CGN_Normal_10000_444_125_20250403202619.csv")
results_1e5_CGNn <- read.csv("simresults_flex_CGN_Uniform_1e+05_444_125_20250403204328.csv")
results_1e3_CGNu <- read.csv("simresults_flex_CGN_Uniform_1000_444_125_20250403203457.csv")
results_1e4_CGNu <- read.csv("simresults_flex_CGN_Uniform_10000_444_125_20250403203757.csv")
results_1e5_CGNu <- read.csv("simresults_flex_CGN_Uniform_1e+05_444_125_20250403204328.csv")
results_1e3_NNCu <- read.csv("simresults_flex_NNC_Uniform_1000_444_125_20250403205722.csv")
results_1e4_NNCu <- read.csv("simresults_flex_NNC_Uniform_10000_444_125_20250403210143.csv")
results_1e5_NNCu <- read.csv("simresults_flex_NNC_Uniform_100000_444_125_20250403210730.csv")
results_1e3_NNCn <- read.csv("simresults_flex_NNC_Normal_1000_444_125_20250403211020.csv")
results_1e4_NNCn <- read.csv("simresults_flex_NNC_Normal_10000_444_125_20250403211322.csv")
results_1e5_NNCn <- read.csv("simresults_flex_NNC_Normal_100000_444_125_20250403211924.csv")

# Create a named list of data frames
dfs <- list(
  "Normal_Skewed6_1e3" = results_1e3_NS,
  "Normal_Skewed6_1e4" = results_1e4_NS,
  "Normal_Skewed6_1e5" = results_1e5_NS,
  "Normal_Uniform_1e3" = results_1e3_Nu,
  "Normal_Uniform_1e4" = results_1e4_Nu,
  "Normal_Uniform_1e5" = results_1e5_Nu,
  "Normal_Normal_1e3" = results_1e3_Nn,
  "Normal_Normal_1e4" = results_1e4_Nn,
  "Normal_Normal_1e5" = results_1e5_Nn,
  "T4_Normal_1e3" = results_1e3_Tn,
  "T4_Normal_1e4" = results_1e4_Tn,
  "T4_Normal_1e5" = results_1e5_Tn,
  "T4_Uniform_1e3" = results_1e3_Tu,
  "T4_Uniform_1e4" = results_1e4_Tu,
  "T4_Uniform_1e5" = results_1e5_Tu,
  "T4_Skewed6_1e3" = results_1e3_Ts,
  "T4_Skewed6_1e4" = results_1e4_Ts,
  "T4_Skewed6_1e5" = results_1e5_Ts,
  "CGN_Normal_1e3" = results_1e3_CGNn,
  "CGN_Normal_1e4" = results_1e4_CGNn,
  "CGN_Normal_1e5" = results_1e5_CGNn,
  "CGN_Uniform_1e3" = results_1e3_CGNu,
  "CGN_Uniform_1e4" = results_1e4_CGNu,
  "CGN_Uniform_1e5" = results_1e5_CGNu,
  "NNC_Uniform_1e3" = results_1e3_NNCu,
  "NNC_Uniform_1e4" = results_1e4_NNCu,
  "NNC_Uniform_1e5" = results_1e5_NNCu,
  "NNC_Normal_1e3" = results_1e3_NNCn,
  "NNC_Normal_1e4" = results_1e4_NNCn,
  "NNC_Normal_1e5" = results_1e5_NNCn
)

# Stack all data frames and add a column for the source name
stacked_results <- bind_rows(dfs, .id = "Simulation")

# Print result
head(stacked_results)

gini_conf_width<-Vectorize(function(auc, brate, n){if(any(is.na(brate), is.na(auc))) {NA} else 2*{presize::prec_auc(auc, brate, n)$conf.width}})
auc_from_gini<-Vectorize(function(x){x/2+.5})
options(scipen=999)
stacked_results %>% group_by(Simulation) %>%
  summarize(
    m_mnv_log  = mean(gini_model-gini_logistic, na.rm=TRUE),
    sd_mnv_log  = sd(gini_model-gini_logistic, na.rm=TRUE),
    m_log_the  = mean(gini_logistic-theor_gini_combined, na.rm=TRUE),
    sd_log_the  = sd(gini_logistic-theor_gini_combined, na.rm=TRUE),
    m_log_logs = mean(gini_logistic-gini_logs, na.rm=TRUE),
    sd_log_logs = sd(gini_logistic-gini_logs, na.rm=TRUE),
    m_mnv_the  = mean(gini_model-theor_gini_combined, na.rm=TRUE),
    sd_mnv_the  = sd(gini_model-theor_gini_combined, na.rm=TRUE),
#    m_mnvs_the = mean(gini_model_rho-theor_gini_combined, na.rm=TRUE),
#    sd_mnvs_the = sd(gini_model_rho-theor_gini_combined, na.rm=TRUE),
    m_mnvs2_the = mean(gini_model_rho2-theor_gini_combined, na.rm=TRUE),
    sd_mnvs2_the = sd(gini_model_rho2-theor_gini_combined, na.rm=TRUE),
    m_logs_the = mean(gini_logs-theor_gini_combined, na.rm=TRUE),
    sd_logs_the = sd(gini_logs-theor_gini_combined, na.rm=TRUE),
    m_mnvs2_logs = mean(gini_model_rho2-gini_logs, na.rm=TRUE),
    sd_mnvs2_logs = sd(gini_model_rho2-gini_logs, na.rm=TRUE),
    m_mnvs2_logs = mean(gini_model_rho2-gini_logs, na.rm=TRUE),
    sd_mnvs2_logs = sd(gini_model_rho2-gini_logs, na.rm=TRUE),
    m_mnvs2_mnv = mean(gini_model_rho2-gini_model, na.rm=TRUE),
    sd_mnvs2_mnv = sd(gini_model_rho2-gini_model, na.rm=TRUE),
gini_se = mean(gini_conf_width(auc_from_gini(gini_logs), actual_bad_rate, as.numeric(substring(Simulation, nchar(Simulation) - 2))), na.rm=TRUE)/(2* qnorm(.975))
  ) %>% print(n=200)

# library(dplyr)
# results_1e3_CGNn %>% filter(round(actual_corr,1)==0.5) %>% filter(round(actual_gini1,1)==0.6)
