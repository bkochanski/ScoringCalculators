library(MASS)

ginic<-function(bc, gc){
  #function for gini when we have cumulative goods and bads vectors
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*
        (bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1
} 

gini_crd<-function(rho=0.5, defrate=0.1, gran=10000) {
  #function translating rho into gini for a given defrate
  drates_i<-pnorm(qnorm(defrate), rho*(qnorm(1:(gran-1)/gran)),sqrt(1-rho^2))
  drates_2i<-(c(1, drates_i)+c(drates_i, 0))/2
  return(ginic(cumsum((drates_2i)/sum(drates_2i)), cumsum((1-drates_2i)/sum(1-drates_2i))))
}

r_from_gini<-function(gini, defaultrate=.1){
  phi_s1<-function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-gini}
  rho_s1<-uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
}

gini_combine_calculator<-function(g1, g2, corr, defaultrate){
  #rho_s1
  phi_s1<-function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-g1}
  rho_s1<-uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
  #rho_s2
  phi_s2<-function(x){gini_crd(rho=x, defrate=defaultrate, gran=100000)-g2}
  rho_s2<-uniroot(phi_s2,lower=0,upper=1,tol = .Machine$double.eps)$root
  
  (a_opt<-(corr*rho_s2-rho_s1)/(corr*rho_s1-rho_s2))
  corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr+1)
  g_result0<-gini_crd(corr_opt, defaultrate, gran=100000)
  g_result1<-if(a_opt<0 | a_opt>1000) {NaN} else {g_result0}
  return(c(new_gini=g_result1, 
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt, gini1=g1, gini2=g2, corr=corr, defrate=defaultrate))
}

#simulation attempt1

# input g1, g2, b, r
# gini_new_from_theor, a_from_theor
# sample1 g1, g2, b, r, gini_new_from_samp1, a_from_samp_1, gini_new_a_from_samp1, a_logit_samp1, gini_a_logit_samp1, amax, a_gini_max, gini_a_gini_max)
# sample2 (...) gini_samp2_a_from_samp1, gini_logit1


#set.seed(500)
m <- 3
n <- 5000
gini1<-.4
gini2<-.6
corrs<-.4
dfrate<-.1
sigma <- matrix(c(1, corrs, r_from_gini(gini1, dfrate),
                  corrs, 1, r_from_gini(gini2, dfrate),
                  r_from_gini(gini1, dfrate), r_from_gini(gini2, dfrate), 1), 
                nrow=3)



tot_results<-list()

for (i in 1:100){
z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma)
df<-z
df[,3]<-(z[,3]<qnorm(dfrate))*1
df<-data.frame(df)
names(df)<-c('s1', 's2', 'default_flag')

rest<-gini_combine_calculator(gini1, gini2, corrs, dfrate)

my_gini(df$default_flag, df$s1)
my_gini(df$default_flag, df$s2)
cor(df$s1, df$s2)
mean(df$default_flag)
res1<-gini_combine_calculator(my_gini(df$default_flag, df$s1), my_gini(df$default_flag, df$s2), cor(df$s2, df$s1), mean(df$default_flag))

gini_a_opt<-my_gini(df$default_flag, res1['a_opt']*df$s1+df$s2)

(mylogit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = "binomial"))
(myprobit <- glm(I(1-default_flag) ~ s1 + s2, data = df, family = binomial(link="probit")))
(logit1_a<-mylogit$coefficients[2]/mylogit$coefficients[3])
(probit1_a<-myprobit$coefficients[2]/myprobit$coefficients[3])

#my_gini(mylogit$y, -mylogit$fitted.values)
gini_logit1<-my_gini(df$default_flag, logit1_a*df$s1+df$s2)

#my_gini(myprobit$y, -myprobit$fitted.values)
gini_probit1<-my_gini(df$default_flag, probit1_a*df$s1+df$s2)

f<-function(x){my_gini(df$default_flag, x*df$s1+df$s2)}
res1_optim<-optimize(f, c(0,1000), maximum = TRUE)
#my_gini(df$default_flag, res1_optim$maximum*df$s1+df$s2)

z2 <- mvrnorm(n,mu=rep(0, m),Sigma=sigma)
#library(psych)
# cor(z)
# cor(z,method='spearman')
df2<-z2
df2[,3]<-(z2[,3]<qnorm(dfrate))*1
df2<-data.frame(df2)
names(df2)<-c('s1', 's2', 'default_flag')

tot_results[[length(tot_results)+1]]<-c(i=i, 
                                    n=n, 
                                    # rest,
                                    res1, 
                                    gini_a_opt=gini_a_opt,
                                    gini_logit1=gini_logit1,
                                    gini_probit1=gini_probit1,
                                    res1_optim, 
                                    sample2_gini_a_opt_t=my_gini(df2$default_flag, rest['a_opt']*df2$s1+df2$s2),
                                    sample2_gini_a_opt_1=my_gini(df2$default_flag, res1['a_opt']*df2$s1+df2$s2),
                                    sample2_gini_logit1=my_gini(df2$default_flag, logit1_a*df2$s1+df2$s2),
                                    sample2_gini_probit1=my_gini(df2$default_flag, probit1_a*df2$s1+df2$s2),
                                    sample2_gini_max1=my_gini(df2$default_flag, res1_optim$maximum*df2$s1+df2$s2))
}

bind_tot_results<-do.call("rbind", tot_results)
df_tot_results<-data.frame(matrix(unlist(bind_tot_results), ncol=length(tot_results[[1]])))
names(df_tot_results)<-names(tot_results[[1]])

hist(df_tot_results$new_gini)
sd(df_tot_results$new_gini)
