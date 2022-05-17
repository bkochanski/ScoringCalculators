#simulation3 - when no solution
#theor_gini_combined<-gini_combine_1(input_gini1, input_gini2, input_corr, input_badrate)
# theor_gini_combined
# i=141
# theor_gini_combined[i]
# gini1_here<-input_gini1[i]
gini1_here<-0.61821019
# gini2_here<-input_gini2[i]
gini2_here<-0.4477561
# corr_here<-input_corr[i]
corr_here<-0.76458430
# badrate_here<-input_badrate[i]
badrate_here<-0.027
phi_s1<-function(x){gini_from_r(rho=x, defrate=badrate_here)-gini1_here}
(rho_s1<-uniroot(phi_s1,lower=0,upper=.9999,tol = .Machine$double.eps)$root)
phi_s2<-function(x){gini_from_r(rho=x, defrate=badrate_here)-gini2_here}
(rho_s2<-uniroot(phi_s2,lower=0,upper=.9999,tol = .Machine$double.eps)$root)
(a_opt<-(corr_here*rho_s2-rho_s1)/(corr_here*rho_s1-rho_s2))
corr_opt<-(a_opt*rho_s1+rho_s2)/sqrt(a_opt^2+2*a_opt*corr_here+1)
g_result0<-gini_from_r(corr_opt, badrate_here)

#sum(is.na(theor_gini_combined)!=is.nan(theor_gini_combined))
