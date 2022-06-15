library(MASS)
#set.seed(500)
m <- 3
n <- 1000
sigma <- matrix(c(1, 0.4, 0.45,
                  0.4, 1, .65,
                  0.45, 0.65, 1), 
                nrow=3)
z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma)
#library(psych)
cor(z)
cor(z,method='spearman')
# pairs.panels(z)
# library(rgl)
# plot3d(z[,1],z[,2],z[,3],pch=20,col='navyblue')


y<-(z[,3]<qnorm(.1))*1
df<-z
df[,3]<-y
df<-data.frame(df)
names(df)<-c('s1', 's2', 'default_flag')

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

my_gini(df$default_flag, df$s1)

InformationValue::somersD(df$default_flag, 1-df$s1)

DescTools::GoodmanKruskalGamma(df$default_flag, -df$s1)


wilcox.test(df$s1~df$default_flag)

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

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

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
           #new_gini_2=g_result0, 
           a_opt=a_opt, score_1_weight=a_opt/(1+a_opt), score_2_weight=1/(1+a_opt), 
           rho1=rho_s1, rho2=rho_s2, new_corr=corr_opt))
}




(res1<-gini_combine_calculator(gini_crd(.45), gini_crd(.65), 0.4, 0.1))
res1[3]/res1[4]

(res2<-gini_combine_calculator(gini_crd(cor(z)[1,3]), gini_crd(cor(z)[2,3]), 0.4, 0.1))
res2[3]/res2[4]

(res3<-gini_combine_calculator(my_gini(df$default_flag, df$s1), my_gini(df$default_flag, df$s2), 0.4, 0.1))
res3[3]/res3[4]


(mylogit <- glm(y ~ s1 + s2, data = df, family = "binomial"))
(myprobit <- glm(y ~ s1 + s2, data = df, family = binomial(link="probit")))
mylogit$coefficients[2]/mylogit$coefficients[3]
myprobit$coefficients[2]/myprobit$coefficients[3]
res1[3]/res1[4]
res2[3]/res2[4]


my_gini(mylogit$y, -mylogit$fitted.values)
my_gini(mylogit$y, -myprobit$fitted.values)
res1[1]
my_gini(df$default_flag, res1[3]*df$s1+res1[4]*df$s2)
res2[1]
my_gini(df$default_flag, res2[3]*df$s1+res2[4]*df$s2)
res3[1]
my_gini(df$default_flag, res3[3]*df$s1+res3[4]*df$s2)

cor(df)
my_gini(df$default_flag, df$s1)
my_gini(df$default_flag, df$s2)
cor(df$default_flag, df$s1)
cor(df$default_flag, df$s2)
cor(df$default_flag, df$s1, method="spearman")
cor(df$default_flag, df$s2, method="spearman")

mean(df$s1>0)
cor(df$s2, df$s1>0)
cor(df$s2, df$s1>0, method="spearman")
my_gini(df$s1>0, df$s2)

mean(df$s1>1)
cor(df$s2, df$s1>1)
cor(df$s2, df$s1>1, method="spearman")

mean(df$s1>2)
cor(df$s2, df$s1>2)
mean(df$s1>3)
cor(df$s2, df$s1>3)


#correlations
cors<-c((-5:9)/10,.99)
ginic4<-function(x){gini_combine_calculator(.6,.6,x,.1)[1]}
ginic4<-Vectorize(ginic4)
ginis<-ginic4(cors)
plot(cors, ginis)
#gini_combine_calculator(.4,.4,-.5,.1)

