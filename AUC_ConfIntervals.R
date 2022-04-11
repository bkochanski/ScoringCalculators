#AUC confidence intervals
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

gini1<-my_gini(df$default_flag, df$s1)
gini2<-my_gini(df$default_flag, df$s2)

auc_from_gini<-function(x){x/2+.5}
auc_from_gini<-Vectorize(auc_from_gini)

auc1<-auc_from_gini(gini1)
auc2<-auc_from_gini(gini2)

library(pROC)
roc1 <- roc(df$default_flag, df$s1)
ci.auc(roc1)
ci.auc(roc1, method="bootstrap")
presize::prec_auc(auc1, mean(df$default_flag), n, conf.level = .95)
