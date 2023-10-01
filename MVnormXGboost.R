library(xgboost)
library(caret)  
library(e1071)

m <- 3
n <- 100000
gini1<-.4
gini2<-.6
corrs<-.4
dfrate<-.1

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


sigma <- matrix(c(1, corrs, r_from_gini(gini1, dfrate),
                  corrs, 1, r_from_gini(gini2, dfrate),
                  r_from_gini(gini1, dfrate), r_from_gini(gini2, dfrate), 1), 
                nrow=3)
z <- MASS::mvrnorm(n,mu=rep(0, m),Sigma=sigma)
simdf<-z
simdf[,3]<-(z[,3]<qnorm(dfrate))*1
simdf<-data.frame(simdf)
names(simdf)<-c('s1', 's2', 'bad')


data<-simdf[,c('s1', 's2', 'bad')]
parts = createDataPartition(data$bad, p = 0.7, list = F)
train = data[parts, ]
test = data[-parts, ]
X_train = data.matrix(train[,-3])                  # independent variables for train
y_train = train[,3]                                # dependent variables for train

X_test = data.matrix(test[,-3])                    # independent variables for test
y_test = test[,3]                                   # dependent variables for test

xgboost_train = xgb.DMatrix(data=X_train, label=y_train)
xgboost_test = xgb.DMatrix(data=X_test, label=y_test)

model <- xgboost(data = xgboost_train,                    # the data   
                 max.depth=3,                          # max depth 
                 nrounds=50)                              # max number of boosting iterations

summary(model)
pred_test = predict(model, xgboost_test)
hist(pred_test)
pred_train = predict(model, xgboost_train)
hist(pred_train)

my_gini_roc<-function(resp, pred){
  #c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  plot(gc, bc)
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 

my_gini_roc(y_test, -pred_test)
my_gini_roc(y_train, -pred_train)


d1<-y_test[order(X_test[,1])]
bc1<-c(0, cumsum(d1)/sum(d1))
gc1<-c(0, cumsum(1-d1)/sum(1-d1))
plot(gc1, bc1, type='l')
d2<-y_test[order(X_test[,2])]
bc2<-c(0, cumsum(d2)/sum(d2))
gc2<-c(0, cumsum(1-d2)/sum(1-d2))
lines(gc2, bc2, col='blue')

w1=0.26
d4<-y_test[order(w1*X_test[,1]+(1-w1)*X_test[,2])]
bc4<-c(0, cumsum(d4)/sum(d4))
gc4<-c(0, cumsum(1-d4)/sum(1-d4))
lines(gc4, bc4, col='orange')

d3<-y_test[order(-pred_test)]
bc3<-c(0, cumsum(d3)/sum(d3))
gc3<-c(0, cumsum(1-d3)/sum(1-d3))
lines(gc3, bc3, col='red')

linmodel<-glm(bad~s1+s2, data=train)
linmodel$coefficients[2]/(linmodel$coefficients[2]+linmodel$coefficients[3])
pred_l<-predict(linmodel, newdata=test)

linmodel_inter<-glm(bad~s1*s2, data=train)
pred_li<-predict(linmodel_inter, newdata=test)

my_gini_roc(test$bad, -pred_test)
my_gini_roc(test$bad, -pred_li)
my_gini_roc(test$bad, -pred_l)
