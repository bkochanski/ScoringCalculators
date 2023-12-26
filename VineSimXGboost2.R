library(xgboost)
library(caret)  
library(e1071)

data<-simdf[,c('s1n', 's2n', 'bad')]
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

w1=0.3358
d4<-y_test[order(w1*X_test[,1]+(1-w1)*X_test[,2])]
bc4<-c(0, cumsum(d4)/sum(d4))
gc4<-c(0, cumsum(1-d4)/sum(1-d4))
lines(gc4, bc4, col='orange')

d3<-y_test[order(-pred_test)]
bc3<-c(0, cumsum(d3)/sum(d3))
gc3<-c(0, cumsum(1-d3)/sum(1-d3))
lines(gc3, bc3, col='red')

linmodel<-glm(bad~s1n+s2n, data=train)
linmodel$coefficients[2]/(linmodel$coefficients[2]+linmodel$coefficients[3])
pred_l<-predict(linmodel, newdata=test)


linmodel_inter<-glm(bad~s1n*s2n, data=train)
pred_li<-predict(linmodel_inter, newdata=test)

logmodel<-glm(bad~s1n+s2n, data=train, family="binomial")
pred_log<-predict(logmodel, newdata=test)

logmodel_inter<-glm(bad~s1n*s2n, data=train, family="binomial")
pred_logi<-predict(logmodel_inter, newdata=test)


my_gini_roc(test$bad, -pred_test)
my_gini_roc(test$bad, -pred_li)
my_gini_roc(test$bad, -pred_l)
my_gini_roc(test$bad, -pred_log)
my_gini_roc(test$bad, -pred_logi)

gini_combine_calculator(my_gini_roc(train$bad, train$s1n), 
                        my_gini_roc(train$bad, train$s2n),
                        cor(train$s1n, train$s2n),
                        mean(train$bad))

