a<-MASS::mvrnorm(n=100, mu=c(0,0,0), Sigma = matrix(c(1,0,0,0,1,0,0,0,1), nrow=3))
a[,3]<-(a[,3]< -0.6)*1
pred1<--a[,1]
targ1<-a[,3]
pred1[2]<-pred1[1]
2*bigstatsr::AUC(pred1, targ1)-1

#Gini
2/length(targ1)*(mean(rank(pred1)[targ1==1])-mean(rank(pred1)[targ1==0]))

#AUC
(mean(rank(pred1)[targ1==1])-mean(rank(pred1)[targ1==0]))/length(targ1)+1/2

