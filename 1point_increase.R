test<-better_gini_calculator()
test$res1$approval_change[2]

profit1point<-Vectorize(function(gini1, B, a0){
  as.numeric(better_gini_calculator(gini1, gini1+.01, B, a0)$res1$profit_change[3])
})
profit1perc<-Vectorize(function(gini1, B, a0){
  as.numeric(better_gini_calculator(gini1, gini1*1.01, B, a0)$res1$profit_change[3])
})
approval1point<-Vectorize(function(gini1, B, a0){
  as.numeric(better_gini_calculator(gini1, gini1+.01, B, a0)$res1$approval_change[2])
})



better_gini_calculator(.4, .41, .05, .6)

ginis<-40:95/100
B=.01
#f=profit1point
f=approval1point
plot(ginis, f(ginis, B, .4)*100, ylim=c(0,6), ylab='')
points(ginis, f(ginis, B, .5)*100, col='red')
points(ginis, f(ginis, B, .6)*100, col='blue')
points(ginis, f(ginis, B, .7)*100, col='green')
points(ginis, f(ginis, B, .8)*100, col='orange')
points(ginis, f(ginis, B, .9)*100, col='purple')

