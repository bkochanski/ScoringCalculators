GINI1<-.8
GINI2<-.6
B<-.1
a0<-.5

FuncMidNormal<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
#FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}
GiniP<-function(f, x){2*(integrate(f,x,1)$value-(1-x)*f(x))/((1-x)*(1-f(x)))-1}
 
a3_in<-function(GINI1, GINI2, B, a0)
{ tryCatch({
y0<-function(x){FuncMidNormal(x,GINI1)}
y1<-function(x){FuncMidNormal(x,GINI2)}

phi0<-function(x){(1-B)*(1-x)+B*(1-y0(x))-a0}
x0<-as.numeric(uniroot(phi0,lower=0,upper=1,tol = .Machine$double.eps)[1])
b0 <- B*(1-y0(x0))/a0

deriv0 <-  numDeriv::grad(y0, x0)
mbr0 <- (1+(1-B)/B/deriv0)^(-1)
ir0 <-  mbr0/(1-mbr0)
profit0 <- a0*(ir0*(1-b0)-b0)
# ginip0 <- GiniP(y0,x0)

phi3<-function(x){numDeriv::grad(y1, x)-deriv0}
x3 <- as.numeric(uniroot(phi3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])
a3 <- ((1-B)*(1-x3)+B*(1-y1(x3)))
b3 <- B*(1-y1(x3))/a3

# deriv3 <- numDeriv::grad(y1, x3)
# mbr3 <- (1+(1-B)/B/deriv3)^(-1)
#ir3 <- ir0
profit3 <- a3*(ir0*(1-b3)-b3)
# ginip3 <-GiniP(y1,x3)

# b3_change <- b3/b0-1
# profit3/profit0-1

a3/a0-1}, error=function(err) NA)
}

options(scipen=99)

a3_inaVV<-Vectorize(function(y){
a3_inaV<-Vectorize(function(x){a3_in(.6,.4,y,x)})
a3_inaV(1:99/100)})

m<-a3_inaVV(1:60/100)
#View(m)

m.df <- reshape2::melt(m, c("x", "y"), value.name = "z")

library(ggplot2)
ggplot(m.df, aes(x = x, y = y, z = z)) +
  stat_contour(aes(colour = ..level..)) +geom_contour_filled() + 
  xlab("initial approval [%]") +
  ylab("population bad rate [%]")

# xlab('Gini 2') + ylab('Correlation') 

# m.df[m.df$x==75 & m.df$y==50,]
# m.df[m.df$x==50 & m.df$y==75,]
