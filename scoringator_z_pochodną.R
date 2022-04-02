GINI1<-.45
GINI2<-.65
B<-0.1
a0<-0.4

# GINI1<-.45
# GINI2<-.46
# B<-0.1
# a0<-0.6


#FuncNormalAlt<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}
# y_e<-function(x){FuncNormalAlt(x,GINI1)}
# y_n<-function(x){FuncNormalAlt(x,GINI2)}
y_e<-function(x){FuncMidFractal(x,GINI1)}
y_n<-function(x){FuncMidFractal(x,GINI2)}

#scenario 1 keep approval rate stable
a1<-a0
phi_a1<-function(x){(1-B)*(1-x)+B*(1-y_e(x))-a1}
x0<-uniroot(phi_a1,lower=0,upper=1,tol = .Machine$double.eps)
x0<-as.numeric(x0[1])
x0
y_e(x0)

#numDeriv::grad(y_e, x0)

#(1+GINI1)/(1-GINI1)*.5*(1-x0)^(2*GINI1/(1-GINI1))+(1-GINI1)/(1+GINI1)*0.5*x0^(-2*GINI1/(1+GINI1))

(deriv0<-numDeriv::grad(y_e, x0))

phi_yn_x1<-function(x){(1-B)*(1-x)+B*(1-y_n(x))-a0}
x1<-uniroot(phi_yn_x1,lower=0,upper=1,tol = .Machine$double.eps)
x1<-as.numeric(x1[1])
x1
y_n(x1)

phi_yn_x3<-function(x){numDeriv::grad(y_n, x)-deriv0}
x2<-uniroot(phi_yn_x3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)
x2<-as.numeric(x2[1])
x2
y_n(x2)
#check the derivative
#delta<-1/10^6
#(y_n(x2+delta)-y_n(x2-delta))/2/delta

interest_rate<-B/deriv0

(b0=B*(1-y_e(x0))/a0)
(profit0<-a0*(interest_rate*(1-b0)-b0))
(b1=B*(1-y_n(x1))/a0)
(b1/b0-1)*100
(profit1<-a0*(interest_rate*(1-b1)-b1))
profit1/profit0-1

(a2<-((1-B)*(1-x2)+B*(1-y_n(x2))))
(b2=B*(1-y_n(x2))/a2)
(b2/b0-1)*100
(profit2<-a2*(interest_rate*(1-b2)-b2))
profit2/profit0-1
a2/a0-1


(y_e(x0)-y_n(x1))/(1-y_e(x0))*100

phi3<-function(x){B*(1-y_n(x))/((1-B)*(1-x)+B*(1-y_n(x)))-b0}
x1prime<-uniroot(phi3,lower=0,upper=.99999999999,tol = .Machine$double.eps)
(x1prime<-as.numeric(x1prime[1]))
(a1<-((1-B)*(1-x1prime)+B*(1-y_n(x1prime))))
(a1/a0-1)*100
(profit3<-a1*(interest_rate*(1-b0)-b0))
profit3/profit0-1
derivative3<-(y_n(x1prime+delta)-y_n(x1prime-delta))/2/delta
(marginal3<-B/derivative3)


plot(c(x0, x1, x1prime, x2),c(y_e(x0), y_n(x1), y_n(x1prime), y_n(x2)), main='', xlab='', ylab='', lwd=4, 
     xlim=c(0,1), ylim=c(0,1))
curve(y_e, add=TRUE, lwd=1)
curve(y_n, add=TRUE, lty=2)

print("Bad rate reduction: ")
(b1/b0-1)*100
print("Approval rate improvement: ")
(a1/a0-1)*100
