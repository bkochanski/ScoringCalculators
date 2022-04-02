GINI1<-.45
GINI2<-.65
B<-0.1
a0<-0.6
sh<-0.83

# GINI1<-.45
# GINI2<-.46
# B<-0.1
# a0<-0.6


FuncNormalAlt<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}
FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}


# y_e<-function(x){FuncNormalAlt(x,GINI1)}
# y_n<-function(x){FuncNormalAlt(x,GINI2)}
# y_e<-function(x){FuncMidFractal(x,GINI1)}
# y_n<-function(x){FuncMidFractal(x,GINI2)}
y_e<-function(x){FuncBiNormal(x,GINI1, sh)}
y_n<-function(x){FuncBiNormal(x,GINI2, sh)}

GiniP<-function(f, x){2*(integrate(f,x,1)$value-(1-x)*f(x))/((1-x)*(1-f(x)))-1
}

# basic scenario
phi_x0<-function(x){(1-B)*(1-x)+B*(1-y_e(x))-a0}
x0<-uniroot(phi_x0,lower=0,upper=1,tol = .Machine$double.eps)
(x0<-as.numeric(x0[1]))
y_e(x0)

(b0=B*(1-y_e(x0))/a0)
(deriv0<-numDeriv::grad(y_e, x0))
(mbr0<-(1+(1-B)/B/deriv0)^(-1))
ir0<-mbr0/(1-mbr0)
(profit0<-a0*(ir0*(1-b0)-b0))

### TEST
plot(c(x0),c(y_e(x0)), main='', xlab='', ylab='', lwd=4, 
     xlim=c(0,1), ylim=c(0,1))
curve(y_e, add=TRUE, lwd=1)


(ginip0<-GiniP(y_e,x0))
#GiniP(y_e,0)


#numDeriv::grad(y_e, x0)

#(1+GINI1)/(1-GINI1)*.5*(1-x0)^(2*GINI1/(1-GINI1))+(1-GINI1)/(1+GINI1)*0.5*x0^(-2*GINI1/(1+GINI1))

#scenario 1 keep approval rate stable
a1<-a0
phi_yn_x1<-function(x){(1-B)*(1-x)+B*(1-y_n(x))-a1}
x1<-uniroot(phi_yn_x1,lower=0,upper=1,tol = .Machine$double.eps)
(x1<-as.numeric(x1[1]))
y_n(x1)
(b1=B*(1-y_n(x1))/a1)
(deriv1<-numDeriv::grad(y_n, x1))
(mbr1<-(1+(1-B)/B/deriv1)^(-1))
ir1<-ir0
(profit1<-a1*(ir1*(1-b1)-b1))
(ginip1<-GiniP(y_n,x1))


(a1_change<-a1/a0-1)
(b1_change<-b1/b0-1)
(profit1_change<-profit1/profit0-1)

# scenario 2 keep bad rate stable
b2<-b0
phi_x2<-function(x){B*(1-y_n(x))/((1-B)*(1-x)+B*(1-y_n(x)))-b2}
x2<-uniroot(phi_x2,lower=0.001,upper=.999,tol = .Machine$double.eps)
(x2<-as.numeric(x2[1]))
y_n(x2)
(a2<-((1-B)*(1-x2)+B*(1-y_n(x2))))
(deriv2<-numDeriv::grad(y_n, x2))
(mbr2<-(1+(1-B)/B/deriv2)^(-1))
ir2<-ir0
(profit2<-a2*(ir2*(1-b2)-b2))
(ginip2<-GiniP(y_n,x2))

(a2_change<-a2/a0-1)
(b2_change<-b2/b0-1)
(profit2_change<-profit2/profit0-1)


#scenario 3: equal marginal rates / profit maximized

phi_yn_x3<-function(x){numDeriv::grad(y_n, x)-deriv0}
x3<-uniroot(phi_yn_x3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)
(x3<-as.numeric(x3[1]))
y_n(x3)
(a3<-((1-B)*(1-x3)+B*(1-y_n(x3))))
(b3=B*(1-y_n(x3))/a3)
(deriv3<-numDeriv::grad(y_n, x3))
(mbr3<-(1+(1-B)/B/deriv3)^(-1))
ir3<-ir0
(profit3<-a3*(ir3*(1-b3)-b3))
(ginip3<-GiniP(y_n,x3))

(a3_change<-a3/a0-1)
(b3_change<-b3/b0-1)
(profit3_change<-profit3/profit0-1)


plot(c(x0, x1, x2, x3),c(y_e(x0), y_n(x1), y_n(x2), y_n(x3)), main='', xlab='', ylab='', lwd=4, 
     xlim=c(0,1), ylim=c(0,1))
curve(y_e, add=TRUE, lwd=1)
curve(y_n, add=TRUE, lty=2)

print("Bad rate reduction: ")
(b1/b0-1)*100
print("Approval rate improvement: ")
(a2/a0-1)*100

data.frame(scenario=1:3, 
          existing_portfolio_gini=rep(ginip0,3),
          new_portfolio_gini=c(ginip1, ginip2, ginip3),
          existing_marginal_bad_rate=rep(mbr0,3),
          new_marginal_bad_rate=c(mbr1, mbr2, mbr3)
          )

data.frame(scenario=1:3, existing_approval=100*rep(a0,3), 
           new_approval=100*c(a1, a2, a3),
           existing_badrate=100*rep(b0,3),
           new_badrate=100*c(b1,b2,b3),
           existing_profit=100*rep(profit0,3),
           new_profit=100*c(profit1, profit2, profit3)
           )

data.frame(no=1:3, scenario=c('keep approval', 'keep bad rate', 'keep marginal bad rate'), approval_change=100*c(a1_change, a2_change, a3_change),
           badrate_change=100*c(b1_change, b2_change, b3_change),
           profit_change=100*c(profit1_change, profit2_change, profit3_change))


