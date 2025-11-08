library(ggplot2)
library(ggrepel)
library(numDeriv)
#library(data.table)


# FUNCTIONS:
# FuncMidNormal<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
# FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}
# FuncBiFractal<-function(x,g, beta){beta*(1-(1-x)^((1+g)/(1-g)))+(1-beta)*x^((1-g)/(1+g))}
# FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}

# In this test we are using proper binormal curve
# FuncBiProper
Gini2Da <- function(g, c_){
  ffind <- function(x){pnorm(x/sqrt(2)) + 2 * VGAM::pbinorm(-x/sqrt(2), 0, cov12 = -(1-c_^2)/(1+c_^2))-(g+1)/2}
  return(uniroot(ffind, lower = 0, upper = 10)$root)
}
H <- function(x){1*(x>=0)}
fpf <- function(v_c, c, d_a) {
  return(pnorm(-(1-c)*v_c-d_a/2*sqrt(1+c^2)) + ( pnorm(-(1-c)*v_c+d_a/2/c*sqrt(1+c^2)) - H(c)))
}
tpf <- function(v_c, c, d_a) {
  return(pnorm(-(1+c)*v_c+d_a/2*sqrt(1+c^2)) + ( pnorm(-(1+c)*v_c+d_a/2/c*sqrt(1+c^2)) - H(c)))
}
v_c_from_fpf <- function(fpfF, fpf_){
  uniroot(function(x){fpfF(x)-fpf_}, lower = -100, upper = 100)$root
}
v_c_from_tpf <- function(tpfF, tpf_){
  uniroot(function(x){tpfF(x)-tpf_}, lower = -100, upper = 100)$root
}

FuncBiProper <- Vectorize(function(x,c,d_a){tpf(v_c_from_fpf(function(x){fpf(x,c,d_a)},x), c, d_a)})
FuncBiProperDeriv <- Vectorize(function(x, c1, d_a1){
#  s = v_c_from_tpf(function(x_){tpf(x_,c1,d_a1)},x)
s = v_c_from_fpf(function(x_){fpf(x_,c1,d_a1)},x)
#  res = ((1-c1)*(dnorm(-(1-c1)*s-d_a1/2*sqrt(1+c1^2)) + dnorm(-(1-c1)*s+d_a1/2/c1*sqrt(1+c1^2))))/((1+c1)*(dnorm(-(1+c1)*s+d_a1/2*sqrt(1+c1^2)) + dnorm(-(1+c1)*s+d_a1/2/c1*sqrt(1+c1^2)) ))
  res = ((1+c1)*(dnorm(-(1+c1)*s+d_a1/2*sqrt(1+c1^2)) + dnorm(-(1+c1)*s+d_a1/2/c1*sqrt(1+c1^2)) ))/((1-c1)*(dnorm(-(1-c1)*s-d_a1/2*sqrt(1+c1^2)) + dnorm(-(1-c1)*s+d_a1/2/c1*sqrt(1+c1^2))))
  return(res)
})

# function GiniP calculates Gini above the cutoff (truncated Gini) calculation
# f - ROC curve model (e.g. binormal with set Gini & shape parameters)
# x - cut-off point representation on the ROC curve (x coordinate: cumulative good proportion below the cutoff)
GiniP<-function(f, x){2*(integrate(f,x,1)$value-(1-x)*f(x))/((1-x)*(1-f(x)))-1}

# Calculator function
# Please note default c1 and c2 is set to one, so the ROC curves default to MidNormal
better_gini_calculatorP<-function(GINI1=.45, GINI2=.65, B=.1, a0=.6, c1=0, c2=0){
  da1 <- Gini2Da(GINI1, c1)
  da2 <- Gini2Da(GINI2, c2)

  # setting binormal functions with GINI1 and GINI2 parameters
  y0<-function(x){FuncBiProper(x,c1,da1)}
  y1<-function(x){FuncBiProper(x,c2,da2)}
  
  # searching for representation of the cutoff point 
  phi0<-function(x){(1-B)*(1-x)+B*(1-y0(x))-a0}
  x0<-as.numeric(uniroot(phi0,lower=0,upper=1,tol = .Machine$double.eps)[1])
  
  # determining portfolio bad rate for scorecard 1
  b0 <- B*(1-y0(x0))/a0
  
  # determining marginal bad rate at the cutoff for scorecard 1 
  #deriv0 <-  numDeriv::grad(y0, x0)
  deriv0 <- FuncBiProperDeriv(x0, c1, da1)
  mbr0 <- (1+(1-B)/B/deriv0)^(-1)
  
  # determining interest rate assuming profit = 0 at marginal bad rate 
  ir0 <-  mbr0/(1-mbr0)
  
  # calculation of the portfolio profit
  profit0 <- a0*(ir0*(1-b0)-b0)
  
  # Gini for scorecard 1 on accepted portfolio
  ginip0 <- GiniP(y0,x0)
  
  ## 1. Bad rate reduction scenario (keep approval)
  
  # same approval
  a1 <- a0
  
  # finding the cutoff representation on ROC curve 2 for scenario 1
  phi1<-function(x){(1-B)*(1-x)+B*(1-y1(x))-a1}
  x1<-as.numeric(uniroot(phi1,lower=0,upper=1,tol = .Machine$double.eps)[1])
  
  # new portfolio bad rate
  b1 <- B*(1-y1(x1))/a0
  
  # new marginal bad rate at the cutoff
  #deriv1 <- numDeriv::grad(y1, x1)
  deriv1 <- FuncBiProperDeriv(x1, c2, da2)
  mbr1 <- (1+(1-B)/B/deriv1)^(-1)
  
  # same interest rate
  ir1 <- ir0
  
  # profit for scorecard 2
  profit1 <- a1*(ir1*(1-b1)-b1)
  
  # calculating Gini on approved portfolio for scorecard 2 
  ginip1 <- GiniP(y1,x1)
  
  # approval change, portfolio bad rate change, profit change
  a1_change <- a1/a0-1
  b1_change <- b1/b0-1
  profit1_change <- profit1/profit0-1
  
  ## 2. Approval rate increase scenario (keep bad rate)
  
  # same portfolio bad rate
  b2 <- b0
  
  # finding the cutoff representation on ROC curve 2 for scenario 2
  phi2<-function(x){B*(1-y1(x))/((1-B)*(1-x)+B*(1-y1(x)))-b2}
  x2 <- as.numeric(uniroot(phi2,lower=0.001,upper=.999,tol = .Machine$double.eps)[1])
  
  # new approval rate
  a2 <- ((1-B)*(1-x2)+B*(1-y1(x2)))
  
  # as above
  # deriv2 <- numDeriv::grad(y1, x2)
  deriv2 <- FuncBiProperDeriv(x2, c2, da2)
  mbr2 <- (1+(1-B)/B/deriv2)^(-1)
  ir2 <- ir0
  profit2 <- a2*(ir2*(1-b2)-b2)
  ginip2 <- GiniP(y1, x2)
  
  a2_change <- a2/a0-1
  b2_change <- b2/b0-1
  profit2_change <- profit2/profit0-1
  
  ## 3. Profit increase scenario (keep marginal bad rate)
  
  # same marginal bad rate - finding the cutoff representation on ROC curve 2 for scenario 3 
  #numDeriv::grad(y1, x)
  # phi3<-function(x){numDeriv::grad(y1, x)-deriv0}
  phi3<-function(x){FuncBiProperDeriv(x, c2, da2)-deriv0}
  x3 <- as.numeric(uniroot(phi3,lower=0+1/1000,upper=1-1/1000,tol = .Machine$double.eps)[1])
  
  # new approval rate and portfolio bad rate in scenario 3
  a3 <- ((1-B)*(1-x3)+B*(1-y1(x3)))
  b3 <- B*(1-y1(x3))/a3
  
  # as above
  #deriv3 <- numDeriv::grad(y1, x3)
  deriv3 <- FuncBiProperDeriv(x3, c2, da2)
  mbr3 <- (1+(1-B)/B/deriv3)^(-1)
  ir3 <- ir0
  profit3 <- a3*(ir3*(1-b3)-b3)
  ginip3 <- GiniP(y1,x3)
  
  a3_change <- a3/a0-1
  b3_change <- b3/b0-1
  profit3_change <- profit3/profit0-1

  
  
  # f3<-function(x){format(round(100*x,3), nsmall=3)}
  # fperc<-function(x){paste(f3(x),"%", sep="")}
  
  return(list(
  inputs = list(GINI1=GINI1, GINI2=GINI2, B=B, a0=a0, c1=c1, c2=c2)
  , curves = list(y0=y0, y1=y1)
  , cutoffs = data.frame(xcord = c(x0, x1, x2, x3), 
                   ycord = c(y0(x0), y1(x1), y1(x2), y1(x3)))
  , res1=data.frame(no=1:3, scenario=c('keep approval (bad rate reduction)', 'keep bad rate (approval rate improvement)', 'keep marginal bad rate (profit increase)'), 
                   approval_change=c(a1_change, a2_change, a3_change),
                   bad_rate_change=c(b1_change, b2_change, b3_change),
                   profit_change=c(profit1_change, profit2_change, profit3_change))
  , res2=data.frame(scenario=1:3, 
                   existing_approval = rep(a0,3), 
                   new_approval = c(a1, a2, a3),
                   existing_bad_rate = rep(b0,3),
                   new_bad_rate = c(b1, b2, b3),
                   existing_profit = rep(profit0, 3),
                   new_profit = c(profit1, profit2, profit3))
  , res3=data.frame(scenario=1:3, 
                   existing_portfolio_gini = rep(ginip0,3),
                   new_portfolio_gini = c(ginip1, ginip2, ginip3),
                   existing_marginal_bad_rate = rep(mbr0,3),
                   new_marginal_bad_rate = c(mbr1, mbr2, mbr3))))
}

# Plot function  
plot_better_gini <- function(res){
df <- data.frame(xcord=res$cutoffs$xcord, ycord=res$cutoffs$ycord, point_id=c("a", "b", "c", "d"))
df$labels <- c('', paste("[1] Bad rate change:\n", ifelse(res$res1$bad_rate_change[1]>0,"+",""), round(100*res$res1$bad_rate_change[1], 3), "%", sep = ""),
               paste("[2] Approval change:\n", ifelse(res$res1$approval_change[2]>0,"+",""), round(100*res$res1$approval_change[2], 3), "%", sep = ""),
               paste("[3] Profit change:\n", ifelse(res$res1$profit_change[3]>0,"+",""), round(100*res$res1$profit_change[3], 3), "%", sep = ""))

ggplot(data.frame(x=c(0,1)), aes(x)) +
  stat_function(fun = res$curves$y0, geom = "line", aes(colour = "y0"), lwd = 1.1) +
  stat_function(fun = res$curves$y1, geom = "line", aes(colour = "y1"), linetype= "dashed", lwd = 1.1) +
  scale_colour_manual("ROC Curve", values = c("dodgerblue2", "olivedrab"), breaks = waiver(), labels=c("Scorecard 1", "Scorecard 2")) +
  geom_point(data = df[1,], aes(xcord, ycord, shape = factor(point_id)),  colour = "black", fill = "dodgerblue2", size = 4, stroke = 1.5) +
  geom_point(data = df[c(2,3, 4),], aes(xcord, ycord, shape = factor(point_id)), colour = "black", fill = "olivedrab", size = 4, stroke = 1.5) +
  scale_shape_manual("Cut-off\npoint",
                     values = c('a' = 21, 'b' = 22, 'c' = 24, 'd' = 23),
                     labels = c("Currently", "[1] Bad rate reduction\n scenario (keep approval)", "[2] Approval rate\n improvement scenario\n (keep bad rate)", "[3] Profit increase\n scenario (keep marginal\n bad rate)")) +
  xlab("cumulative good proportion") + ylab("cumulative bad proportion") + ggtitle("") + 
  theme(legend.position="bottom", 
        legend.box = "vertical", 
        axis.text=element_text(size=14), 
        axis.title=element_text(size=16,face="bold"),
        legend.text=element_text(size=11)) +
  guides(color = guide_legend(order = 1)) +
  coord_fixed() +
  geom_text_repel(data = df[2:4,],
                  mapping = aes(x=xcord, y=ycord, label = labels),
                  size=5,
                  min.segment.length = 0,
                  hjust = 1, vjust = -1, show.legend = FALSE)
}

# Example
(resP<-better_gini_calculatorP(.45, .65, .1, .6))
plot_better_gini(resP)
 
(res<-better_gini_calculator(.45, .65, .1, .6))
plot_better_gini(res)
# 
resP$res1
res$res1

#CHECK ##########
# FuncBiNormal(0:10/10, .7, 1)-FuncBiProper(0:10/10, 0, Gini2Da(.7, 0))
# 
numDeriv::grad(function(x){FuncBiNormal(x, .7, 1)}, .5)
FuncBiProperDeriv(.5, 0, Gini2Da(.7, 0))
# numDeriv::grad(function(x){FuncBiProper(x, 0, Gini2Da(.7, 0))}, .5)
# 
# ProperBiNormalDeriv1 <- Vectorize(function(x, c1, d_a1){
#   s = v_c_from_tpf(function(x_){tpf(x_,c1,d_a1)},x)
#   res = ((1-c1)*(dnorm(-(1-c1)*s-d_a1/2*sqrt(1+c1^2)) + dnorm(-(1-c1)*s+d_a1/2/c1*sqrt(1+c1^2))))/((1+c1)*(dnorm(-(1+c1)*s+d_a1/2*sqrt(1+c1^2)) + dnorm(-(1+c1)*s+d_a1/2/c1*sqrt(1+c1^2)) ))
#   return(res)
# })
# 
# x0 <- .5
# numDeriv::grad(function(x){FuncBiNormal(x, .7, 1)}, x0)
# #numDeriv::grad(function(x){FuncBiProper(x, 0, Gini2Da(.7, 0))}, x0)
# ProperBiNormalDeriv1(x0, 0, Gini2Da(.7,0))
# 
