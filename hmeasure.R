#install.packages("hmeasure")
library(hmeasure)

my_gini<-function(resp, pred){
  c<-pred[order(pred)]
  d<-resp[order(pred)]
  bc<-c(0, cumsum(d)/sum(d))
  gc<-c(0, cumsum(1-d)/sum(1-d))
  sum((gc[2:(length(gc))]-gc[1:(length(gc)-1)])*(bc[2:(length(bc))]+bc[1:(length(bc)-1)]))-1} 


gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf, subdivisions=500, rel.tol = .Machine$double.eps^.24)$value/defrate/(1-defrate)-1
}

r_from_gini<-function(gini, defaultrate=.1){
  phi_s1<-function(x){gini_from_r(rho=x, defrate=defaultrate)-gini}
  uniroot(phi_s1,lower=0,upper=1,tol = .Machine$double.eps)$root
}

m <- 3
n <- 5000
gini1<-.4
gini2<-.6
corrs<-.4
dfrate<-.1
sigma <- matrix(c(1, corrs, r_from_gini(gini1, dfrate),
                  corrs, 1, r_from_gini(gini2, dfrate),
                  r_from_gini(gini1, dfrate), r_from_gini(gini2, dfrate), 1), 
                nrow=3)
z <- MASS::mvrnorm(n,mu=rep(0, m),Sigma=sigma)
df<-z
df[,3]<-(z[,3]<qnorm(dfrate))*1
df<-data.frame(df)
names(df)<-c('s1', 's2', 'default_flag')

my_gini(df$default_flag, df$s1)
hmeasure::HMeasure(1-df$default_flag, df$s1)$metrics


#generate scores and defaults from a ROC curve - is it possible?
# złych losujemy z unif (0,1)
# dobrych losujemy z unif(0,1) i przekształcamy przez krzywą

# FUNCTIONS:
# FuncMidNormal<-function(x,g){pnorm(qnorm((g+1)/2)*sqrt(2)+qnorm(x))}
# FuncMidNormal<-Vectorize(FuncMidNormal)
# 
# FuncMidFractal<-function(x,g){0.5*(1-(1-x)^((1+g)/(1-g)))+0.5*x^((1-g)/(1+g))}
# 

FuncBiFractal<-function(x,g, beta){beta*(1-(1-x)^((1+g)/(1-g)))+(1-beta)*x^((1-g)/(1+g))}
FuncBiFractal<-Vectorize(FuncBiFractal)

FuncBiNormal<-function(x,g,shape){pnorm(qnorm((g+1)/2)*sqrt(1+shape^2)+shape*qnorm(x))}
FuncBiNormal<-Vectorize(FuncBiNormal)

plot(1:99/100, FuncBiNormal(1:99/100, .5, 1/1.4))

defrate1<-0.1
size<-100000
badscores<-runif(defrate1*size)
goodscores<-FuncBiNormal(runif(size-length(badscores)), .5, 1.4)
df2<-data.frame(score=c(badscores, goodscores), bad=c(rep(1, length(badscores)), rep(0, length(goodscores))))

my_gini(df2$bad, df2$score)
hmeasure::HMeasure(1-df2$bad, df2$score)$metrics[c("H")]


scorecards<-1000
fginiin<-runif(scorecards, .2, .9)
fshapein<-runif(scorecards, .7, 1.4)
fbetain<-runif(scorecards, 0, 1)
fdefratein<-rep(.1, scorecards)
#fdefratein<-runif(scorecards, .01, .30)
fginiout<-c()
fh<-c()
size<-1000

for (i in 1:length(fginiin)) {
  print(i)
  defrate<-fdefratein[i]
  badscores<-runif(defrate*size)
  # goodscores<-FuncBiNormal(runif(size-length(badscores)), fginiin[i], fshapein[i])
  goodscores<-FuncBiFractal(runif(size-length(badscores)), fginiin[i], fbetain[i])
  df2<-data.frame(score=c(badscores, goodscores), bad=c(rep(1, length(badscores)), rep(0, length(goodscores))))
  fginiout[i]<-my_gini(df2$bad, df2$score)
  fh[i]<-as.numeric(hmeasure::HMeasure(1-df2$bad, df2$score)$metrics[c("H")])
}


library(ggplot2)
ggplot(data.frame(fbetain, fginiin, fginiout, fshapein, fh), aes(fginiin,fbetain, col=fh)) +
geom_point(size=4) +
  scale_colour_gradientn(colours = terrain.colors(10))
