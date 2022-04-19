gini_crd<-function(rho=0.5, defrate=0.1, gran=10000) {
  #function translating rho into gini for a given defrate
  drates_i<-pnorm(qnorm(defrate), rho*(qnorm(1:(gran-1)/gran)),sqrt(1-rho^2))
  drates_2i<-(c(1, drates_i)+c(drates_i, 0))/2
  return(ginic(cumsum((drates_2i)/sum(drates_2i)), cumsum((1-drates_2i)/sum(1-drates_2i))))
}

gini_crd()

F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}

a<-c(0,F1(0, .1, .4))
b<-c(0,F2(0, .1, .4))

test<-function(x){F1(x, .1, .5)*F2(x, .1, .5)}
test(5)

gini_from_r<-function(rho=.5, defrate=.1){
  
  integrate(function(x){F1(x, defrate, rho)*F2(x, defrate, rho)}, 
            lower=-.5, upper=5)$value
}

gini_crd()
gini_from_r(.5, .1)

F2(-Inf,.1,.5)

length(Inf)==1

integrate(function(x){F1(x, .1, .5)*F2(x, .1, .5)}, lower=-.5, upper=.5)
x=-.5
F1(x, .1, .5)
F2(x, .1, .5)
integrate(function(x){F2(x,.1,.5)}, lower=-Inf, upper=Inf)$value
integrate(function(x){F1(x,.1,.5)}, lower=-Inf, upper=Inf)$value

F1(.4,.4,.4)
F1(.400000001,.4,.4)-F1(.4,.4,.4)

ftest1<-function(x){F1(x,.4,.4)}
ftest1(.4)
ftest1<-Vectorize(ftest1)
curve(ftest1, from=-1, to=1)

ftest2<-Vectorize(function(x){F1(x,.4,.4)})
curve(ftest2, from=-1, to=1)

 
gini_from_r<-function(rho=.5, defrate=.1){
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
            lower=-Inf, upper=Inf)$value/defrate/(1-defrate)-1
}

gini_crd(.4,.1)
gini_crd(.4,.1, 1000)
gini_crd(.4,.1, 10000)
gini_crd(.4,.1, 100000)
gini_crd(.4,.1, 1000000)
gini_crd(.4,.1, 10000000)
gini_crd(.4,.1, 100000000)

gini_from_r(.4,.1)

F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}

#======================== Final:

gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf)$value/defrate/(1-defrate)-1
}

gini_from_r(.4,.1)



# Error in integrate(function(x) { : maximum number of subdivisions reached 
gini_combine_calculator(my_gini(df$default_flag, df$s1), my_gini(df$default_flag, 
                                                                 df$s2), cor(df$s2, df$s1), mean(df$default_flag))

gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf)$value/defrate/(1-defrate)-1
}

gini_combine_calculator(0.3241679, 0.292137, 0.01201619, 0.04160888)
#?integrate

gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf, subdivisions=200)$value/defrate/(1-defrate)-1
}

gini_combine_calculator(0.3241679, 0.292137, 0.01201619, 0.04160888)


gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf, subdivisions=500)$value/defrate/(1-defrate)-1
}

r_from_gini(.284, .1677)

gini_from_r<-function(rho=.5, defrate=.1){
  F1<-function(s, d, rho){integrate(function(x){pnorm((qnorm(d)-rho*x)/sqrt(1-rho^2))*dnorm(x)}, lower=-Inf, upper=s)$value}
  F2<-function(s,d,rho){pnorm((-qnorm(d)+rho*s)/sqrt(1-rho^2))*dnorm(s)}
  F1_<-Vectorize(function(x){F1(x, defrate, rho)})
  F2_<-Vectorize(function(x){F2(x, defrate, rho)})
  2*integrate(function(x){F1_(x)*F2_(x)}, 
              lower=-Inf, upper=Inf, subdivisions=500, rel.tol = .Machine$double.eps^.25)$value/defrate/(1-defrate)-1
}

gini_crd_test<-function(rho=0.5, defrate=0.1, gran=10000) {
  #function translating rho into gini for a given defrate
  drates_i<-pnorm(qnorm(defrate), rho*(qnorm(1:(gran-1)/gran)),sqrt(1-rho^2))
  drates_2i<-(c(1, drates_i)+c(drates_i, 0))/2
  return(ginic(cumsum((drates_2i)/sum(drates_2i)), cumsum((1-drates_2i)/sum(1-drates_2i))))
}

gini_from_r()
gini_crd_test(gran=100000000)

gini_combine_calculator(0.4517029, 0.450631, 0.3589866, 0.2633053) 

