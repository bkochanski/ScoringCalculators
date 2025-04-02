# Wygenerować dane z Claytona i BB8
# Tutaj kod
summary(simdf)

#Znaleźć maksymalne Gini za pomocą modelu uczenia maszynowego XGBoost
# około 0.65

# Obliczyć maksymalne Gini z modelu MVN
(gini1=2*bigstatsr::AUC(-simdf$s1, simdf$bad)-1)
(gini2=2*bigstatsr::AUC(-simdf$s2, simdf$bad)-1)
(cors = cor(simdf$s1, simdf$s2))
cor(simdf$s1n, simdf$s2n)
(rhos = cor(simdf$s1, simdf$s2, method="spearman"))
(rhos2 = 2*sin(rhos*pi/6))
(badrate=mean(simdf$bad))

source('Gini_combined_master.R')

gini_combine_calculator(gini1, gini2, cors, badrate)
gini_combine_calculator(gini1, gini2, rhos, badrate)
gini_combine_calculator(gini1, gini2, rhos2, badrate)

library(splines)
mylogs1 <- glm(I(1-bad) ~ ns(s1, df=2) + ns(s2, df=2), data = simdf, family = binomial)
(gini_logs1<-my_gini(mylogs1$y, -mylogs1$fitted.values))

mylogs2 <- glm(I(1-bad) ~ ns(s1n, df=2) + ns(s2n, df=2), data = simdf, family = binomial)
(gini_logs2<-my_gini(mylogs2$y, -mylogs2$fitted.values))

mylogs3 <- glm(I(1-bad) ~ ns(qnorm(pobs(s1)), df=2) + ns(qnorm(pobs(s2)), df=2), data = simdf, family = binomial)
(gini_logs3<-my_gini(mylogs3$y, -mylogs3$fitted.values))


a<-ns(simdf$s1n, df=5)

mylogs1_list <- list()
gini_logs1_list <- list()

for (i in 1:5) {
mylogs1_list[[i]] <- glm(I(1-bad) ~ ns(s1, df=i) + ns(s2, df=i), data = simdf, family = binomial)
gini_logs1_list[i]<-my_gini(mylogs1_list[[i]]$y, -mylogs1_list[[i]]$fitted.values)
}

mylogs2_list <- list()
gini_logs2_list <- list()

for (i in 1:8) {
  print(i)
  mylogs2_list[[i]] <- glm(I(1-bad) ~ ns(s1n, df=i) + ns(s2n, df=i), data = simdf, family = binomial)
  gini_logs2_list[i]<-my_gini(mylogs2_list[[i]]$y, -mylogs2_list[[i]]$fitted.values)
}

plot(as.numeric(gini_logs1_list))
plot(as.numeric(gini_logs2_list))

# mylogs2 <- glm(I(1-bad) ~ ns(s1n, df=5) + ns(s2n, df=5), data = simdf, family = binomial)
# (gini_logs2<-my_gini(mylogs2$y, -mylogs2$fitted.values))

#normalizacja z danych
library(VineCopula)
hist(qnorm(pobs(simdf$s1)))

mylogs3_list <- list()
gini_logs3_list <- list()

for (i in 1:8) {
  print(i)
  mylogs3_list[[i]] <- glm(I(1-bad) ~ ns(qnorm(pobs(s1)), df=i) + ns(qnorm(pobs(s2)), df=i), data = simdf, family = binomial)
  gini_logs3_list[i]<-my_gini(mylogs3_list[[i]]$y, -mylogs3_list[[i]]$fitted.values)
}

plot(as.numeric(gini_logs3_list))

mylogs4_list <- list()
gini_logs4_list <- list()

for (i in 1:8) {
  print(i)
  mylogs4_list[[i]] <- glm(I(1-bad) ~ I(qnorm(pobs(s1))) + ns(qnorm(pobs(s1)), df=i) + I(qnorm(pobs(s1))) + ns(qnorm(pobs(s2)), df=i), data = simdf, family = binomial)
  gini_logs4_list[i]<-my_gini(mylogs4_list[[i]]$y, -mylogs4_list[[i]]$fitted.values)
}

plot(as.numeric(gini_logs4_list))

summary(mylogs3_list[[2]])

# Wizualizacja splines
nfirst<-2000
splinestmp <- ns(simdf$s1[1:nfirst], df=4)
#plot(splinestmp[,1], splinestmp[,2])
plot(simdf$s1[1:nfirst], splinestmp[,1], ylim=c(-.5,1), pch='.')
points(simdf$s1[1:nfirst], splinestmp[,2], col='blue', pch='.')
points(simdf$s1[1:nfirst], splinestmp[,3], col='red', pch='.')
points(simdf$s1[1:nfirst], splinestmp[,4], col='green', pch='.')


