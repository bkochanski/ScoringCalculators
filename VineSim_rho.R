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

gini_combine_calculator(gini1, gini2, rhos2, badrate)

library(splines)
mylogs1 <- glm(I(1-bad) ~ ns(s1, df=5) + ns(s2, df=5), data = simdf, family = binomial)
(gini_logs1<-my_gini(mylogs1$y, -mylogs1$fitted.values))

mylogs2 <- glm(I(1-bad) ~ ns(s1n, df=5) + ns(s2n, df=5), data = simdf, family = binomial)
(gini_logs2<-my_gini(mylogs2$y, -mylogs2$fitted.values))


a<-ns(simdf$s1n, df=5)
plot(a[,1])

# Obliczyć maksymalne Gini z modelu

