# define R-vine tree structure matrix
Matrix <- c(3, 2, 1,
            0, 2, 1,
            0, 0, 1)
Matrix <- matrix(Matrix, 3, 3)
# define R-vine pair-copula family matrix
family <- c(0, 1, 3,
            0, 0, 10,
            0, 0, 0)
family <- matrix(family, 3, 3)
# define R-vine pair-copula parameter matrix

par <- c(0, .5, 0.7,
         0, 0, 3.32,
         0, 0, 0)
par <- matrix(par, 3, 3)
# define second R-vine pair-copula parameter matrix
par2 <- c(0, 0, 0,
          0, 0, 0.88,
          0, 0, 0)
par2 <- matrix(par2, 3, 3)
# define RVineMatrix object
library(VineCopula)
RVM <- RVineMatrix(Matrix = Matrix, family = family,
                   par = par, par2 = par2,
                   names = c("V1", "V2", "V3"))

# plot(RVM,1)
# plot(RVM,2)
# contour(RVM,1)
# contour(RVM,2)

# simulate a sample of size 300 from the R-vine copula model
#set.seed(123)
sizesim<-100000
simdata <- RVineSim(sizesim, RVM)
#psych::pairs.panels(simdata)

simdf<-data.frame(simdata)
names(simdf)<-c("latentunif", "s1", "s2")
simdf$latent<-(qnorm(simdf$latentunif))
mean(simdf$latent)
sd(simdf$latent)
mean(simdf$latentunif<.1)
simdf$bad<-1*(simdf$latentunif<.1)

rocs1_x <- cumsum((1-simdf$bad[order(simdf$s1)]))/sum(1-simdf$bad)
rocs1_y <- cumsum((simdf$bad[order(simdf$s1)]))/sum(simdf$bad)
plot(rocs1_x, rocs1_y, type="l")
(ginis1 <- 2*bigstatsr::AUC(-simdf$s1, simdf$bad)-1)
(ginis2 <- 2*bigstatsr::AUC(-simdf$s2, simdf$bad)-1)
(corr_B_spearman <- cor(simdf$s1, simdf$s2, method="spearman"))
(badrate_ <- mean(simdf$bad))

rocs2_x <- cumsum((1-simdf$bad[order(simdf$s2)]))/sum(1-simdf$bad)
rocs2_y <- cumsum((simdf$bad[order(simdf$s2)]))/sum(simdf$bad)
points(rocs2_x, rocs2_y, type="l", col='blue')

source('Gini_combined_master.R')

(mvn_res<-gini_combine_calculator(ginis1, ginis2, corr_B_spearman, badrate_, method="spearman"))
simdf$smvn <- mvn_res["score_1_weight"]*simdf$s1 + mvn_res["score_2_weight"]*simdf$s2
rocs_mvn_x <- cumsum((1-simdf$bad[order(simdf$smvn)]))/sum(1-simdf$bad)
rocs_mvn_y <- cumsum((simdf$bad[order(simdf$smvn)]))/sum(simdf$bad)
points(rocs_mvn_x, rocs_mvn_y, type="l", col='red')
(ginismvn <- 2*bigstatsr::AUC(-simdf$smvn, simdf$bad)-1)
mvn_res["new_gini"]

library(splines)
flex_logistic <- glm(I(1-bad) ~ ns(qnorm(pobs(s1)), df=2) + ns(qnorm(pobs(s2)), df=2), data = simdf, family = binomial)
simdf$slog <- predict(flex_logistic)
rocs_log_x <- cumsum((1-simdf$bad[order(simdf$slog)]))/sum(1-simdf$bad)
rocs_log_y <- cumsum((simdf$bad[order(simdf$slog)]))/sum(simdf$bad)
points(rocs_log_x, rocs_log_y, type="l", col='green')
(ginismvn <- 2*bigstatsr::AUC(-simdf$slog, simdf$bad)-1)

rocsdf <- data.frame(rocs1_x, rocs1_y, rocs2_x, rocs2_y, rocs_mvn_x, rocs_mvn_y, rocs_log_x, rocs_log_y)

library(ggplot2)
plot1<-ggplot(rocsdf) +
  geom_line(aes(rocs1_x, rocs1_y), lty=1, col="grey", linewidth=1) + xlab("False positive rate") + ylab("False negative rate") +
  geom_line(aes(rocs2_x, rocs2_y), lty=3, col="grey", linewidth=1) +
  geom_line(aes(rocs_mvn_x, rocs_mvn_y), lty=1, linewidth=1) +
  geom_line(aes(rocs_log_x, rocs_log_y), lty=2)
plot1 + theme_bw()

