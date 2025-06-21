source('Gini_combined_master.R')

#data<-readxl::read_excel('C:/Users/blaze/Downloads/data.xlsx')
#data<-readxl::read_excel('C:/Users/Błażej/Downloads/data.xlsx')
data<-readxl::read_excel('/Users/blakocha/Downloads/data.xlsx')

data$Region <- c("Mexico", "Brazil", "Ukraine", "Ecuador")

gini_from_auc<-function(x){(x-.5)*2}
gini_from_auc<-Vectorize(gini_from_auc)

auc_from_gini<-function(x){x/2+.5}
auc_from_gini<-Vectorize(auc_from_gini)

#gini_from_auc(.7)
#auc_from_gini(.4)

data$Gini_Bureau<-gini_from_auc(data$AUC_Bureau)
data$Gini_Psych<-gini_from_auc(data$AUC_Psych)
data$Gini_Combined<-gini_from_auc(data$AUC_Combined)

#View(data)

gini_predict<-function(g1, g2, corr, defaultrate){unname(gini_combine_calculator(g1, g2, corr, defaultrate)[1])}
gini_predict<-Vectorize(gini_predict)

data$Gini_Predict_r<-gini_predict(data$Gini_Bureau, data$Gini_Psych, data$r_Bureau_Psych, data$Bad_Rate)
data$Gini_Predict_rho<-gini_predict(data$Gini_Bureau, data$Gini_Psych, data$rho_Bureau_Psych, data$Bad_Rate)

data$AUC_predict_r<-auc_from_gini(data$Gini_Predict_r)
data$AUC_predict_rho<-auc_from_gini(data$Gini_Predict_rho)

AUC_lwr<-function(auc, brate, n){presize::prec_auc(auc, brate, n)$lwr}
AUC_lwr<-Vectorize(AUC_lwr)
AUC_upr<-function(auc, brate, n){presize::prec_auc(auc, brate, n)$upr}
AUC_upr<-Vectorize(AUC_upr)

data$AUC_predict_r_lwr<-AUC_lwr(data$AUC_predict_r, data$Bad_Rate, data$N)
data$AUC_predict_r_upr<-AUC_upr(data$AUC_predict_r, data$Bad_Rate, data$N)

data$Gini_predict_r_lwrc<-gini_from_auc(data$AUC_predict_r_lwr)
data$Gini_predict_r_uprc<-gini_from_auc(data$AUC_predict_r_upr)

data$Gini_predict_r_lwr <- data$Gini_Predict_r - (data$Gini_Predict_r - data$Gini_predict_r_lwrc)/1.96
data$Gini_predict_r_upr <- data$Gini_Predict_r + (data$Gini_Predict_r - data$Gini_predict_r_lwrc)/1.96

#presize::prec_auc(auc_from_gini(.376), .1677, 1306)$lwr

#data$AUC_predict_rlwr<-auc_from_gini(data$Gini_Predict_r)


gini_from_auc(presize::prec_auc(auc_from_gini(.376), .1677, 1306)$lwr)
gini_from_auc(presize::prec_auc(auc_from_gini(.376), .1677, 1306)$upr)


library(ggplot2)
p5 <- ggplot(data, aes(y=paste0(Sample, ': ', Region, '\nN=', N, '\nbad rate=', Bad_Rate*100, '%\ncorr=', r_Bureau_Psych))) +
  geom_point(aes(x=Gini_Bureau), col='gray28', shape=5, size=3) +
  geom_text(aes(x=Gini_Bureau, label=paste('Bureau:', Gini_Bureau)), hjust=0, vjust=1.5, col='gray28', angle = -20) +
  geom_point(aes(x=Gini_Psych), col='gray28', shape=9, size=3)+
  geom_text(aes(x=Gini_Psych, label=paste('Psych:', Gini_Psych)), hjust=0, vjust=-1, col='gray28', angle = 20) +
  geom_point(aes(x=Gini_Combined), col='black', shape=2, size=3)+
  geom_text(aes(x=Gini_Combined, label=paste('Log. reg.:', Gini_Combined)), hjust=0, vjust=-1, angle = 20) +
  geom_point(aes(x=Gini_Predict_r), col='black', shape=17, size=3)+
  geom_errorbar(aes(xmin = Gini_predict_r_lwr, xmax = Gini_predict_r_upr), col='black', linetype=5, width=.1)+
  geom_text(aes(x=Gini_Predict_r, label=paste('MVN calc.:', round(Gini_Predict_r,3))), hjust=0, vjust=1.5, col='black', angle=-20) +
  #  geom_point(aes(x=Gini_Predict_rho), col='dark blue')+
  #  geom_text(aes(x=Gini_Predict_rho, label=paste('Predicted (rho):', round(Gini_Predict_rho,3))), hjust=0, vjust=1.5, col='dark blue') +
  scale_y_discrete(limits=rev)+
  xlab('')+ylab('')+xlim(c(0.2, .75)) +
  theme_bw() + theme(text = element_text(size = 15))

print(p5)

# Save as EPS
ggsave("fig05.eps", plot = p5, device = "eps", width = 7, height = 5)
