cors<-c((0:9)/10, .98)
source('Gini_combined_master.R')
ginic4<-function(x){gini_combine_calculator(.6,.6,x,.1)[1]}
ginic4<-Vectorize(ginic4)
ginis<-ginic4(cors)

library(ggplot2)
p2 <- ggplot(data.frame(cors, ginis), aes(x=cors, y=ginis)) + 
  geom_point() +
  #  geom_text_repel(label=paste0('r=', round(cors,2), ' Gini=', round(ginis,3)), hjust=0.4, vjust=1) +
  xlab("Correlation") + ylab("Gini of the combinded scorecard") + ggtitle("") + ylim(c(.6, .85)) + theme_bw()

plot(p2)

# Save as EPS
ggsave("fig02.eps", plot = p2, device = "eps", width = 7, height = 5)
