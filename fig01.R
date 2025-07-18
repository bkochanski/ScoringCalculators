#res40000<-read.csv("ScoringCalculators/simresults4000020220428233229.csv")
#res10000<-read.csv("ScoringCalculators/simresults1000020220428223014.csv")
#res1000<-read.csv("ScoringCalculators/simresults100020220428212247.csv")
#res100000<-read.csv("ScoringCalculators/simresults10000020220429000403.csv")

res40000<-read.csv("simresults4000020220428233229.csv")
res10000<-read.csv("simresults1000020220428223014.csv")
res1000<-read.csv("simresults100020220428212247.csv")
res100000<-read.csv("simresults10000020220429000403.csv")

# sum(!is.na(res1000$diff.probit))
# sum(!is.na(res10000$diff.probit))
# sum(!is.na(res100000$diff.probit))

# boxplot(list(res1000$diff.logistic[!is.na(res1000$diff.logistic)],
#              res10000$diff.logistic[!is.na(res10000$diff.logistic)],
#              res100000$diff.logistic[!is.na(res100000$diff.logistic)]),
#         horizontal=TRUE
# )

library(ggplot2)

combined_data <- data.frame(
  Value = c(res1000$diff.logistic[!is.na(res1000$diff.logistic)],
            res10000$diff.logistic[!is.na(res10000$diff.logistic)],
            res100000$diff.logistic[!is.na(res100000$diff.logistic)]),
  Group = rep(c("n = 1000", "n = 10000", "n = 100000"), 
              c(
                length(res1000$diff.logistic[!is.na(res1000$diff.logistic)]),
                length(res10000$diff.logistic[!is.na(res10000$diff.logistic)]),
                length(res100000$diff.logistic[!is.na(res100000$diff.logistic)])
              ))
)

# Create a horizontal boxplot
library(dplyr)
p <- ggplot(combined_data, aes(x = Value, y = Group)) +
  geom_boxplot() +
  labs(x = "", y = "") + theme_bw()

print(p)

# Save as EPS
ggsave("fig01.eps", plot = p, device = "eps")
