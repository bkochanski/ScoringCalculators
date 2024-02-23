res40000<-read.csv("ScoringCalculators/simresults4000020220428233229.csv")
res10000<-read.csv("ScoringCalculators/simresults1000020220428223014.csv")
res1000<-read.csv("ScoringCalculators/simresults100020220428212247.csv")
res100000<-read.csv("ScoringCalculators/simresults10000020220429000403.csv")

boxplot(list(res1000$diff.logistic[!is.na(res1000$diff.logistic)],
             res10000$diff.logistic[!is.na(res10000$diff.logistic)],
             res100000$diff.logistic[!is.na(res100000$diff.logistic)]),
        horizontal=TRUE
)
