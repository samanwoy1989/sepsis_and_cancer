# Check how the direction of these pathways 
# look in SS survivors 
############################################
library(ssgeosurv)
data("ss.list")
data("ss.surv.list")

result.df <- data.frame()
pathids = sel

for (j in 1:length(ss.surv.list)){
  for (i in 1: length(pathids)) {
    path <- pathids[i]
    which.surv <- which(ss.surv.list[[j]]$Outcome=="survivor")
    which.nonsurv <- which(ss.surv.list[[j]]$Outcome=="nonsurvivor")
    path.eset <- ss.surv.list[[j]][featureNames(ss.surv.list[[j]])%in%genes.by.pathway[[path]],]
    rttStat = rowttests(path.eset, "Outcome")$statistic
    pathScore = sum(rttStat)/sqrt(length(rttStat))
    survScore = (-1)*pathScore # as rowttests takes alphabetically, i.e., nonsurv - surv
    result.df[i,j] = survScore
    }
}
rownames(result.df) = pathids
colnames(result.df) = names(ss.surv.list)

is.up.in.surv = apply(result.df, 1, function(x) sum(x>0)>0.5*length(x))
print.noquote(paste0(sum(is.up.in.surv), " out of ", length(pathids), 
      " are up-regulated in survivors compared to non-survivors in septic shock."))

