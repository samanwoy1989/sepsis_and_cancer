# This code performs analysis of sepsis and cancer datasets
# To accompany the manuscript BMC Cancer
# Samanwoy Mukhopadhyay and Saroj Kant Mohapatra
# Revised for re-submission to BMC Cancer
# Dec-Feb 2020

# clean up work area
rm(list=ls()); graphics.off()

###############################
# preliminaries
source("Rcode/prelims.R")

###############################
# get pathway score data 
source("Rcode/getData.R")
###############################

# There are 23 studies (17 cancers + 6 SS)
length(tc.eset) # 17 cancer studies
length(ss.eset) # 6 SS studies

# 687 cancer cases, 445 SS cases
sum(sapply(tc.eset, function(eset) ncol(eset[,eset$Target=="Case"])))
sum(sapply(ss.eset, function(eset) ncol(eset[,eset$Group=="Septic_Shock"])))

# 272 pathways, each pathway with > 10 genes
sum(sapply(genes.by.pathway, length)>10)

############################################################
# Hierarchical clustering of diseases using 90 pathways
# Fig S1: dendrogram_90_pathways
# Get the data for selected pathways
# Legend Test:
# Supplementary Figure S1\
# For the 90 pathways significantly associated with SS, pathway scores
# were computed for the cancer studies. Scores for SS (6 studies) and cancer
# (17 studies) were subjected to hierarchical clustering (details in main text).
# As shown in this dendrogram, there is clear segregation of the cancer studies
# into two groups: sepsis-like (SLC) and other cancer (CA).
############################################################
sstcga = gsea.list$Zobs[sel.pathways,]
d = dist(t(sstcga))
hc = hclust(d)
plot(hc, hang=-0.5, axes=F, xlab="", ylab="", main="Hierarchical clustering of diseases\n (with 90 pathways)")

# 90 pathways significantly enriched (i.e., deviated in the same direction)
# in each of the SS study
is.signif = rowSums(gsea.list$pvals[,1:6]<0.01)==6
is.same.direction = abs(rowSums(sign(gsea.list$Zobs[,1:6])))==6
sum(is.signif & is.same.direction) # 90 pathways

#############################################################
# Find the groups of cancer in the sepsis cluster
# Fig 2A: heatmap (of 66 pathways by all studies)
#############################################################
source("Rcode/drawPathwayHeatmap.R")
rm(hc, d, is.signif, is.same.direction)

# 66 pathways are selected from 90 pathways
# on basis of significant enrichment in at least 80% of studies 
# of one group of cancer (SLC or CA)
is.slc = (rowSums(gsea.list$pvals[sel.pathways,ids.sepsis.like]<0.01)>=0.8*6)
is.ca = (rowSums(gsea.list$pvals[sel.pathways,ids.cancer]<0.01)>=0.8*11)
sum(is.ca | is.slc)
# 66
rm(is.slc, is.ca)

###############################
# Fig 2B
# Draw a box plot of mean pathway score for the three disease groups
# for the following pathways
# (1) all 66, (2) CA-only, (3) SLC-only, (4) both CA and SLC 
###############################
# get the pathway scores for 66 pathways
# This shows up-regulation in SLC as (with greater magnitude) SS
# and down-regulation in CA
drawBoxplot3DiseaseGroup()

##########################################################################
# Average number of pathways up-regulated in disease group
# SS: 66, SLC: 59.7, CA: 35.4
##########################################################################
options(digits=4)
cat("Average number of pathways significantly dysregulated in a disease group:\n")
sapply(list("SS"=ids.sepsis, "SLC"=ids.sepsis.like, "CA"=ids.cancer), function(sids) {
  p1 = gsea.list$pvals[sel,sids]
  z1 =  gsea.list$Zobs[sel,sids]
  nPathsUp = mean(colSums(p1<0.01 & z1>0))
  #nPathsUp = apply(gsea.list$pvals[sel,]<0.01, 2, sum)[sids]
  return(nPathsUp)
})

############################################################################################################
# Network analysis
# Fig S3   Degree distribution and box plots
# Fig 3    Network diagram
############################################################################################################
source("Rcode/net.R")

# The network consists of 244 nodes and 5304 edges
length(V(gp)) # 244 vertices
length(E(gp)) # 5304 edges
rm(gp)

############################################################################################################
# Draw individual Pathway gene heatmaps
# Fig: S2
############################################################################################################
figs2 = "Results/Fig_S2.pdf"
if(!file.exists(figs2)) {
  #fcounter = 1
  p66ann = read.table(file="Metadata/pathways66annot.csv", header=TRUE, sep=",")
  rownames(p66ann) = p66ann$Pathway
  pathway.groups = as.character(sort(unique(p66ann[,3])))
  pdf(file=figs2, onefile=TRUE)
  op = par()
  par(mar=c(2,2,10,2))
  plot(1:2, type="n", xlab="", ylab="", axes=F,
           onefile=TRUE,
           main="Gene-level heatmap for 66 pathways",
           cex.main=2)
  par(op)
  for(grp in pathway.groups) {
    pathids = as.character(p66ann[grp==p66ann[,3],1])
    for(pathid in pathids) {
      print(pathid)
      drawPathwayGeneHeatmap(pathid, groupid=grp)
    }
  }
  dev.off()
  rm(grp, pathid, p66ann)
} else {
  print.noquote(paste0("Figure exists - ", figs2))
}
rm(figs2)
graphics.off()
#####################################################################
# Table 2
# Machine learning (SVM and NN) - based Classification of samples
# based on 66 pathway expression
#####################################################################
# Prepare validation data for machine learning
source("Rcode/prepareData4MLvalidation.R")

# 542 cancer cases and 180 control samples in validation
print(rowSums(sapply(vset.list, function(eset) table(eset$Group))))

###############################################
# function to get misclassification rate
getMiss = function(cmatrix) {
  100*sum(cmatrix[2]+cmatrix[3])/sum(cmatrix)
}

# Support Vector Machine: 5-fold cross-validation
library(MLInterfaces)
set.seed(1234)
xvsv = MLearn(label~., pathES, svmI, xvalSpec("LOG", 5, balKfold.xvspec(5)))
cmat = confuMat(xvsv)
print(cmat)
nCorrect = sum(cmat[1]+cmat[4])
nTotal = sum(cmat)
accu = 100*(nCorrect)/nTotal
cat(paste0(nCorrect, " out of ", nTotal, " (SVM Accuracy: ", formatC(accu),
           "%) correctly classified.\n"))
print(paste0("Misclassification rate for SVM: ", formatC(getMiss(cmat), 2), "%"))
rm(cmat, nCorrect, nTotal, accu)

# Neural Net: : 5-fold cross-validation
#nn1 = MLearn(label~., pathES, nnetI, which.train, size=3, decay=0.1)
load("Results/mlearn.models/nn_miss-3-3.rda")
cmat = confuMat(xvnn)
print(cmat)
nCorrect = sum(cmat[1]+cmat[4])
nTotal = sum(cmat)
accu = 100*(nCorrect)/nTotal
cat(paste0(nCorrect, " out of ", nTotal, " (NN Accuracy: ", formatC(accu),
           "%) correctly classified.\n"))
print(paste0("Misclassification rate for NN: ", formatC(getMiss(confuMat(xvnn)), 2), "%"))
rm(cmat, nCorrect, nTotal, accu)
rm(pathES, getMiss) # This needs to be done because the input of 

###############################################
# Preapare the TCGA data for surviaval analysis
###############################################
source("Rcode/prepData4SurvivalAnalysis.R")

# check the percentage of pathways significantly (p < 0.1)
# associated with survival
path66SLC = path66ByProj[, c("HNSC","LIHC","KIRC")]
nSignif = colSums(path66SLC<0.1)
nSignifPerc = 100*nSignif/66
dframe = data.frame(numPaths = nSignif, Total = rep(66,3), Percentage = nSignifPerc)
knitr::kable(dframe)

rm(path66SLC, nSignif, nSignifPerc, dframe)

# Figure 4: survival analysis and K-M plot
source("Rcode/drawSurvPlots.R")
gc(reset=TRUE)

#####################################################################################
# Uncomment the following line for checking direction of gene expression in survivors
# Warning: Takes a long time !!
source("Rcode/check_gexp_survivors_SS.R")
#####################################################################################

#####################################################################################
# Checking Viral percentage of TCGA samples from Different Reports
#####################################################################################
source("Rcode/checking_viral_percentage_TCGA_data.R")
