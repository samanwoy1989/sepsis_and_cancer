# getData.R

################################ 
# For each pathway run permutation testing for pathway enrichment
source("Rcode/runSSgsea.R")

##################################################################################################################################################################
# Find out the significant paths in all 6 SS datasets
##################################################################################################################################################################
#ss.gsea.list$pvals =  apply(ss.gsea.list$pvals, 2, function(x) p.adjust(x, "fdr"))
is.significant = apply(ss.gsea.list$pvals, 1, function(x) all(x<sig.threshold))
is.same.direction = apply(ss.gsea.list$Zobs, 1, function(x) abs(sum(sign(x)))==length(x))
sel.pathways = names(which(is.significant & is.same.direction))
rm(is.significant, is.same.direction)

######################################################################################################################
# From here the TCGA data analysis starts
######################################################################################################################
################################
# Performs Pathway GSEA analysis for 17 TCGA DATA 
source("Rcode/runTCGAgsea.R")

######################################################################################################################
# Merge SS and TCGA pathway result
# zobs (Enrichment score) and p-value are retained
######################################################################################################################
gsea.list = list("Zobs"=data.frame(ss.gsea.list$Zobs, tcga.gsea.list$Zobs),
                 "pvals"=data.frame(ss.gsea.list$pvals, tcga.gsea.list$pvals))

colnames(gsea.list$Zobs)[match(params1$Tissue, colnames(gsea.list$Zobs))] = as.character(params1$Code)
colnames(gsea.list$pvals)[match(params1$Tissue, colnames(gsea.list$pvals))] = as.character(params1$Code)

rm(ss.gsea.list, tcga.gsea.list)
