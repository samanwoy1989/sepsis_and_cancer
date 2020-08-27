# get pathway gene lfc heatmap
#drawPathwayGeneHeatmap.R
# Input: KEGG id
# Output: data.frame of genes by studies lfc, heatmap (side effect)

drawPathwayGeneHeatmap = function(kegg="hsa05032", sort.by.var=TRUE, groupid="", ...) {
  # Get the common genes of the pathway present in all studies
  egs.tc = names(which(table(unlist(sapply(tc.eset, featureNames)))==length(tc.eset)))
  egs.ss = names(which(table(unlist(sapply(ss.eset, featureNames)))==length(ss.eset)))
  egs = intersect(egs.tc, egs.ss)
  pathegs = intersect(egs, genes.by.pathway[[kegg]])

  # get the log-fold change for each gene in the pathway
  lfc.tc = sapply(tc.eset, function(eset) {rowttests(eset[pathegs,], "Target")[,"dm"]})
  lfc.ss = sapply(ss.eset, function(eset) {rowttests(eset[pathegs,], "Group")[,"dm"]})
  lfc = cbind(lfc.tc, lfc.ss)
  rownames(lfc) = pathegs
  # Arrange the data according to disease groups
  ids = c(ids.cancer, ids.sepsis.like, ids.sepsis)
  lfc = lfc[, ids]
  sidecols = c(rep("blue", length(ids.cancer)), 
               rep("cyan", length(ids.sepsis.like)), 
               rep("purple", length(ids.sepsis)))
  gsyms = links(org.Hs.egSYMBOL[rownames(lfc)])[,"symbol"]
  gnames = links(org.Hs.egGENENAME[rownames(lfc)])[,"gene_name"]
  rownames(lfc) = apply(data.frame(gsyms, gnames), 1, paste, collapse=" | ")
  
  # get median gene expression across disease groups
  lfc.ss = rowMedians(lfc[,ids.sepsis])
  lfc.slc = rowMedians(lfc[,ids.sepsis.like])
  lfc.ca = rowMedians(lfc[,ids.cancer])
  
  # Format data for plotting
  lfc1 = as.matrix(data.frame(lfc.ca, lfc.slc, lfc.ss))
  rownames(lfc1) = rownames(lfc)
  colnames(lfc1) = c("CA","SLC","SS")
  lfc=lfc1
  rm(lfc1)
  # keep those genes with same direction of transcriptional change
  # between sepsis and SLC
  lfc = lfc[sign(lfc[,"SS"])==sign(lfc[,"SLC"]),]
  
  # sort the genes
  o = order(lfc[,"CA"]-lfc[,"SLC"])
  lfc = lfc[o,]
  
  titlestr = paste(strsplit(pathways.list[paste("path:",kegg,sep="")], split=" - Homo")[[1]][1], "\n", kegg, sep="")

  par(mar=c(6,4,3,2))
  #par(mar = rep(2, 4))
  heatmap.2(lfc, trace="none",
            col= colorRampPalette(c("darkgreen","green4","white","red","darkred"))(79), 
            sepcolor="white",
            sepwidth=c(0.5, 0.5), dendrogram = "none", 
            Colv = FALSE, Rowv=FALSE, key = T, keysize = 1,
            margins = c(10, 25), 
            cexCol=1.5,
            ...)
  title(main=titlestr, sub=groupid, col.sub="blue", cex.sub=1.5, outer=F)
  return(lfc)
}