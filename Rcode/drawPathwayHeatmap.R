# This code Draws heat map of pathway 
# enrichment scores of 23 datasets
##################################
ids = cutree(hc, k=2)
ids.sepsis = colnames(gsea.list$Zobs)[substr(colnames(gsea.list$Zobs), 1, 3)=="GSE"]
ids.sepsis.like = setdiff(names(which(ids==1)), ids.sepsis)
ids.cancer = names(which(ids==2))
ids = c(ids.sepsis, ids.sepsis.like, ids.cancer)
sidecols = c(rep("purple", length(ids.sepsis)), rep("cyan", length(ids.sepsis.like)), 
             rep("blue", length(ids.cancer)))

################################################
# Find the pathways that are significant in 80% 
################################################
thresh <- 0.8
pmatrix = gsea.list$pvals[sel.pathways,]
pathways.sl = names(which(apply(pmatrix[,ids.sepsis.like], 1, function(x) sum(x<0.01)>= thresh*length(x))))
pathways.c = names(which(apply(pmatrix[,ids.cancer], 1, function(x) sum(x<0.01)>= thresh*length(x))))
pathways.sl.only = setdiff(pathways.sl, pathways.c)
pathways.c.only = setdiff(pathways.c, pathways.sl)
pathways.both = intersect(pathways.sl, pathways.c)
sel = c(pathways.c.only, pathways.both, pathways.sl.only)
keggnames = as.character(pathways.list[paste0("path:",sel)])
keggnames = as.character(sapply(keggnames, function(x) strsplit(x, split=" - Homo")[[1]][1]))

# changes the names to TCGA Codes
plotdat = as.matrix(sstcga[sel,ids])
heatmap.2(plotdat, trace="none", 
          col= colorRampPalette(c("darkgreen","green4","white","red","darkred"))(79), 
          sepcolor="white",
          #RowSideColors = rowSideCols,
          #cellnote=cellnotemat, notecex=0.7,
          #notecol="black",
          sepwidth=c(0.5, 0.5), 
          Rowv = FALSE, dendrogram = "column", key = T, 
          margins = c(6, 19), ColSideColors = sidecols, 
          labRow=keggnames, cexCol = 0.8, cexRow = 0.6)
legend("topright", 
       legend = paste( c("CA", "SLC", "SS")), pch=15,
       col = c("blue", "cyan", "purple"),
         cex=1, bty="n")
rm(thresh, keggnames, pathways.c, pathways.sl, plotdat, pathways.both, 
   pathways.c.only, pathways.sl.only, pmatrix, sidecols, ids)

