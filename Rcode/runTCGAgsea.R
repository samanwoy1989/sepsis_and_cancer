#########################################################################
# This code runs GSEA Permutaion analysis to calculate Pathway score for
# 17 TCGA  datasets 
# Samanwoy Sat Mar  2 08:01:38 IST 2019
#########################################################################
fn <- "Results/pathway_score_TCGA_nperm_10k.rda"
if(file.exists(file= fn)) {
    cat("Loading TCGA GSEA file ...")
    load(file= fn )
    cat(" done!\n")
} else {
    start.time <- Sys.time()
    #########################################################################
    # Creating incidence Matrix for GSEA analysis
    #########################################################################
    # Create the incidence matrix
    pathway.names <- names(genes.by.pathway)
    pathway.genes <- unique(unlist(genes.by.pathway))
    numPathways <- length(pathway.names)
    numGenes <- length(pathway.genes)
    Am <- matrix(0, nrow=numPathways, ncol=numGenes)
    rownames(Am) <- pathway.names
    colnames(Am) <- pathway.genes
    for(i in 1:length(genes.by.pathway)) {
        Am[i,genes.by.pathway[[i]]] <- 1
    }
    # Reduce the incidence matrix by removing all gene sets that have 10 genes
    # or fewer
    selectedRows = (rowSums(Am)>10)
    inc.mat = Am[selectedRows, ]
    rm(i)
    ##########################################################
    tcga.gsea.res <- list()
    for ( i in 1: length(tc.eset)) {
        selg <- intersect(colnames(inc.mat), featureNames(tc.eset[[i]]))
        inc.mat <- inc.mat[,selg]
        gtable <- exprs(tc.eset[[i]])[selg,]
        set.seed(123)
        NPERM = 10000
        # Performing differntial expression
        fac <- pData(tc.eset[[i]])$Target
        rtt <- rowttests(as.matrix(gtable), factor(fac))
        rttStat <- rtt[["statistic"]]
        obs <- inc.mat %*% rttStat
        obs <- as.vector(obs)
        names(obs) <- rownames(inc.mat)
        ## performing pathway permutaion GSEA
        mylist <- list()
        drawPermutHist <- function(pathid=pathid) {
          zobs <- as.numeric(obs[pathid])
          permMat <- matrix(0, nrow = nrow(gtable), ncol = NPERM)
          j <- 1L
          while (j < (NPERM + 1)) {
            p1 <- sample(fac)
            permMat[, j] <- rowttests(as.matrix(gtable), p1, tstatOnly = TRUE)[["statistic"]]
            j <- j + 1L
        }
        ##permVec <- as.vector(inc.mat["hsa05202",] %*% permMat) ##^^original for null distribution
        permVec <- as.vector(inc.mat[pathid,] %*% permMat)
        
        #  hist(permVec, xlab="Simulated Pathway Score", ylim=c(0,300),
        #       xlim=range(c(permVec, zobs)),
        #       main="Pathway hsa05202 [Transcriptional misregulation in cancer - Human]")
        #  abline(v=zobs, lwd=3, col="red")
        #  text(x=c(30, 0), y=c(200,200), labels=c("Observed Score", "Null distribution"), 
        #       col=c("red","black"), cex=1.3 )
        permut_pval = ifelse(zobs<0,sum(permVec<zobs),sum(permVec>zobs))
        permut_pval = permut_pval/NPERM
        mylist[i]<- permut_pval
        res <- c("p.val"=permut_pval,"zobs"=zobs)
        return(res)
        }
    # The following line of code may take long time to execute ...
    pval.score <- sapply(rownames(inc.mat), drawPermutHist)
    tcga.gsea.res[[i]] <-pval.score
    }
    # Naming the list with study IDs
    names(tcga.gsea.res) <- names(tc.eset)
    # Transposing the result matrices within each studies and combinig them to a big matrix
    tcga.gsea.pathway.score.list <- list()
    for(i in 1: length(tcga.gsea.res)){
        temp.mat <- t(tcga.gsea.res[[i]])
        colnames(temp.mat) <- paste(colnames(temp.mat), names(tc.eset)[i], sep="_")
        tcga.gsea.pathway.score.list[[i]] <- temp.mat
    }
    tcga.mat <- do.call(cbind, tcga.gsea.pathway.score.list)
save(tcga.mat, file = fn)

# measuring time takes to run
end.time <- Sys.time()
time.taken <- end.time - start.time
write.table(c(start.time, end.time, time.taken), "Results/TCGA_RES/sys.time.10K.txt")
}
##################
# Create empty list
tcga.gsea.list <- list()
tcga.zobs <-tcga.mat[,seq(2, ncol(tcga.mat), 2)]
colnames(tcga.zobs) <- sapply(strsplit(colnames(tcga.zobs), "zobs_"), "[[", 2)
###############
tcga.pvals <-tcga.mat[,seq(1, ncol(tcga.mat), 2)]
colnames(tcga.pvals) <- sapply(strsplit(colnames(tcga.pvals), "p.val_"), "[[", 2)
tcga.gsea.list[[1]] <- tcga.zobs 
tcga.gsea.list[[2]] <- tcga.pvals
names(tcga.gsea.list) = c("Zobs","pvals")
rm(fn, tcga.pvals, tcga.zobs, tcga.mat)