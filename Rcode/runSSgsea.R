#########################################################################
# This code runs GSEA Permutaion analysis to calculate Pathway score for
# 6 Septic Shock datasets 
# Samanwoy Fri Mar  1 12:53:42 IST 2019
#########################################################################
fn <- "Results/pathway_score_gsea_nperm_SS6_10k.rda"
if(file.exists(file= fn)) {
    cat("Loading SS GSEA file ...")
    load(file= fn )
    cat(" done!\n")
} else {
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
    ##########################################################
    ss.gsea.res <- list()
    for ( i in 1: length(ss.eset)) {
        selg <- intersect(colnames(inc.mat), featureNames(ss.eset[[i]]))
        inc.mat <- inc.mat[,selg]
        gtable <- exprs(ss.eset[[i]])[selg,]
        set.seed(123)
        NPERM = 10000
        # Performing differntial expression
        fac <- ss.eset[[i]]$Group
        require(genefilter)
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
    ss.gsea.res[[i]] <-pval.score
    }
    # Naming the list with study IDs
    names(ss.gsea.res) <- names(ss.eset)
    # Transposing the result matrices within each studies and combinig them to a big matrix
    ss.gsea.pathway.score.list <- list()
    for(i in 1: length(ss.gsea.res)){
        temp.mat <- t(ss.gsea.res[[i]])
        colnames(temp.mat) <- paste(colnames(temp.mat), names(ss.eset)[i], sep="_")
        ss.gsea.pathway.score.list[[i]] <- temp.mat
    }
    ss.mat <- do.call(cbind, ss.gsea.pathway.score.list)
    save(ss.mat, file = "Results/pathway_score_gsea_nperm_SS6_10k.rda")
}
##################
# Create empty list
ss.gsea.list <- list()
ss.zobs <-ss.mat[,seq(2, ncol(ss.mat), 2)]
colnames(ss.zobs) <- sapply(strsplit(colnames(ss.zobs), "zobs_"), "[[", 2)
###############
ss.pvals <-ss.mat[,seq(1, ncol(ss.mat), 2)]
colnames(ss.pvals) <- sapply(strsplit(colnames(ss.pvals), "p.val_"), "[[", 2)
ss.gsea.list[[1]] <- ss.zobs 
ss.gsea.list[[2]] <- ss.pvals
names(ss.gsea.list) = c("Zobs","pvals")
rm(ss.mat, ss.pvals, ss.zobs, fn)