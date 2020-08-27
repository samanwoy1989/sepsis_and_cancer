## Prepare data for validation of machine learning (SVM, NN)  analysis 
## (Samplewise tumor-normal data) 
## Each sample has 66 pathway scores
## Data are from GEO
## First run code in main.R up to 
# ...# Prepare validation data for machine learning
     # source("Rcode/prepareData4MLvalidation.R")
#####################################################
# Get pathway expression data for selected 90 pathways
pathids = sel
#pathids = rownames(sstcga)

# Function to get normalixed gene expression matrix
# for a given expression set
getValSampleGexp = function(eset, paired=TRUE) {
  if(paired==TRUE) {   # for paired data, i.e., tumour, adjacent normal (for GSE32665 and GSE103236)
    eset$PTID = factor(eset$PTID)
    #which.paired = eset$PTID%in%eset$PTID[which.paired]
    #which.paired = which(duplicated(eset$PTID))
    #eset = eset[,which.paired]
    to.keep = which(eset$PTID%in%names(which(table(eset$PTID)==2)))
    eset = eset[,to.keep]
    uniptid = unique(eset$PTID)
    which.case = which(eset$Group!="Zcontrol")
    x.case = exprs(eset[,which.case])
    colnames(x.case) = eset$PTID[which.case]
    x.case = x.case[,uniptid]
    which.con = which(eset$Group=="Zcontrol")
    x.con = exprs(eset[,which.con])
    colnames(x.con) = eset$PTID[which.con]
    x.con = x.con[,uniptid]
    gexp = x.case-x.con
  } else {   # for unpaired data, i.e., cancer and control subjects (BRCA COLRECT NSCLC PRAD GSE112790 GSE33479 GSE84984)
    require(genefilter)
    which.con =  which(eset$Group=="Zcontrol")
    which.case = which(eset$Group!="Zcontrol")
    gexp = sapply(which.case, function(caseid) {
      ids = c(caseid, which.con)
      lfc = rowttests(eset[,ids], "Group")[,"dm"]
    })
    colnames(gexp) = sampleNames(eset)[which.case]
    rownames(gexp) = featureNames(eset)
  }
  return(gexp)
}

################################################################
# the validation expression sets are loaded from tcnibmgML package
# the name vset.list is not changed because it makes more sense 
# as the data structure (a list of expression sets) is guessed 
# from the name vset.list
vset.list = tcnibmgML

# Get normalized gene expression for each case
is.paired=sapply(vset.list, function(eset) "PTID"%in%varLabels(eset)) %>% as.logical
list.val.gexp = list()
for(i in 1:length(vset.list)) {
  eset = vset.list[[i]]
  paired = is.paired[i]
  list.val.gexp[[i]] = getValSampleGexp(eset=eset, paired=paired)
  # if names of eset includes GSE then change to tissue code
  if(substr(names(vset.list)[i],1,3)=="GSE") {
    names(list.val.gexp)[i] = setdiff(unique(eset$Group), "Zcontrol")
  } else {
    names(list.val.gexp)[i] = names(vset.list)[i]
  }
}

# Function to get pathway score (zobs) given a pathway id and 
# a gene expression matrix
getPathScore = function(pathid=pathid, gexp) {
  curr.egs = intersect(genes.by.pathway[[pathid]], rownames(gexp))
  x = gexp[curr.egs,]
  apply(x, 2, function(y) {
    y = sum(y)/sqrt(length(y))
  })
}

# Function to get the pathway expression matrix for a gene expression matrix
getZobs = function(gexp) {
  sapply(pathids, function(pathid) {
    getPathScore(pathid, gexp)
  }) %>% t  
}

# get a list of pathway score matrix for all studies
list.val.zobs = sapply(list.val.gexp, getZobs)

# assay data
x = list.val.zobs[[1]]
for(i in 2:length(list.val.zobs)) {
  x = cbind(x, list.val.zobs[[i]])
}
# cancer type, hard-coding
caType = rep("CA", length(list.val.zobs))
names(caType) = names(list.val.zobs)
caType[c("LIHC","STAD")] = "SLC"
# Number of samples for each tissue
nSamples = sapply(list.val.zobs, ncol)
tissue = rep(names(nSamples), times=nSamples)
label = rep("CA", length(tissue))
label[tissue%in%c("LIHC", "STAD")] = "SLC"
pheno.dat = as.data.frame(cbind(tissue, label))
rownames(pheno.dat) = colnames(x)

# make an expression set
assaydata = x
phenodata = new("AnnotatedDataFrame", pheno.dat)
pathES = new("ExpressionSet", exprs=assaydata, phenoData=phenodata)