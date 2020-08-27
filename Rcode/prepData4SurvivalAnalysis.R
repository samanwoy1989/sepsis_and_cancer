## Prepare data for Survival analysis 
## (Samplewise tumor-normal data) 
## Each sample has 66 pathway scores
##
## First run code in main.R up to 
# ...
# source("Rcode/drawPathwayHeatmap.R")
# before calling this script
#################################################
# Get pathway expression data for selected 90 pathways
pathids = sel

# Function to get normalixed gene expression matrix
# for a given expression set
getSampleGexp = function(eset) {
  # use all samples with "A" label (exclude B, C, etc.)
  to.keep = which((substr(colnames(eset),16,16))=="A")
  eset = eset[,to.keep]
  # get sample id
  all.ids = substr(colnames(eset),1,12)
  smpl.ids = which(table(all.ids)==2) %>% names
  t.id = paste0(smpl.ids, ".01A_mR")
  n.id = paste0(smpl.ids, ".11A_mR")
  gexp = exprs(eset[,t.id]) - exprs(eset[,n.id])
  colnames(gexp) = smpl.ids
  return(gexp)
}
# Get normalized gene expression for each case
list.gexp = sapply(tc.eset, getSampleGexp)

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
list.zobs = sapply(list.gexp, getZobs)

# Format data for classification
# Transpose data such that each row is a subject
list.zobs.trans = sapply(list.zobs, t)
list.input.slc = list.zobs.trans[ids.sepsis.like]
list.input.ca = list.zobs.trans[ids.cancer]

# format data into a single matrix
# first, create the matrix for SLC
id = ids.sepsis.like[1]
nsamples = nrow(list.input.slc[[id]])
tissue = rep(id, nsamples)
label = rep("SLC", nsamples)
x = data.frame(list.input.slc[[id]], Tissue = tissue, label = label)

for(i in 2:length(ids.sepsis.like)) {
  id = ids.sepsis.like[i]
  nsamples = nrow(list.input.slc[[id]])
  tissue = rep(id, nsamples)
  label = rep("SLC", nsamples)
  x = rbind(x, data.frame(list.input.slc[[id]], Tissue = tissue, label = label))
}

# add the CA samples
for(i in 1:length(ids.cancer)) {
  id = ids.cancer[i]
  nsamples = nrow(list.input.ca[[id]])
  tissue = rep(id, nsamples)
  label = rep("CA", nsamples)
  x = rbind(x, data.frame(list.input.ca[[id]], Tissue = tissue, label = label))
}

# make an expression set
assaydata = t(x[,1:66])
phenodata = new("AnnotatedDataFrame", x[,67:68])
pathES = new("ExpressionSet", exprs=assaydata, phenoData=phenodata)

# remove unnecessary variables
rm(list=c("assaydata","getPathScore","getZobs","id","label","list.input.ca","list.input.slc","list.zobs",
     "list.zobs.trans","nsamples","pathids","phenodata","tissue","x"))

###################################################
# Supplementary Table S3
# Survival analysis
# For each of the selected pathways
# do survival analysis and save the p-values 
# in a matrix
###################################################

fn_tableS3 = "Results/pathbycancertypeSurvPval.csv"
if(file.exists(file=fn_tableS3)) {
  path66ByProj = read.table(file=fn_tableS3, header=TRUE, sep="\t")
} else {
  proj.ids = c("HNSC","KIRC","LIHC","BRCA","LUAD","LUSC","PRAD","THCA")
  path66ByProj = matrix(NA, nrow=length(sel), ncol=length(proj.ids))
  colnames(path66ByProj) = proj.ids
  rownames(path66ByProj) = sel
  
  for(pathid in rownames(path66ByProj)) {
    pathstr = pathways.list[paste0("path:",pathid)]
    pathstr = strsplit(pathstr, split=" -")[[1]][1]
    pathstr = paste0("Pathway: ", pathstr, " [", pathid, "]")
    #graphics.off(); rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)
    for(proj in colnames(path66ByProj)) {
      graphics.off(); rm(sfit, dat, info, pathX, projES, to.keep, sel.ptids)
      # get the clinical metadata from an external file and add to the expression set
      info = read.table(file=paste0("Metadata/",proj,"_surv_clinical.tsv"), header=T, sep="\t")
      to.keep = (apply(info[,2:3], 1, function(x) !all(x=="--")))
      info = info[to.keep,]
      is.na(info) <- info=="--"
      
      rownames(info) = gsub("-",".", info$submitter_id)
      
      sel.ptids = intersect(rownames(info), sampleNames(pathES))
      #print.noquote(length(sel.ptids))
      projES = pathES[, sel.ptids]
      info = info[sel.ptids,]
      
      pathX = as.vector(exprs(projES[pathid, ]))
      names(pathX) = sampleNames(projES)
      dat = data.frame(info, pathX)
      dat$day = as.integer(as.character(dat$days_to_death))
      dat$day[dat$outcome==0] = as.integer(as.character(dat$days_to_last_follow_up[dat$outcome==0]))
      dat$high_expressing = as.integer(dat$pathX > median(dat$pathX))
      dat = dat[!is.na(dat$day),]
      #print.noquote(nrow(dat))
      sfit <- survfit(Surv(day, outcome)~high_expressing, data=dat)
      res = ggsurvplot(sfit, conf.int=FALSE, pval=TRUE, risk.table=TRUE, 
                       legend.labs=c("Low", "High"), 
                       legend.title=paste0(pathid, " Score"),  
                       palette=c("dodgerblue2", "orchid2"), 
                       title=paste0("Kaplan-Meier Curve for ", proj, " Survival\n",pathstr), 
                       risk.table.height=.25)
      pval = surv_pvalue(sfit, data = NULL, method = "survdiff",
                         test.for.trend = FALSE, combine = FALSE)$pval
      cat(paste0(pathid, " ", proj, " ", formatC(pval, digits=2), "\n"))
      print(res)
      path66ByProj[pathid, proj] = pval
    }
  }
  
  selnames = sapply(strsplit(pathways.list[paste0("path:",sel)], split=" -"), function(x) x[[1]][1]) %>% as.character
  rownames(path66ByProj) = paste0(selnames, " | ", sel)
  write.table(path66ByProj, file="Results/pathbycancertypeSurvPval.csv", row.names=TRUE, col.names = TRUE, quote=F, sep="\t")
}

ids.slc = intersect(colnames(path66ByProj), ids.sepsis.like)
ids.ca  = intersect(colnames(path66ByProj), ids.cancer)
percSignifSurv = apply(path66ByProj, 2, function(x) 100*sum(x<0.1)/length(x))
mycols = c(rep("blue",length(ids.ca)), rep("cyan",length(ids.slc)))
names(mycols) = c(ids.ca, ids.slc)
percSignifSurv = percSignifSurv[c(ids.ca, ids.slc)]
