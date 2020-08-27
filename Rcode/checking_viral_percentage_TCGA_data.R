# this code calculates the amount of 
# viral samples in our analysis in TCGA cancer data sets
# from few published articles

####################################
# Use the packages for loading data
###################################
#require("Biobase")
#require("xlsx")
#require("ssnibmg")
#require("tcnibmg")
#require("ssgeosurv")
# Load the data
################
#data("ss.eset")
#data("tc.eset")
#data("ss.surv.list")
###############################
# preliminaries
#source("../analysis/Rcode/prelims.R")
#######################################
#               Tang 2013
#######################################
our.metadata = pData(tc.eset[[1]])
for(i in 2:length(tc.eset)) {
  our.metadata = rbind(our.metadata, pData(tc.eset[[i]]))
}
our.smpls = as.character(unlist(sapply(tc.eset, function(eset) 
  {sampleNames(eset)[which(eset$Target=="Case")]
  })))
our.metadata = pData(tc.eset[[1]])
for(i in 2:length(tc.eset)) {
  our.metadata = rbind(our.metadata, pData(tc.eset[[i]]))
}
our.metadata = our.metadata[our.smpls,]

# Read data from Tang 2013
############################
#tang.dat = read.csv("tang_2013.csv")
tang.dat = read.xlsx2(file="Metadata/41467_2013_BFncomms3513_MOESM350_ESM.xlsx",
                      sheetIndex = 1, startRow = 4)
smpl.ids = as.character(tang.dat$Sample.barcode)
smpl.ids = unlist(lapply(smpl.ids, function(x) {
  (paste0(substr(gsub("-", ".", x), 1, 16), "_mR"
  ))}))
which.no.dups = which(!duplicated(smpl.ids))
tang.dat = tang.dat[which.no.dups,]
rownames(tang.dat) = smpl.ids[which.no.dups]
rm(which.no.dups, smpl.ids)

# Finding the overlapping index and subsetting the tang data 
# for overlapping data only
################################
slc.ids = c("HNSC","LIHC","KIRC","CHOL", "ESCA","STAD")
nSamples = length(our.smpls)
overlp.ind = which(rownames(tang.dat)%in%our.smpls)
nOverlap = length(overlp.ind)
nOverlapPerc = 100*nOverlap/nSamples
overlp.dat = tang.dat[overlp.ind,]
tissue.dist.overlap = table(our.metadata[rownames(overlp.dat),"Code"])
nSLC = sum(tissue.dist.overlap[slc.ids])
nCA = sum(tissue.dist.overlap)-nSLC
nSLCPerc = 100*nSLC/nOverlap
nCAPerc = 100*nCA/nOverlap
nVirus = sum(as.character(overlp.dat$Virus.sequence.ID)!="N/A")

sample.ids.virus = rownames(overlp.dat)[which(overlp.dat$Virus.sequence.ID!="N/A")]
tissue.dist.overlap.with.virus = 
  table(our.metadata[sample.ids.virus,"Code"])
nVirusSLC = sum(tissue.dist.overlap.with.virus[slc.ids])
nVirusCA = nVirus-nVirusSLC
nVirusSLCPerc = 100*nVirusSLC/nSLC
nVirusCAPerc = 100*nVirusCA/nCA

# for chi-squared test
count = matrix(c(nVirusSLC, nSLC-nVirusSLC,
                 nVirusCA, nCA-nVirusCA), 
               nrow=2)
colnames(count) = c("SLC","CA")
rownames(count) = c("Virus","noVirus")
pval = formatC(chisq.test(count)$p.value, digits=1)

# Create and Print the data table
Tang2013= c("Samples analysed in current manuscript", "Overlap with published study", "SLC samples", "CA samples", "Viral integration in SLC", "Viral integration in CA")
Numbers = c(nSamples, nOverlap, nSLC, nCA, nVirusSLC,  nVirusCA)
Percentage =c(NA, nOverlapPerc, nSLCPerc, nCAPerc, nVirusSLCPerc, nVirusCAPerc)

tab.dat = as.data.frame(cbind(Tang2013, Numbers, Percentage=formatC(Percentage,3)), row.names = NULL)
print(knitr::kable(tab.dat))
cat("P-value = ", pval, "\n")

rm(nCA, count, nCAPerc, nOverlap, nOverlapPerc, nSLC, nSLCPerc, Numbers, nVirus,
   nVirusCA, nVirusCAPerc, nVirusSLC, nVirusSLCPerc, overlp.dat, overlp.ind,
   Percentage, sample.ids.virus, tab.dat, tang.dat, Tang2013, tissue.dist.overlap,
   tissue.dist.overlap.with.virus, pval)

#######################################
#               Cao 2016
#######################################
cat("\n\n")
cao.dat = read.xlsx2(file="Metadata/41598_2016_BFsrep28294_MOESM2_ESM.xls",
                     sheetIndex = 1, startRow = 2)
smpl.ids = as.character(cao.dat$Sample)
smpl.ids = unlist(lapply(smpl.ids, function(x) {
  paste0(substr(gsub("-", ".", x), 1, 16), "_mR")
}))
which.no.dups = which(!duplicated(smpl.ids))
cao.dat = cao.dat[which.no.dups,]
rownames(cao.dat) = smpl.ids[which.no.dups]
rm(which.no.dups, smpl.ids)

# Format the data for one sum score across all the data
#######################################################
temp.mat = as.matrix(cao.dat[,-c(1,2)])>0
temp.vec = rowSums(temp.mat)
cao.dat = cbind(cao.dat[,c(1,2)], virInt= temp.vec)

# Finding the overlapping index and 
# subsetting the Cao data 
# for overlapping data only
################################
nSamples = length(our.smpls)
overlp.ind = which(rownames(cao.dat)%in%our.smpls)
nOverlap = length(overlp.ind)
nOverlapPerc = 100*nOverlap/nSamples
overlp.dat = cao.dat[overlp.ind,]
tissue.dist.overlap = table(our.metadata[rownames(overlp.dat),"Code"])
nSLC = sum(tissue.dist.overlap[slc.ids])
nCA = sum(tissue.dist.overlap)-nSLC
nSLCPerc = 100*nSLC/nOverlap
nCAPerc = 100*nCA/nOverlap
nVirus = sum(as.character(overlp.dat$virInt)!="0")

sample.ids.virus = rownames(overlp.dat)[which(overlp.dat$virInt!="0")]
tissue.dist.overlap.with.virus = 
  table(our.metadata[sample.ids.virus,"Code"])
nVirusSLC = sum(tissue.dist.overlap.with.virus[slc.ids])
nVirusCA = nVirus-nVirusSLC
nVirusSLCPerc = 100*nVirusSLC/nSLC
nVirusCAPerc = 100*nVirusCA/nCA

# for chi-squared test
count = matrix(c(nVirusSLC, nSLC-nVirusSLC,
                 nVirusCA, nCA-nVirusCA), 
               nrow=2)
colnames(count) = c("SLC","CA")
rownames(count) = c("Virus","noVirus")
pval = formatC(chisq.test(count)$p.value, digits=1)

# Create and Print the data table
Cao2016= c("Samples analysed in current manuscript", "Overlap with published study", "SLC samples", "CA samples", "Viral integration in SLC", "Viral integration in CA")
Numbers = c(nSamples, nOverlap, nSLC, nCA, nVirusSLC,  nVirusCA)
Percentage =c(NA, nOverlapPerc, nSLCPerc, nCAPerc, nVirusSLCPerc, nVirusCAPerc)

tab.dat = as.data.frame(cbind(Cao2016, Numbers, Percentage=formatC(Percentage,3)), row.names = NULL)
print(knitr::kable(tab.dat))
cat("P-value = ", pval, "\n")

rm(nCA, count, nCAPerc, nOverlap, nOverlapPerc, nSLC, nSLCPerc, Numbers, nVirus,
   nVirusCA, nVirusCAPerc, nVirusSLC, nVirusSLCPerc, overlp.dat, overlp.ind,
   Percentage, sample.ids.virus, tab.dat, cao.dat, Cao2016, tissue.dist.overlap,
   tissue.dist.overlap.with.virus, pval)

#######################################################
#               made from Table 1 of kazemian2015
#######################################################
cat("\n\n")
data("tc.eset")
codes = names(tc.eset)
kaze2015 = t(data.frame(
  "BLCA"=c(119, 9),
  "BRCA"=c(125, 1),
  "CHOL"=c(0, 0),
  "COAD"=c(44, 10),
  "ESCA"=c(0, 0),
  "HNSC"=c(123, 23),
  "KICH"=c(66, 0),
  "KIRC"=c(67, 2),
  "KIRP"=c(120, 0),
  "LIHC"=c(115, 32),
  "LUAD"=c(125, 0),
  "LUSC"=c(125, 5),
  "PRAD"=c(124, 2),
  "READ"=c(36, 9),
  "STAD"=c(0, 0),
  "THCA"=c(123, 0), 
  "UCEC"=c(168, 30)))

colnames(kaze2015) = c("NumSamples","ViralInt")
slc.ids = intersect(rownames(kaze2015),c("HNSC","LIHC","KIRC","CHOL", "ESCA","STAD"))
ca.ids = setdiff(rownames(kaze2015), slc.ids)
nSLC = sum(kaze2015[slc.ids, "NumSamples"])
nCA = sum(kaze2015[, "NumSamples"])-nSLC

nSLC.V = sum(kaze2015[slc.ids, "ViralInt"])
nCA.V = sum(kaze2015[ca.ids, "ViralInt"])

percSLC.V = 100*nSLC.V/nSLC
percCA.V = 100*nCA.V/nCA

# Create and Print the data table
Kazemian2015= c("SLC samples", "CA samples", "Viral integration in SLC", "Viral integration in CA")
Numbers = c(nSLC, nCA, nSLC.V,  nCA.V)
Percentage = c((nSLC/sum(kaze2015[,1]))*100, (nCA/sum(kaze2015[,1]))*100, percSLC.V, percCA.V)

tab.dat = as.data.frame(cbind(Kazemian2015, Numbers, Percentage=formatC(Percentage, 3)), row.names = NULL)
print(knitr::kable(tab.dat))

# for chi-squared test
count = matrix(c(nSLC.V, nSLC-nSLC.V, nCA.V, nCA-nCA.V), nrow=2)
colnames(count) = c("SLC","CA")
rownames(count) = c("Virus","noVirus")
pval = formatC(chisq.test(count)$p.value, digits=1)
cat("P-value = ", pval, "\n")





