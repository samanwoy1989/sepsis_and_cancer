# This code performs preliminary tasks, loading libraries, metadata, etc.
# Samanwoy Mukhopadhyay

# Significance threshold
sig.threshold <- 0.01

# Loads libraries
####################################################################
require("org.Hs.eg.db")
require("genefilter")
require("gplots")
require("rgl")
require("pca3d")
require("igraph")
require("stats")
require("survival")
require("survminer")
require("xlsx")

# Use the packages for data
###################################
library("ssnibmg") # Mukhopadhyay et.al 2017
library("ssgeosurv") # Septic shock data from GEO with survival information 
library("tcnibmg") # TCGA data for cancer
library("tcnibmgML") # Cancer data from GEO for ML-Validatation of models developed with TCGA data

# Get the SS data "ssnibmg" 
data("ss.eset")
# Get the TCGA data 
data("tc.eset")
# Get the data of ss survivors from GEO
data("ss.surv.list")
# Get the data for SLC/CA data from GEO for validatation of signature
# from TCGA data analysis
data("tcnibmgML")

# Params1: TCGA Code-to-Tissue mapping
####################################################################
params1 = read.table("Metadata/params.txt", header=TRUE) 

# Function definitions
####################################################################
source("Rcode/funcdefs/drawPathwayGeneHeatmap.R")
source("Rcode/funcdefs/drawBoxplot3DiseaseGroup.R")

# Load the gene by pathway list (KEGG Database)
#########################################################
load("Metadata/gbyp.rda")
