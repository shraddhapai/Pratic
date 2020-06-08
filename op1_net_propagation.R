#Clean workspace and memory ----
rm(list=ls())
gc()

#Set working directory----
gps0=getwd()
gps0=paste(gps0,"/%s",sep="")
rootDir=gps0
setwd(gsub("%s","",rootDir))
set.seed(8)

#Load libraries----
library("data.table")
library("matrixStats")
library("plyr")
library("doParallel")
library("parallel")
library("netSmooth")
library("clusterExperiment")
library("scater")
library("netDx")
source("funs/funs_propagation.R")

#Set input and output paths ----
#Input
nets_path="input/nets_l.rda"
datasets_path="input/datasets.rda"
pathway_list="input/Human_AllPathways_March_01_2020_symbol.gmt"
#Output
outDir="output/OV"

#Load data
load(nets_path)
load(datasets_path)

#Set data
no_cores=8L
cat(names(geno_l),"\n");cancer_name="OV";
cat(names(nets_l),"\n");net_name="CANCER";
datList=list()
groupList=list()

#Start use case of propagation integrated with prediction ----

#Select cancer somatic mutation dataset
i_geno=grep(cancer_name,names(geno_l),fixed=T)
#Select gene-gene interaction network to use for the propagation
j_net=grep(net_name,names(nets_l),fixed=T)

#Retrieve dataset and info  
geno_name=names(geno_l)[i_geno]
geno=geno_l[[i_geno]]
pheno=pheno_l[[i_geno]]
#Retrieve adjancy matrix  
net_name=names(nets_l)[j_net]
net=nets_l[[j_net]]$net_adj

#Remove genes in the dataset which don't appear in the network
geno=complete_m(geno,net)
#Propagation
cl <- makeCluster(no_cores)
registerDoParallel(cl)
  prop_net=prop_m(geno,net,cl,no_cores=no_cores)
stopCluster(cl)
#Set the name of the project and future resulting directory
project_title=paste(geno_name,"_x",net_name,sep="")
#Apply binarization
genoP=discr_prop_l(prop_net,geno,project_title)

#Setup to build the predictor
pathwayList <- readPathways(pathway_list,MIN_SIZE=5L,MAX_SIZE=300L)
exprdat <- SummarizedExperiment(genoP, colData=pheno)
groupList[["genetic"]] <- pathwayList

#Run the predictor
buildPredictor(dataList=exprdat,groupList=groupList,
  makeNetFunc=makeNets, ### custom network creation function
  outDir=outDir, ## absolute path
  numCores=no_cores, numSplits=2L
)
print("END")
#nFoldCV=10L, CVcutoff=9L,numSplits=30L, CVmemory=8L
