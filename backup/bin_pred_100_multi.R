
################################


rm(list=ls())
rootData="/home/wr/%s"
#Load libraries ----
require(netDx)
require(netDx.examples)
library(matrixStats)
require(data.table)

#Build the path of the input and output files ----
setwd(sprintf(rootData,"WJ/env/data_bin/"))
rda_files = list.files(pattern="*.rda$")
pred_dir=sprintf(rootData,"WJ/predictions/BIN_30_%s")


run_time_files=c("/home/wr/WJ/env/data_bin/BIN_GBM_time.txt",
"/home/wr/WJ/env/data_bin/BIN_HNSC_time.txt",
"/home/wr/WJ/env/data_bin/BIN_KIRC2_time.txt",
"/home/wr/WJ/env/data_bin/BIN_LUAD_time.txt",
"/home/wr/WJ/env/data_bin/BIN_LUSC_time.txt",
"/home/wr/WJ/env/data_bin/BIN_OV_time.txt",
"/home/wr/WJ/env/data_bin/BIN_SKCM_time.txt")

for(file in rda_files){
start.time <- Sys.time()
  envFile=file
  outDir=sprintf(pred_dir,gsub(".rda","",file))
  print("I am running with this envFile and outDir:")
  print(envFile)
  print(outDir)
  #Load data ----
  load(file=envFile)
  pheno[]=lapply(pheno,as.character)

#PSN function passed to netDx ----
makeNets <- function(dataList,groupList,netDir,numCores,...) {
  netList <- c(); netList2 <- c()
  
  # create genetic nets
  if (!is.null(groupList[["genetic"]])) {netList <- makeMutNets(dataList[["genetic"]],
                                                                groupList[["genetic"]],
                                                                netDir,numC=numCores)
  }
  cat(sprintf("\t%i genetic-pathway nets\n", length(netList)))
  
  cat(sprintf("Total of %i nets\n", length(netList)))
  
  return(netList)
}

# g geno matrix, genes by patients (columns) - binary
# pList list of genesets
# outDir - dir where nets are to be written
makeMutNets <- function(g,pList,oDir,numC) {
  g <- t(g) # transpose to have genes as columns
  cl	<- makeCluster(numC)
  registerDoParallel(cl)
  
  numPat <- c()
  netList <- foreach(k=1:length(pList)) %do% {
    idx <- which(colnames(g) %in% pList[[k]])
    
    if (length(idx)>0) {
      has_mut <- rowSums(g[,idx,drop=FALSE])
      has_mutp <- names(has_mut)[which(has_mut>0)]
      
      if (length(has_mutp)>=6) {
        cat(sprintf("%s: %i patients\n", names(pList)[k],
                    length(has_mutp)))
        #numPat <- c(numPat, length(has_mutp))
        pat_pairs <- t(combinat::combn(has_mutp,2));
        pat_pairs <- cbind(pat_pairs,1);
        outFile <- sprintf("%s/%s_cont.txt",oDir,names(pList)[k])
        write.table(pat_pairs, file=outFile,sep="\t",
                    col=FALSE,row=FALSE,quote=FALSE)
        basename(outFile)
      } else NULL
    } else {
      NULL
    }
  }
  stopCluster(cl)
  unlist(netList)
}

#Setup to run predictor -----
datList 	<- list()
groupList	<- list()

#pathFile	<- sprintf("%s/extdata/Human_160124_AllPathways.gmt",path.package("netDx.examples"))
pathFile="/home/wr/Human_April_2018.gmt"
pathwayList <- readPathways(pathFile,MIN_SIZE=5L,MAX_SIZE=300L)

datList[["genetic"]] <- geno;
groupList[["genetic"]] <- pathwayList

#Run the predictor ----
if (file.exists(outDir)){
  cat("the directory of output already exists")
  x <- readline("Do you want to overwrite the directory: y|n?")
  if(x=="y"){unlink("outDir", recursive=TRUE)}
  else {
    cat("Stop the execution")
    stop(call. = FALSE)
  }
}

runPredictor_nestedCV(pheno,
                      dataList=datList,groupList=groupList,
                      makeNetFunc=makeNets, ### custom network creation function
                      outDir=outDir, ## absolute path
                      numCores=18L,nFoldCV=10L, CVcutoff=9L,numSplits=30L, CVmemory=16L
)

Sys.sleep(120)

  print("I finished a prediction")
end.time <- Sys.time()
time.taken <- end.time - start.time
file_run=paste(envFile,"time.txt",sep="")
write.table(time.taken,file=run_time_files[1],quote = FALSE)
run_time_files=run_time_files[-1]
}
