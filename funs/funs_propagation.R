#' Remove genes in the dataset which don't appear in the network
#'
#' @details This function is included in the netDx use case which involves
#'   propagating the sparse matrix of patient's profiles to reduce its sparsity.
#'   This function removes the genes which appear in the patient's profiles and
#'   that are not represented as nodes in the network selected for the
#'   propagation. In fact, the network-based propagation computes a score for
#'   each gene/node depending on how much is connected to all the other ones. In
#'   case an expressed or mutated gene in the sparse matrix is not represented,
#'   it cannot get the new score of propagation.
#' @param m_profiles (data.frame) sparse matrix of patient profiles. Rownames
#'   are unique genes. Colnames are unique patients. A cell is a numeric value.
#' @param net (data.frame) adjancy matrix format of a network
#' @return (data.frame) sparse matrix of patient profiles with a number of
#'   rows/genes lower or equal to its starting value
#' @examples
#' load(nets_path)
#' load(datasets_path)
#' geno=geno_l[[1]]
#' net=nets_l[[1]]$net_adj
#' geno=complete_m(geno,net)
#' 
#' @export
complete_m = function(m_profiles,net){
  #bind the genes which are not expressed in the input matrix of patient profiles but they exist in the network
  #Recover patient genes and network genes
  pat_genes=rownames(m_profiles)
  network_genes=rownames(net)
  
  #Find genes which exist in the patients but do no exist in the network
  network_genes_2rem=setdiff(pat_genes,network_genes)
  #Remove
  m_profiles=m_profiles[-match(network_genes_2rem,pat_genes),]
  pat_genes=rownames(m_profiles)
  
  #Find genes which exist in the network but do no exist in the patients
  network_genes_2add=setdiff(network_genes,pat_genes)
  #Normalize
  m2bind=matrix(0,length(network_genes_2add),length(colnames(m_profiles)))
  m_profiles=rbind(m_profiles,m2bind)
  
  #Finish
  rownames(m_profiles)=c(pat_genes,network_genes_2add)
  return(m_profiles)
}


#' This function applies the random walk with restart propagation algorithm to a
#' matrix of patients profiles
#'
#' @details A network is an undirected graph G defined by a set of nodes
#'   corresponding to genes, and edges connecting nodes with an experimental
#'   evidence of interaction. A priori nodes are genes for which an information
#'   is known. A novel node is a candidate for being associated to the nodes
#'   above based on their information. A node prediction task leads to detect
#'   novel nodes and propagation techniques are largely applied for the purpose.
#'   Network-based propagation algorithms for node prediction transfer the
#'   information from a priori nodes to any other node in a network. Each node
#'   gets an imputation value which assesses how much information got. The
#'   prediction is based on the guilty-by-association principle. A node with a
#'   high imputation value has a high probability to be associated to a priori
#'   nodes. E.g. in a house where room A has one heater, if room B is the second
#'   hottest room it means that B is close to A and that there is a high
#'   probability that they share a door or wall. These algorithms exploit the
#'   global topology of the network. However, when they are applied to detect if
#'   unknown nodes are functionally associated to known ones, they may suffer of
#'   a drawback depending by the context. In biology, two functionally related
#'   fragments interact physically (direct interaction) or interact indirectly
#'   thanks to one or very few mediators. Therefore, exploring too far
#'   similarities between nodes can introduce noise in the prediction. We apply
#'   a random walk with restart propagation algorithm which resolution is set to
#'   0.2 for giving high values only to the close neighbours of the a priori
#'   nodes.
#' @param mat (data.frame) sparse matrix of patient profiles. Rownames
#'   are unique genes. Colnames are unique patients. A cell is a numeric value.
#' @param net (data.frame) adjancy matrix format of a network
#' @param cl (SOCKcluster) cluster object created with makeCluster function from parallel
#' @param no_cores (numeric) number of cores used to create the cluster object
#' @return (data.frame) continuous matrix of patient profiles in which each gene
#'   has the final propagation score
#' @examples
#' load(nets_path)
#' load(datasets_path)
#' geno=geno_l[[1]]
#' net=nets_l[[1]]$net_adj
#' geno=complete_m(geno,net)
#' cl <- makeCluster(no_cores)
#' registerDoParallel(cl)
#' prop_net=prop_m(geno,net,cl,no_cores=no_cores)
#' stopCluster(cl)
#' @export
prop_m=function(mat,net,cl,no_cores){
  #Split the matrix into sections, each one will be processed by one core
  inds <- split(seq_len(ncol(mat)), sort(rep_len(seq_len(no_cores), ncol(mat))))
  res.l=list()
  #Apply parallelized propagation
  res.l=foreach(k = 1:length(inds),.packages=c("netSmooth","scater","clusterExperiment")) %dopar% {
    nS.res=netSmooth(mat[,inds[[k]]], net , alpha=0.2, verbose = 'auto', normalizeAdjMatrix = c("columns")) 
    return(nS.res)
  }
  #Merge the results
  nS.res=do.call(cbind, res.l)
  #Return the final propagated matrix
  return(nS.res)
}

#' Apply discretization to the matrix resulted from the propagation on the
#' sparse patient matrix
#'
#' @details This function is included in the netDx use case which involves
#'   propagating the sparse matrix of patient's profiles to reduce its sparsity.
#'   This function applies discretization on the propagated matrix of patient
#'   profiles. It sets to 1 the genes which got the highest propagation value.
#'   While, the remaining genes are set to 0. This discretization is driven by
#'   the fact that higher is the propagation value and higher is the chance that
#'   the gene is involved in the patient condition and expression/mutation
#'   profile. On the contrary, genes which got either a medium or a low value
#'   are not trustable.
#' @param m_profiles (data.frame) continous matrix of patient profiles resulted from
#' applying the network-based propagation algorithm on a binary somatic mutation sparse 
#' matrix.
#' @param m_bin (data.frame) binary somatic mutation sparse matrix. Rownames
#'   are unique genes. Colnames are unique patients. A cell contains a zero or a one.
#' @return (data.frame) binary somatic mutation matrix which sparsity has been decreased
#' @examples
#' load(nets_path)
#' load(datasets_path)
#' geno=geno_l[[1]]
#' net=nets_l[[1]]$net_adj
#' geno=complete_m(geno,net)
#' cl <- makeCluster(no_cores)
#' registerDoParallel(cl)
#' prop_net=prop_m(geno,net,cl,no_cores=no_cores)
#' stopCluster(cl)
#' genoP=discr_prop_l(prop_net,geno,"project_title")
#' @export
discr_prop_l=function(m_prop,m_bin,name_dataset){
  m_prop=apply(-m_prop,2,rank)
  n_muts=colSums2(m_bin)
  n_topXmuts=c(3)
  
  m_props_l=list()
  for(k_top in 1:length(n_topXmuts)){
    name_prop=paste(name_dataset,"_x",n_topXmuts[k_top],sep="")
    n_new_muts=n_muts*n_topXmuts[k_top]
    for(i_col in 1:length(n_new_muts)){
      m_prop[m_prop[,i_col]<=n_new_muts[i_col],i_col]=1
      m_prop[m_prop[,i_col]>n_new_muts[i_col],i_col]=0
    }
    m_props_l[[name_prop]]=m_prop
  }
  
  if(length(m_props_l)!=1){
    return(m_props_l)
  }
  if(length(m_props_l)==1){
    return(m_prop)
  }
}


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