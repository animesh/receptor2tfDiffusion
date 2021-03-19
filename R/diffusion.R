#' Diffusion model
#'
#' To calculate Network diffusion, there are mainly two steps: creating Signalling Network (SigNet) and
#' use SigNet output in calculate single-sample network connectivity using diffusion model.
#'
#' Signalling Network
#'
#' Signalling etwork quantify network connectivity by calculating signal strength between receptors and TFs.
#'
#' @param network a dataframe of two columns "from" and "to" with strings representing gene IDs
#' @param nodeW
#' @param outputNode output node
#' @param inputNode input node containing receptors in the signalling network
#' @param inputSignal minimum amount of signal placed as input for signal propagation
#' @param n number of iteration
#' @keywords Network Signalling
#' @export
#' @return
#' @examples
#' toydata <- data(diffusionData)
#' signalNet<- signalOnNetwork(diffusionData$network, nodeW=M[,1], outputNode = 'FLI1',inputNode = 'TNFRSF1B', inputSignal = 0.99,n = 2000)
#' { ... }

signalOnNetwork <- function(network,nodeW,outputNode = 'FLI1',inputNode = 'TNFRSF1B', inputSignal = 0.99,n = 2000){

  ##Test that output and input exist in the network
  if(!(all(outputNode %in% network[,1] | outputNode %in% network[,2]))){
    print('output node not in network')
  }
  if(!(inputNode %in% network[,1] | inputNode %in% network[,2])){
    print('input node not in network')
  }

  connectivityMatrix <- matrix(0,nrow = length(nodeW),ncol = length(nodeW))

  rownames(connectivityMatrix) <- names(nodeW)

  colnames(connectivityMatrix) <- names(nodeW)

  deltaMatrix <- connectivityMatrix

  for(ii in 1:nrow(network)){

    network$edgeW[ii] <- nodeW[network$from[ii]]*nodeW[network$to[ii]]

  }

  network$edgeW <- network$edgeW/(500*(max(network$edgeW)))

  signal <- rep(0,length(nodeW))

  names(signal) <- names(nodeW)

  Enode <- rep(0,n)
  if(length(outputNode) > 1){
    Enode <- matrix(0,length(outputNode),n)
    rownames(Enode) <- outputNode
  }

  signal[inputNode] <- signal[inputNode] + inputSignal

  for(ii in 1:n){
    #signal[inputNode] <- signal[inputNode] + inputSignal
    delta <- signal[network$from] - signal[network$to]
    network$edgeFlux <- delta*network$edgeW

    for(jj in names(signal)){
      signal[jj] <- signal[jj] - sum(network$edgeFlux[which(network$from == jj)])
      signal[jj] <- signal[jj] + sum(network$edgeFlux[which(network$to == jj)])
    }

    if(length(outputNode) == 1){
      Enode[ii] <- signal[outputNode]
    }
    if(length(outputNode) > 1){
      Enode[,ii] <- signal[outputNode]
    }
  }

  return(list(Enode = Enode,signal = signal))
}

#' Diffusion model
#'
#' This function calculates single-sample network connectivity by using aggregate accumulation of signal per node as a function of time.
#'
#' @param receptors vector of string with ids.
#' @param TFs vector of string with gene ids of Transcription Factors
#' @param network a dataframe of two columns "from" and "to" with strings representing gene IDs
#' @param nCores number of cores for parallel processing
#' @param nTicks maximum number of time steps. Defaults to 2000.
#' @keywords Diffusion Model
#' @export
#' @return list of two components: Enode which is end node and Signal i.e. final state of network
#' @examples
#' data(diffusionData)
#' diffusionMap <- function(receptors, TFs, M, network, nCores=2, nTicks=2000)
#' { ... }
#' signalNet



diffusionMap <- function(receptors, TFs, M, network, nCores=2, nTicks=2000){
  if(!all(receptors %in% rownames(M))) {
     stop ("Please include receptors in the gene expression data (M) ")
  }
  if(!all(TFs %in% rownames(M))){
     stop("Please include TFs in the gene expression data (M) ")
  }
  if(!all(c(network[,1], network[,2]) %in% rownames(M))) {
     stop("Please include all network nodesin the gene expression data (M) ")
  }

  library(doParallel)
  registerDoParallel(cores=nCores)

  testf <- function(receptor){
     aa <-   signalOnNetwork(network,nodeW = nodeW,outputNode = selectedTFs,inputNode = receptor, inputSignal = 0.99,n = nTicks)
    diff_time <- apply(aa[[1]],1,function(x){min(which(x > 0.5*max(x)))})
  }

  sample_dtime <- vector('list',ncol(M))
  names(sample_dtime) <- colnames(M)

  for(jj in colnames(M)){
    nodeW <- M[unique(c(network[,1],network[,2])),jj]
    diff_time <- foreach(receptor = receptors, .combine = rbind) %dopar% testf(receptor)
    rownames(diff_time) <- receptors
    diff_time[diff_time > nTicks] <- nTicks + 1
    sample_dtime[[jj]] <- diff_time
  }

  X <- matrix(0,nrow = ncol(M), ncol = length(sample_dtime[[1]]))
  colnames(X)<- as.character(1:ncol(X))

  count <- 0

  for (receptor in receptors){
    for (tf in colnames(sample_dtime[[1]])){
      count <- count +1
      colnames(X)[count] <- paste(receptor,tf,sep = '_')
      for(ii in 1:54){
        X[ii,count] <- sample_dtime[[ii]][receptor,tf]
      }
    }
  }

  return(X)
}
