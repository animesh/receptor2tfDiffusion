#' receptor2tfDiffusion
#'
#' To calculate Network diffusion, there are mainly two steps: creating Signalling Network (SigNet) and
#' use SigNet output in calculate single-sample network connectivity using diffusion model.
#'
#'
#' @param network a dataframe of two columns "from" and "to" with strings representing gene IDs
#' @param nodeW weight on the network node
#' @return
#' @examples
#' #network <- initializeSignallingNetwork(network = network, nodeW = nodeW)
#' #{ ... }
#' @export

initializeSignallingNetwork <- function(network,nodeW){

  if(!all(c(network$from,network$to) %in% names(nodeW))){
    stop('network genes missing in expression (nodeW)')
  }

  for(ii in 1:nrow(network)){
    network$edgeW[ii] <- nodeW[network$from[ii]]*nodeW[network$to[ii]]
  }

  ##Consider this one
  network$edgeW <- network$edgeW/(500*(max(network$edgeW)))
  network$flux <- 0


  signal <- rep(0,length(nodeW))
  names(signal) <- names(nodeW)
  nodes <- data.frame(nodeW= nodeW,signal = signal)
  network <- data.table(network)
  return(list(network= network,signal = signal))
}


#' Quantify network connectivity by calculating signal strength between receptors and TFs.
#'
#' @param network output from initilizeNetwork list with 'network' and 'signal' components
#' @param nodeW weight on the network node
#' @param outputNode output node
#' @param inputNode input node containing receptors in the signalling network
#' @param inputSignal minimum amount of signal placed as input for signal propagation
#' @param n number of iteration
#' @keywords Network Signalling
#' @return
#' @examples
#' #data(diffusionData)
#' #network <- initializeSignallingNetwork(network = network, nodeW = nodeW)
#' #signalNet<- signalOnNetwork(diffusionData$network, nodeW=M[,1], outputNode = 'FLI1',inputNode = 'TNFRSF1B', inputSignal = 0.99,n = 2000)
#' #{ ... }
#' @export
#'
#'
signalOnNetwork <-function(network,outputNode,inputNode,inputSignal = 0.99, n = 2000){

  signal <- network$signal
  network <- network$network
  signal[inputNode] <- inputSignal
  setkey(network,'from')
  Enode <- matrix(0,nrow = n,ncol = length(outputNode))
  colnames(Enode) <- outputNode

  for(ii in 1:n){
    network$delta <- signal[network$from]-signal[network$to]
    network$flux <- network$delta*network$edgeW
    f1 <- network[,.(flux = sum(flux)),by = from]
    r1 <- network[,.(flux = sum(flux)),by = to]
    signal[f1$from] <- signal[f1$from] - f1$flux
    signal[r1$to] <- signal[r1$to]+r1$flux
    Enode[ii,outputNode] <- signal[outputNode]
  }
  return(list(Enode=t(Enode),signal = signal))
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
#' @param M gene expression data
#' @keywords Diffusion Model
#' @return list of two components: Enode which is end node and Signal i.e. final state of network
#' @examples
#' data(diffusionData)
#' # TFs <- diffusionData$TFs
#' # diffusionMap <- function(receptors, TFs, M, network, nCores=2, nTicks=2000)
#' #{ ... }
#' @export

diffusionMap <- function(receptors, TFs, M, network, nCores=2, nTicks=2000){
  setDTthreads(1)
  if(!all(receptors %in% rownames(M))) {
    stop ("Please include receptors in the gene expression data (M) ")
  }
  if(!all(TFs %in% rownames(M))){
    stop("Please include TFs in the gene expression data (M) ")
  }
  if(!all(c(network[,1], network[,2]) %in% rownames(M))) {
    stop("Please include all network nodes in the gene expression data (M) ")
  }


  library(doParallel)
  registerDoParallel(cores=nCores)

  testf <- function(receptor){
    # print(fastSignalOnNetwork)
    aa <-   fastSignalOnNetwork(startingNet,outputNode = TFs,
                                inputNode = receptor, inputSignal = 0.99,n = nTicks)
    diff_time <- apply(aa[[1]],1,function(x){
      min(which(x > 0.5 * max(x)))
    })
  }

  sample_dtime <- vector('list',ncol(M))
  names(sample_dtime) <- colnames(M)

  for(jj in colnames(M)){
    nodeW <- M[unique(c(network[,1],network[,2])),jj]
    startingNet <- initializeSignallingNetwork(network,nodeW)
    diff_time <- foreach(receptor = receptors, .export=c('fastSignalOnNetwork','setkey'),.combine = rbind) %dopar% testf(receptor)
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
      for(ii in 1:ncol(M)){
        X[ii,count] <- sample_dtime[[ii]][receptor,tf]
      }
    }
  }

  return(X)
}



#' Network Graph
#'
#' This function generates network graph from receptor to TF by using shortest path with neighboring network nodes.
#'
#' @param network a dataframe of two columns "from" and "to" with strings representing gene IDs
#' @param receptors receptors in the network node
#' @param TFs transcription factors in the network node
#' @keywords getSubnetwork
#' @return list of two components: network subgraph from receptor i.e. start node to TF i.e. end node
#' @examples
#' # g <- graph_from_data_frame(network, directed = F, vertices = NULL)
#' #{ ... }
#' @export

getSubnetwork <- function(network,receptors,TFs){

  g <- graph_from_data_frame(network, directed = F, vertices = NULL)

  receptor2TFsubgraph <- function(g,from = 'OSMR',to = 'RELA'){
    gg <- all_shortest_paths(g, from = from, to = to)
    vs <- unique(unlist(gg$res))
    nb <- neighbors(g,vs)
    sg <- subgraph(g, unique(c(vs,nb)))
    return(sg)
  }

  nodes <- NULL

  for (tf in TFs){
    for (R in receptors){
      gg <- receptor2TFsubgraph(g,from = tf,to = R)
      gg <- igraph::simplify(gg,remove.multiple = F, remove.loops = T)
      nodes <- unique(c(nodes,names(V(gg))))
    }
  }

  newNetwork <- network[network[,1] %in% nodes & network[,2] %in% nodes,]
  return(newNetwork)
}



