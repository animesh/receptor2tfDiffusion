# receptor2tfDiffusion: 

**Diffusion based network model**
install_github('Amy3100/receptor2tfDiffusion',force = T)
library(receptor2tfDiffusion)
data(diffusionData)
network <- diffusionData$network 
nodeW <- diffusionData$M[unique(c(network[,1],network[,2])),24]

##Extract subnetworks
##rundiffusion
##extract subnetwork for matrix creation 

diffusion_curve <- signalOnNetwork(network,nodeW,outputNode = 'FLI1',inputNode = 'TNFRSF1B',n = 2000)
plot(diffusion_curve$Enode)

R <- c('TNFRSF11B','OSMR','IL6R')
TFs <- c('NFKB1','HNF4A','MAF','IRF7')
library(igraph)
netLight <- getSubnetwork(network,receptors = R, TFs = TFs)

t50 <- diffusionMap(receptors = R, TFs = TFs, M, network = netLight, nCores = 2, nTicks = 1000) 

library(limma)
##limma stuff
dMat <- model.matrix(~0 + as.factor(sampleSheet$group))
colnames(dMat) <- c('ctrl','postResistant','postResponse','preResistant','preResponse')
fit <- lmFit(t50,dMat)
contrastMatrix <- makeContrasts(preResistant-ctrl,preResistant-preResponse,levels = dMat)
fit2 <- contrasts.fit(fit,contrasts = contrastMatrix)
decideTests(eBayes(fit2),p.value = 0.01)
