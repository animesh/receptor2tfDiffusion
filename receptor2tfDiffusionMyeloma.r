#setup#####
sudo apt update -qq
sudo apt install --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt update 
sudo apt -y upgrade 
sudo apt -y install --no-install-recommends r-base
sudo apt -y install libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev 
install.packages("pkgdown")
install.packages("devtools")
install.packages("BiocManager")
sudo apt -y install gfortran libblas-dev liblapack-dev
BiocManager::install("igraph")
library(igraph)
BiocManager::install("limma")
library(limma)
library(devtools)
install.packages("doParallel")
library(doParallel)
install.packages("data.table")
library(data.table)
#receptor2tfDiffusion#####
devtools::install_github('Amy3100/receptor2tfDiffusion',force = T)
library(receptor2tfDiffusion)
#data####
data(diffusionData)
names(diffusionData)
network <- diffusionData$network 
nodeW <- diffusionData$M[unique(c(network[,1],network[,2])),24]
diffusion_curve <- signalOnNetwork(network,nodeW,outputNode = 'FLI1',inputNode = 'TNFRSF1B',n = 2000)
hist(diffusion_curve$signal, main = "Final signal intensity", xlab = "Signal" )
#par(mar=c(0,0,0,0))
plot(diffusion_curve$Enode, main= "Diffusion curve" , ylab= "signal at FLI1", xlab= "time ticks")
R <- c('TNFRSF11B','OSMR','IL6R')
TFs <- c('NFKB1','HNF4A','MAF','IRF7')
netLight <- getSubnetwork(network,receptors = R, TFs = TFs)
t50 <- diffusionMap(receptors = R, TFs = TFs, M= diffusionData$M, network = netLight, nCores = 2, nTicks = 1000) 
sampleSheet <- diffusionData$sampleSheet
dMat <- model.matrix(~0 + as.factor(sampleSheet$group))
colnames(dMat) <- c('ctrl','postResistant','postResponse','preResistant','preResponse')
contrastMatrix <- makeContrasts(preResistant-ctrl,preResistant-preResponse,levels = dMat)
fit <- lmFit(t(t50),dMat)
fit2 <- contrasts.fit(fit,contrasts = contrastMatrix)
decideTests(eBayes(fit2),p.value = 0.01)
