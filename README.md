# receptor2tfDiffusion - An R package for patient-specific diffusion modelling

## Diffusion based signalling model in R

Diffusion modelling is a conceptually simple way to quantify network connectivity. Network connectivity between cell surface receptors and transcription factors may be interesting for studying individual differences in disease progression or drug response especially in autoimmune diseases. Here we look at an example data set from inflammatory bowel disease (IBD) and use gene expression to make individualized statements about network connectivity.

### Example: patient-specific network diffusion analysis in Ulcerative colitis 

Install the 'receptor2tfDiffusion' package using the 'devtools' library. The igraph package will also be helpful.

```
install.packages("devtools")
library(devtools)
devtools::install_github('Amy3100/receptor2tfDiffusion',force = T)
```

For the analysis, following packages are required:

```
library(igraph)
library(limma)
library(devtools)
library(receptor2tfDiffusion)
library(doParallel)
```
The package come with a pre normalized dataset as well as a network model and a list of transcription factors (TFs) and receptors. This contains an object `M` with gene expression data for UC patients and controls and an object `sampleSheet` contains metadata about patients pre and post treatment along with their treatment outcomes.

We need to load a toy dataset in R:

```
data(diffusionData)
names(diffusionData)
```

The 'diffusionData' object consists of a network, a 'data frame' with a 'from' and 'to' list of gene symbols which describes the connectivity of the signalling network believed to be relevant for IBD. The package simulates the spread of signal from a receptor. By monitoring a transcription factor over time the speed degree if connectivity to the receptor can be quantified. 

```
network <- diffusionData$network 
```


The function 'signalOnNetwork' takes in a network, gene expression values for the genes making up the network and names of the input (receptor) and output (TF) nodes. Other kinds of genes or a proteins can be used but conceptually going from receptors to TFs makes the most sense. Let us look at how signal travels from the receptor TNFRSF1B to the inflammation controlling TF FLI1. To do this, we run the function 'diffusion_curve'. Node weights are given from gene expression data. In this case we choose a sample at random e.g. sample No. 24.  

```
nodeW <- diffusionData$M[unique(c(network[,1],network[,2])),24]
diffusion_curve <- signalOnNetwork(network,nodeW,outputNode = 'FLI1',inputNode = 'TNFRSF1B',n = 2000)
```

The diffusion curve is a R 'list' object with two components 'Enode' is the signal increasing over time at the output node, and signal.

```
hist(diffusion_curve$signal, main ="Final signal intensity", xlab="Signal" )
```
We visualize diffusion plot calculated for given receptor 'TNFRSF1B' and TF 'FLI1' with 2000 iterations :

```
#par(mar=c(0,0,0,0))
plot(diffusion_curve$Enode, main= "Diffusion curve" , ylab="signal at FLI1", xlab="time ticks")
```

The 'signal' gives the amount of estimated signal that has accumulated in every node in the network. This is useful during debugging start of a project to see if the signal has stabilized globally or if large amounts of the network remains without signal due to e.g. a network connectivity issue.

The diffusion calculation is time consuming. For mapping large sets of receptors and TFs we recommend using a server with multiple cores. We use the library "doParallel" (mentioned in list of required libraries) to facilitate the large inputs in the function "diffusionMap".

```
R <- c('TNFRSF11B','OSMR','IL6R')
TFs <- c('NFKB1','HNF4A','MAF','IRF7')
netLight <- getSubnetwork(network,receptors = R, TFs = TFs)
```
The "netLight" dataframe now contains a smaller network relevant for the three receptors in R and the four TFs in 'TFs'

We can the run a full network connectivity mapping between all these candidate receptor and TFs. Since keeping the full curves is tedious and memory intensive. And also hard for further analysis. We represent each connectivity with one number: 't50'. This is the number of time ticks it takes for a curve to reach 1/2 its maximum. 

We can calculate t50 using the "diffusionMap" function. Feel free to change the number of cores seed to reflect the system you are on. May be grab a cup of coffee if it turns tedious for your machine. 

```
t50 <- diffusionMap(receptors = R, TFs = TFs, M= diffusionData$M, network = netLight, nCores = 2, nTicks = 1000) 
```
The t50s can the be used for data analysis or machine learning. The matrix is samples by TF receptor pairs so you can use it directly in PCA (e.g. "prcomp") but it needs to be transposed before use in tool like "limma" which expects variables in the rows.

```
sampleSheet <- diffusionData$sampleSheet
dMat <- model.matrix(~0 + as.factor(sampleSheet$group))
colnames(dMat) <- c('ctrl','postResistant','postResponse','preResistant','preResponse')
contrastMatrix <- makeContrasts(preResistant-ctrl,preResistant-preResponse,levels = dMat)
fit <- lmFit(t(t50),dMat)
fit2 <- contrasts.fit(fit,contrasts = contrastMatrix)
decideTests(eBayes(fit2),p.value = 0.01)
```
As we see there are receptor-TF pairs that are strongly related to both drug resistance and inflammation.

## Session Info
sessionInfo()

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)
