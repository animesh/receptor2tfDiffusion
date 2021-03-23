#' diffusionData
#'
#' Network connectivity between receptors and TFs using diffusion
#'
#' Gene expression data from IBD patients with matching
#' signalling network and candidate receptors-transcription factors.
#'
#' @docType data
#'
#' @format An object of class \code{"list"};
#'
#' @format Variable "M": Data frame containing expression data for 24501 genes and 54 samples. Variable "sampleSheet":
#' Data frame containing information on whether patients responded or resistant to the anti-TNF treatment. 54 samples and 12 phenotype columns.
#'
#' @keywords datasets
#' @name diffusionData
#' @usage data(diffusionData)
#' @references Singh et al. (2021) ....
#'
#' @examples
#' data(diffusionData)
#' receptors <- diffusionData$receptors[1:4]
#' TF <- diffusionData$TFs[1]
#' nodeW <- diffusionData$M[,34]
#' network <- diffusionData$network
#' \donttest{diffusion_curve <- signalOnNetwork(network,nodeW,outputNode = TF,inputNode = receptors)}
NULL
