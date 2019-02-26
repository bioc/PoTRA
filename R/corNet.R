#' Create the results object for PoTRA
#' More detailed description.
#'
#' @param mydata Data Frame
#' @param genelist Data Frame
#' @param Num.sample.normal Integer
#' @param Num.sample.tumor Integer
#' @param Pathway.database s4
#'
#' @examples
#' library(BiocGenerics)
#' library(graph)
#' library(graphite)
#' library(igraph)
#' humanKEGG <- pathways("hsapiens", "kegg")
#' load('./data/PoTRA-vignette.RData')
#' results <- PoTRA.corN(mydata=mydata,genelist=genelist,Num.sample.normal=Num.sample.normal,Num.sample.case=Num.sample.case,Pathway.database=Pathway.database,PR.quantile=PR.quantile)
#'

#' @export

