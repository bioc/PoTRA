#' PoTRA Package
#'
#' \description{The PoTRA analysis is based on topological ranks of genes in biological pathways. PoTRA can be used to detect pathways involved in disease. The PoTRA package contains one function for creating the PoTRA results object.}
#'
#'
#' \usage{Results <- PoTRA.corN(mydata=mydata, genelist=genelist, Num.sample.normal=Num.sample.normal, Num.sample.case=Num.sample.case, Pathway.database=Pathway.database, PR.quantile=PR.quantile)}
#'
#
#' \arguments{
#'   \item{mydata}{Dataframe that contains rownames consisting of entrez gene identifiers and columns representing normals and samples from gene expression data.}
#'   \item{genelist}{Dataframe that consists of a single column of entrez gene identifiers (the same as those found in the rownames(mydata)).}
#'   \item{Num.sample.normal}{Represents normal samples in the mydata dataframe.}
#'   \item{Num.sample.case}{Represents case samples in the mydata dataframe.}
#'   \item{Pathway.database}{Object contains gene lists from KEGG, Reactome, PathBank, Panther and PharmGKB databases, and these are made available through the devel graphite library.}
#'   \item{PR.quantile}{Object contains the percentile of PageRank scores as a cutoff for hub genes, the recommended percentile is 0.95.}
#' }
#'
#' \details{ 
#'
#' Pathway.database = humanKEGG, etc (see below for additional database definitions.)
#'
#'  \item{KEGG}{humanKEGG = pathways('hsapiens','kegg')}
#'  \item{Reactome}{humanReactome = pathways('hsapiens','reactome')}
#'  \item{Biocarta}{humanBiocarta = pathways('hsapiens','pathbank')}
#'  \item{PharmGKB}{humanPharmGKB = pathways('hsapiens','pharmgkb')}
#' }  

#' @docType package
#' @author Chaoxing Li, Li Liu, Valentin Dinu \email{valentin.dinu@asu.edu}
#'
#' @name PoTRA
#' @citation Li C, Liu L, Dinu V. 2018. Pathways of topological rank analysis (PoTRA): a novel method to detect pathways involved in hepatocellular carcinoma. PeerJ 6:e4571 https://doi.org/10.7717/peerj.4571
NULL
