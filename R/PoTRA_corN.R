#' The PoTRA analysis is based on topological ranks of genes in biological pathways. PoTRA can be used to detect pathways involved in disease. The PoTRA package contains one function for creating the PoTRA results object.
#'
#' @param mydata Dataframe that contains rownames consisting of entrez gene identifiers and columns representing normals and samples from gene expression data Row names must represent gene identifiers (entrez). A minimum of 18,000 genes are recommended.
#' @param genelist The object genelist is a dataframe that consists of a single column of entrez gene identifiers (the same as those found in the rownames(mydata)).
#' @param Num.sample.normal The number of normal samples in the mydata dataframe. 
#' @param Num.sample.case The number of case samples in the mydata dataframe.
#' @param Pathway.database The pathway database, such as KEGG, Reactome, Biocarta and PharmGKB. 
#' @param PR.quantile Contains the percentile of PageRank scores as a cutoff for hub genes, the recommended percentile is 0.95.
#'
#' 
#' @return None
#'
#' @examples
#' Results <- PoTRA.corN(mydata=mydata, genelist=genelist, Num.sample.normal=Num.sample.normal, Num.sample.case=Num.sample.case, Pathway.database=Pathway.database, PR.quantile=PR.quantile)
#'
#' Pathway.database (options): 
#' humanKEGG=pathways('hsapiens','kegg')
#' humanReactome=pathways('hsapiens','reactome')
#' humanBiocarta=pathways('hsapiens','biocarta')
#' humanPharmGKB=pathways('hsapiens','pharmgkb')
#'
#' Pathway.database=humanKEGG
#' Pathway.database=humanReactome
#' Pathway.database=humanBiocarta
#' Pathway.database=humanPharmGKB
#'
#' @export
PoTRA.corN <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database,PR.quantile) {
    
    require(BiocGenerics)
    require(graph)
    require(graphite)
    require(igraph)  
    require(org.Hs.eg.db) 
    Fishertest<-c()
    TheNumOfHubGene.normal<-c()
    TheNumOfHubGene.case<-c()
    E.normal<-c()
    E.case<-c()
    length.pathway<-c()
    kstest<-c()
    pathwaynames <- c()
    
    #humanPharmGKB <- pathways("hsapiens", "pharmgkb")
    #humanBiocarta <- pathways("hsapiens", "biocarta")
    #humanKEGG <- pathways("hsapiens", "kegg")
    
    for (x in 1:length(Pathway.database[1:length(Pathway.database)])){
        
        print(x)
        p0 <-Pathway.database[[x]]
        pathwaynames[x] <- p0@title
        p <- convertIdentifiers(p0, "entrez")
        g<-pathwayGraph(p) 
        nodelist<-nodes(g)
        graph.path<-igraph.from.graphNEL(g)    
        ## Plot graph
        #plot(graph.path, edge.arrow.size=.5, vertex.color="gold", vertex.size=5, 
        #vertex.frame.color="gray", vertex.label.color="black", 
        #vertex.label.cex=0.8, vertex.label.dist=2, edge.curved=0.2) 
        
        
        genelist_reformatted <- sprintf('ENTREZID:%s', genelist) 
        length.intersect<-length(intersect(unlist(nodelist),unlist(genelist_reformatted)))
        length.pathway[x]<-length.intersect
        graph.path<-induced_subgraph(graph.path, as.character(intersect(unlist(nodelist),unlist(genelist_reformatted))))
        #plot(graph.path)   
        
        
        if (length.intersect<5){
            next
        }else{
            
            #collect expression data of genes for a specific pathway across normal and tumor samples.
            
            path<-data.frame(matrix(0,length.intersect,(Num.sample.normal+Num.sample.case)))
            a<- c()
            
            for (j in 1:length.intersect){
                a[j]<-intersect(unlist(nodelist),unlist(genelist_reformatted))[j]  
                
                path[j,]<-mydata[which(genelist_reformatted==a[j]),]  #collect expression data of genes for a specific pathway across normal and tumor samples.
            }
            
            ##Construct a gene-gene network for normal samples and calculate PageRank values for each gene in this network.
            
            cor.normal <- apply(path[,1:Num.sample.normal], 1, function(x) { apply(path[,1:Num.sample.normal], 1, function(y) { cor.test(x,y)[[3]] })})
            
            cor.normal<-as.matrix(cor.normal) 
            
            cor.normal.adj<-matrix(p.adjust(cor.normal,method="fdr"),length.intersect,length.intersect)
            
            cor.normal.adj[ cor.normal.adj > 0.05 ] <- 0
            cor.normal.adj[ is.na(cor.normal.adj)] <- 0
            diag(cor.normal.adj) <- 0
            graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
            E.normal[x]<-length(E(graph.normal))
            
            PR.normal<-page.rank(graph.normal,directed=FALSE)$vector
            
            ##Construct a gene-gene network for tumor samples and calculate PageRank values for each gene in this network.
            
            cor.case <- apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(x) { apply(path[,(Num.sample.normal+1):(Num.sample.normal+Num.sample.case)], 1, function(y) { cor.test(x,y)[[3]] })})
            
            cor.case<-as.matrix(cor.case) 
            cor.case.adj<-matrix(p.adjust(cor.case,method="fdr"),length.intersect,length.intersect)
            
            cor.case.adj[ cor.case.adj > 0.05 ] <- 0
            cor.case.adj[ is.na(cor.case.adj)] <- 0
            diag(cor.case.adj) <- 0
            graph.case<-graph.adjacency(cor.case.adj,weighted=TRUE,mode="undirected")
            E.case[x]<-length(E(graph.case))
            
            PR.case<-page.rank(graph.case,directed=FALSE)$vector
            
            
            matrix.HCC<- matrix("NA",2*length.intersect,2)
            rownames(matrix.HCC)<-as.character(c(PR.normal,PR.case))
            colnames(matrix.HCC)<-c("Disease_status","PageRank")
            
            matrix.HCC[,1]<-c(rep("Normal",length.intersect), rep("Cancer",length.intersect))
            
            loc.largePR<-which(as.numeric(rownames(matrix.HCC))>=quantile(PR.normal,PR.quantile))
            loc.smallPR<-which(as.numeric(rownames(matrix.HCC))<quantile(PR.normal,PR.quantile))
            
            matrix.HCC[loc.largePR,2]<-"large_PageRank"
            matrix.HCC[loc.smallPR,2]<-"small_PageRank"
            
            table.HCC<-list(1,2)
            names(table.HCC)<-c("Disease_status","PageRank")
            
            table.HCC$Disease_status<-matrix("NA",2*length.intersect,2)
            table.HCC$PageRank<-matrix("NA",2*length.intersect,2)
            
            table.HCC$Disease_status<-matrix.HCC[,1]
            table.HCC$PageRank<-matrix.HCC[,2]
            
            cont.HCC<-table(table.HCC$Disease_status,table.HCC$PageRank)
            TheNumOfHubGene.normal[x]<-cont.HCC[2]
            TheNumOfHubGene.case[x]<-cont.HCC[1]
            
            if (dim(cont.HCC)[1]!=dim(cont.HCC)[2]){
                Fishertest[x]<-1
            }else{
                Fishertest[x]<-fisher.test(cont.HCC)$p.value
            }
            tryCatch(
                kstest[x]<-ks.test(PR.normal,PR.case)$p.value, 
                error=function(e) print("."))  
            
            
            if (E.normal[x]<E.case[x]){
                Fishertest[x]<-1
                kstest[x]<-1
            }else{
                Fishertest[x]<-Fishertest[x]
                kstest[x]<-kstest[x]
            }
            
        }
    }
    return(list(Fishertest.p.value=Fishertest,KStest.p.value=kstest,LengthOfPathway=length.pathway,TheNumberOfHubGenes.normal=TheNumOfHubGene.normal,TheNumberOfHubGenes.case=TheNumOfHubGene.case,TheNumberOfEdges.normal=E.normal,TheNumberOfEdges.case=E.case,PathwayName=pathwaynames))
}
