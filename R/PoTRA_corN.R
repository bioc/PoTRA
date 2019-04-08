#' @export
#'

PoTRA.corN <- function(mydata,genelist,Num.sample.normal,Num.sample.case,Pathway.database, PR.quantile) {
    
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

for (x in seq_along(Pathway.database[seq_along(Pathway.database)])){       

    print(x)

    p0 <-Pathway.database[[x]]
    pathwaynames[x] <- p0@title
    p <- convertIdentifiers(p0, "entrez")
    g<-pathwayGraph(p) 
    nodelist<-nodes(g)
    graph.path<-igraph.from.graphNEL(g)          
    genelist_reformatted <- sprintf('ENTREZID:%s', genelist) 

    length.intersect<-length(intersect(unlist(nodelist),
        unlist(genelist_reformatted)))

    length.pathway[x]<-length.intersect
    graph.path<-induced_subgraph(graph.path, 
        as.character(intersect(unlist(nodelist),
            unlist(genelist_reformatted))))

    if (length.intersect<5){
        next
        }
    else{

        #########################################################    
        #collect expression data of genes for a specific pathway# 
        #across normal and tumor samples.                       #  
        #########################################################    
            sum_samples = Num.sample.normal+Num.sample.case 
            path<-data.frame(matrix(0, length.intersect, sum_samples))

            a<- c()

            for (j in seq_len(length.intersect)){
                a[j]<-intersect(unlist(nodelist),unlist(genelist_reformatted))[j]  
		#################################################
		#collect expression data of genes for a specific# 
		#pathway across normal and tumor samples.       #
                #################################################
		path[j,]<-mydata[which(genelist_reformatted==a[j]),]  
            }

	    ################################################################
            #Construct a gene-gene network for normal samples and calculate# 
	    #PageRank values for each gene in this network.	           #
            ################################################################ 
            len_samp_norm = seq_len(Num.sample.normal)
	   
	    cor.normal_1<-apply(path[,len_samp_norm], 1, function(y) {
	        cor.test(x,y)[[3]] 
            })

	    cor.normal<-apply(path[,len_samp_norm],1, function(x) {
	        cor.normal_1
	    })

            cor.normal<-as.matrix(cor.normal) 

	    stat_adj_norm = p.adjust(cor.normal,method="fdr")
            cor.normal.adj<-matrix(stat_adj_norm,length.intersect,length.intersect)

            cor.normal.adj[ cor.normal.adj > 0.05 ] <- 0
            cor.normal.adj[ is.na(cor.normal.adj)] <- 0
            diag(cor.normal.adj) <- 0
            graph.normal<-graph.adjacency(cor.normal.adj,weighted=TRUE,mode="undirected")
            E.normal[x]<-length(E(graph.normal))

            PR.normal<-page.rank(graph.normal,directed=FALSE)$vector

	    ###############################################################
            #Construct a gene-gene network for tumor samples and calculate# 
            #PageRank values for each gene in this network.		  #
            ###############################################################
	    range_start = Num.sample.normal+1
	    range_end = Num.sample.normal+Num.sample.case
	    range = range_start:range_end

	    cor.case_1<-apply(path[,range], 1, function(y) {
                cor.test(x,y)[[3]] 
	    })	

	    cor.case <- apply(path[,range], 1, function(x) {
		cor.case_1
            })

	    cor.case<-as.matrix(cor.case) 

	    stat_adj_case = p.adjust(cor.case,method="fdr")
	    cor.case.adj<-matrix(stat_adj_case,length.intersect,length.intersect)

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

	    quant=quantile(PR.normal,PR.quantile)

            loc.largePR<-which(as.numeric(rownames(matrix.HCC))>=quant)
            loc.smallPR<-which(as.numeric(rownames(matrix.HCC))<quant)

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
return(list(Fishertest.p.value=Fishertest,
    KStest.p.value=kstest,
    LengthOfPathway=length.pathway,
    TheNumberOfHubGenes.normal=TheNumOfHubGene.normal,
    TheNumberOfHubGenes.case=TheNumOfHubGene.case,
    TheNumberOfEdges.normal=E.normal,
    TheNumberOfEdges.case=E.case,
    PathwayName=pathwaynames))
}
