#BiocManager::install("STRINGdb")
analyzePPINetwork<-function(geneList,plotName="",ppiEdgeThreshold=400,ppiPValue=0.01,minGenes=5,get_largest_connected_subnet = TRUE)
{
  require(STRINGdb)
  require(igraph)
  require(influential)

  string_db <- STRINGdb$new(version="11", species=9606,score_threshold=ppiEdgeThreshold, input_directory="")
  geneMapped <- string_db$map(data.frame(gene=geneList), "gene", removeUnmappedRows = TRUE)
  subnetwork <- string_db$get_subnetwork(geneMapped$STRING_id)
  g<-components(subnetwork)
  largestCluster<-names(sort(table(g$membership),decreasing = T)[1])
  genesInCluster<-igraph::groups(g)[[largestCluster]]
  enrichment<-string_db$get_ppi_enrichment(genesInCluster)
  #print(length(genesInCluster))
  #print(enrichment$enrichment)

  if(length(genesInCluster) >= minGenes)
  {
    if(is.null(enrichment$enrichment) || enrichment$enrichment <= ppiPValue)
    {
      if(get_largest_connected_subnet)
      {
        subnetwork1 <- string_db$get_subnetwork(genesInCluster)
      }else
      {
        subnetwork1 <- subnetwork
      }
      b<-influential::betweenness(subnetwork1)
      d<-influential::degree(subnetwork1)
      n1<-influential::neighborhood.connectivity(subnetwork1)
      lh.index <- influential::lh_index(subnetwork1)
      cr<-influential::clusterRank(subnetwork1)
      ci<-influential::collective.influence(subnetwork1)
      # Hub genes based on network centralities
      MyData <- data.frame(DC=d,NC=n1,BC=b,CR=cr,LH_index=lh.index,CI=ci)
      tryCatch({
        My.vertices.IVI <- ivi.from.indices(DC = MyData$DC,
                                            NC = MyData$NC,
                                            BC = MyData$BC,CR=MyData$CR,LH_index = MyData$LH_index,CI=MyData$CI)
        IVIScores<-data.frame(STRING_id=row.names(MyData),IVIScores=My.vertices.IVI)
        topGenes<-merge(x=IVIScores,y=geneMapped,by="STRING_id",all.x=TRUE,all.y=FALSE)
        colnames(MyData)<-c("Degree","Neighborhood.Connectivity","Betweenness","ClusterRank","LH_index","Collective.Influence")
        MyData$STRING_id<-row.names(MyData)
        topGenes1<-merge(x=topGenes,y=MyData,by="STRING_id",all.x=TRUE,all.y=FALSE)
        topGenes1<-topGenes1[order(topGenes1$IVIScores,decreasing = T),]
        print(dim(topGenes1))
        #ppi_cluster_pngs[i] <- file.path(out_dir,paste0("ppi_cluster_",i,".png"))
        png(plotName,width=8,height=8,units = "in", res=600)
        string_db$plot_network(genesInCluster, add_link=F, add_summary=F)
        dev.off()
        return(topGenes1)
      },
      error=function(error)
      {
        message("Original error message:")
        message(error)
        return(NULL)
      },finally = function()
      {
        png(plotName,width=8,height=8,units = "in", res=600)
        string_db$plot_network(genesInCluster, add_link=F, add_summary=F)
        dev.off()
      })
    }
  }else
    {
      print("No Enrichment of PPI network")
      return(NULL)
    }
}

parent<-"/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/WGCNA/datKME/"
module<-c("brown","yellow","red","salmon")
for(color in module)
{
  genes<-scan(paste0(parent,"module_MM.",color,".txt"),character())
  hubGenes<-analyzePPINetwork(genes,paste0(parent,"Module-datKME",color,".STRING.PPI.700.png"),ppiEdgeThreshold = 700,get_largest_connected_subnet = FALSE)
  write.table(hubGenes,paste0(parent,"Module-DatKME",color,".hubGenes.PPI.700.txt"),quote = F,row.names = F,sep = "\t")
}
