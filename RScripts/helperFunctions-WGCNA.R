getGOEnrichmentWebGestaltRWGCNA<-function(geneList,enrichDatabase,moduleName,parent)
{
  require(enrichR)
  require(ggplot2)
  require(stringr)
  projectName<-paste0(moduleName,"_",enrichDatabase)
  enrichResult<-NULL
  tryCatch({
    enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                                enrichDatabase=enrichDatabase,#c("pathway_KEGG","geneontology_Biological_Process_noRedundant","pathway_Reactome"),
                                interestGene=geneList,#interestGeneFile=geneFile,
                                interestGeneType="genesymbol", referenceSet="genome", isOutput=TRUE,
                                outputDirectory=paste0(parent,"EnrichmentAnalysis"),
                                projectName=projectName,saveRawGseaResult=TRUE,
                                minNum = 5,maxNum = 2000,sigMethod = "fdr")
    projectName<-stringr::str_replace_all(projectName,pattern = "\\.",replacement = "_")
    #print(projectName)
    title<-paste0("Enriched terms from ",enrichDatabase)
    if(!is.null(enrichResult)){
      g<-getOntologyHistogram(enrichResult,title = title)
      ggsave(paste0(parent,"EnrichmentAnalysis/","Project_",projectName,"/OntologyTerms-Hist.pdf"),width = 6,height = 3,units = "in",scale = 2,g)
      g<-getOntologyDotplot(enrichResult,title = title)
      ggsave(paste0(parent,"EnrichmentAnalysis/","Project_",projectName,"/OntologyTerms-Dotplot.pdf"),width = 6,height = 3,units = "in",scale = 2,g)
    }},
    error=function(error)
    {
      message("Original error message:")
      message(error)
    },finally = function()
    {
      return(enrichResult)
    }
  )
}

extractGeneNames<-function(geneNames)
{
  l<-sapply(geneNames,FUN=function(x){unlist(base::strsplit(x,split = "\\|"))})
  if(length(class(l)) == 2)##For matrix object return
  {
    genesDF<-t(l)
  }else
  {
    genesDF<-do.call(rbind,l)
  }
  l1<-sapply(geneNames,FUN=function(x){unlist(strsplit(x,split = "\\."))})
  if(length(class(l1)) == 2)
  {
    ensemblID<-t(l1)
  }else
  {
    ensemblID<-do.call(rbind,l1)
  }

  return(data.frame(EnsemblID=ensemblID[,1],GeneName=genesDF[,2],row.names = geneNames))
}
