###Calculate Cancer-specific genes using the TPM values
library(TissueEnrich)
library(plyr)
parent<-"~/myPART/AllSamplesPipeliner/"
fileName<-"RSEM.genes.TPM.all_samples.txt"

tpmValues<-read.table(paste0(parent,fileName),header = TRUE,sep = "\t",row.names = 1)
row.names(tpmValues)<-paste0(row.names(tpmValues),"|",tpmValues$GeneName)
#tpmValues<-tpmValues[row.names(allcounts),]
tpmValues1<-as.data.frame(t(sapply(unique(tpmValues$GeneName),FUN = function(col)
  colMeans(tpmValues[tpmValues$GeneName == col,c(seq(2,(ncol(tpmValues))))]))))
row.names(tpmValues1)<-unique(tpmValues$GeneName)
#colnames(tpmValues)<-gsub("^X", "",  colnames(tpmValues))
colnames(tpmValues1) <- gsub("^Sample_","",colnames(tpmValues1))
colnames(tpmValues1) <- unlist(lapply(colnames(tpmValues1),function(x){i<-unlist(gregexpr(pattern ="_S",x));return(substr(x,0,i[length(i)]-1));}))
tpmValuesFull<-tpmValues1
#tpmValues<-tpmValues1[apply(tpmValues1,MARGIN = 1, function(x) any(x >= 1)),]

features<-read.table(paste0(parent,"../features.txt"),header = TRUE,sep = "\t",row.names = 1)
##Filter Cancer with 1 sample, normals, NET, and ACC
#tableTypes<-table(features$New.Diagnosis)
#filtNames<-c(names(tableTypes[tableTypes == 1]),"Neuroendocrine","Adrenocortical carcinoma","Normal")
#samplesToCompare<-row.names(features[features$Subtype %in% c("Endocrine") | features$New.Diagnosis %in% c("Normal"),])
features<-features[features$Subtype %in% c("Endocrine") | features$New.Diagnosis %in% c("Normal"),]
###########################################################
NETSamples<-features[features$New.Diagnosis == "Neuroendocrine",]
features[features$New.Diagnosis == "Neuroendocrine",]$New.Diagnosis<-paste0(NETSamples$New.Diagnosis,"-",NETSamples$PrimaryCancerTissue)
NormalSamples<-features[features$New.Diagnosis == "Normal",]
features[features$New.Diagnosis == "Normal",]$New.Diagnosis<-paste0(NormalSamples$New.Diagnosis,"-",NormalSamples$Tissue)
features<-features[features$New.Diagnosis != "Normal-Human Universal Reference Total RNA",]
###Removing the samples with bad quality (low TIN) from the analysis
features<-features[!(row.names(features) %in% c("47_RT00085_X007R_resent","76_RT00105_X336R")),]
features<-features[order(features$New.Diagnosis),]

######Check the cancer names in latest file
# f<-read.xlsx("~/myPART/NGS clinical data for Ashish 4-26-21v2-no_name.xlsx")
# RTCancerMapping<-f[,c("Subject.Code","cancer/tumor.name.-.MyPART.Tumor.Pathologies")]
tpmValuesFull1<-tpmValuesFull[,row.names(features)]
tpmValues<-tpmValuesFull[,row.names(features)]
tpmValues<-tpmValues[apply(tpmValues,MARGIN = 1, function(x) any(x >= 1)),]

#tpmValues<-vsdNormalizedCounts
tTpmValues<-data.frame(t(tpmValues))
#tTpmValues<-tTpmValues[,c(1:10)]
tTpmValues$group<-features$New.Diagnosis
tTpmValuesMean <- as.data.frame(t(sapply(unique(tTpmValues$group),FUN = function(col)
  colMeans(tTpmValues[tTpmValues$group == col,c(seq(1,(ncol(tTpmValues) - 1)))]))))
se<-SummarizedExperiment(assays = SimpleList(as.matrix(t(tTpmValuesMean))),
                         rowData = colnames(tTpmValuesMean),colData = row.names(tTpmValuesMean))
output<-teGeneRetrieval(se,expressedGeneThreshold = 1)
cellSpecificGenesOutput<-data.frame(assay(output))
cellSpecificGenesOutput %>% filter(Group %in% c("Tissue-Enriched")) %>% group_by(Tissue) %>% tally()
save(cellSpecificGenesOutput,tTpmValuesMean,tpmValues,tpmValuesFull1,features,file= paste0("~/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/","cancerSpecificGenes-TPM.rda"))
#save(cellSpecificGenesOutput,tTpmValuesMean,tpmValues,features,file= paste0(parent,"cancerSpecificGenes-TPM.rda"))


###Mapping the Cancer-Specific genes to the drug data from DrugBank and NCI-60
library("rDGIdb")
parent<-"/Users/jaina13/myPART/AllSamplesPipeliner/"
#o<-load(file= paste0(parent,"cancerSpecificGenes-TPM.rda"))
o<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/","cancerSpecificGenes-TPM.rda"))
drugTargets<-read.table(paste0(parent,"CancerSpecificGenes/DrugBankTargets-Format.tsv"),header = TRUE,sep = "\t")
drugBankGeneMapping<-read.table(paste0(parent,"CancerSpecificGenes/genesDrugBank.tsv"),header = TRUE,sep = "\t")
drugTargetsMap<-merge(x=drugTargets,y=drugBankGeneMapping,by.x="targetID",by.y="gene_claim_name",all.x=TRUE,all.y=FALSE)
#teGenes<-getTissueSpecificGenes(rdaFilePath = "~/myPART/TissueEnrichCombineExpression.rda")
proteinCodingGenes<-read.table(paste0("~/myPART/FusionAnalysis/hg38.proteinCodingGenes.Ensembl.txt"),sep = "\t",header = T)
l<-list()
for(tissue in unique(features$New.Diagnosis))
{
  if(!grepl("Normal",tissue))
  {
    #c<-cellSpecificGenesOutput %>% dplyr::filter(Group %in% c("Tissue-Enriched","Group-Enriched")) #%>% dplyr::filter(Tissue == tissue)
    c<-cellSpecificGenesOutput %>% dplyr::filter(Group %in% c("Tissue-Enriched")) %>% dplyr::filter(Tissue == tissue)
    #normalGenes<-unique(c[grepl("Normal",c$Tissue),]$Gene)
    #c<-c[!(c$Gene %in% normalGenes),] %>% dplyr::filter(Tissue == tissue)

  #resultFilter <- queryDGIdb(c$Gene,sourceDatabases = NULL,geneCategories = c("CLINICALLY ACTIONABLE","TUMOR SUPPRESSOR"),interactionTypes = c("inhibitor","antagonist","suppressor"))
  #geneSummary<-byGene(resultFilter)
  #geneWithDrugs<-geneSummary[geneSummary$DistinctDrugCount > 0,]
  c2<-sapply(c$Gene,FUN = function(x){
    return(paste0(drugTargetsMap[drugTargetsMap$gene_name %in% x,2],collapse = ","))
  })
  c$drugBank<-c2
  c<-c[c$Gene %in% proteinCodingGenes$Gene.name,]
  l[[str_replace_all(tissue,pattern = "[ /]",replacement = ".")]]<-c
  #g<-getGOEnrichmentWebGestaltRWGCNA(c$Gene,"geneontology_Biological_Process_noRedundant",str_replace_all(tissue,pattern = "[ /]",replacement = "."),paste0(parent,"EndocrineSubgroupResults/CancerSpecificGenes/"))
  #g<- getGOEnrichmentWebGestaltRWGCNA(c$Gene,"pathway_Reactome",str_replace_all(tissue,pattern = "[ /]",replacement = "."),paste0(parent,"EndocrineSubgroupResults/CancerSpecificGenes/"))
  }
}
hs <- createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF", fontSize = 12,fontName = "Arial Narrow", fgFill = "#4F80BD")
write.xlsx(l, paste0(parent,"/EndocrineSubgroupResults/CancerSpecificGenes/CancerSpecificGenes.xlsx"), colWidths = c(NA, "auto", "auto"),colNames = TRUE, borders = "rows", headerStyle = hs)

##Compare drugs targetting csGenes across cancers
o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
csGenes<-do.call(rbind,l)
csGenes$isDrug<-"No"
csGenes[csGenes$drugBank != "",]$isDrug<-"Yes"
csGenes<-csGenes[csGenes$isDrug =="Yes",]
csTargetDrugs<-lapply(unique(csGenes$Tissue), function(x){
  drugs<-unlist(lapply(csGenes[csGenes$Tissue == x,]$drugBank, function(x1){return(unlist(stringr::str_split(x1,",")))}))
  return(unique(drugs))
})
names(csTargetDrugs)<-unique(csGenes$Tissue)
Reduce(intersect,csTargetDrugs)

##BarPlot for the total Number of Tissue-specific genes
o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
csGenes<-do.call(rbind,l)
csGenes$isDrug<-"No"
csGenes[csGenes$drugBank != "",]$isDrug<-"Yes"
csGenes<-csGenes[csGenes$isDrug =="Yes",]
g<-ggplot(csGenes, aes(fill=isDrug, x=(Tissue))) +
#g<-ggplot(csGenes, aes(x=(Tissue))) +
  geom_bar(position="stack") +
  scale_fill_manual(values = c("red","blue")) +
  #labs(x="", y = "xCell Enrichment Score") +
  labs(x="", y = "# cancer-specific genes") +
  #scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  #ggtitle("xCell Enrichment of cells in ACC RNA-Seq data") +
  #theme_ipsum() +
  xlab("")
ggsave(paste0(parent,"cancerSpecificGenesBarplot-DrugTarget-AllSamples.pdf"),plot = g,height = 10,width = 8)

##Compare cs genes
csGenesFilt<-csGenes[csGenes$Tissue %in% unique(csGenesEndocrine$Tissue),]
commonCSGenes<-intersect(csGenesFilt[csGenesFilt$drugBank != "",]$Gene,csGenesEndocrine[csGenesEndocrine$drugBank != "",]$Gene)
csEndocrineSpecific<-setdiff(csGenesEndocrine[csGenesEndocrine$drugBank != "",]$Gene,csGenesFilt[csGenesFilt$drugBank != "",]$Gene)
csAllSpecific<-setdiff(csGenesFilt[csGenesFilt$drugBank != "",]$Gene,csGenesEndocrine[csGenesEndocrine$drugBank != "",]$Gene)

d<-rbind(csGenesFilt[csGenesFilt$drugBank != "",],csGenesEndocrine[csGenesEndocrine$drugBank != "",])
##Filtering the csAllSpecific Genes as they are calculated with the bad NET-SI sample
d<-d[!(d$Gene %in% csAllSpecific), ]
d$Common<-"F"#d$uniqVar %in% commonVars
d$Common[d$Gene %in% commonCSGenes] <- "Common"
d$Common[d$Gene %in% csEndocrineSpecific] <- "Endocrine"
d$Common[d$Gene %in% csAllSpecific] <- "All"

g<-ggplot(d, aes(fill=Common,x=Tissue)) +
  geom_histogram(stat="count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1))
ggsave("~/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGeneCountComp-DrugTargets.pdf",width = 12,height = 6,g)

###Expression heatmap of the Cancer-specific genes
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(randomcoloR)

o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
featuresFilt<-features
hm_data<-log2(as.matrix(tpmValues[unique(csGenes$Gene),])+1)
column_annotations = HeatmapAnnotation(df = featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")])
col_age = colorRamp2(c(0,max(featuresFilt$Age,na.rm = T)), c("white","green"))
#col_rin = colorRamp2(c(0,max(featuresFilt$RIN)), c("white","yellow"))
n <- distinctColorPalette(length(levels(as.factor(featuresFilt$Tissue))))
names(n)<-(levels(as.factor(featuresFilt$Tissue)))
col_tissue <-n[featuresFilt$Tissue]

n1 <- distinctColorPalette(length(levels(as.factor(featuresFilt$New.Diagnosis))))
names(n1)<-(levels(as.factor(featuresFilt$New.Diagnosis)))
col_disease <-n1[featuresFilt$New.Diagnosis]

df<-featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")]
colnames(df)<-c("Tissue","Diagnosis","Age","Sex","Race","TumorSite")

hm_data1<-hm_data
colnames(hm_data1)<-featuresFilt[colnames(hm_data),]$New.Diagnosis##row.names(features)#
dend <- as.dendrogram(hclust(dist(t(hm_data1),method = "euclidean"),method = "complete"))#"ward.D2"))
groupCodes <-col_disease[featuresFilt[colnames(hm_data),]$New.Diagnosis]
labels_colors(dend) <- groupCodes[order.dendrogram(dend)]
pdf(paste0(parent,"/PCsSampleClustering-CancerSpecificGenes-complete.pdf"),width = 14,height = 14)
par(mar = c(2, 2, 2, 25))
plot(dend,horiz = TRUE)
dev.off()

#col_tissue <- distinctColorPalette()
ha = HeatmapAnnotation(
  #df = features[,c(5,6,7,13,15)],
  #df = featuresFilt[,c(3,4,5,6,13,7)],
  df = df,
  col = list(Tissue = col_tissue,
             Diagnosis = col_disease,#c("Normal" = "black","Adrenocortical carcinoma" = "red"),#col_cluster,#c("Adrenocortical carcinoma" = "red", "Normal" = "black", "Carcinoid tumor" = "blue","Gastrointestinal stromal tumor"="brown"),
             Age = col_age,
             Sex = c("Male"="green","Female"="brown"),
             #RIN = col_rin,
             Race =c("White"="#B35806","Unknown"="#FDBC6B","Other"="white","Black or African American"="#E3E4EF","Native Hawaiian or Other Pacific Islander"="#8D81B6","Asian"="blue"),
             TumorSite=c("metastasis"="green","primary site"="brown","recurrence"="white")
             #cluster = col_cluster
  )
)

cheatmap <- ComplexHeatmap::Heatmap(hm_data,
                                    col=colorRamp2(seq(-max(abs(hm_data), na.rm = T), max(abs(hm_data), na.rm = T), length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    bottom_annotation = ha,#column_annotations,
                                    show_column_names=F,
                                    show_row_names=F,
                                    column_names_rot = 45,
                                    cluster_rows = FALSE,
                                    show_heatmap_legend = F) # Turning off to control the placement
pdf(paste0(parent,"PCsSampleClustering-Subtypes-All-L.pdf"))
draw(cheatmap, show_annotation_legend = TRUE)
dev.off()

##Make expression boxplots
# o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/SampleTPMValues.rda"))
o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
tpmValuesFull<-tpmValues
genes<-c("GNRHR","SPINK6","TPH1","CELA2A","KCNJ6","F5","GHRHR")
genes<-c("RET")
for(gene in genes)
{
# row.names(mapping)<-mapping$gene_id
#row.names(tpmValues)<-paste0(row.names(tpmValues),"_",mapping$GeneName)
# tpmValues$GeneName<-mapping$GeneName
tpmvaluesGene<-tpmValuesFull1[gene,]
d<-reshape2::melt(tpmvaluesGene)
f<-features[,c("RTNo","PrimaryCancerTissue","New.Diagnosis","TumorSite")]
f$New.Diagnosis[f$New.Diagnosis == "Adrenocortical carcinoma-recurrence"]<-"Adrenocortical carcinoma-metastasis"
f$samples<-row.names(f)
# d<-d[d$GeneName == gene,]
final_data <- sqldf::sqldf("SELECT * FROM d LEFT OUTER JOIN f where d.variable == samples")

g<-ggplot(final_data,aes(x=New.Diagnosis,y=value))+
  theme_classic(base_size = 10) +
  theme(text=element_text(face = "bold"),axis.text = element_text(size = 15),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.2)) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
  geom_point()+
  labs(x="", y = "TPM Values") +
  ggtitle(paste0("Expression of ",gene))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  coord_flip()
#ggsave(paste0(parent,gene,"-expBoxplot.pdf"),plot = g,height = 10,width = 8)
}
#finalT<-final_data1[,c(7,4,5,6,3)]
#colnames(finalT)<-c(colnames(finalT)[1:4],"TPM-Values")
#write.table(finalT,paste0(parent,gene,"-expression.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)

tpmGenes<-tpmValuesFull[c("SMARCB1","NF2"),]
tpmGenes<-tpmGenes[,(row.names(features[features$RTNo %in% mutRTMap[mutRTMap$Gene == "SMARCB1" | mutRTMap$Gene == "NF2",]$RTNo,]))]
cor(unlist(tpmGenes[1,]),unlist(tpmGenes[2,]))
tpmGenes<-tpmGenes[,(row.names(features[grepl("Adrenocortical carcinoma",features$New.Diagnosis),]))]
cor(unlist(tpmGenes[1,]),unlist(tpmGenes[2,]))


####################Plot Gene expression Dashboard##########
genes<-c("GNRHR","SPINK6","KCNJ8","MMP9","CALCA","CALCB","CALCR","GAST","GHRHR","AVPR1B","TPH1","DAPK2","IL5","RET","SLC6A2","ASCL1","NOTCH1","NEUROD1")
#genes<-c("GNRHR","ASCL1")
geneExp<-getGeneExpressionBoxplot(genes=genes)
depmapCrispr<-getGeneEffectsInDepMap(genes=genes)
# i2dashboard(
#   title = "Cancer-specific genes",
#   author = "",
#   interactive = TRUE,
#   theme = "yeti")->dashboard
# datadir(dashboard) <- parent
plotsList<-list()
for(gene in genes)
{
  p<-plotly::subplot(ggplotly(geneExp[[gene]]),depmapCrispr[[gene]]$combined,nrows = 1,widths = c(0.4,0.6),heights = c(1),titleX=TRUE,titleY = TRUE,margin = 0.08)
  plotsList[[gene]]<-p
  # dashboard %<>% i2dash::add_page(
  #   page = gene,
  #   title = gene,
  #   layout="2x2_grid",
  #   menu = NULL)
  # dashboard %<>%
  #   i2dash::add_component(ggplotly(geneExp[[gene]]),
  #                         page = gene,
  #                         title = "GeneExpressionData") %>%
  #   i2dash::add_component((depmapCrispr[[gene]]$combined),
  #                         page = gene,
  #                         title = "Dependenacy Map using DepMap CRISPR-Cas9 data")
}
getMAFDashboard(MAFfilePath = NULL,plotList = plotsList,outputFileTitle = "Cancer Specific Genes",outputFilePath=paste0(parent,"/EndocrineSubgroupResults/"),outputFileName = "CancerSpecificGenes.html")

# dashboard %>% assemble(file = paste0(parent,"MyDashboard.Rmd"), pages = genes)
# rmarkdown::render(paste0(parent,"MyDashboard.Rmd"),
#                   output_format="all", output_file=paste0("CancerSpecificGenes-Dashboard.html"),
#                   output_dir = parent,
#                   intermediates_dir = parent)

getGeneExpressionBoxplot<-function(genes=NULL)
{

  o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
  tpmValuesFull<-tpmValues
  #genes<-c("GNRHR","SPINK6","TPH1","CELA2A","KCNJ6","F5","GHRHR")
  #genes<-c("RET")
  print(setdiff(genes,row.names(tpmValuesFull)))
  plotsLists<-list()
  for(gene in genes)
  {
    # row.names(mapping)<-mapping$gene_id
    #row.names(tpmValues)<-paste0(row.names(tpmValues),"_",mapping$GeneName)
    # tpmValues$GeneName<-mapping$GeneName
    tpmvaluesGene<-tpmValuesFull1[gene,]
    d<-reshape2::melt(tpmvaluesGene)
    f<-features[,c("RTNo","PrimaryCancerTissue","New.Diagnosis","TumorSite")]
    f$New.Diagnosis[f$New.Diagnosis == "Adrenocortical carcinoma-recurrence"]<-"Adrenocortical carcinoma-metastasis"
    f$samples<-row.names(f)
    # d<-d[d$GeneName == gene,]
    final_data <- sqldf::sqldf("SELECT * FROM d LEFT OUTER JOIN f where d.variable == samples")
    #print(dim(final_data))
    g<-ggplot(final_data,aes(x=New.Diagnosis,y=value))+
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 15),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
      ggplot2::theme(plot.title = element_text(hjust = 0.5))+
      geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
      ggplot2::geom_point(ggplot2::aes(label=TumorSite,label2=samples))+
      labs(x="", y = "TPM Values") +
      ggtitle(paste0("Expression of ",gene))+
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
      coord_flip()
    plotsLists[[gene]]<- g
    #ggsave(paste0(parent,gene,"-expBoxplot.pdf"),plot = g,height = 10,width = 8)
  }
  return(plotsLists)
}

getGeneEffectsInDepMap<-function(genes=NULL,depMapVersion=NULL,groupBySubtypes=FALSE){

  #genes<-c("CHEK1","MMP9","DCTN1","KLC4","HSP90AB1","CDKN1B","SOX9","CLN6","SYT8","NOS2")
  if (is.null(genes)) {
    stop("Need to enter the gene or gene list")
  }
  require(depmap)
  require(ExperimentHub)
  require(tidyverse)
  require(gridExtra)
  require(cowplot)
  eh <- ExperimentHub::ExperimentHub()
  dep<-AnnotationHub::query(eh, "depmap")
  ##Get Crispr dependency Scores
  ##Version 21Q1
  # pedDepMapDataInfo<-read.table("~/myPART/DepMapPedSampleInfo.txt",sep = "\t",header = T)
  # pedDepMapDataInfo <- pedDepMapDataInfo[pedDepMapDataInfo$Class_for_Manuscript == "Pediatric",]
  crisprData<-eh[["EH5358"]]#depmap::depmap_crispr()
  metaData<-eh[["EH5362"]]#depmap::depmap_metadata()
  geneExpDataDepMap<-eh[["EH5360"]]
  cnvDataDepMap<-eh[["EH5359"]]
  mutationCallsDepMap<-eh[["EH5361"]]

  #t<-data.frame(crisprData %>% dplyr::filter(depmap_id %in% pedDepMapDataInfo$DepMap_ID))
  t<-data.frame(crisprData)
  plotsLists<-list()
  for(gene in genes)
  {
    #gene<-"CHEK1"
    print(gene)
    #d$gene<-unlist(lapply(d$gene,FUN=function(x){return(unlist(str_split(x," "))[1])}))
    d <- t[t$gene_name %in% gene,]
    d1 <- merge(d,metaData,by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name,primary_disease,subtype_disease,dependency,cell_line.x,cell_line_name)
    d1$subtype_disease[is.na(d1$subtype_disease)]<-d1$primary_disease[is.na(d1$subtype_disease)]
    medianDep<-d1 %>% group_by(primary_disease) %>% summarise(Mean=mean(dependency),Median=median(dependency))

    d2 <- merge(d1,geneExpDataDepMap[geneExpDataDepMap$gene_name %in% gene,],by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name.x,primary_disease,subtype_disease,dependency,rna_expression,cell_line.x,cell_line_name)
    #print(paste0(genes,"-",round(min(medianDep$Median),digits = 2)))
    d3<- merge(d2,mutationCallsDepMap[mutationCallsDepMap$gene_name %in% gene,],by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name.x,primary_disease,subtype_disease,dependency,rna_expression,cell_line.x,cell_line_name,is_deleterious,is_tcga_hotspot,is_cosmic_hotspot)
    d3$color[d3$is_deleterious==TRUE]<-"blue"
    d3$color[d3$is_tcga_hotspot == TRUE]<-"red"
    d3$color[d3$is_cosmic_hotspot == TRUE]<-"red"
    d3$color[is.na(d3$color)]<-"black"
    d3$size<-ifelse(d3$rna_expression >=log(5+1),"high","low")
    colorName<-unique(d3$color)
    g<-ggplot2::ggplot(d3,ggplot2::aes(x=primary_disease,y=dependency))+
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(text=ggplot2::element_text(face = "bold"),axis.text = ggplot2::element_text(size = 10),
                     axis.title = ggplot2::element_text(size = 15,face = "bold"),legend.background = ggplot2::element_rect(colour = "black")) +
      ggplot2::geom_boxplot()+
      ggplot2::geom_point(ggplot2::aes(text=sprintf("CellLine: %s<br>Subtype: %s<br>Primary: %s<br>Log2(TPM+1): %s<br>Dependency: %s", cell_line_name, subtype_disease,primary_disease,rna_expression,dependency),
                                       color=color,size=size))+
      ggplot2::scale_size_manual(values=c(2.2,1)) +
      ggplot2::scale_color_manual(values=colorName) +
      ggplot2::guides(color="none") +
      ggplot2::theme(legend.position='none') +
      #ggplot2::geom_point(ggplot2::aes())+
      #ggplot2::geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=0.5, notch=FALSE)+
      ggplot2::labs(x="", y = "Normalized Dependency Score") +
      ggplot2::geom_hline(yintercept=-1, linetype="dashed", color = "red") +
      #geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.5) +
      ggplot2::coord_flip()
    ggly <- plotly::ggplotly(g,tooltip="text")
    # add hover info
    # hoverinfo <- with(d3, paste0("CellLine: ", cell_line_name, "</br></br>",
    #                              "Subtype: ", subtype_disease, "</br>",
    #                              "Primary: ", primary_disease, "</br>",
    #                              "Log2(TPM+1): ", rna_expression, "</br>",
    #                              "Dependency: ", dependency))
    # ggly$x$data[[1]]$text <- hoverinfo
    # ggly$x$data[[1]]$hoverinfo <- c("text", "boxes")

    g1<-ggplot2::ggplot(d1,aes(x=dependency)) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 12),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
      ggplot2::labs(x="", y = "") +
      ggplot2::geom_histogram() +
      ggplot2::geom_vline(xintercept=-1, linetype="dashed", color = "red")
    #p <- cowplot::plot_grid(g1, g, align = "v",ncol = 1,rel_heights=c(1,2))
    p<-plotly::subplot(plotly::ggplotly(g1),ggly,nrows = 2,heights = c(0.25,0.75),titleX=TRUE)
    #save_plot(paste0("/Users/jaina13/myPART/WGSData/tumor-only-somatic-mafs/DepMap/Demap-AllCellLines-",round(min(medianDep$Median),digits = 2),"-CRISPR-",gene,".pdf"), p,base_height = 10,base_width = 8)
    plotsLists[[gene]]<- list(ggHistogram=g1,ggBoxplot=ggly,combined=p)
  }
  return(plotsLists)
}

