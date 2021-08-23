library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(randomcoloR)
library(ggplot2)
library(dplyr)

###########Figure 1A###################
source("~/myPART/MyPART-Analysis/RScripts/helperFunctions.R")
parent<-"~/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/"
l1<-load(file = paste0(parent,"vsdNormalizedCounts.rda"))
features<-featuresFilt
vsdNormalizedCounts <- vsdNormalizedCounts[,row.names(features)]
features$Tissue <- tools::toTitleCase(features$Tissue)
features$TumorSite <- tools::toTitleCase(features$TumorSite)
pca <- prcomp(t(as.matrix(vsdNormalizedCounts)),scale. = T)
eigs <- pca$sdev^2
var_explained <- paste(sprintf(eigs/sum(eigs)*100,fmt = '%.2f'),"%",sep = "")
pcaData<-pca$x
colnames(pcaData)<-paste(colnames(pcaData),rep("(",ncol(pcaData)),var_explained,rep(")",ncol(pcaData)),sep = "")

pca_exp<-pcaData[,1:5]
featuresFilt<-features#features[grepl("Normal|Gastrointestinal",features$New.Diagnosis),]
hm_data<-as.matrix(t(pca_exp[row.names(featuresFilt),]))
column_annotations = HeatmapAnnotation(df = featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")])
col_age = colorRamp2(c(0,max(featuresFilt$Age,na.rm = T)), c("white","blue"))
#col_rin = colorRamp2(c(0,max(featuresFilt$RIN)), c("white","yellow"))
n <- distinctColorPalette(length(levels(as.factor(featuresFilt$Tissue))))
names(n)<-(levels(as.factor(featuresFilt$Tissue)))
col_tissue <-n[featuresFilt$Tissue]

n1 <- distinctColorPalette(length(levels(as.factor(featuresFilt$New.Diagnosis))))
names(n1)<-(levels(as.factor(featuresFilt$New.Diagnosis)))
col_disease <-n1[featuresFilt$New.Diagnosis]

df<-featuresFilt[,c("Tissue","New.Diagnosis","Age","Sex","Race","TumorSite")]
colnames(df)<-c("Tissue","Diagnosis","Age","Sex","Race","TumorSite")
#col_tissue <- distinctColorPalette()
ha = HeatmapAnnotation(
  #df = features[,c(5,6,7,13,15)],
  #df = featuresFilt[,c(3,4,5,6,13,7)],
  df = df,
  col = list(Tissue = col_tissue,
             Diagnosis = col_disease,#c("Normal" = "black","Adrenocortical carcinoma" = "red"),#col_cluster,#c("Adrenocortical carcinoma" = "red", "Normal" = "black", "Carcinoid tumor" = "blue","Gastrointestinal stromal tumor"="brown"),
             Age = col_age,
             Sex = c("Male"="turquoise","Female"="brown"),
             #RIN = col_rin,
             Race =c("White"="#B35806","Unknown"="#FDBC6B","Other"="black","Black or African American"="#E3E4EF","Native Hawaiian or Other Pacific Islander"="#8D81B6","Asian"="blue"),
             TumorSite=c("Metastasis"="turquoise","Primary Site"="brown","Recurrence"="black")
             #cluster = col_cluster
  ),
  annotation_name_gp =  gpar(fontsize = 12),
  annotation_legend_param = list(title_gp = gpar(fontsize = 15, fontface = 2),labels_gp = gpar(fontsize = 12))
)

cheatmap <- ComplexHeatmap::Heatmap(hm_data,
                                    col=colorRamp2(seq(-max(abs(pca_exp), na.rm = T), max(abs(pca_exp), na.rm = T), length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    bottom_annotation = ha,#column_annotations,
                                    show_column_names=FALSE,
                                    #column_names_rot = 45,
                                    cluster_rows = FALSE,
                                    cluster_columns = TRUE,
                                    clustering_distance_columns = "euclidean",
                                    clustering_method_columns = "complete",
                                    show_heatmap_legend = F
                                    ) # Turning off to control the placement
pdf(paste0(parent,"PCsSampleClustering-Subtypes-All.pdf"),width = 14,height = 7)
ComplexHeatmap::draw(cheatmap, show_annotation_legend = TRUE)
dev.off()


###########Figure 1B###################
parent<-"/Users/jaina13/myPART/AllSamplesPipeliner/"
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
    c<-cellSpecificGenesOutput %>% dplyr::filter(Group %in% c("Tissue-Enriched")) %>% dplyr::filter(Tissue == tissue)
    c2<-sapply(c$Gene,FUN = function(x){
      return(paste0(drugTargetsMap[drugTargetsMap$gene_name %in% x,2],collapse = ","))
    })
    c$drugBank<-c2
    c<-c[c$Gene %in% proteinCodingGenes$Gene.name,]
    l[[str_replace_all(tissue,pattern = "[ /]",replacement = ".")]]<-c
  }
}
#o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
csGenes<-do.call(rbind,l)
csGenes$isDrug<-"No"
csGenes[csGenes$drugBank != "",]$isDrug<-"Yes"
#csGenes<-csGenes[csGenes$isDrug =="Yes",]
g<-ggplot(csGenes, aes(fill=isDrug, x=(Tissue))) +
  #g<-ggplot(csGenes, aes(x=(Tissue))) +
  geom_bar(position="stack") +
  scale_fill_manual(values = c("grey","blue")) +
  guides(fill=guide_legend(title="isDrugTarget")) +
  labs(x="", y = "# Cancer-Specific Genes") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.title = element_text(size=14,face = "bold"),legend.text = element_text(size=14),axis.title = element_text(size = 14,face = "bold"),axis.text = element_text(size = 16,face = "bold",colour = "black")) +
  theme(plot.margin = unit(c(1,1,1,2), "cm")) +
  #margin(t, r, l, b)
  xlab("")
ggsave(paste0(parent,"EndocrineSubgroupResults/cancerSpecificGenesBarplot-DrugTarget-AllSamples.pdf"),plot = g,height = 8,width = 10)

###########Figure 1C [Expression Boxplots]###################
o2<-load(file= paste0("/Users/jaina13/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/CancerSpecificGenes/cancerSpecificGenes-TPM.rda"))
#genes<-c("GNRHR","SPINK6","TPH1","CELA2A","KCNJ6","F5","GHRHR")
gene<-"TPH1"
tpmvaluesGene<-tpmValuesFull1[gene,]
d<-reshape2::melt(tpmvaluesGene)
f<-features[,c("RTNo","PrimaryCancerTissue","New.Diagnosis","TumorSite")]
f$New.Diagnosis[f$New.Diagnosis == "Adrenocortical carcinoma-recurrence"]<-"Adrenocortical carcinoma-metastasis"
f$samples<-row.names(f)
final_data <- sqldf::sqldf("SELECT * FROM d LEFT OUTER JOIN f where d.variable == samples")
g<-ggplot(final_data,aes(x=New.Diagnosis,y=value))+
  theme_classic(base_size = 10) +
  theme(text=element_text(face = "bold",colour = "black"),title = element_text(size = 16),axis.text = element_text(size = 15,colour = "black"),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.2)) +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
  geom_point(size=1.5)+
  guides(size="none") +
  labs(x="", y = "Transcripts Per Million (TPM)") +
  #ggtitle(paste0("Expression of ",gene))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 30))+
  theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
  coord_flip()
ggsave(paste0(parent,"EndocrineSubgroupResults/",gene,"-expBoxplot.pdf"),plot = g,height = 10,width = 12)

###########Figure 1D [DepMap Boxplots]###################
require(depmap)
require("ExperimentHub")
require(tidyverse)
require(gridExtra)
require(cowplot)
eh <- ExperimentHub()
dep<-query(eh, "depmap")
##Get Crispr dependency Scores
##Version 21Q1
#pedDepMapDataInfo<-read.table("~/myPART/DepMapPedSampleInfo.txt",sep = "\t",header = T)
#pedDepMapDataInfo <- pedDepMapDataInfo[pedDepMapDataInfo$Class_for_Manuscript == "Pediatric",]
crisprData<-eh[["EH5358"]]#depmap::depmap_crispr()
metaData<-eh[["EH5362"]]#depmap::depmap_metadata()

#genes<-hubGenes$gene[1:10]#c("STAR","CYP17A1","CYP21A2","FDXR","CYP11A1","POR","CYP11B1","NR5A1","FDX1","NR0B1")#c("AURKB","SMARCA4")
#genes<-c("FGF4","PLAT","CCKBR","KCNA1","KIT","GRIN3A")
#genes<-c("PLXNA2","TNRC6C" ,"LTBP1" , "PJA2" ,  "PTCH1"  ,"USO1")
gene<-"TPH1"
#t<-data.frame(crisprData %>% dplyr::filter(depmap_id %in% pedDepMapDataInfo$DepMap_ID))
t<-data.frame(crisprData)
plotsLists<-list()
#for(gene in genes)
{
  #gene<-"CHEK1"
  #d$gene<-unlist(lapply(d$gene,FUN=function(x){return(unlist(str_split(x," "))[1])}))
  d <- t[t$gene_name %in% gene,]
  d1 <- merge(d,metaData,by="depmap_id",all.x=TRUE,all.y=FALSE) %>% dplyr::select(depmap_id,gene_name,primary_disease,subtype_disease,dependency,cell_line.x)
  d1$subtype_disease[is.na(d1$subtype_disease)]<-d1$primary_disease[is.na(d1$subtype_disease)]
  d1$subtype_disease[d1$subtype_disease == "Undifferentiated"]<-paste0(d1$primary_disease[d1$subtype_disease == "Undifferentiated"],"(",d1$subtype_disease[d1$subtype_disease == "Undifferentiated"],")")
  medianDep<-d1 %>% group_by(subtype_disease) %>% summarise(Mean=mean(dependency),Median=median(dependency))
  #print(paste0(genes,"-",round(min(medianDep$Median),digits = 2)))
  g<-ggplot2::ggplot(d1,aes(x=primary_disease,y=dependency))+
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 15,colour = "black"),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
    ggplot2::geom_boxplot(outlier.colour="black", outlier.shape=16,outlier.size=1, notch=FALSE)+
    #ggplot2::geom_point(size=1.5)+
    #ggplot2::guides(size="none") +
    ggplot2::labs(x="", y = "Normalized Dependency Score") +
    ggplot2::geom_hline(yintercept=-1, linetype="dashed", color = "red") +
    theme(plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
    #geom_dotplot(binaxis='y', stackdir='center',dotsize = 0.5) +
    coord_flip()
  # g1<-ggplot2::ggplot(d1,aes(x=dependency)) +
  #   ggplot2::theme_classic(base_size = 10) +
  #   ggplot2::theme(text=element_text(face = "bold"),axis.text = element_text(size = 12),axis.title = element_text(size = 15,face = "bold"),legend.background = element_rect(colour = "black")) +
  #   ggplot2::labs(x="", y = "") +
  #   ggplot2::geom_histogram() +
  #   ggplot2::geom_vline(xintercept=-1, linetype="dashed", color = "red")
  # p <- cowplot::plot_grid(g1, g, align = "v",ncol = 1,rel_heights=c(1,2))
  #save_plot(paste0(parent,"Demap-PedCellLines-",round(min(medianDep$Median),digits = 2),"-CRISPR-",gene,".pdf"), p,base_height = 10,base_width = 8)
  #plotsLists[[gene]]<- list(ggHistogram=g1,ggBoxplot=g,combined=p)
  plotsLists[[gene]]<- list(ggBoxplot=g)
}
ggsave(paste0(parent,"EndocrineSubgroupResults/",gene,"-DepMapBoxplot.pdf"),plot = plotsLists[[gene]]$ggBoxplot,height = 10,width = 12)

