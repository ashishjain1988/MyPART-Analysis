library(DESeq2)
library(apeglm)
library(edgeR)
library(WebGestaltR)

source("~/myPART/MyPART-Analysis/RScripts/helperFunctions.R")

parent<-"~/myPART/AllSamplesPipeliner/EndocrineSubgroupResults/"
l1<-load(file = paste0(parent,"vsdNormalizedCounts.rda"))
features<-featuresFilt

pca <- prcomp(t(as.matrix(vsdNormalizedCounts)),scale. = T)
eigs <- pca$sdev^2
var_explained <- paste(sprintf(eigs/sum(eigs)*100,fmt = '%.2f'),"%",sep = "")
pcaData<-pca$x
colnames(pcaData)<-paste(colnames(pcaData),rep("(",ncol(pcaData)),var_explained,rep(")",ncol(pcaData)),sep = "")
pcaData<-pcaData[,c(1:5)]

#palette <- distinctColorPalette(n)
###Dendrogram
library("dendextend")
library(randomcoloR)
pcaData<-pca$x[,c(1:5)]
group<-as.factor(features$New.Diagnosis)
n <- length(levels(group))
palette <- distinctColorPalette(n)
names(palette)<-levels(group)
groupCodes <-palette[features$New.Diagnosis]

clusters<-c(1:length(levels(group)))
names(clusters)<-levels(group)
clusterName<-clusters[features$New.Diagnosis]

row.names(pcaData)<-row.names(features)#features$New.Diagnosis##
dend <- as.dendrogram(hclust(dist(pcaData,method = "euclidean"),method = "complete"))#"ward.D2"))
labels_colors(dend) <- groupCodes[order.dendrogram(dend)]
#dend_iris <- color_branches(dend,col=groupCodes, clusters=clusterName)

pdf(paste0(parent,"/PCsSampleClustering-All-Dendrogram-Sample-complete-sample.pdf"),width = 14,height = 14)
#pdf(paste0(parent,"/PCsSampleClustering-complete.pdf"),width = 14,height = 14)
par(mar = c(2, 2, 2, 25))
plot(dend,horiz = TRUE)
dev.off()


# Principal Components Analysis: 5 PCs as heatmap input
#l<-load(file = paste0("/Users/jaina13/myPART/AllSamplesPipeliner/vsdNormalizedCounts.rda"))
#features<-read.table(paste0("~/myPART/features.txt"),header = TRUE,sep = "\t",row.names = 1)
NETSamples<-features[features$New.Diagnosis == "Neuroendocrine",]
features[features$New.Diagnosis == "Neuroendocrine",]$New.Diagnosis<-paste0(NETSamples$New.Diagnosis,"-",NETSamples$PrimaryCancerTissue)
#vsdNormalizedCounts<-vsdNormalizedCounts[,row.names(features)]
tableTypes<-table(features$New.Diagnosis)
filtNames<-c(names(tableTypes[tableTypes == 1]),"Neuroendocrine-small intestine","Neuroendocrine-lung","Neuroendocrine-pancreas","Adrenocortical carcinoma","Normal")
features<-features[features$New.Diagnosis %in% filtNames,]
vsdNormalizedCounts <- vsdNormalizedCounts[,row.names(features)]
rnaSamples<-row.names(features[features$Tissue == "Human Universal Reference Total RNA",])
lowQualitySamples<-intersect(row.names(features),c("47_RT00085_X007R_resent","76_RT00105_X336R"))
vsdNormalizedCounts <- vsdNormalizedCounts[,!(colnames(vsdNormalizedCounts) %in% c(rnaSamples,lowQualitySamples))]
features<-features[colnames(vsdNormalizedCounts),]
features$Tissue <- tools::toTitleCase(features$Tissue)
pca <- prcomp(t(as.matrix(vsdNormalizedCounts)),scale. = T)
eigs <- pca$sdev^2
var_explained <- paste(sprintf(eigs/sum(eigs)*100,fmt = '%.2f'),"%",sep = "")
pcaData<-pca$x
colnames(pcaData)<-paste(colnames(pcaData),rep("(",ncol(pcaData)),var_explained,rep(")",ncol(pcaData)),sep = "")

pca_exp<-pcaData[,1:5]
featuresFilt<-features#features[grepl("Normal|Gastrointestinal",features$New.Diagnosis),]
hm_data<-as.matrix(t(pca_exp[row.names(featuresFilt),]))
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(randomcoloR)
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
                                      col=colorRamp2(seq(-max(abs(pca_exp), na.rm = T), max(abs(pca_exp), na.rm = T), length.out = 20),
                                                     rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                      bottom_annotation = ha,#column_annotations,
                                      show_column_names=F,
                                      column_names_rot = 45,
                                      cluster_rows = FALSE,
                                      show_heatmap_legend = F) # Turning off to control the placement
pdf(paste0(parent,"EndocrineSubgroupResults/PCsSampleClustering-Subtypes-All-L.pdf"))
draw(cheatmap, show_annotation_legend = TRUE)
dev.off()

###Clinical Features Correlation Plot
require(ComplexHeatmap)
require(circlize)
f<-data.frame(Diagnosis=as.numeric(as.factor(features$New.Diagnosis)))
#f$Tissue<-as.numeric(as.factor(features$Tissue))
f$Sex<-as.numeric(as.factor(features$Sex))
f$Race<-as.numeric(as.factor(features$Race))
f$Age<-(features$Age)
#f$RIN<-features$RIN
#f<-features
cormatrix<-round(cor(pcaData[,1:5],f, method = "spearman",use = "na.or.complete"),2)
cheatmap <- ComplexHeatmap::Heatmap(t(cormatrix),
                                    col=colorRamp2(seq(-1, 1, length.out = 20),
                                                   rev(colorRampPalette(brewer.pal(9, "PuOr"))(20))),
                                    cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                      grid.text(t(cormatrix)[i, j], x, y)
                                    },
                                    bottom_annotation = NULL,
                                    show_column_names=T,
                                    column_names_rot = 45,
                                    cluster_rows = T,
                                    cluster_columns = F,
                                    show_heatmap_legend = T)
pdf(paste0(parent,"PCsFeaturesCorrelation.pdf"))
#pdf(outputFile)
draw(cheatmap, show_annotation_legend = TRUE)
dev.off()
